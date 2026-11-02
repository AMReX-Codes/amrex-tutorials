#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>

// A unique number for each cell of the periodic domain: ghost cells
// outside the domain get the number of their periodic image.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cell_id (int i, int j, int k, int n)
{
    auto wrap = [n] (int m) { return ((m % n) + n) % n; };
    return amrex::Real(wrap(i) + n * (wrap(j) + n * wrap(k)));
}

// Max over valid and the first ng ghost cells of |mf - f(i,j,k)|
template <typename F>
amrex::Real max_error (amrex::MultiFab const& mf, int ng, F const& f)
{
    amrex::MultiFab diff(mf.boxArray(), mf.DistributionMap(), 1, ng);
    auto const& da = diff.arrays();
    auto const& ma = mf.const_arrays();
    amrex::ParallelFor(diff, amrex::IntVect(ng),
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        da[b](i,j,k) = ma[b](i,j,k) - f(i,j,k);
    });
    return diff.norm0(0, ng);
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        using amrex::Real;
        using namespace amrex::literals;

        int n_cell = 32;
        int max_grid_size = 16;
        int dst_max_grid_size = 12;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("dst_max_grid_size", dst_max_grid_size);
        }

        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        amrex::DistributionMapping dm(ba);
        amrex::Geometry geom(domain,
                             amrex::RealBox({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                            {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}),
                             amrex::CoordSys::cartesian, {AMREX_D_DECL(1, 1, 1)});

        int const n = n_cell;
        auto id = [=] AMREX_GPU_DEVICE (int i, int j, int k) { return cell_id(i,j,k,n); };
        auto twice_id = [=] AMREX_GPU_DEVICE (int i, int j, int k) { return 2*cell_id(i,j,k,n); };
        bool passed = true;
        auto report = [&] (char const* what, Real err) {
            amrex::Print() << what << ": max error = " << err << "\n";
            passed = passed && (err == 0.0_rt);
        };

        // FillBoundary: copy valid data into ghost cells, from neighboring
        // boxes and, with geom.periodicity(), across periodic boundaries.
        amrex::MultiFab a(ba, dm, 1, 2);
        a.setVal(-1.0_rt); // valid and ghost cells
        auto const& aa = a.arrays();
        amrex::ParallelFor(a, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            aa[b](i,j,k) = cell_id(i,j,k,n);
        });

        // Fill only the first of the two ghost layers.
        a.FillBoundary(amrex::IntVect(1), geom.periodicity());
        report("FillBoundary, first ghost layer", max_error(a, 1, id));
        bool const untouched = a.min(0, 2) == -1.0_rt; // the second layer still has -1
        amrex::Print() << "FillBoundary, second ghost layer untouched: "
                       << (untouched ? "yes" : "no") << "\n";
        passed = passed && untouched;

        // The non-blocking version starts the communication and returns.
        // Work that does not need the ghost cells of a can be done before
        // FillBoundary_finish, overlapping it with the communication. Here,
        // the destination of the next section is set up.
        a.FillBoundary_nowait(geom.periodicity());

        // A BoxArray shifted by n_cell/4 and chopped differently, so that it
        // partly lies outside the domain, with its own DistributionMapping.
        amrex::BoxArray ba_dst(amrex::shift(domain, amrex::IntVect(n_cell/4)));
        ba_dst.maxSize(dst_max_grid_size);
        amrex::DistributionMapping dm_dst(ba_dst);
        amrex::MultiFab dst(ba_dst, dm_dst, 1, 1);
        dst.setVal(-1.0_rt);

        a.FillBoundary_finish();
        report("FillBoundary, both ghost layers", max_error(a, 2, id));

        // ParallelCopy: copy between MultiFabs with different BoxArrays and
        // DistributionMappings, e.g., after regridding. Here the source is
        // the valid cells of a (0 ghost cells), and the destination is the
        // valid and 1 ghost cell of dst. With geom.periodicity(), cells
        // outside the domain are filled from their periodic images.
        dst.ParallelCopy(a, 0, 0, 1, amrex::IntVect(0), amrex::IntVect(1), geom.periodicity());
        report("ParallelCopy", max_error(dst, 1, id));

        // ParallelAdd adds instead of copying.
        dst.ParallelAdd(a, 0, 0, 1, amrex::IntVect(0), amrex::IntVect(1), geom.periodicity());
        report("ParallelAdd", max_error(dst, 1, twice_id));

        // SumBoundary is the reverse of FillBoundary: values in ghost cells
        // are added to the valid cells they overlap. This is what deposition
        // needs, e.g., of particles near a box boundary. Here every valid
        // cell deposits 1 into itself and its 3^dim neighbors, including
        // ghost cells. Each cell of a box computes how much it receives from
        // the valid cells of the same box, which avoids races.
        amrex::MultiFab c(ba, dm, 1, 1);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(c); mfi.isValid(); ++mfi)
        {
            amrex::Box const& vbx = mfi.validbox();
            amrex::Array4<Real> const& ca = c.array(mfi);
            amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                int count = 0;
                for (int kk = k-AMREX_D_PICK(0,0,1); kk <= k+AMREX_D_PICK(0,0,1); ++kk) {
                for (int jj = j-AMREX_D_PICK(0,1,1); jj <= j+AMREX_D_PICK(0,1,1); ++jj) {
                for (int ii = i-1;                    ii <= i+1;                    ++ii) {
                    if (vbx.contains(amrex::IntVect(AMREX_D_DECL(ii,jj,kk)))) { ++count; }
                }}}
                ca(i,j,k) = Real(count);
            });
        }

        c.SumBoundary(geom.periodicity());
        // Now every valid cell has received 1 from each of its neighbors.
        // The ghost cells are not updated; call FillBoundary if needed.
        Real const expected = AMREX_D_TERM(3.0_rt, *3.0_rt, *3.0_rt);
        report("SumBoundary", std::max(c.max(0) - expected, expected - c.min(0)));

        if (passed) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
