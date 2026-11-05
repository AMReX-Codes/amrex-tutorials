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
    // ((m % n) + n) % n maps m into [0, n) also for negative m
    i = ((i % n) + n) % n;
    j = ((j % n) + n) % n;
    k = ((k % n) + n) % n;
    return amrex::Real(i + n * (j + n * k));
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

        // FillBoundary: copy valid data into ghost cells, from neighboring
        // boxes and, with geom.periodicity(), across periodic boundaries.
        amrex::MultiFab a(ba, dm, 1, 2);
        a.setVal(-1.0_rt); // valid and ghost cells
        auto const& aa = a.arrays();
        amrex::ParallelFor(a, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            aa[b](i,j,k) = cell_id(i,j,k,n);
        });

        // The expected result: cell_id in the valid and all ghost cells.
        amrex::MultiFab ref(ba, dm, 1, 2);
        auto const& ra = ref.arrays();
        amrex::ParallelFor(ref, ref.nGrowVect(), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            ra[b](i,j,k) = cell_id(i,j,k,n);
        });
        amrex::MultiFab err(ba, dm, 1, 2);

        // Fill only the first of the two ghost layers.
        a.FillBoundary(amrex::IntVect(1), geom.periodicity());
        amrex::MultiFab::Copy(err, a, 0, 0, 1, 2);
        amrex::MultiFab::Subtract(err, ref, 0, 0, 1, 2);
        amrex::Print() << "FillBoundary, first ghost layer: max error = " << err.norm0(0, 1) << "\n";

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
        amrex::MultiFab::Copy(err, a, 0, 0, 1, 2);
        amrex::MultiFab::Subtract(err, ref, 0, 0, 1, 2);
        Real const err_fillboundary = err.norm0(0, 2);
        amrex::Print() << "FillBoundary, both ghost layers: max error = " << err_fillboundary << "\n";

        // ParallelCopy: copy between MultiFabs with different BoxArrays and
        // DistributionMappings, e.g., after regridding. Here the source is
        // the valid cells of a (0 ghost cells), and the destination is the
        // valid and 1 ghost cell of dst. With geom.periodicity(), cells
        // outside the domain are filled from their periodic images.
        dst.ParallelCopy(a, 0, 0, 1, amrex::IntVect(0), amrex::IntVect(1), geom.periodicity());

        amrex::MultiFab ref_dst(ba_dst, dm_dst, 1, 1);
        auto const& rda = ref_dst.arrays();
        amrex::ParallelFor(ref_dst, ref_dst.nGrowVect(), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            rda[b](i,j,k) = cell_id(i,j,k,n);
        });
        amrex::MultiFab err_dst(ba_dst, dm_dst, 1, 1);
        amrex::MultiFab::Copy(err_dst, dst, 0, 0, 1, 1);
        amrex::MultiFab::Subtract(err_dst, ref_dst, 0, 0, 1, 1);
        amrex::Print() << "ParallelCopy: max error = " << err_dst.norm0(0, 1) << "\n";

        // ParallelAdd adds instead of copying, so dst is now 2*cell_id.
        dst.ParallelAdd(a, 0, 0, 1, amrex::IntVect(0), amrex::IntVect(1), geom.periodicity());
        amrex::MultiFab::Copy(err_dst, dst, 0, 0, 1, 1);
        amrex::MultiFab::Subtract(err_dst, ref_dst, 0, 0, 1, 1);
        amrex::MultiFab::Subtract(err_dst, ref_dst, 0, 0, 1, 1);
        amrex::Print() << "ParallelAdd: max error = " << err_dst.norm0(0, 1) << "\n";

        // SumBoundary is the reverse of FillBoundary: values in ghost cells
        // are added to the valid cells they overlap. This is what deposition
        // needs, e.g., of particles near a box boundary. Here every valid
        // cell deposits 1 into the 3^dim cells around it, itself included,
        // some of them ghost cells. Each cell of a box computes how much it
        // receives from the valid cells of the same box, which avoids races.
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
        Real const err_sumboundary = std::max(c.max(0) - expected, expected - c.min(0));
        amrex::Print() << "SumBoundary: max error = " << err_sumboundary << "\n";

        if (err_fillboundary == 0.0_rt && err_sumboundary == 0.0_rt) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
