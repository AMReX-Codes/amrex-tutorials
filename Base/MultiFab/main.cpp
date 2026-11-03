#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <type_traits>

// 1 + prod_d sin(2 pi x_d) at the center of cell iv; periodic on [0,1]^d.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real phi_value (amrex::IntVect const& iv,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    using namespace amrex::literals;
    amrex::Real r = 1.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        amrex::Real x = prob_lo[d] + (iv[d] + 0.5_rt) * dx[d];
        r *= std::sin(2.0_rt * amrex::Math::pi<amrex::Real>() * x);
    }
    return 1.0_rt + r;
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
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
        }

        // Box: the index space of the domain, cells 0 to n_cell-1 in each direction.
        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));

        // BoxArray: the domain chopped into boxes of at most max_grid_size cells
        // in each direction.
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);

        // DistributionMapping: which MPI rank owns which box.
        amrex::DistributionMapping dm(ba);

        // Geometry: the physical domain [0,1]^d, periodic in all directions.
        amrex::RealBox real_box({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
        amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian,
                             {AMREX_D_DECL(1, 1, 1)});

        amrex::Print() << "Domain " << domain << " split into " << ba.size()
                       << " boxes on " << amrex::ParallelDescriptor::NProcs()
                       << " MPI rank(s)\n";

        // A MultiFab with 2 components and 1 layer of ghost cells. Each rank
        // allocates data (an FArrayBox) only for the boxes it owns.
        int const ncomp = 2;
        int const nghost = 1;
        amrex::MultiFab mf(ba, dm, ncomp, nghost);
        mf.setVal(0.0_rt);

        auto const prob_lo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();

        // Method 1: ParallelFor over the whole MultiFab. On GPU, all boxes
        // are handled by a single kernel launch, which is best for many
        // small boxes. On CPU, it is tiled and threaded with OpenMP.
        auto const& ma = mf.arrays();
        amrex::ParallelFor(mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            ma[b](i,j,k,0) = phi_value(amrex::IntVect(AMREX_D_DECL(i,j,k)), prob_lo, dx);
        });

        // Method 2: MFIter loop over the boxes (tiles on CPU) with an Array4
        // per box. Use it when per-box work is needed, e.g., temporary
        // FArrayBoxes or several kernels on the same box.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box const& bx = mfi.tilebox();
            amrex::Array4<Real> const& a = mf.array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                a(i,j,k,1) = phi_value(amrex::IntVect(AMREX_D_DECL(i,j,k)), prob_lo, dx);
            });
        }

        // Fill ghost cells from neighboring boxes, including periodic images.
        mf.FillBoundary(geom.periodicity());

        // Reference: the function evaluated directly in valid and ghost cells.
        amrex::MultiFab exact(ba, dm, 1, nghost);
        auto const& ea = exact.arrays();
        amrex::ParallelFor(exact, exact.nGrowVect(),
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            ea[b](i,j,k) = phi_value(amrex::IntVect(AMREX_D_DECL(i,j,k)), prob_lo, dx);
        });

        // Built-in reductions over all ranks.
        Real const cell_volume = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
        Real const integral = mf.sum(0) * cell_volume;
        amrex::Print() << "integral = " << integral << " (exact: 1)\n"
                       << "max      = " << mf.max(0) << "\n"
                       << "min      = " << mf.min(0) << "\n"
                       << "norm0    = " << mf.norm0(0) << "\n"
                       << "norm1    = " << mf.norm1(0) << "\n"
                       << "norm2    = " << mf.norm2(0) << "\n";

        // Both components, ghost cells included, must match the reference.
        amrex::MultiFab diff(ba, dm, 1, nghost);
        Real err = 0.0_rt;
        for (int n = 0; n < ncomp; ++n) {
            amrex::MultiFab::Copy(diff, mf, n, 0, 1, nghost);
            amrex::MultiFab::Subtract(diff, exact, 0, 0, 1, nghost);
            err = std::max(err, diff.norm0(0, nghost));
        }
        amrex::Print() << "max error including ghost cells = " << err << "\n";

        amrex::WriteSingleLevelPlotfile("plt_multifab", mf, {"phi_parallelfor", "phi_mfiter"},
                                        geom, 0.0_rt, 0);

        constexpr bool is_float = std::is_same_v<Real,float>;
        if (err > (is_float ? 1.e-5_rt : 1.e-12_rt) ||
            std::abs(integral - 1.0_rt) > (is_float ? 1.e-3_rt : 1.e-12_rt)) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
