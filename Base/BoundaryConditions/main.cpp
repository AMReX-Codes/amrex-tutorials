#include <AMReX.H>
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <utility>

using amrex::Real;
using namespace amrex::literals;

// Fills ghost cells outside the domain on BCType::ext_dir faces with the
// given Dirichlet values. AMReX fills the other BC types (foextrap,
// reflect_even, ...) itself; periodic ghost cells are filled by FillBoundary.
struct DirichletFill
{
    amrex::GpuArray<Real,AMREX_SPACEDIM> phi_lo;
    amrex::GpuArray<Real,AMREX_SPACEDIM> phi_hi;

    AMREX_GPU_DEVICE
    void operator() (amrex::IntVect const& iv, amrex::Array4<Real> const& dest,
                     int dcomp, int numcomp, amrex::GeometryData const& geom,
                     Real /*time*/, amrex::BCRec const* bcr, int bcomp,
                     int /*orig_comp*/) const
    {
        amrex::Box const& domain = geom.Domain();
        for (int n = 0; n < numcomp; ++n) {
            amrex::BCRec const& bc = bcr[bcomp+n];
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (bc.lo(d) == amrex::BCType::ext_dir && iv[d] < domain.smallEnd(d)) {
                    dest(iv, dcomp+n) = phi_lo[d];
                } else if (bc.hi(d) == amrex::BCType::ext_dir && iv[d] > domain.bigEnd(d)) {
                    dest(iv, dcomp+n) = phi_hi[d];
                }
            }
        }
    }
};

// Laplacian of phi at cell iv. For ext_dir, the ghost cell holds the value
// on the boundary face, which is only dx/2 away from the cell center.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real laplacian (amrex::IntVect const& iv, amrex::Array4<Real const> const& phi,
                amrex::BCRec const& bc, amrex::Box const& domain,
                amrex::GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
    Real lap = 0.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        auto const e = amrex::IntVect::TheDimensionVector(d);
        Real const flo = (bc.lo(d) == amrex::BCType::ext_dir && iv[d] == domain.smallEnd(d))
            ? 2.0_rt : 1.0_rt;
        Real const fhi = (bc.hi(d) == amrex::BCType::ext_dir && iv[d] == domain.bigEnd(d))
            ? 2.0_rt : 1.0_rt;
        lap += (fhi * (phi(iv+e) - phi(iv)) - flo * (phi(iv) - phi(iv-e))) * dxinv[d] * dxinv[d];
    }
    return lap;
}

// Exact solution: linear steady state plus a decaying sine mode in x.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_solution (Real x, Real t, Real phi_lo, Real phi_hi, Real diffusion_coef)
{
    Real const pi = amrex::Math::pi<Real>();
    return phi_lo + (phi_hi - phi_lo) * x
        + std::exp(-pi * pi * diffusion_coef * t) * std::sin(pi * x);
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        Real diffusion_coef = 1.0_rt;
        Real cfl = 0.9_rt;
        Real stop_time = 0.05_rt;
        Real phi_lo = 1.0_rt; // Dirichlet value at x = 0
        Real phi_hi = 0.0_rt; // Dirichlet value at x = 1
        int plot_int = -1;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("diffusion_coef", diffusion_coef);
            pp.query("cfl", cfl);
            pp.query("stop_time", stop_time);
            pp.query("phi_lo", phi_lo);
            pp.query("phi_hi", phi_hi);
            pp.query("plot_int", plot_int);
        }

        // Boundary conditions, one BCRec per component:
        //   x: Dirichlet (ext_dir) on both sides
        //   y: zero gradient, as even reflection (lo) and extrapolation (hi)
        //   z: periodic (int_dir), which must match the Geometry
        amrex::Vector<amrex::BCRec> bcs(1);
        bcs[0].setLo(0, amrex::BCType::ext_dir);
        bcs[0].setHi(0, amrex::BCType::ext_dir);
#if (AMREX_SPACEDIM >= 2)
        bcs[0].setLo(1, amrex::BCType::reflect_even);
        bcs[0].setHi(1, amrex::BCType::foextrap);
#endif
#if (AMREX_SPACEDIM == 3)
        bcs[0].setLo(2, amrex::BCType::int_dir);
        bcs[0].setHi(2, amrex::BCType::int_dir);
#endif
        amrex::Array<int,AMREX_SPACEDIM> is_periodic{};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            is_periodic[d] = (bcs[0].lo(d) == amrex::BCType::int_dir) ? 1 : 0;
        }

        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        amrex::DistributionMapping dm(ba);
        amrex::Geometry geom(domain,
                             amrex::RealBox({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                            {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}),
                             amrex::CoordSys::cartesian, is_periodic);
        auto const prob_lo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();
        auto const dxinv = geom.InvCellSizeArray();

        // Physical boundary filler: AMReX handles foextrap, reflect_even, etc.
        // and calls DirichletFill for every ghost cell outside the domain.
        DirichletFill dirichlet_fill{};
        dirichlet_fill.phi_lo[0] = phi_lo;
        dirichlet_fill.phi_hi[0] = phi_hi;
        amrex::PhysBCFunct<amrex::GpuBndryFuncFab<DirichletFill>> physbc
            (geom, bcs, amrex::GpuBndryFuncFab<DirichletFill>(dirichlet_fill));

        amrex::MultiFab phi_old(ba, dm, 1, 1);
        amrex::MultiFab phi_new(ba, dm, 1, 1);

        auto const& init = phi_new.arrays();
        amrex::ParallelFor(phi_new, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            init[b](i,j,k) = exact_solution(x, 0.0_rt, phi_lo, phi_hi, diffusion_coef);
        });

        Real const dt = cfl * dx[0] * dx[0] / (2.0_rt * AMREX_SPACEDIM * diffusion_coef);
        Real time = 0.0_rt;
        int step = 0;
        amrex::BCRec const bc = bcs[0];

        while (time < stop_time * (1.0_rt - 1.e-6_rt))
        {
            std::swap(phi_old, phi_new);
            Real const dt_step = std::min(dt, stop_time - time);

            // Ghost cells: first those inside the (periodic) domain, then
            // those outside the domain from the physical BCs.
            phi_old.FillBoundary(geom.periodicity());
            physbc(phi_old, 0, 1, phi_old.nGrowVect(), time, 0);

            Real const fac = dt_step * diffusion_coef;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(phi_new, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box const& bx = mfi.tilebox();
                auto const& po = phi_old.const_array(mfi);
                auto const& pn = phi_new.array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (amrex::IntVect const& iv)
                {
                    pn(iv) = po(iv) + fac * laplacian(iv, po, bc, domain, dxinv);
                });
            }

            time += dt_step;
            ++step;
            if (plot_int > 0 && step % plot_int == 0) {
                amrex::WriteSingleLevelPlotfile(amrex::Concatenate("plt", step, 5), phi_new,
                                                {"phi"}, geom, time, step);
            }
        }

        amrex::MultiFab error(ba, dm, 1, 0);
        auto const& pa = phi_new.const_arrays();
        auto const& ea = error.arrays();
        amrex::ParallelFor(error, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            ea[b](i,j,k) = pa[b](i,j,k) - exact_solution(x, time, phi_lo, phi_hi, diffusion_coef);
        });
        Real const err = error.norm0();
        amrex::Print() << "Finished " << step << " steps at time " << time
                       << ", max error = " << err << "\n";

        if (err > 1.e-3_rt) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
