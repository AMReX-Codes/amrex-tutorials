#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <utility>

// prod_d sin(2 pi x_d) at the center of cell iv
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real sine_mode (amrex::IntVect const& iv,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& prob_lo,
                       amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
{
    using namespace amrex::literals;
    amrex::Real r = 1.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        amrex::Real x = prob_lo[d] + (iv[d] + 0.5_rt) * dx[d];
        r *= std::sin(2.0_rt * amrex::Math::pi<amrex::Real>() * x);
    }
    return r;
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
        Real diffusion_coef = 1.0_rt;
        Real cfl = 0.9_rt;
        Real stop_time = 0.01_rt;
        int plot_int = -1; // plotfile interval in steps; < 0 for no plotfiles
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("diffusion_coef", diffusion_coef);
            pp.query("cfl", cfl);
            pp.query("stop_time", stop_time);
            pp.query("plot_int", plot_int);
        }

        // Periodic unit domain
        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        amrex::DistributionMapping dm(ba);
        amrex::Geometry geom(domain,
                             amrex::RealBox({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                            {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}),
                             amrex::CoordSys::cartesian, {AMREX_D_DECL(1, 1, 1)});
        auto const prob_lo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();
        auto const dxinv = geom.InvCellSizeArray();

        // One ghost cell for the 2*AMREX_SPACEDIM+1 point stencil
        amrex::MultiFab phi_old(ba, dm, 1, 1);
        amrex::MultiFab phi_new(ba, dm, 1, 1);

        // Initial condition: phi = 1 + prod_d sin(2 pi x_d)
        auto const& init = phi_new.arrays();
        amrex::ParallelFor(phi_new, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            init[b](i,j,k) = 1.0_rt + sine_mode(amrex::IntVect(AMREX_D_DECL(i,j,k)), prob_lo, dx);
        });

        // Explicit time step limit: dt <= dx^2 / (2 * dim * D)
        Real const dt = cfl * dx[0] * dx[0] / (2.0_rt * AMREX_SPACEDIM * diffusion_coef);

        Real time = 0.0_rt;
        int step = 0;
        if (plot_int > 0) {
            amrex::WriteSingleLevelPlotfile(amrex::Concatenate("plt", step, 5), phi_new,
                                            {"phi"}, geom, time, step);
        }

        while (time < stop_time * (1.0_rt - 1.e-6_rt))
        {
            std::swap(phi_old, phi_new);
            Real const dt_step = std::min(dt, stop_time - time);

            // Fill ghost cells, here all from periodic images and neighbors.
            phi_old.FillBoundary(geom.periodicity());

            // phi_new = phi_old + dt * D * Laplacian(phi_old)
            auto const& pold = phi_old.const_arrays();
            auto const& pnew = phi_new.arrays();
            Real const fac = dt_step * diffusion_coef;
            amrex::ParallelFor(phi_new, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto const& po = pold[b];
                pnew[b](i,j,k) = po(i,j,k) + fac * (
                      AMREX_D_TERM((po(i+1,j,k) - 2.0_rt*po(i,j,k) + po(i-1,j,k)) * dxinv[0]*dxinv[0],
                                 + (po(i,j+1,k) - 2.0_rt*po(i,j,k) + po(i,j-1,k)) * dxinv[1]*dxinv[1],
                                 + (po(i,j,k+1) - 2.0_rt*po(i,j,k) + po(i,j,k-1)) * dxinv[2]*dxinv[2]));
            });

            time += dt_step;
            ++step;
            if (step % 10 == 0) {
                amrex::Print() << "Step " << step << ", time = " << time << "\n";
            }
            if (plot_int > 0 && step % plot_int == 0) {
                amrex::WriteSingleLevelPlotfile(amrex::Concatenate("plt", step, 5), phi_new,
                                                {"phi"}, geom, time, step);
            }
        }

        // Exact solution: each sine mode decays as exp(-D |k|^2 t), |k|^2 = dim (2 pi)^2.
        Real const decay = std::exp(-diffusion_coef * AMREX_SPACEDIM * 4.0_rt
                                    * amrex::Math::pi<Real>() * amrex::Math::pi<Real>() * time);
        amrex::MultiFab error(ba, dm, 1, 0);
        auto const& pa = phi_new.const_arrays();
        auto const& ea = error.arrays();
        amrex::ParallelFor(error, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            Real exact = 1.0_rt + decay * sine_mode(amrex::IntVect(AMREX_D_DECL(i,j,k)), prob_lo, dx);
            ea[b](i,j,k) = pa[b](i,j,k) - exact;
        });
        Real const err = error.norm0() / decay; // relative to the mode amplitude
        amrex::Print() << "Finished " << step << " steps at time " << time
                       << ", mode amplitude = " << decay
                       << ", relative max error = " << err << "\n";

        if (err > 0.02_rt) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
