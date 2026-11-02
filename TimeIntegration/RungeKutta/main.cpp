#include <AMReX.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_TimeIntegrator.H>

#include <cmath>
#include <iomanip>
#include <limits>
#include <string>

using namespace amrex;

namespace {

struct Problem
{
    Real kappa = 0.01_rt;  // diffusion coefficient
    GpuArray<Real,AMREX_SPACEDIM> vel{AMREX_D_DECL(1.0_rt, 0.5_rt, 0.25_rt)};
};

// rhs = -vel.grad(phi) + kappa*lap(phi) with second-order centered differences
void compute_rhs (MultiFab& rhs, MultiFab& phi, Geometry const& geom, Problem const& prob)
{
    // The integrator passes its stage data here; fill its ghost cells first.
    phi.FillBoundary(geom.periodicity());

    auto const dxi = geom.InvCellSizeArray();
    auto const vel = prob.vel;
    auto const kappa = prob.kappa;
    auto const& r = rhs.arrays();
    auto const& p = phi.const_arrays();
    ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        auto const& a = p[b];
        Real const adv = AMREX_D_TERM(
              vel[0]*(a(i+1,j,k) - a(i-1,j,k))*(0.5_rt*dxi[0]),
            + vel[1]*(a(i,j+1,k) - a(i,j-1,k))*(0.5_rt*dxi[1]),
            + vel[2]*(a(i,j,k+1) - a(i,j,k-1))*(0.5_rt*dxi[2]));
        Real const lap = AMREX_D_TERM(
              (a(i+1,j,k) - 2._rt*a(i,j,k) + a(i-1,j,k))*(dxi[0]*dxi[0]),
            + (a(i,j+1,k) - 2._rt*a(i,j,k) + a(i,j-1,k))*(dxi[1]*dxi[1]),
            + (a(i,j,k+1) - 2._rt*a(i,j,k) + a(i,j,k-1))*(dxi[2]*dxi[2]));
        r[b](i,j,k) = -adv + kappa*lap;
    });
}

// Exact solution of the spatially discretized system for phi(t=0) = prod_d sin(2 pi x_d).
// Each 1D factor is advected with a modified speed and decays at a modified
// rate, so the error measured against it is the time integration error only.
void exact_solution (MultiFab& phi, Geometry const& geom, Problem const& prob, Real time)
{
    constexpr Real kw = 2._rt*Math::pi<Real>();
    auto const dx = geom.CellSizeArray();
    auto const problo = geom.ProbLoArray();
    GpuArray<Real,AMREX_SPACEDIM> shift{};
    Real decay = 1._rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real const s = std::sin(0.5_rt*kw*dx[d]);
        shift[d] = prob.vel[d]*std::sin(kw*dx[d])/(kw*dx[d])*time;
        decay *= std::exp(-4._rt*prob.kappa*s*s/(dx[d]*dx[d])*time);
    }
    auto const& a = phi.arrays();
    ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        IntVect const iv(AMREX_D_DECL(i,j,k));
        Real v = decay;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real const x = problo[d] + (iv[d] + 0.5_rt)*dx[d];
            v *= std::sin(kw*(x - shift[d]));
        }
        a[b](i,j,k) = v;
    });
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 16;
        int max_grid_size = 8;
        Real stop_time = 0.25_rt;
        Real dt = 0.01_rt;  // time step of the coarsest run
        int nruns = 3;      // each run halves the time step
        Problem prob;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("stop_time", stop_time);
            pp.query("dt", dt);
            pp.query("nruns", nruns);
            pp.query("kappa", prob.kappa);
            Vector<Real> vel;
            if (pp.queryarr("velocity", vel)) {
                AMREX_ALWAYS_ASSERT(vel.size() >= AMREX_SPACEDIM);
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { prob.vel[d] = vel[d]; }
            }
        }

        // The TimeIntegrator reads its method from ParmParse:
        //   integration.type    = ForwardEuler | RungeKutta | SUNDIALS
        //   integration.rk.type = 1 (forward Euler), 2 (trapezoid), 3 (SSPRK3), 4 (RK4),
        //                         or 0 (user Butcher tableau)
        // queryAdd stores our defaults when they are not given in the inputs.
        std::string int_type = "RungeKutta";
        int rk_type = 4;
        ParmParse("integration").queryAdd("type", int_type);
        ParmParse("integration.rk").queryAdd("type", rk_type);

        // Formal order of accuracy of the selected method (0: unknown)
        int order = 0;
        if (int_type == "ForwardEuler" || int_type == "0") {
            order = 1;
        } else if ((int_type == "RungeKutta" || int_type == "1") && rk_type >= 1 && rk_type <= 4) {
            order = rk_type;
        }

        Box const domain(IntVect(0), IntVect(n_cell-1));
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);
        Geometry geom(domain, RealBox(AMREX_D_DECL(0._rt,0._rt,0._rt),
                                      AMREX_D_DECL(1._rt,1._rt,1._rt)),
                      CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});

        MultiFab phi(ba, dm, 1, 1);  // one ghost cell for the stencil
        MultiFab exact(ba, dm, 1, 0);
        exact_solution(exact, geom, prob, stop_time);

        amrex::Print() << "integration.type = " << int_type;
        if (int_type == "RungeKutta" || int_type == "1") {
            amrex::Print() << ", integration.rk.type = " << rk_type;
        }
        amrex::Print() << "\n\n        dt     steps  RHS evals      max error   rate\n";

        Real err_prev = 0;
        for (int irun = 0; irun < nruns; ++irun)
        {
            Real const dt_run = dt / Real(1 << irun);

            Real time = 0;
            exact_solution(phi, geom, prob, time);

            // The integrator allocates its stage data like phi.
            TimeIntegrator<MultiFab> integrator(phi, time);

            int nrhs = 0;
            integrator.set_rhs([&] (MultiFab& rhs, MultiFab& state, Real /*t*/)
            {
                ++nrhs;
                compute_rhs(rhs, state, geom, prob);
            });

            // Called after each step (set_post_stage_action: after each stage),
            // e.g., for applying limiters.  Here we count the steps.
            int nsteps = 0;
            integrator.set_post_step_action([&] (MultiFab& /*state*/, Real /*t*/) { ++nsteps; });

            integrator.set_time_step(dt_run);
            integrator.set_max_steps(static_cast<int>(std::ceil(stop_time/dt_run)) + 1);

            // Take as many steps as needed to reach stop_time.  The last step is
            // shortened if needed.
            integrator.evolve(phi, stop_time);

            MultiFab::Subtract(phi, exact, 0, 0, 1, 0);
            Real const err = phi.norm0();

            Real rate = 0;
            if (irun > 0) {
                rate = std::log2(err_prev/err);
                // Skip the check once round-off dominates the error.
                if (order > 0 && err > 1000*std::numeric_limits<Real>::epsilon()
                    && rate < order - 0.2_rt) {
                    status = 1;
                }
            }

            amrex::Print() << std::scientific << std::setprecision(3)
                           << "  " << dt_run << "  " << std::setw(6) << nsteps
                           << "  " << std::setw(9) << nrhs << "  " << err;
            if (irun > 0) {
                amrex::Print() << std::fixed << std::setprecision(2) << "  " << std::setw(5) << rate;
            }
            amrex::Print() << "\n";
            err_prev = err;
        }

        if (order > 0) {
            amrex::Print() << "\nExpected order of accuracy: " << order << "\n";
        }
        if (!(err_prev < 0.1_rt)) { status = 1; }  // also catches NaN
        if (status != 0) {
            amrex::Print() << "FAILED\n";
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
