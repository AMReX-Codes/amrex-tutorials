#include <AMReX.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_TimeIntegrator.H>

#include <cmath>
#include <iomanip>
#include <string>

using namespace amrex;

namespace {

struct Problem
{
    Real kappa = 0.1_rt;  // diffusion coefficient
    GpuArray<Real,AMREX_SPACEDIM> vel{AMREX_D_DECL(1.0_rt, 0.5_rt, 0.25_rt)};
};

// rhs = -vel.grad(phi) + kappa*lap(phi), or only one of the two terms, with
// second-order centered differences
void compute_rhs (MultiFab& rhs, MultiFab& phi, Geometry const& geom, Problem const& prob,
                  bool advection, bool diffusion)
{
    // The integrator passes its stage data here; fill its ghost cells first.
    phi.FillBoundary(geom.periodicity());

    auto const dxi = geom.InvCellSizeArray();
    auto const vel = prob.vel;
    Real const cadv = advection ? 1._rt : 0._rt;
    Real const kappa = diffusion ? prob.kappa : 0._rt;
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
        r[b](i,j,k) = -cadv*adv + kappa*lap;
    });
}

// Exact solution of the spatially discretized system for phi(t=0) = prod_d sin(2 pi x_d).
// Returns the amplitude of the solution.
Real exact_solution (MultiFab& phi, Geometry const& geom, Problem const& prob, Real time)
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
    return decay;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 16;
        int max_grid_size = 8;
        Real stop_time = 0.1_rt;
        Real dt = 0.02_rt;        // fixed step of the implicit, IMEX and MRI methods
        Real fast_dt = 0.002_rt;  // fast step of the MRI methods
        Real max_error = 1.e-3_rt;
        Vector<std::string> methods;
        Problem prob;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("stop_time", stop_time);
            pp.query("dt", dt);
            pp.query("fast_dt", fast_dt);
            pp.query("max_error", max_error);
            if (!pp.queryarr("methods", methods)) {
                methods = {"ERK", "DIRK", "IMEX-RK", "EX-MRI"};
            }
            pp.query("kappa", prob.kappa);
            Vector<Real> vel;
            if (pp.queryarr("velocity", vel)) {
                AMREX_ALWAYS_ASSERT(vel.size() >= AMREX_SPACEDIM);
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { prob.vel[d] = vel[d]; }
            }
        }

        // Integrator tolerances, unless given in the inputs file
        {
            ParmParse pp("integration");
            Real rel_tol = 1.e-6_rt;  // for adaptive steps and implicit solves
            Real abs_tol = 1.e-10_rt;
            pp.queryAdd("rel_tol", rel_tol);
            pp.queryAdd("abs_tol", abs_tol);
            // Stop at stop_time instead of stepping past it and interpolating.
            ParmParse("integration.sundials").add("stop_time", stop_time);
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
        Real const amplitude = exact_solution(exact, geom, prob, stop_time);

        struct Result { std::string method; int steps, fast_steps, nrhs; Real error; };
        Vector<Result> results;

        for (auto const& method : methods)
        {
            amrex::Print() << "\n=== integration.sundials.type = " << method << "\n";

            // Same as "integration.sundials.type = <method>" in the inputs file.
            // Other integration.sundials.* parameters (e.g., method names) are
            // read from the inputs as usual.
            ParmParse("integration.sundials").add("type", method);

            Real time = 0;
            exact_solution(phi, geom, prob, time);

            // Choose the SUNDIALS backend in the code instead of with integration.type.
            TimeIntegrator<MultiFab> integrator(IntegratorTypes::Sundials, phi, time);

            // Right-hand side functions that count their calls
            int nrhs = 0;
            auto make_rhs = [&] (bool advection, bool diffusion)
            {
                return [&, advection, diffusion] (MultiFab& rhs, MultiFab& state, Real /*t*/)
                {
                    ++nrhs;
                    compute_rhs(rhs, state, geom, prob, advection, diffusion);
                };
            };

            if (method == "ERK" || method == "DIRK") {
                // A single right-hand side, treated explicitly (ERK) or implicitly (DIRK)
                integrator.set_rhs(make_rhs(true, true));
            } else if (method == "IMEX-RK") {
                // Stiff diffusion implicit, advection explicit
                integrator.set_imex_rhs(make_rhs(false, true), make_rhs(true, false));
            } else if (method == "EX-MRI") {
                // Slow advection with step dt; fast diffusion with substeps of
                // fast_dt, explicit (integration.sundials.fast_type = ERK) or
                // implicit (DIRK).
                integrator.set_rhs(make_rhs(true, false));
                integrator.set_fast_rhs(make_rhs(false, true));
                integrator.set_fast_time_step(fast_dt);
            } else {
                amrex::Abort("This tutorial supports ERK, DIRK, IMEX-RK and EX-MRI");
            }

            if (method == "ERK") {
                // The stable explicit step is limited by diffusion; let SUNDIALS
                // choose the step from the error estimate.
                integrator.set_adaptive_step();
            } else {
                integrator.set_time_step(dt);
            }

            int steps = 0;
            int fast_steps = 0;
            integrator.set_post_step_action([&] (MultiFab& /*state*/, Real /*t*/) { ++steps; });
            integrator.set_post_fast_step_action([&] (MultiFab& /*state*/, Real /*t*/) { ++fast_steps; });

            integrator.evolve(phi, stop_time);

            MultiFab::Subtract(phi, exact, 0, 0, 1, 0);
            Real const error = phi.norm0() / amplitude;
            results.push_back({method, steps, fast_steps, nrhs, error});
            if (!(error < max_error)) { status = 1; }  // also catches NaN
        }

        amrex::Print() << "\nmethod      steps  fast steps  RHS evals  relative error\n";
        for (auto const& r : results) {
            amrex::Print() << std::left << std::setw(10) << r.method << std::right
                           << std::setw(7) << r.steps << std::setw(12) << r.fast_steps
                           << std::setw(11) << r.nrhs << std::scientific << std::setprecision(3)
                           << std::setw(16) << r.error << "\n";
        }
        if (status != 0) {
            amrex::Print() << "FAILED: error exceeds " << max_error << "\n";
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
