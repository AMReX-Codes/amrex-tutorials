#include "ConvDiffOp.H"

#include <AMReX.H>
#include <AMReX_GMRES.H>
#include <AMReX_GMRES_MLMG.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <type_traits>

using namespace amrex;

namespace {

using RealVect_t = GpuArray<Real,AMREX_SPACEDIM>;

// Exact solution u = prod_d sin(pi x_d), zero on the boundary of the unit
// square/cube.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_u (RealVect_t const& x)
{
    Real r = 1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { r *= std::sin(Math::pi<Real>()*x[d]); }
    return r;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_du (RealVect_t const& x, int dir)
{
    constexpr Real pi = Math::pi<Real>();
    Real r = 1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        r *= (d == dir) ? pi*std::cos(pi*x[d]) : std::sin(pi*x[d]);
    }
    return r;
}

// Fill exact with u and rhs with lap(u) - c.grad(u).
void init_problem (MultiFab& exact, MultiFab& rhs, Geometry const& geom,
                   RealVect_t const& vel)
{
    auto const problo = geom.ProbLoArray();
    auto const dx = geom.CellSizeArray();
    auto const& exact_a = exact.arrays();
    auto const& rhs_a = rhs.arrays();
    ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        IntVect const iv(AMREX_D_DECL(i,j,k));
        RealVect_t x;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            x[d] = problo[d] + (Real(iv[d])+Real(0.5))*dx[d];
        }
        Real const u = exact_u(x);
        Real r = -Real(AMREX_SPACEDIM)*Math::pi<Real>()*Math::pi<Real>()*u;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { r -= vel[d]*exact_du(x,d); }
        exact_a[b](i,j,k) = u;
        rhs_a[b](i,j,k) = r;
    });
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        Vector<Real> velocity{AMREX_D_DECL(Real(20),Real(10),Real(5))};
        bool use_precond = true;
        int verbose = 1;
        Real tol_rel = std::is_same_v<Real,float> ? Real(1.e-5) : Real(1.e-10);
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
            pp.query("use_precond", use_precond);
            pp.query("verbose", verbose);
            pp.query("tol_rel", tol_rel);
        }

        Box domain(IntVect(0), IntVect(n_cell-1));
        RealBox real_box({AMREX_D_DECL(0._rt,0._rt,0._rt)}, {AMREX_D_DECL(1._rt,1._rt,1._rt)});
        Geometry geom(domain, real_box, CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        RealVect_t vel;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { vel[d] = velocity[d]; }

        MultiFab exact(ba, dm, 1, 0);
        MultiFab rhs(ba, dm, 1, 0);
        init_problem(exact, rhs, geom, vel);

        Real const h = geom.CellSize(0);
        Real const err_tol = Real(10)*h*h;

        // Laplacian with homogeneous Dirichlet boundaries, used by MLMG as
        // the preconditioner.  maxorder 2 matches the stencil in ConvDiffOp.
        MLPoisson mlpoisson({geom}, {ba}, {dm});
        mlpoisson.setMaxOrder(2);
        mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)},
                              {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)});
        mlpoisson.setLevelBC(0, nullptr); // homogeneous

        MLMG mlmg(mlpoisson);
        mlmg.setVerbose(0);
        // A bottom solver with a convergence test (the default BiCGStab)
        // would make the preconditioner nonlinear.
        mlmg.setBottomSolver(MLMG::BottomSolver::smoother);

        // Part 1: GMRES for the user-defined operator.
        {
            ConvDiffOp op(geom, ba, dm, vel, mlmg);
            op.usePrecond(use_precond);

            GMRES<MultiFab,ConvDiffOp> gmres;
            gmres.define(op);
            gmres.setVerbose(verbose);
            gmres.setRestartLength(30);
            gmres.setMaxIters(1000);

            MultiFab u = op.makeVecLHS();
            u.setVal(0);
            gmres.solve(u, rhs, tol_rel, Real(0));

            MultiFab::Subtract(u, exact, 0, 0, 1, 0);
            Real const err = u.norm0();
            amrex::Print() << "Convection-diffusion with GMRES: "
                           << gmres.getNumIters() << " iterations, status "
                           << gmres.getStatus() << ", max error = " << err << "\n";
            if (gmres.getStatus() != 0 || err > err_tol) { status = 1; }
        }

        // Part 2: GMRESMLMG solves an MLMG problem, here lap(u) = f, with
        // GMRES preconditioned by MLMG.  This can help when MLMG alone
        // converges slowly.
        {
            MultiFab lap_rhs(ba, dm, 1, 0);
            init_problem(exact, lap_rhs, geom, RealVect_t{});

            MultiFab u(ba, dm, 1, 1);
            u.setVal(0); // GMRESMLMG uses u as the initial guess

            GMRESMLMG gmsolver(mlmg);
            gmsolver.getGMRES().setVerbose(verbose);
            gmsolver.solve(u, lap_rhs, tol_rel, Real(0));

            MultiFab::Subtract(u, exact, 0, 0, 1, 0);
            Real const err = u.norm0();
            amrex::Print() << "Poisson with GMRESMLMG: "
                           << gmsolver.getNumIters() << " iterations, max error = "
                           << err << "\n";
            if (err > err_tol) { status = 1; }
        }

        if (status) {
            amrex::Print() << "FAILED: error larger than " << err_tol << "\n";
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
