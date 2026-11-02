#include <AMReX.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <type_traits>

using namespace amrex;

// Solve the heat equation d(phi)/dt = kappa*lap(phi) on a periodic domain
// with backward Euler:
//
//     (1 - dt*kappa*lap) phi^{n+1} = phi^n,
//
// which is MLABecLaplacian's  a*alpha*phi - b*div(beta*grad(phi)) = rhs
// with a = 1, alpha = 1, b = dt and beta = kappa.

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int nsteps = 20;
        Real kappa = 1;
        Real dt_factor = 10; // dt in units of the explicit stability limit
        int verbose = 1;
        Real tol_rel = std::is_same_v<Real,float> ? Real(1.e-5) : Real(1.e-10);
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nsteps", nsteps);
            pp.query("kappa", kappa);
            pp.query("dt_factor", dt_factor);

            ParmParse ppm("mlmg");
            ppm.query("verbose", verbose);
            ppm.query("tol_rel", tol_rel);
        }

        Box domain(IntVect(0), IntVect(n_cell-1));
        RealBox real_box({AMREX_D_DECL(0._rt,0._rt,0._rt)}, {AMREX_D_DECL(1._rt,1._rt,1._rt)});
        Geometry geom(domain, real_box, CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        auto const problo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();

        Real const dx2 = dx[0]*dx[0];
        Real const dt = dt_factor * dx2 / (Real(2*AMREX_SPACEDIM)*kappa);

        // Initial condition: a single Fourier mode, phi = prod_d sin(2 pi x_d).
        MultiFab phi_new(ba, dm, 1, 1);
        MultiFab phi_old(ba, dm, 1, 0);
        MultiFab phi0(ba, dm, 1, 0);
        auto const& phi0_a = phi0.arrays();
        ParallelFor(phi0, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            IntVect const iv(AMREX_D_DECL(i,j,k));
            Real r = 1;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Real const x = problo[d] + (Real(iv[d])+Real(0.5))*dx[d];
                r *= std::sin(Real(2)*Math::pi<Real>()*x);
            }
            phi0_a[b](i,j,k) = r;
        });
        MultiFab::Copy(phi_new, phi0, 0, 0, 1, 0);

        // Build the operator and the solver once and reuse them every step.
        MLABecLaplacian mlabec({geom}, {ba}, {dm});
        mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                         LinOpBCType::Periodic,
                                         LinOpBCType::Periodic)},
                           {AMREX_D_DECL(LinOpBCType::Periodic,
                                         LinOpBCType::Periodic,
                                         LinOpBCType::Periodic)});
        mlabec.setLevelBC(0, nullptr); // no boundary data needed for periodic
        mlabec.setScalars(Real(1), dt);
        mlabec.setACoeffs(0, Real(1));
        mlabec.setBCoeffs(0, kappa);

        MLMG mlmg(mlabec);
        mlmg.setVerbose(verbose);

        Real time = 0;
        for (int step = 1; step <= nsteps; ++step) {
            MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);
            // phi_new holds the previous solution, a good initial guess.
            mlmg.solve({&phi_new}, {&phi_old}, tol_rel, Real(0));
            time += dt;
            amrex::Print() << "Step " << step << ": time = " << time
                           << ", MLMG iterations = " << mlmg.getNumIters() << "\n";
        }

        // The mode is an eigenfunction of the discrete Laplacian with
        // eigenvalue -lam, so each step multiplies it by 1/(1+dt*kappa*lam).
        // The PDE solution decays as exp(-kappa*D*(2 pi)^2*t) instead.
        Real lam = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real const s = std::sin(Math::pi<Real>()*dx[d]);
            lam += Real(4)*s*s/(dx[d]*dx[d]);
        }
        Real const g_discrete = std::pow(Real(1)/(Real(1)+dt*kappa*lam), nsteps);
        Real const g_pde = std::exp(-kappa*Real(AMREX_SPACEDIM*4)*Math::pi<Real>()*Math::pi<Real>()*time);

        MultiFab diff(ba, dm, 1, 0);
        MultiFab::LinComb(diff, Real(1), phi_new, 0, -g_discrete, phi0, 0, 0, 1, 0);
        Real const err_discrete = diff.norm0() / g_discrete;
        MultiFab::LinComb(diff, Real(1), phi_new, 0, -g_pde, phi0, 0, 0, 1, 0);
        Real const err_pde = diff.norm0() / g_pde;

        amrex::Print() << "Decay factor: numerical " << phi_new.norm0()/phi0.norm0()
                       << ", discrete exact " << g_discrete
                       << ", PDE exact " << g_pde << "\n"
                       << "Relative error vs. discrete exact solution: " << err_discrete << "\n"
                       << "Relative error vs. PDE exact solution:      " << err_pde << "\n";

        // The first error only measures how well the linear systems are
        // solved.  The second also includes the O(dt) time error.
        Real const tol = std::is_same_v<Real,float> ? Real(1.e-3) : Real(1.e-7);
        if (err_discrete > tol) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
