#include <AMReX.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <type_traits>

using namespace amrex;

AMREX_ENUM(BCKind, dirichlet, neumann, periodic);

namespace {

// Exact solution phi = prod_d cos(2 pi x_d) on the unit square/cube.  Its
// normal derivative vanishes on all faces and it is periodic, so it works
// with all three boundary types.  lap(phi) = -D (2 pi)^2 phi.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_phi (GpuArray<Real,AMREX_SPACEDIM> const& x)
{
    Real r = 1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        r *= std::cos(Real(2)*Math::pi<Real>()*x[d]);
    }
    return r;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        Array<BCKind,AMREX_SPACEDIM> bc{AMREX_D_DECL(BCKind::dirichlet,
                                                     BCKind::dirichlet,
                                                     BCKind::dirichlet)};
        int verbose = 1;
        int max_iter = 100;
        Real tol_rel = std::is_same_v<Real,float> ? Real(1.e-5) : Real(1.e-10);
        Real tol_abs = 0;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            // e.g., bc = dirichlet neumann periodic (one per direction)
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                pp.query_enum_case_insensitive("bc", bc[d], d);
            }

            ParmParse ppm("mlmg");
            ppm.query("verbose", verbose);
            ppm.query("max_iter", max_iter);
            ppm.query("tol_rel", tol_rel);
            ppm.query("tol_abs", tol_abs);
        }

        amrex::Print() << "Boundary conditions:";
        for (auto b : bc) { amrex::Print() << " " << getEnumNameString(b); }
        amrex::Print() << "\n";

        // Domain boundary types for the operator
        Array<LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
        Array<int,AMREX_SPACEDIM> is_periodic{};
        bool has_dirichlet = false;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (bc[d] == BCKind::dirichlet) {
                lobc[d] = hibc[d] = LinOpBCType::Dirichlet;
                has_dirichlet = true;
            } else if (bc[d] == BCKind::neumann) {
                lobc[d] = hibc[d] = LinOpBCType::Neumann;
            } else {
                lobc[d] = hibc[d] = LinOpBCType::Periodic;
                is_periodic[d] = 1;
            }
        }

        // Solve on two grids.  For a second-order method, the error drops
        // by a factor of four when the grid spacing is halved.
        Real err[2];
        for (int igrid = 0; igrid < 2; ++igrid)
        {
            int const n = (igrid == 0) ? n_cell : 2*n_cell;

            Box domain(IntVect(0), IntVect(n-1));
            RealBox real_box({AMREX_D_DECL(0._rt,0._rt,0._rt)}, {AMREX_D_DECL(1._rt,1._rt,1._rt)});
            Geometry geom(domain, real_box, CoordSys::cartesian, is_periodic);

            BoxArray ba(domain);
            ba.maxSize(max_grid_size);
            DistributionMapping dm(ba);

            auto const problo = geom.ProbLoArray();
            auto const probhi = geom.ProbHiArray();
            auto const dx = geom.CellSizeArray();

            MultiFab rhs(ba, dm, 1, 0);
            MultiFab exact(ba, dm, 1, 0);
            auto const& rhs_a = rhs.arrays();
            auto const& exact_a = exact.arrays();
            ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                IntVect const iv(AMREX_D_DECL(i,j,k));
                GpuArray<Real,AMREX_SPACEDIM> x;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    x[d] = problo[d] + (Real(iv[d]) + Real(0.5)) * dx[d];
                }
                exact_a[b](i,j,k) = exact_phi(x);
                rhs_a[b](i,j,k) = -Real(AMREX_SPACEDIM*4) * Math::pi<Real>() * Math::pi<Real>()
                    * exact_phi(x);
            });

            // phi needs one ghost cell.  Outside the domain, it holds the
            // Dirichlet boundary values, which live on the domain faces.  Its
            // valid cells hold the initial guess.
            MultiFab phi(ba, dm, 1, 1);
            auto const& phi_a = phi.arrays();
            ParallelFor(phi, IntVect(1), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                IntVect const iv(AMREX_D_DECL(i,j,k));
                GpuArray<Real,AMREX_SPACEDIM> x;
                bool inside = true;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    if (iv[d] < domain.smallEnd(d)) {
                        x[d] = problo[d];
                        inside = false;
                    } else if (iv[d] > domain.bigEnd(d)) {
                        x[d] = probhi[d];
                        inside = false;
                    } else {
                        x[d] = problo[d] + (Real(iv[d]) + Real(0.5)) * dx[d];
                    }
                }
                phi_a[b](i,j,k) = inside ? Real(0) : exact_phi(x);
            });

            MLPoisson mlpoisson({geom}, {ba}, {dm});

            // Domain boundary types, then the boundary values (in phi's ghost cells).
            mlpoisson.setDomainBC(lobc, hibc);
            mlpoisson.setLevelBC(0, &phi);

            MLMG mlmg(mlpoisson);
            mlmg.setVerbose(verbose);
            mlmg.setMaxIter(max_iter);
            // Stop when max|rhs - L(phi)| <= max(tol_abs, tol_rel*R), with R the
            // larger of max|rhs| and the initial residual.
            mlmg.solve({&phi}, {&rhs}, tol_rel, tol_abs);

            if (!has_dirichlet) {
                // The solution is unique only up to a constant.  The exact solution
                // has zero mean, so remove the mean of the numerical one.
                Real const mean = phi.sum(0) / Real(domain.d_numPts());
                phi.plus(-mean, 0, 1, 0);
            }

            MultiFab::Subtract(exact, phi, 0, 0, 1, 0);
            err[igrid] = exact.norm0();
            amrex::Print() << "n_cell = " << n << ": max error = " << err[igrid] << "\n";
        }

        Real const ratio = err[0] / err[1];
        amrex::Print() << "Convergence rate = " << std::log2(ratio) << "\n";

        if (ratio < Real(3.5)) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
