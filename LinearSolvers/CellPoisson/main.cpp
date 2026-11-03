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

struct Params
{
    int max_grid_size = 16;
    Array<BCKind,AMREX_SPACEDIM> bc{AMREX_D_DECL(BCKind::dirichlet,
                                                 BCKind::dirichlet,
                                                 BCKind::dirichlet)};
    int verbose = 1;
    int max_iter = 100;
    Real tol_rel = std::is_same_v<Real,float> ? Real(1.e-5) : Real(1.e-10);
    Real tol_abs = 0;
};

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

// Solve lap(phi) = rhs on an n_cell^D grid and return the max-norm error.
Real solve_poisson (int n_cell, Params const& p)
{
    Box domain(IntVect(0), IntVect(n_cell-1));
    RealBox real_box({AMREX_D_DECL(0._rt,0._rt,0._rt)}, {AMREX_D_DECL(1._rt,1._rt,1._rt)});
    Array<int,AMREX_SPACEDIM> is_periodic{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        is_periodic[d] = (p.bc[d] == BCKind::periodic);
    }
    Geometry geom(domain, real_box, CoordSys::cartesian, is_periodic);

    BoxArray ba(domain);
    ba.maxSize(p.max_grid_size);
    DistributionMapping dm(ba);

    MultiFab rhs(ba, dm, 1, 0);
    MultiFab exact(ba, dm, 1, 0);
    // phi needs one ghost cell.  Outside the domain, it holds the Dirichlet
    // boundary values.  Its valid cells hold the initial guess.
    MultiFab phi(ba, dm, 1, 1);

    auto const problo = geom.ProbLoArray();
    auto const probhi = geom.ProbHiArray();
    auto const dx = geom.CellSizeArray();
    auto const dlo = lbound(domain);
    auto const dhi = ubound(domain);

    auto const& rhs_a = rhs.arrays();
    auto const& exact_a = exact.arrays();
    auto const& phi_a = phi.arrays();
    ParallelFor(phi, IntVect(1), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        IntVect const iv(AMREX_D_DECL(i,j,k));
        IntVect const lo(AMREX_D_DECL(dlo.x,dlo.y,dlo.z));
        IntVect const hi(AMREX_D_DECL(dhi.x,dhi.y,dhi.z));
        GpuArray<Real,AMREX_SPACEDIM> x;
        bool inside = true;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (iv[d] < lo[d]) {
                x[d] = problo[d]; // boundary value lives on the domain face
                inside = false;
            } else if (iv[d] > hi[d]) {
                x[d] = probhi[d];
                inside = false;
            } else {
                x[d] = problo[d] + (Real(iv[d]) + Real(0.5)) * dx[d];
            }
        }
        Real const phi_x = exact_phi(x);
        phi_a[b](i,j,k) = inside ? Real(0) : phi_x;
        if (rhs_a[b].contains(i,j,k)) {
            exact_a[b](i,j,k) = phi_x;
            rhs_a[b](i,j,k) = -Real(AMREX_SPACEDIM*4) * Math::pi<Real>() * Math::pi<Real>() * phi_x;
        }
    });

    MLPoisson mlpoisson({geom}, {ba}, {dm});

    // Domain boundary types, then the boundary values (in phi's ghost cells).
    Array<LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
    bool has_dirichlet = false;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (p.bc[d] == BCKind::dirichlet) {
            lobc[d] = hibc[d] = LinOpBCType::Dirichlet;
            has_dirichlet = true;
        } else if (p.bc[d] == BCKind::neumann) {
            lobc[d] = hibc[d] = LinOpBCType::Neumann;
        } else {
            lobc[d] = hibc[d] = LinOpBCType::Periodic;
        }
    }
    mlpoisson.setDomainBC(lobc, hibc);
    mlpoisson.setLevelBC(0, &phi);

    MLMG mlmg(mlpoisson);
    mlmg.setVerbose(p.verbose);
    mlmg.setMaxIter(p.max_iter);
    // Stop when max|rhs - L(phi)| <= max(tol_abs, tol_rel*R), with R the
    // larger of max|rhs| and the initial residual.
    mlmg.solve({&phi}, {&rhs}, p.tol_rel, p.tol_abs);

    if (!has_dirichlet) {
        // The solution is unique only up to a constant.  The exact solution
        // has zero mean, so remove the mean of the numerical one.
        Real const mean = phi.sum(0) / Real(domain.d_numPts());
        phi.plus(-mean, 0, 1, 0);
    }

    MultiFab::Subtract(exact, phi, 0, 0, 1, 0);
    return exact.norm0();
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        Params p;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", p.max_grid_size);
            // e.g., bc = dirichlet neumann periodic (one per direction)
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                pp.query_enum_case_insensitive("bc", p.bc[d], d);
            }

            ParmParse ppm("mlmg");
            ppm.query("verbose", p.verbose);
            ppm.query("max_iter", p.max_iter);
            ppm.query("tol_rel", p.tol_rel);
            ppm.query("tol_abs", p.tol_abs);
        }

        amrex::Print() << "Boundary conditions:";
        for (auto bc : p.bc) { amrex::Print() << " " << getEnumNameString(bc); }
        amrex::Print() << "\n";

        // Solve on two grids.  For a second-order method, the error drops
        // by a factor of four when the grid spacing is halved.
        Real const err_c = solve_poisson(n_cell, p);
        Real const err_f = solve_poisson(2*n_cell, p);
        Real const ratio = err_c / err_f;

        amrex::Print() << "n_cell = " << n_cell << ": max error = " << err_c << "\n"
                       << "n_cell = " << 2*n_cell << ": max error = " << err_f << "\n"
                       << "Convergence rate = " << std::log2(ratio) << "\n";

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
