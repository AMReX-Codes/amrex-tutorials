#include <AMReX.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <type_traits>

using namespace amrex;

// Project a cell-centered velocity u onto a divergence-free field:
//
//     div(sigma grad(phi)) = div(u),   u_new = u - sigma*grad(phi),
//
// with phi on the nodes.  The domain is periodic in x and has walls
// (u.n = 0, i.e., Neumann for phi) in the other directions.
//
// This is an approximate projection: the nodal Laplacian is not exactly
// the product of the discrete divergence and gradient, so div(u_new) is
// O(dx^2) rather than zero.

namespace {

using RealVect_t = GpuArray<Real,AMREX_SPACEDIM>;

// Divergence-free part, from the stream function sin(2 pi x) sin^2(pi y).
// It is periodic in x and has zero normal component at the walls.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real u_div_free (RealVect_t const& x, int dir)
{
    constexpr Real pi = Math::pi<Real>();
    if (dir == 0) {
        return pi*std::sin(2*pi*x[0])*std::sin(2*pi*x[1]);
    } else if (dir == 1) {
        Real const s = std::sin(pi*x[1]);
        return -2*pi*std::cos(2*pi*x[0])*s*s;
    } else {
        return 0;
    }
}

// Gradient part, grad(q) with q = cos(2 pi x) prod_{d>0} cos(pi x_d).
// Its normal component also vanishes at the walls.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real u_gradient (RealVect_t const& x, int dir)
{
    constexpr Real pi = Math::pi<Real>();
    Real r = 1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real const kd = (d == 0) ? 2*pi : pi;
        r *= (d == dir) ? -kd*std::sin(kd*x[d]) : std::cos(kd*x[d]);
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
        int verbose = 1;
        Real tol_rel = std::is_same_v<Real,float> ? Real(1.e-4) : Real(1.e-10);
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);

            ParmParse ppm("mlmg");
            ppm.query("verbose", verbose);
            ppm.query("tol_rel", tol_rel);
        }

        Box domain(IntVect(0), IntVect(n_cell-1));
        RealBox real_box({AMREX_D_DECL(0._rt,0._rt,0._rt)}, {AMREX_D_DECL(1._rt,1._rt,1._rt)});
        Geometry geom(domain, real_box, CoordSys::cartesian, {AMREX_D_DECL(1,0,0)});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        // Cell-centered velocity with one ghost cell, needed by the nodal
        // divergence.  Ghost cells outside the walls are not used, but must
        // not hold NaNs.
        MultiFab vel(ba, dm, AMREX_SPACEDIM, 1);
        MultiFab vel_exact(ba, dm, AMREX_SPACEDIM, 0);
        vel.setVal(0);
        {
            auto const problo = geom.ProbLoArray();
            auto const dx = geom.CellSizeArray();
            auto const& vel_a = vel.arrays();
            auto const& vex_a = vel_exact.arrays();
            ParallelFor(vel_exact, IntVect(0), AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
            {
                IntVect const iv(AMREX_D_DECL(i,j,k));
                RealVect_t x;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    x[d] = problo[d] + (Real(iv[d])+Real(0.5))*dx[d];
                }
                vex_a[b](i,j,k,n) = u_div_free(x,n);
                vel_a[b](i,j,k,n) = u_div_free(x,n) + u_gradient(x,n);
            });
        }

        // Nodal Poisson operator with constant sigma = 1.  For variable
        // density, use sigma = 1/rho with setSigma instead.  The empty list
        // of factories means no embedded boundaries.
        Vector<FabFactory<FArrayBox> const*> const no_eb{};
        MLNodeLaplacian linop({geom}, {ba}, {dm}, LPInfo(), no_eb, Real(1));
        linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                        LinOpBCType::Neumann,
                                        LinOpBCType::Neumann)},
                          {AMREX_D_DECL(LinOpBCType::Periodic,
                                        LinOpBCType::Neumann,
                                        LinOpBCType::Neumann)});

        BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
        MultiFab rhs(nba, dm, 1, 0);
        MultiFab phi(nba, dm, 1, 1);
        phi.setVal(0);

        // rhs = div(u) on the nodes
        linop.compDivergence({&rhs}, {&vel});
        Real const div0 = rhs.norm0();

        // With walls and periodic boundaries only, phi is determined up to a
        // constant.  MLMG handles this singular system.
        MLMG mlmg(linop);
        mlmg.setVerbose(verbose);
        mlmg.solve({&phi}, {&rhs}, tol_rel, Real(0));

        // u -= sigma*grad(phi)
        linop.updateVelocity({&vel}, {&phi});

        linop.compDivergence({&rhs}, {&vel});
        Real const div1 = rhs.norm0();

        MultiFab::Subtract(vel_exact, vel, 0, 0, AMREX_SPACEDIM, 0);
        Real err = 0;
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            err = std::max(err, vel_exact.norm0(n));
        }

        amrex::Print() << "max |div(u)| before projection: " << div0 << "\n"
                       << "max |div(u)| after projection:  " << div1 << "\n"
                       << "max error in the divergence-free part: " << err << "\n";

        // Both the remaining divergence and the error are O(dx^2).
        Real const h = geom.CellSize(0);
        if (div1 > Real(5)*h*h*div0 || err > Real(50)*h*h) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
