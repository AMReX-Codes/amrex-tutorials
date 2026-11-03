#include <AMReX.H>
#include <AMReX_MLCurlCurl.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <type_traits>

using namespace amrex;

// Solve  curl(alpha curl E) + beta E = f  for the edge-centered field E on
// the unit square/cube with perfect-conductor walls (tangential E = 0).
//
// E always has three components.  E_d lives on the edges along direction d:
// it is cell-centered in direction d and nodal in the others.  In 2D, E_z
// is nodal in x and y.
//
// Exact solution: E_d = prod_{e != d} sin(pi x_e) (product over the AMREX_SPACEDIM
// coordinates).  It vanishes on the walls tangential to it and div(E) = 0,
// so curl(curl E) = -lap(E) = n_d pi^2 E_d with n_d factors in the product.

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_E (int dir, GpuArray<Real,AMREX_SPACEDIM> const& x)
{
    Real r = 1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        if (d != dir) { r *= std::sin(Math::pi<Real>()*x[d]); }
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
        Real alpha = 1;
        Real beta = 1;
        int verbose = 1;
        Real tol_rel = std::is_same_v<Real,float> ? Real(1.e-4) : Real(1.e-10);
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("alpha", alpha);
            pp.query("beta", beta);

            ParmParse ppm("mlmg");
            ppm.query("verbose", verbose);
            ppm.query("tol_rel", tol_rel);
        }

        Box domain(IntVect(0), IntVect(n_cell-1));
        RealBox real_box({AMREX_D_DECL(0._rt,0._rt,0._rt)}, {AMREX_D_DECL(1._rt,1._rt,1._rt)});
        Geometry geom(domain, real_box, CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        // Edge-centered fields: an Array of three MultiFabs.
        Array<MultiFab,3> E, rhs, E_exact;
        for (int dir = 0; dir < 3; ++dir) {
            IntVect itype(1);
            if (dir < AMREX_SPACEDIM) { itype[dir] = 0; }
            BoxArray const eba = amrex::convert(ba, itype);
            E      [dir].define(eba, dm, 1, 1);
            rhs    [dir].define(eba, dm, 1, 0);
            E_exact[dir].define(eba, dm, 1, 0);

            // Initial guess.  It also supplies the values of the tangential
            // components on Dirichlet boundaries, here zero.
            E[dir].setVal(0);

            int nfactors = 0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { if (d != dir) { ++nfactors; } }
            Real const fac = alpha*Real(nfactors)*Math::pi<Real>()*Math::pi<Real>() + beta;

            auto const problo = geom.ProbLoArray();
            auto const dx = geom.CellSizeArray();
            auto const ixt = rhs[dir].ixType();
            auto const& rhs_a = rhs[dir].arrays();
            auto const& ex_a = E_exact[dir].arrays();
            ParallelFor(rhs[dir], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                IntVect const iv(AMREX_D_DECL(i,j,k));
                GpuArray<Real,AMREX_SPACEDIM> x;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    Real const offset = ixt.cellCentered(d) ? Real(0.5) : Real(0);
                    x[d] = problo[d] + (Real(iv[d]) + offset) * dx[d];
                }
                Real const e = exact_E(dir, x);
                ex_a[b](i,j,k) = e;
                rhs_a[b](i,j,k) = fac * e;
            });
        }

        MLCurlCurl mlcc({geom}, {ba}, {dm});
        mlcc.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                       LinOpBCType::Dirichlet,
                                       LinOpBCType::Dirichlet)},
                         {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                       LinOpBCType::Dirichlet,
                                       LinOpBCType::Dirichlet)});
        mlcc.setScalars(alpha, beta);

        // Nodal points shared by neighboring boxes must hold the same rhs.
        // Ours is computed from a formula, so this is only for illustration.
        mlcc.prepareRHS({&rhs});

        // MLMG templated on the field type Array<MultiFab,3>
        MLMGT<Array<MultiFab,3>> mlmg(mlcc);
        mlmg.setVerbose(verbose);
        mlmg.solve({&E}, {&rhs}, tol_rel, Real(0));

        Real const h = geom.CellSize(0);
        Real const tol = Real(2)*h*h;
        Array<char const*,3> const names{"Ex", "Ey", "Ez"};
        for (int dir = 0; dir < 3; ++dir) {
            MultiFab::Subtract(E_exact[dir], E[dir], 0, 0, 1, 0);
            Real const err = E_exact[dir].norm0();
            amrex::Print() << names[dir] << ": max error = " << err << "\n";
            if (err > tol) { status = 1; }
        }

        if (status) {
            amrex::Print() << "FAILED: error larger than " << tol << "\n";
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
