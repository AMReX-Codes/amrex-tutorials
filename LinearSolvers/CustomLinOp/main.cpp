#include "MultiCompNodeLap.H"

#include <AMReX.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <type_traits>

using namespace amrex;

// Solve  sum_m C_nm lap(phi_m) = f_n  on the nodes of the unit square/cube
// with phi = 0 on the boundary, using the user-defined MultiCompNodeLap.
// Exact solution: phi_n = prod_d sin((n+1) pi x_d).

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        Vector<Real> coeff{2.0_rt, 1.0_rt,
                           0.5_rt, 1.0_rt};
        int verbose = 1;
        Real tol_rel = std::is_same_v<Real,float> ? Real(1.e-5) : Real(1.e-10);
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            // queryarr does not shrink a non-empty vector, so read into an
            // empty one.
            Vector<Real> coeff_in;
            if (pp.queryarr("coeff", coeff_in)) { coeff = coeff_in; }

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

        // The operator takes cell-centered grids; the unknowns are nodal.
        MultiCompNodeLap linop({geom}, {ba}, {dm}, coeff);
        linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet)},
                          {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet)});
        int const ncomp = linop.getNComp();

        BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
        MultiFab phi(nba, dm, ncomp, 1);
        MultiFab rhs(nba, dm, ncomp, 0);
        MultiFab exact(nba, dm, ncomp, 0);

        // Initial guess.  For a nodal solver, it also provides the values
        // on Dirichlet boundary nodes, here zero.
        phi.setVal(0);

        GpuArray<Real,MultiCompNodeLap::max_ncomp*MultiCompNodeLap::max_ncomp> c{};
        for (int i = 0; i < ncomp*ncomp; ++i) { c[i] = coeff[i]; }
        auto const problo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();
        auto const& rhs_a = rhs.arrays();
        auto const& exact_a = exact.arrays();
        ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            IntVect const iv(AMREX_D_DECL(i,j,k));
            GpuArray<Real,MultiCompNodeLap::max_ncomp> lap{};
            for (int m = 0; m < ncomp; ++m) {
                Real const km = Real(m+1)*Math::pi<Real>();
                Real phi_m = 1;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    phi_m *= std::sin(km*(problo[d] + Real(iv[d])*dx[d]));
                }
                exact_a[b](i,j,k,m) = phi_m;
                lap[m] = -Real(AMREX_SPACEDIM)*km*km*phi_m;
            }
            for (int n = 0; n < ncomp; ++n) {
                Real r = 0;
                for (int m = 0; m < ncomp; ++m) { r += c[n*ncomp+m]*lap[m]; }
                rhs_a[b](i,j,k,n) = r;
            }
        });

        amrex::Print() << "Solving for " << ncomp << " components with C =";
        for (auto x : coeff) { amrex::Print() << " " << x; }
        amrex::Print() << "\n";

        MLMG mlmg(linop);
        mlmg.setVerbose(verbose);
        mlmg.solve({&phi}, {&rhs}, tol_rel, Real(0));

        Real const h = geom.CellSize(0);
        // The error is O(dx^2) and grows with the wave numbers, i.e., ncomp.
        Real const tol = Real(5*ncomp*ncomp)*h*h;
        MultiFab::Subtract(exact, phi, 0, 0, ncomp, 0);
        for (int n = 0; n < ncomp; ++n) {
            Real const err = exact.norm0(n);
            amrex::Print() << "Component " << n << ": max error = " << err << "\n";
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
