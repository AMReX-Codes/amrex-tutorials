#include <AMReX.H>
#include <AMReX_FFT_Poisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>
#include <type_traits>
#include <utility>

using namespace amrex;

namespace {

using BCPair = std::pair<FFT::Boundary,FFT::Boundary>;

// Manufactured solution phi = prod_d f_d(pi*c_d*x_d), with f_d = sin or cos.
// Then lap(phi) = -pi^2 * sum_d c_d^2 * phi.
struct Exact
{
    GpuArray<Real,AMREX_SPACEDIM> c;
    GpuArray<int,AMREX_SPACEDIM> use_cos;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Real operator() (GpuArray<Real,AMREX_SPACEDIM> const& x) const
    {
        Real r = 1;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real a = Math::pi<Real>() * c[d] * x[d];
            r *= use_cos[d] ? std::cos(a) : std::sin(a);
        }
        return r;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Real laplacian (GpuArray<Real,AMREX_SPACEDIM> const& x) const
    {
        Real c2 = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { c2 += c[d]*c[d]; }
        return -Math::pi<Real>()*Math::pi<Real>() * c2 * (*this)(x);
    }
};

// Solve lap(phi) = rhs with the given solver and check the result.
template <typename Solver>
bool run_case (std::string const& name, Solver& solver, Geometry const& geom,
               BoxArray const& ba, DistributionMapping const& dm, Exact const& exact)
{
    auto const problo = geom.ProbLoArray();
    auto const dx = geom.CellSizeArray();

    MultiFab rhs(ba, dm, 1, 0);
    MultiFab phi_exact(ba, dm, 1, 0);
    {
        auto const& rhs_a = rhs.arrays();
        auto const& exact_a = phi_exact.arrays();
        ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            IntVect iv(AMREX_D_DECL(i,j,k));
            GpuArray<Real,AMREX_SPACEDIM> x;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                x[d] = problo[d] + (Real(iv[d])+Real(0.5))*dx[d];
            }
            rhs_a[b](i,j,k) = exact.laplacian(x);
            exact_a[b](i,j,k) = exact(x);
        });
    }

    // One ghost cell: the solver fills it from the boundary conditions,
    // so that we can apply the discrete Laplacian below.
    MultiFab phi(ba, dm, 1, 1);
    solver.solve(phi, rhs);

    // The FFT solvers invert the standard second-order (2*D+1)-point
    // Laplacian exactly, so the discrete residual is at round-off level.
    MultiFab resid(ba, dm, 1, 0);
    {
        auto const& resid_a = resid.arrays();
        auto const& phi_a = phi.const_arrays();
        auto const& rhs_a = rhs.const_arrays();
        ParallelFor(resid, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            auto const& p = phi_a[b];
            Real lap = AMREX_D_TERM(
                (p(i-1,j,k) - Real(2)*p(i,j,k) + p(i+1,j,k)) / (dx[0]*dx[0]),
              + (p(i,j-1,k) - Real(2)*p(i,j,k) + p(i,j+1,k)) / (dx[1]*dx[1]),
              + (p(i,j,k-1) - Real(2)*p(i,j,k) + p(i,j,k+1)) / (dx[2]*dx[2]));
            resid_a[b](i,j,k) = lap - rhs_a[b](i,j,k);
        });
    }
    Real res = resid.norm0() / rhs.norm0();

    // Compared with the continuous solution, the error is O(dx^2).
    MultiFab::Subtract(phi_exact, phi, 0, 0, 1, 0);
    Real err = phi_exact.norm0();

    Real h = geom.CellSize(0);
    for (int d = 1; d < AMREX_SPACEDIM; ++d) { h = std::max(h, Real(geom.CellSize(d))); }
    Real const res_tol = std::is_same_v<Real,float> ? Real(1.e-3) : Real(1.e-10);
    Real const err_tol = Real(10) * h * h;

    amrex::Print() << name << "\n"
                   << "    relative residual " << res << "\n"
                   << "    max error         " << err << "\n";
    return res <= res_tol && err <= err_tol;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
        }

        // Unit domain [0,1]^D
        Box domain(IntVect(0), IntVect(n_cell-1));
        RealBox rb({AMREX_D_DECL(Real(0),Real(0),Real(0))},
                   {AMREX_D_DECL(Real(1),Real(1),Real(1))});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        bool pass = true;

        // 1. Periodic in all directions. The solution is unique up to a constant,
        //    and the right-hand side must have zero mean. The solver sets the
        //    mean of phi to that of the rhs (zero here).
        {
            Geometry geom(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
            FFT::Poisson<MultiFab> solver(geom);
            Exact exact{{AMREX_D_DECL(Real(2),Real(2),Real(2))},
                        {AMREX_D_DECL(0,0,0)}};
            pass = run_case("Periodic (FFT::Poisson)", solver, geom, ba, dm, exact) && pass;
        }

        // 2. Homogeneous Dirichlet (odd) and Neumann (even) boundaries, using
        //    sine and cosine transforms: x Dirichlet on both sides, y Neumann on
        //    both sides, z Dirichlet at zlo and Neumann at zhi.
        {
            Array<BCPair,AMREX_SPACEDIM> bc{AMREX_D_DECL(
                std::make_pair(FFT::Boundary::odd,  FFT::Boundary::odd),
                std::make_pair(FFT::Boundary::even, FFT::Boundary::even),
                std::make_pair(FFT::Boundary::odd,  FFT::Boundary::even))};
            Geometry geom(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
            FFT::Poisson<MultiFab> solver(geom, bc);
            Exact exact{{AMREX_D_DECL(Real(1),Real(1),Real(0.5))},
                        {AMREX_D_DECL(0,1,0)}};
            pass = run_case("Dirichlet/Neumann (FFT::Poisson)", solver, geom, ba, dm, exact)
                && pass;
        }

#if (AMREX_SPACEDIM == 3)
        // 3. FFT in x and y (periodic here), tridiagonal solve in z (Neumann at
        //    zlo, Dirichlet at zhi). PoissonHybrid also accepts a non-uniform
        //    dz: solver.solve(phi, rhs, dz).
        {
            Array<BCPair,AMREX_SPACEDIM> bc{
                BCPair{FFT::Boundary::periodic, FFT::Boundary::periodic},
                BCPair{FFT::Boundary::periodic, FFT::Boundary::periodic},
                BCPair{FFT::Boundary::even,     FFT::Boundary::odd}};
            Geometry geom(domain, rb, CoordSys::cartesian, {1,1,0});
            FFT::PoissonHybrid<MultiFab> solver(geom, bc);
            Exact exact{{Real(2),Real(2),Real(0.5)}, {0,0,1}};
            pass = run_case("Periodic x,y + Neumann/Dirichlet z (FFT::PoissonHybrid)",
                            solver, geom, ba, dm, exact) && pass;
        }
#endif

        if (pass) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
