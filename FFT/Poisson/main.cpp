#include <AMReX.H>
#include <AMReX_FFT_Poisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <type_traits>
#include <utility>

using namespace amrex;

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

        Real const pi = Math::pi<Real>();

        // The FFT solvers invert the standard second-order (2*D+1)-point
        // Laplacian exactly, so the discrete residual is at round-off level.
        // Compared with the continuous solution, the error is O(dx^2).
        Real const res_tol = std::is_same_v<Real,float> ? Real(1.e-3) : Real(1.e-10);

        // 1. Periodic in all directions. The solution is unique up to a constant,
        //    and the right-hand side must have zero mean. The solver sets the
        //    mean of phi to that of the rhs (zero here).
        {
            Geometry geom(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
            FFT::Poisson<MultiFab> solver(geom);
            auto const problo = geom.ProbLoArray();
            auto const dx = geom.CellSizeArray();

            // phi = prod_d sin(2 pi x_d), lap(phi) = -D (2 pi)^2 phi
            MultiFab rhs(ba, dm, 1, 0);
            MultiFab phi_exact(ba, dm, 1, 0);
            auto const& rhs_a = rhs.arrays();
            auto const& exact_a = phi_exact.arrays();
            ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real const AMREX_D_DECL(x = problo[0] + (Real(i)+Real(0.5))*dx[0],
                                        y = problo[1] + (Real(j)+Real(0.5))*dx[1],
                                        z = problo[2] + (Real(k)+Real(0.5))*dx[2]);
                Real const p = AMREX_D_TERM(std::sin(2*pi*x), * std::sin(2*pi*y), * std::sin(2*pi*z));
                exact_a[b](i,j,k) = p;
                rhs_a[b](i,j,k) = -AMREX_SPACEDIM * 4*pi*pi * p;
            });

            // One ghost cell: the solver fills it from the boundary conditions,
            // so that we can apply the discrete Laplacian below.
            MultiFab phi(ba, dm, 1, 1);
            solver.solve(phi, rhs);

            MultiFab resid(ba, dm, 1, 0);
            auto const& resid_a = resid.arrays();
            auto const& phi_a = phi.const_arrays();
            ParallelFor(resid, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto const& p = phi_a[b];
                Real lap = AMREX_D_TERM(
                    (p(i-1,j,k) - Real(2)*p(i,j,k) + p(i+1,j,k)) / (dx[0]*dx[0]),
                  + (p(i,j-1,k) - Real(2)*p(i,j,k) + p(i,j+1,k)) / (dx[1]*dx[1]),
                  + (p(i,j,k-1) - Real(2)*p(i,j,k) + p(i,j,k+1)) / (dx[2]*dx[2]));
                resid_a[b](i,j,k) = lap - rhs_a[b](i,j,k);
            });
            Real const res = resid.norm0() / rhs.norm0();
            MultiFab::Subtract(phi_exact, phi, 0, 0, 1, 0);
            Real const err = phi_exact.norm0();

            amrex::Print() << "Periodic (FFT::Poisson)\n"
                           << "    relative residual " << res << "\n"
                           << "    max error         " << err << "\n";
            if (res > res_tol) { status = 1; }
        }

        // 2. Homogeneous Dirichlet (odd) and Neumann (even) boundaries, using
        //    sine and cosine transforms.
        {
            Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> bc;
            bc[0] = {FFT::Boundary::odd,  FFT::Boundary::odd};  // Dirichlet on both sides
            bc[1] = {FFT::Boundary::even, FFT::Boundary::even}; // Neumann on both sides
#if (AMREX_SPACEDIM == 3)
            bc[2] = {FFT::Boundary::odd,  FFT::Boundary::even}; // Dirichlet at zlo, Neumann at zhi
#endif
            Geometry geom(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
            FFT::Poisson<MultiFab> solver(geom, bc);
            auto const problo = geom.ProbLoArray();
            auto const dx = geom.CellSizeArray();

            // phi = sin(pi x) cos(pi y) sin(pi z / 2), lap(phi) = -(1 + 1 + 1/4) pi^2 phi
            MultiFab rhs(ba, dm, 1, 0);
            MultiFab phi_exact(ba, dm, 1, 0);
            auto const& rhs_a = rhs.arrays();
            auto const& exact_a = phi_exact.arrays();
            ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real const AMREX_D_DECL(x = problo[0] + (Real(i)+Real(0.5))*dx[0],
                                        y = problo[1] + (Real(j)+Real(0.5))*dx[1],
                                        z = problo[2] + (Real(k)+Real(0.5))*dx[2]);
                Real const p = AMREX_D_TERM(std::sin(pi*x), * std::cos(pi*y), * std::sin(Real(0.5)*pi*z));
                exact_a[b](i,j,k) = p;
                rhs_a[b](i,j,k) = -(AMREX_D_TERM(Real(1), + Real(1), + Real(0.25))) * pi*pi * p;
            });

            MultiFab phi(ba, dm, 1, 1);
            solver.solve(phi, rhs);

            MultiFab resid(ba, dm, 1, 0);
            auto const& resid_a = resid.arrays();
            auto const& phi_a = phi.const_arrays();
            ParallelFor(resid, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto const& p = phi_a[b];
                Real lap = AMREX_D_TERM(
                    (p(i-1,j,k) - Real(2)*p(i,j,k) + p(i+1,j,k)) / (dx[0]*dx[0]),
                  + (p(i,j-1,k) - Real(2)*p(i,j,k) + p(i,j+1,k)) / (dx[1]*dx[1]),
                  + (p(i,j,k-1) - Real(2)*p(i,j,k) + p(i,j,k+1)) / (dx[2]*dx[2]));
                resid_a[b](i,j,k) = lap - rhs_a[b](i,j,k);
            });
            Real const res = resid.norm0() / rhs.norm0();
            MultiFab::Subtract(phi_exact, phi, 0, 0, 1, 0);
            Real const err = phi_exact.norm0();

            amrex::Print() << "Dirichlet/Neumann (FFT::Poisson)\n"
                           << "    relative residual " << res << "\n"
                           << "    max error         " << err << "\n";
            if (res > res_tol) { status = 1; }
        }

#if (AMREX_SPACEDIM == 3)
        // 3. FFT in x and y (periodic here), tridiagonal solve in z (Neumann at
        //    zlo, Dirichlet at zhi). PoissonHybrid also accepts a non-uniform
        //    dz: solver.solve(phi, rhs, dz).
        {
            Array<std::pair<FFT::Boundary,FFT::Boundary>,AMREX_SPACEDIM> bc;
            bc[0] = {FFT::Boundary::periodic, FFT::Boundary::periodic};
            bc[1] = {FFT::Boundary::periodic, FFT::Boundary::periodic};
            bc[2] = {FFT::Boundary::even,     FFT::Boundary::odd}; // Neumann at zlo, Dirichlet at zhi
            Geometry geom(domain, rb, CoordSys::cartesian, {1,1,0});
            FFT::PoissonHybrid<MultiFab> solver(geom, bc);
            auto const problo = geom.ProbLoArray();
            auto const dx = geom.CellSizeArray();

            // phi = sin(2 pi x) sin(2 pi y) cos(pi z / 2), lap(phi) = -(4 + 4 + 1/4) pi^2 phi
            MultiFab rhs(ba, dm, 1, 0);
            MultiFab phi_exact(ba, dm, 1, 0);
            auto const& rhs_a = rhs.arrays();
            auto const& exact_a = phi_exact.arrays();
            ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real const x = problo[0] + (Real(i)+Real(0.5))*dx[0];
                Real const y = problo[1] + (Real(j)+Real(0.5))*dx[1];
                Real const z = problo[2] + (Real(k)+Real(0.5))*dx[2];
                Real const p = std::sin(2*pi*x) * std::sin(2*pi*y) * std::cos(Real(0.5)*pi*z);
                exact_a[b](i,j,k) = p;
                rhs_a[b](i,j,k) = -Real(8.25) * pi*pi * p;
            });

            MultiFab phi(ba, dm, 1, 1);
            solver.solve(phi, rhs);

            MultiFab resid(ba, dm, 1, 0);
            auto const& resid_a = resid.arrays();
            auto const& phi_a = phi.const_arrays();
            ParallelFor(resid, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto const& p = phi_a[b];
                Real lap = (p(i-1,j,k) - Real(2)*p(i,j,k) + p(i+1,j,k)) / (dx[0]*dx[0])
                         + (p(i,j-1,k) - Real(2)*p(i,j,k) + p(i,j+1,k)) / (dx[1]*dx[1])
                         + (p(i,j,k-1) - Real(2)*p(i,j,k) + p(i,j,k+1)) / (dx[2]*dx[2]);
                resid_a[b](i,j,k) = lap - rhs_a[b](i,j,k);
            });
            Real const res = resid.norm0() / rhs.norm0();
            MultiFab::Subtract(phi_exact, phi, 0, 0, 1, 0);
            Real const err = phi_exact.norm0();

            amrex::Print() << "Periodic x,y + Neumann/Dirichlet z (FFT::PoissonHybrid)\n"
                           << "    relative residual " << res << "\n"
                           << "    max error         " << err << "\n";
            if (res > res_tol) { status = 1; }
        }
#endif

        if (status == 0) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
        }
    }
    amrex::Finalize();
    return status;
}
