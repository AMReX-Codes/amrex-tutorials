#include <AMReX.H>
#include <AMReX_FFT.H>
#include <AMReX_FFT_Poisson.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace {

#if (AMREX_SPACEDIM == 2)
// E1(x) + ln(x), where E1 is the exponential integral. It is finite at x = 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real E1_plus_log (Real x)
{
    constexpr Real euler_gamma = Real(0.5772156649015329);
    if (x <= Real(1)) {
        // Power series: -gamma - sum_n (-x)^n / (n n!)
        Real sum = 0;
        Real term = 1;
        for (int n = 1; n <= 20; ++n) {
            term *= -x / Real(n);
            sum += term / Real(n);
        }
        return -euler_gamma - sum;
    } else if (x < Real(50)) {
        // Continued fraction, evaluated from the tail
        Real t = x + Real(81);
        for (int n = 40; n >= 1; --n) {
            t = x + Real(2*n-1) - Real(n*n)/t;
        }
        return std::exp(-x)/t + std::log(x);
    } else {
        return std::log(x);
    }
}
#endif

// Free-space potential of a Gaussian charge with total charge 1 and width
// sigma, solving lap(phi) = rho.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_potential (Real r, Real sigma)
{
#if (AMREX_SPACEDIM == 3)
    // phi = -erf(r/(sqrt(2) sigma)) / (4 pi r)
    constexpr Real fourpi = Real(4)*Math::pi<Real>();
    if (r < Real(1.e-6)*sigma) {
        return -std::sqrt(Real(2)/Math::pi<Real>()) / (fourpi*sigma);
    } else {
        return -std::erf(r/(std::sqrt(Real(2))*sigma)) / (fourpi*r);
    }
#else
    // phi = (ln(r) + E1(r^2/(2 sigma^2))/2) / (2 pi)
    Real two_sigma2 = Real(2)*sigma*sigma;
    Real x = r*r/two_sigma2;
    return (std::log(two_sigma2) + E1_plus_log(x)) / (Real(4)*Math::pi<Real>());
#endif
}

Real rel_error (MultiFab const& phi, MultiFab const& phi_exact)
{
    MultiFab d(phi.boxArray(), phi.DistributionMap(), 1, 0);
    MultiFab::LinComb(d, Real(1), phi, 0, Real(-1), phi_exact, 0, 0, 1, 0);
    return d.norm0() / phi_exact.norm0();
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        Real sigma = Real(0.1); // width of the Gaussian charge
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("sigma", sigma);
        }

        // Domain [-1/2,1/2]^D. The charge is (almost) zero at the domain
        // boundary; outside the domain there is nothing but free space.
        Box domain(IntVect(0), IntVect(n_cell-1));
        RealBox rb({AMREX_D_DECL(Real(-0.5),Real(-0.5),Real(-0.5))},
                   {AMREX_D_DECL(Real( 0.5),Real( 0.5),Real( 0.5))});
        Geometry geom(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        auto const problo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();

        MultiFab rho(ba, dm, 1, 0);
        MultiFab phi_exact(ba, dm, 1, 0);
        {
            // Normalized Gaussian: rho = exp(-r^2/(2 sigma^2)) / (2 pi sigma^2)^(D/2)
            Real norm = std::pow(Real(2)*Math::pi<Real>()*sigma*sigma,
                                 Real(-0.5)*Real(AMREX_SPACEDIM));
            auto const& rho_a = rho.arrays();
            auto const& exact_a = phi_exact.arrays();
            ParallelFor(rho, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                AMREX_D_TERM(Real x = problo[0] + (Real(i)+Real(0.5))*dx[0];,
                             Real y = problo[1] + (Real(j)+Real(0.5))*dx[1];,
                             Real z = problo[2] + (Real(k)+Real(0.5))*dx[2];)
                Real r2 = AMREX_D_TERM(x*x, +y*y, +z*z);
                rho_a[b](i,j,k) = norm * std::exp(-r2/(Real(2)*sigma*sigma));
                exact_a[b](i,j,k) = exact_potential(std::sqrt(r2), sigma);
            });
        }

        // The open-BC solution is the convolution phi = G * rho with the
        // free-space Green's function G. OpenBCSolver computes it with FFTs on
        // a domain doubled in each direction (Hockney's method). The user
        // provides G integrated over a cell, here with 2^D-point Gauss
        // quadrature, which also avoids the singularity at r = 0.
        FFT::OpenBCSolver<Real> solver(geom.Domain());
        {
            auto const lo = geom.Domain().smallEnd();
            Real const gfac = Real(1)/std::sqrt(Real(12));
#if (AMREX_SPACEDIM == 3)
            Real const fac = -dx[0]*dx[1]*dx[2] / (Real(8) * Real(4)*Math::pi<Real>());
#else
            Real const fac = dx[0]*dx[1] / (Real(4) * Real(2)*Math::pi<Real>());
#endif
            // (i,j,k) - lo is the separation between two cells in cell units.
            solver.setGreensFunction([=] AMREX_GPU_DEVICE (int i, int j, int k) -> Real
            {
                amrex::ignore_unused(k); // unused in 2D
                IntVect const iv(AMREX_D_DECL(i,j,k));
                Real r = 0;
#if (AMREX_SPACEDIM == 3)
                for (int gz = -1; gz <= 1; gz += 2) {
#endif
                for (int gy = -1; gy <= 1; gy += 2) {
                for (int gx = -1; gx <= 1; gx += 2) {
                    IntVect const g(AMREX_D_DECL(gx,gy,gz));
                    Real d2 = 0;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        Real xg = (Real(iv[d]-lo[d]) + Real(g[d])*gfac) * dx[d];
                        d2 += xg*xg;
                    }
#if (AMREX_SPACEDIM == 3)
                    r += Real(1)/std::sqrt(d2);   // G = -1/(4 pi r)
#else
                    r += Real(0.5)*std::log(d2);  // G = ln(r)/(2 pi)
#endif
                }}
#if (AMREX_SPACEDIM == 3)
                }
#endif
                return fac * r;
            });
        }

        MultiFab phi(ba, dm, 1, 0);
        solver.solve(phi, rho);

        // Second-order accurate: the error scales with (dx/sigma)^2.
        Real const tol = Real(0.15) * (dx[0]/sigma) * (dx[0]/sigma);
        Real err = rel_error(phi, phi_exact);
        amrex::Print() << "FFT::OpenBCSolver   relative error: " << err << "\n";
        bool pass = err < tol;

#if (AMREX_SPACEDIM == 3)
        // In 3D, FFT::PoissonOpenBC wraps the steps above with the same
        // Green's function.
        FFT::PoissonOpenBC<MultiFab> poisson(geom);
        MultiFab phi2(ba, dm, 1, 0);
        poisson.solve(phi2, rho);
        Real err2 = rel_error(phi2, phi_exact);
        amrex::Print() << "FFT::PoissonOpenBC  relative error: " << err2 << "\n";
        pass = pass && err2 < tol;
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
