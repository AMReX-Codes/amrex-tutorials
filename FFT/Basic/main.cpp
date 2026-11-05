#include <AMReX.H>
#include <AMReX_FFT.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <type_traits>

using namespace amrex;

namespace {

// Test field: a smooth low-wavenumber part plus a high-wavenumber "noise" part
// that the low-pass filter removes.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real low_part (AMREX_D_DECL(Real x, Real y, Real z))
{
    constexpr Real tpi = Real(2)*Math::pi<Real>();
    return AMREX_D_TERM(std::sin(tpi*x), *std::cos(tpi*Real(2)*y), *std::cos(tpi*z));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real dlow_dx (AMREX_D_DECL(Real x, Real y, Real z))
{
    constexpr Real tpi = Real(2)*Math::pi<Real>();
    return tpi * AMREX_D_TERM(std::cos(tpi*x), *std::cos(tpi*Real(2)*y), *std::cos(tpi*z));
}

constexpr int khigh = 8;
constexpr Real ahigh = Real(0.25);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real high_part (AMREX_D_DECL(Real x, Real y, Real z))
{
    constexpr Real tk = Real(2)*Math::pi<Real>()*Real(khigh);
    return ahigh * AMREX_D_TERM(std::cos(tk*x), *std::sin(tk*y), *std::cos(tk*z));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real dhigh_dx (AMREX_D_DECL(Real x, Real y, Real z))
{
    constexpr Real tk = Real(2)*Math::pi<Real>()*Real(khigh);
    return -ahigh * tk * AMREX_D_TERM(std::sin(tk*x), *std::sin(tk*y), *std::cos(tk*z));
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int kcut = 4; // low-pass filter cutoff wavenumber
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("kcut", kcut);
        }

        // Periodic unit domain [0,1)^D
        Box domain(IntVect(0), IntVect(n_cell-1));
        RealBox rb({AMREX_D_DECL(Real(0),Real(0),Real(0))},
                   {AMREX_D_DECL(Real(1),Real(1),Real(1))});
        Geometry geom(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        // The test field, its x-derivative and its low-wavenumber part at
        // the cell centers.
        MultiFab phi(ba, dm, 1, 0);
        MultiFab dphi_exact(ba, dm, 1, 0);
        MultiFab phi_low_exact(ba, dm, 1, 0);
        {
            auto const problo = geom.ProbLoArray();
            auto const dx = geom.CellSizeArray();
            auto const& phi_a = phi.arrays();
            auto const& dphi_a = dphi_exact.arrays();
            auto const& low_a = phi_low_exact.arrays();
            ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                AMREX_D_TERM(Real x = problo[0] + (Real(i)+Real(0.5))*dx[0];,
                             Real y = problo[1] + (Real(j)+Real(0.5))*dx[1];,
                             Real z = problo[2] + (Real(k)+Real(0.5))*dx[2];)
                phi_a[b](i,j,k) = low_part(AMREX_D_DECL(x,y,z)) + high_part(AMREX_D_DECL(x,y,z));
                dphi_a[b](i,j,k) = dlow_dx(AMREX_D_DECL(x,y,z)) + dhigh_dx(AMREX_D_DECL(x,y,z));
                low_a[b](i,j,k) = low_part(AMREX_D_DECL(x,y,z));
            });
        }

        // R2C does the expensive setup (plans, communication patterns) once;
        // keep it around and reuse it.
        FFT::R2C<Real,FFT::Direction::both> r2c(geom.Domain());

        // The real-to-complex transform stores only half of the spectrum in x:
        // (nx/2+1, ny, nz) complex values, distributed over the ranks.
        auto const& [cba, cdm] = r2c.getSpectralDataLayout();
        cMultiFab phi_hat(cba, cdm, 1, 0);
        amrex::Print() << "Real domain     " << geom.Domain() << "\n"
                       << "Spectral domain " << cba.minimalBox() << "\n";

        // 1. Round trip. FFTs are unnormalized, so forward followed by
        //    backward multiplies the data by the number of points.
        MultiFab phi2(ba, dm, 1, 0);
        r2c.forward(phi, phi_hat);
        r2c.backward(phi_hat, phi2);
        phi2.mult(r2c.scalingFactor()); // 1/(nx*ny*nz)
        MultiFab::Subtract(phi2, phi, 0, 0, 1, 0);
        Real err_roundtrip = phi2.norm0();
        amrex::Print() << "Round trip error:        " << err_roundtrip << "\n";

        int const nx = n_cell;
        int const ny = (AMREX_SPACEDIM >= 2) ? n_cell : 1;
        int const nz = (AMREX_SPACEDIM == 3) ? n_cell : 1;
        Real const scale = r2c.scalingFactor();

        // 2. Spectral derivative d/dx: multiply by i*kx. forwardThenBackward
        //    calls the lambda at every spectral point (i,j,k); it avoids the
        //    communication of copying the spectral data to a cMultiFab.
        MultiFab dphi(ba, dm, 1, 0);
        r2c.forwardThenBackward(phi, dphi,
            [=] AMREX_GPU_DEVICE (int i, int, int, GpuComplex<Real>& sp)
            {
                // i is in [0,nx/2]. Zero the Nyquist mode, it has no sign.
                Real kx = (2*i == nx) ? Real(0) : Real(2)*Math::pi<Real>()*Real(i);
                sp *= GpuComplex<Real>(Real(0), kx*scale);
            });
        Real const dphi_max = dphi_exact.norm0();
        MultiFab::Subtract(dphi_exact, dphi, 0, 0, 1, 0);
        Real err_deriv = dphi_exact.norm0() / dphi_max;
        amrex::Print() << "Spectral d/dx rel. error: " << err_deriv << "\n";

        // 3. Low-pass filter: zero all modes with |k| > kcut. In y and z, indices
        //    above n/2 are negative wavenumbers.
        MultiFab phi_low(ba, dm, 1, 0);
        r2c.forwardThenBackward(phi, phi_low,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, GpuComplex<Real>& sp)
            {
                int ky = (j <= ny/2) ? j : j-ny;
                int kz = (k <= nz/2) ? k : k-nz;
                int k2 = i*i + ky*ky + kz*kz;
                sp *= (k2 > kcut*kcut) ? Real(0) : scale;
            });
        MultiFab::Subtract(phi_low_exact, phi_low, 0, 0, 1, 0);
        Real err_filter = phi_low_exact.norm0();
        amrex::Print() << "Low-pass filter error:   " << err_filter << "\n";

        Real const tol = std::is_same_v<Real,float> ? Real(1.e-4) : Real(1.e-10);
        if (err_roundtrip > tol || err_deriv > tol || err_filter > tol) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
