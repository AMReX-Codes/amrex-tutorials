#include <AMReX.H>
#include <AMReX_FFT.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        std::string outfile = "spectrum.txt";
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("outfile", outfile);
        }
        AMREX_ALWAYS_ASSERT(n_cell > 16);

        // Periodic unit domain [0,1)^D. Wavenumbers are integers in units of 2*pi.
        Box domain(IntVect(0), IntVect(n_cell-1));
        RealBox rb({AMREX_D_DECL(Real(0),Real(0),Real(0))},
                   {AMREX_D_DECL(Real(1),Real(1),Real(1))});
        Geometry geom(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        // Synthetic velocity field made of three Fourier modes with integer
        // |k| = 2, 5 and 8. Each mode a*cos(2 pi k.x) has kinetic energy a^2/4.
        //   u_x  = A cos(2 pi 2y)
        //   u_y  = B sin(2 pi (3x + 4y))     (3D: 3x + 4z)
        //   u_D += C cos(2 pi 8x)            (last component)
        constexpr Real A = Real(1.0), B = Real(0.5), C = Real(0.25);
        constexpr int ncomp = AMREX_SPACEDIM;
        MultiFab vel(ba, dm, ncomp, 0);
        {
            auto const dx = geom.CellSizeArray();
            auto const& va = vel.arrays();
            ParallelFor(vel, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                constexpr Real tpi = Real(2)*Math::pi<Real>();
                AMREX_D_TERM(Real x = (Real(i)+Real(0.5))*dx[0];,
                             Real y = (Real(j)+Real(0.5))*dx[1];,
                             Real z = (Real(k)+Real(0.5))*dx[2];)
#if (AMREX_SPACEDIM == 3)
                va[b](i,j,k,0) = A*std::cos(tpi*Real(2)*y);
                va[b](i,j,k,1) = B*std::sin(tpi*(Real(3)*x + Real(4)*z));
                va[b](i,j,k,2) = C*std::cos(tpi*Real(8)*x);
#else
                va[b](i,j,k,0) = A*std::cos(tpi*Real(2)*y);
                va[b](i,j,k,1) = B*std::sin(tpi*(Real(3)*x + Real(4)*y))
                    +            C*std::cos(tpi*Real(8)*x);
#endif
            });
        }

        // Batched forward transform: all velocity components in one call.
        FFT::Info info{};
        info.setBatchSize(ncomp);
        FFT::R2C<Real,FFT::Direction::forward> r2c(geom.Domain(), info);
        auto const& [cba, cdm] = r2c.getSpectralDataLayout();
        cMultiFab vel_hat(cba, cdm, ncomp, 0);
        r2c.forward(vel, vel_hat);

        // Shell-averaged spectrum E(k) = sum over |k| in [k-1/2,k+1/2) of
        // |u_hat|^2/(2 N^2), so that sum_k E(k) = <|u|^2>/2 (Parseval).
        int const nx = n_cell;
        int const ny = (AMREX_SPACEDIM >= 2) ? n_cell : 1;
        int const nz = (AMREX_SPACEDIM == 3) ? n_cell : 1;
        int const nshells = int(std::sqrt(Real(AMREX_SPACEDIM))*Real(n_cell/2)) + 2;
        Real const npts = geom.Domain().d_numPts();
        Real const fac = Real(0.5) / (npts*npts);

        // Histogram on the device: each spectral point adds its energy to
        // its shell with an atomic add.
        Gpu::DeviceVector<Real> ek_d(nshells, Real(0));
        Real* pek = ek_d.data();
        auto const& vh = vel_hat.const_arrays();
        ParallelFor(vel_hat, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            // Only kx >= 0 is stored. Indices above n/2 are negative wavenumbers.
            int ky = (j <= ny/2) ? j : j-ny;
            int kz = (k <= nz/2) ? k : k-nz;
            int shell = int(std::round(std::sqrt(Real(i*i + ky*ky + kz*kz))));
            Real e = 0;
            for (int n = 0; n < ncomp; ++n) {
                e += amrex::norm(vh[b](i,j,k,n)); // |u_hat|^2
            }
            // Hermitian symmetry: u_hat(-k) = conj(u_hat(k)). Except for kx = 0
            // and the Nyquist mode, each stored point stands for two modes.
            if (i > 0 && 2*i != nx) { e *= Real(2); }
            HostDevice::Atomic::Add(pek+shell, e*fac);
        });

        // Bring the local sums to the host and add them across ranks.
        std::vector<Real> ek(nshells);
        Gpu::copy(Gpu::deviceToHost, ek_d.begin(), ek_d.end(), ek.begin());
        ParallelDescriptor::ReduceRealSum(ek.data(), nshells);

        if (ParallelDescriptor::IOProcessor()) {
            std::ofstream ofs(outfile);
            ofs << "# k E(k)\n";
            for (int s = 0; s < nshells; ++s) {
                ofs << s << " " << ek[s] << "\n";
            }
        }

        // Checks: energy in the expected shells, and Parseval's theorem.
        Real ke_phys = Real(0.5) * MultiFab::Dot(vel, 0, vel, 0, ncomp, 0) / npts;
        Real ke_spec = 0;
        for (auto e : ek) { ke_spec += e; }

        std::vector<Real> ek_exact(nshells, Real(0));
        ek_exact[2] = A*A/Real(4);
        ek_exact[5] = B*B/Real(4);
        ek_exact[8] = C*C/Real(4);
        Real ke_exact = (A*A + B*B + C*C) / Real(4);

        Real err_shell = 0;
        amrex::Print() << "Nonzero shells:\n";
        for (int s = 0; s < nshells; ++s) {
            err_shell = std::max(err_shell, std::abs(ek[s]-ek_exact[s]) / ke_exact);
            if (ek[s] > Real(1.e-3)*ke_exact) {
                amrex::Print() << "  E(" << s << ") = " << ek[s]
                               << "  (exact " << ek_exact[s] << ")\n";
            }
        }
        Real err_parseval = std::abs(ke_spec-ke_phys) / ke_phys;
        amrex::Print() << "Kinetic energy: physical space " << ke_phys
                       << ", spectral space " << ke_spec << "\n"
                       << "Parseval relative error: " << err_parseval << "\n"
                       << "Max shell error:         " << err_shell << "\n"
                       << "Spectrum written to " << outfile << "\n";

        Real const tol = std::is_same_v<Real,float> ? Real(1.e-4) : Real(1.e-10);
        if (err_parseval > tol || err_shell > tol) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
