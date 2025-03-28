#include <AMReX.H>
#include <AMReX_FFT.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        Geometry geom;
        MultiFab vx, vy, vz;
        {
            PlotFileData plot("plot");
            geom.define(plot.probDomain(0),
                        RealBox(plot.probLo(), plot.probHi()),
                        plot.coordSys(), {1,1,1});
            vx = plot.get(0, "velx");
            vy = plot.get(0, "vely");
            vz = plot.get(0, "velz");
        }

        cMultiFab cx, cy, cz;
        {
            // Note that the complex Hermitian output array Y has (nx/2+1,ny,nz) elements.
	    // Y[nx-i,j,k] = Y[i,j,k]*
            FFT::R2C<Real,FFT::Direction::forward> fft(geom.Domain());
            auto const& [ba, dm] = fft.getSpectralDataLayout();
            cx.define(ba,dm,1,0);
            cy.define(ba,dm,1,0);
            cz.define(ba,dm,1,0);
            fft.forward(vx, cx);
            fft.forward(vy, cy);
            fft.forward(vz, cz);
        }

        // For simplicity, we assume the domain is a cube, and we are not
        // going to worry about scaling.

        int nx = geom.Domain().length(0);
        int ny = geom.Domain().length(1);
        int nz = geom.Domain().length(2);
        int nk = nx;
        Gpu::DeviceVector<Real> ke(nk,Real(0.0));
        Real* pke = ke.data();
        auto const& cxa = cx.const_arrays();
        auto const& cya = cy.const_arrays();
        auto const& cza = cz.const_arrays();
        ParallelFor(cx, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            int ki = i;
            int kj = (j <= ny/2) ? j : ny-j;
            int kk = (k <= nz/2) ? k : nz-k;
            Real d = std::sqrt(Real(ki*ki+kj*kj+kk*kk));
            int di = int(std::round(d));
            if (di < nk) {
                Real value = amrex::norm(cxa[b](i,j,k))
                    +        amrex::norm(cya[b](i,j,k))
                    +        amrex::norm(cza[b](i,j,k));
		// Account for Hermitian symmetry in x-direction
	        // Hermitian symmetry Y[nx-i,j,k] = Y[i,j,k]*
                if ((i > 0) && (2*i != nx)) {
		    // Multiply by 2 because we have +ki and -ki
                    value *= Real(2.0);
                }
                HostDevice::Atomic::Add(pke+di, value);
            }
        });

#ifdef AMREX_USE_GPU
        Gpu::HostVector<Real> ke_h(ke.size());
        Gpu::copyAsync(Gpu::deviceToHost, ke.begin(), ke.end(), ke_h.begin());
        Gpu::streamSynchronize();
        Real* pke_h = ke_h.data();
#else
        Real* pke_h = pke;
#endif

        ParallelDescriptor::ReduceRealSum(pke_h, nk);

        if (ParallelDescriptor::IOProcessor()) {
            std::ofstream ofs("spectrum.txt");
            for (int i = 0; i < nk; ++i) {
                ofs << i << " " << pke_h[i] << "\n";
            }
        }
    }
    amrex::Finalize();
}
