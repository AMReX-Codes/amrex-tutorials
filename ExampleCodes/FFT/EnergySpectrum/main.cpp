#include <AMReX.H>
#include <AMReX_FFT.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

int main(int argc, char *argv[]) {
  amrex::Initialize(argc, argv);
  {
    std::string infile("plot");
    std::string outfile("spectrum.txt");
    Vector<std::string> vel_labels = {AMREX_D_DECL("velx", "vely", "velz")};
    ParmParse pp;
    pp.query("infile", infile);
    pp.query("outfile", outfile);
    pp.queryarr("vel_labels", vel_labels, 0, 3);
    Geometry geom;
    MultiFab vx, vy, vz;
    Vector<Real> scaleL(3);
    {
      PlotFileData plot(infile);
      geom.define(plot.probDomain(0), RealBox(plot.probLo(), plot.probHi()),
                  plot.coordSys(), {1, 1, 1});
      // Read data, remove mean field
      vx = plot.get(0, vel_labels[0]);
      vx.mult(1. / vx.sum(0, 0));
      vy = plot.get(0, vel_labels[1]);
      vy.mult(1. / vy.sum(0, 0));
      vz = plot.get(0, vel_labels[2]);
      vz.mult(1. / vz.sum(0, 0));
      Real L = std::min(
          {geom.ProbLength(0), geom.ProbLength(1), geom.ProbLength(2)});
      scaleL = {AMREX_D_DECL(1. / std::round(geom.ProbLength(0) / L),
                             1. / std::round(geom.ProbLength(1) / L),
                             1. / std::round(geom.ProbLength(2) / L))};
    }
    cMultiFab cx, cy, cz;
    {
      // Note that the complex Hermitian output array Y has (nx/2+1,ny,nz)
      // elements. Y[nx-i,j,k] = Y[i,j,k]*
      FFT::R2C<Real, FFT::Direction::forward> fft(geom.Domain());
      auto const &[ba, dm] = fft.getSpectralDataLayout();
      cx.define(ba, dm, 1, 0);
      cy.define(ba, dm, 1, 0);
      cz.define(ba, dm, 1, 0);
      fft.forward(vx, cx);
      fft.forward(vy, cy);
      fft.forward(vz, cz);
    }

    // For simplicity, we are not going to worry about scaling.

    int nx = geom.Domain().length(0);
    int ny = geom.Domain().length(1);
    int nz = geom.Domain().length(2);
    int nk = nx;
    Gpu::DeviceVector<Real> ke(nk, Real(0.0));
    Real *pke = ke.data();
    auto const &cxa = cx.const_arrays();
    auto const &cya = cy.const_arrays();
    auto const &cza = cz.const_arrays();
    ParallelFor(cx, [=] AMREX_GPU_DEVICE(int b, int i, int j, int k) {
      int ki = i;
      int kj = (j <= ny / 2) ? j : ny - j;
      int kk = (k <= nz / 2) ? k : nz - k;
      Real d = std::sqrt(Real(ki * ki * scaleL[0] + kj * kj * scaleL[1] +
                              kk * kk * scaleL[2]));
      int di = int(std::round(d));
      if (di < nk) {
        Real value = amrex::norm(cxa[b](i, j, k)) +
                     amrex::norm(cya[b](i, j, k)) +
                     amrex::norm(cza[b](i, j, k));
        // Account for Hermitian symmetry in x-direction
        // Hermitian symmetry Y[nx-i,j,k] = Y[i,j,k]*
        if ((i > 0) && (2 * i != nx)) {
          // Multiply by 2 because we have +ki and -ki
          value *= Real(2.0);
        }
        HostDevice::Atomic::Add(pke + di, value);
      }
    });

#ifdef AMREX_USE_GPU
    Gpu::HostVector<Real> ke_h(ke.size());
    Gpu::copyAsync(Gpu::deviceToHost, ke.begin(), ke.end(), ke_h.begin());
    Gpu::streamSynchronize();
    Real *pke_h = ke_h.data();
#else
    Real *pke_h = pke;
#endif

    ParallelDescriptor::ReduceRealSum(pke_h, nk);

    if (ParallelDescriptor::IOProcessor()) {
      std::ofstream ofs(outfile.c_str());
      for (int i = 0; i < nk; ++i) {
        ofs << i << " " << pke_h[i] << "\n";
      }
    }
  }
  amrex::Finalize();
}
