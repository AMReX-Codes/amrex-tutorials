#include "FieldSolver.H"
#include "PICParticleContainer.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <cmath>
#include <string>
#include <vector>

using namespace amrex::literals;

// Fourier coefficient of E_x for sin(k x): (2/V) * integral of E_x sin(k x).
amrex::Real modeAmplitude (const amrex::MultiFab& efield, const amrex::Geometry& geom,
                           amrex::Real k)
{
    const auto plo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    auto const& e = efield.const_arrays();
    amrex::Real sum = amrex::ParReduce(amrex::TypeList<amrex::ReduceOpSum>{},
                                       amrex::TypeList<amrex::Real>{},
                                       efield, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k3) noexcept -> amrex::GpuTuple<amrex::Real>
        {
            amrex::ignore_unused(j, k3);
            const amrex::Real x = plo[0] + (i + 0.5_rt) * dx[0];
            return { e[b](i,j,k3,0) * std::sin(k * x) };
        });
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    return 2.0_rt * sum / geom.Domain().d_numPts();
}

void writePlotfile (int step, amrex::Real time, const amrex::Geometry& geom,
                    const amrex::MultiFab& rho, const amrex::MultiFab& phi,
                    const amrex::MultiFab& efield, const PICParticleContainer& pc)
{
    amrex::MultiFab plotmf(rho.boxArray(), rho.DistributionMap(), 2+AMREX_SPACEDIM, 0);
    amrex::MultiFab::Copy(plotmf, rho, 0, 0, 1, 0);
    amrex::MultiFab::Copy(plotmf, phi, 0, 1, 1, 0);
    amrex::MultiFab::Copy(plotmf, efield, 0, 2, AMREX_SPACEDIM, 0);
    const std::string pltfile = amrex::Concatenate("plt", step, 5);
    const amrex::Vector<std::string> mesh_names{"rho", "phi", AMREX_D_DECL("Ex", "Ey", "Ez")};
    amrex::WriteSingleLevelPlotfile(pltfile, plotmf, mesh_names, geom, time, step);

    // Particle component names, without the positions.
    const auto& names = pc.GetRealSoANames();
    const amrex::Vector<std::string> particle_names(names.begin() + AMREX_SPACEDIM, names.end());
    pc.WritePlotFile(pltfile, "electrons", particle_names, {});
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        // Normalized units: epsilon0 = electron mass = background density = 1
        // and electron charge = -1, so the plasma frequency is 1.
        const amrex::Real epsilon0 = 1.0_rt;
        const amrex::Real q = -1.0_rt;
        const amrex::Real m = 1.0_rt;
        const amrex::Real n0 = 1.0_rt;
        const amrex::Real omega_p = std::sqrt(n0*q*q / (epsilon0*m));

        amrex::IntVect n_cell(AMREX_D_DECL(32, 8, 8));
        amrex::IntVect nppc(AMREX_D_DECL(4, 1, 1));
        int max_grid_size = 16;
        int nsteps = 100;
        int plot_int = -1;
        int print_int = 10;
        amrex::Real dt = 0.1_rt;           // in units of 1/omega_p
        amrex::Real perturbation = 0.01_rt; // relative density perturbation
        {
            amrex::ParmParse pp;
            // Read 3 values; a 2D build uses the first two.
            std::vector<int> v;
            if (pp.queryarr("n_cell", v)) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { n_cell[d] = v.at(d); }
            }
            if (pp.queryarr("nppc", v)) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { nppc[d] = v.at(d); }
            }
            pp.query("max_grid_size", max_grid_size);
            pp.query("nsteps", nsteps);
            pp.query("plot_int", plot_int);
            pp.query("print_int", print_int);
            pp.query("dt", dt);
            pp.query("perturbation", perturbation);
        }

        // One wavelength fits in x. The cells are cubes.
        const amrex::Real Lx = 2.0_rt * amrex::Math::pi<amrex::Real>();
        const amrex::Real k = 2.0_rt * amrex::Math::pi<amrex::Real>() / Lx;
        amrex::RealBox real_box;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            real_box.setLo(d, 0.0_rt);
            real_box.setHi(d, Lx * amrex::Real(n_cell[d]) / amrex::Real(n_cell[0]));
        }
        const amrex::Box domain(amrex::IntVect(0), n_cell - 1);
        const amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
        const amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        const amrex::DistributionMapping dm(ba);

        // Electrons, displaced by xi0*sin(k x), which gives the density
        // perturbation n0*perturbation*cos(k x) to first order.
        const amrex::Real xi0 = perturbation / k;
        PICParticleContainer electrons(geom, dm, ba, q, m);
        electrons.initLangmuirWave(nppc, n0, k, xi0);
        amrex::Print() << "Number of particles: " << electrons.TotalNumberOfParticles() << "\n";

        // Fixed ions provide a uniform neutralizing background.
        const amrex::Real rho_ion = -electrons.totalCharge() / geom.ProbDomain().volume();

        // Particles must not move more than one cell per step.
        AMREX_ALWAYS_ASSERT(omega_p * xi0 * dt < geom.CellSize(0));

        // Mesh data. rho and efield need ghost cells for deposition and gather,
        // and phi for the gradient.
        amrex::MultiFab rho(ba, dm, 1, 1);
        amrex::MultiFab phi(ba, dm, 1, 1);
        amrex::MultiFab efield(ba, dm, AMREX_SPACEDIM, 1);
        phi.setVal(0.0);
        FieldSolver solver(geom, ba, dm, epsilon0);

        std::vector<amrex::Real> mode_history;
        amrex::Real energy0 = 0;
        amrex::Real max_energy_error = 0;

        for (int step = 0; step <= nsteps; ++step)
        {
            const amrex::Real time = step * dt;

            // Charge deposition, field solve and gather at time t^n.
            electrons.depositCharge(rho);
            rho.plus(rho_ion, 0, 1, 0);
            solver.solvePhi(phi, rho);
            solver.computeE(efield, phi);
            electrons.gatherField(efield);

            // Leapfrog: velocities live at half steps. Start with v^{-1/2}.
            if (step == 0) { electrons.pushVelocity(-0.5_rt*dt); }

            const amrex::Real field_energy = solver.fieldEnergy(efield);
            const amrex::Real ke_old = electrons.kineticEnergy();
            electrons.pushVelocity(dt);
            const amrex::Real ke_new = electrons.kineticEnergy();
            const amrex::Real energy = field_energy + 0.5_rt*(ke_old + ke_new);
            if (step == 0) { energy0 = energy; }
            max_energy_error = amrex::max(max_energy_error, std::abs(energy - energy0) / energy0);

            mode_history.push_back(modeAmplitude(efield, geom, k));

            if (print_int > 0 && step % print_int == 0) {
                amrex::Print() << "Step " << step << ", t = " << time
                               << ": E mode " << mode_history.back()
                               << ", field energy " << field_energy
                               << ", total energy " << energy << "\n";
            }
            if (plot_int > 0 && step % plot_int == 0) {
                writePlotfile(step, time, geom, rho, phi, efield, electrons);
            }

            if (step < nsteps) { electrons.pushPosition(dt); }
        }

        // The mode amplitude oscillates as cos(omega_p t). Measure the frequency
        // from its zero crossings, and compare the initial amplitude with the
        // linear theory value -q n0 xi0 / epsilon0.
        std::vector<amrex::Real> crossings;
        for (int n = 1; n < static_cast<int>(mode_history.size()); ++n) {
            const amrex::Real a = mode_history[n-1];
            const amrex::Real b = mode_history[n];
            if (a*b < 0) { crossings.push_back((n - 1 + a/(a-b)) * dt); }
        }

        const amrex::Real expected_amplitude = -q * n0 * xi0 / epsilon0;
        const amrex::Real amplitude_error = std::abs(mode_history[0] / expected_amplitude - 1.0_rt);
        amrex::Print() << "Initial E mode amplitude: " << mode_history[0]
                       << " (linear theory " << expected_amplitude << ")\n";

        amrex::Real frequency_error = 1.0_rt;
        if (crossings.size() >= 2) {
            const amrex::Real half_period = (crossings.back() - crossings.front())
                / amrex::Real(crossings.size() - 1);
            const amrex::Real omega = amrex::Math::pi<amrex::Real>() / half_period;
            frequency_error = std::abs(omega / omega_p - 1.0_rt);
            amrex::Print() << "Measured frequency: " << omega
                           << " (plasma frequency " << omega_p << ")\n";
        } else {
            amrex::Print() << "Not enough zero crossings; increase nsteps\n";
        }
        amrex::Print() << "Relative errors: frequency " << frequency_error
                       << ", amplitude " << amplitude_error
                       << ", max energy change " << max_energy_error << "\n";

        // Finite grid and time step effects are below one percent here.
        const amrex::Real tol = 0.02_rt;
        if (frequency_error > tol || amplitude_error > tol || max_energy_error > tol) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
