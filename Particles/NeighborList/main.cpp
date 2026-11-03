#include "MDParticleContainer.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <type_traits>

using namespace amrex::literals;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 8;
        int max_grid_size = 4;
        int nppc = 2;            // particles per cell in each direction
        int nsteps = 100;
        int num_rebuild = 10;    // steps between neighbor list rebuilds
        int print_int = 20;
        amrex::Real v_thermal = 1.0_rt;
        amrex::Real cfl = 0.05_rt;
        MDParams params;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nppc", nppc);
            pp.query("nsteps", nsteps);
            pp.query("num_rebuild", num_rebuild);
            pp.query("print_int", print_int);
            pp.query("v_thermal", v_thermal);
            pp.query("cfl", cfl);
            pp.query("sigma", params.sigma);
            pp.query("stiffness", params.stiffness);
        }

        // Periodic domain with unit cells. The neighbor list radius is one
        // cell, so one layer of neighbor cells is enough.
        const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        const amrex::RealBox real_box({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                      {AMREX_D_DECL(amrex::Real(n_cell), amrex::Real(n_cell),
                                                    amrex::Real(n_cell))});
        const amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
        const amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);
        params.r_list = geom.CellSize(0);
        AMREX_ALWAYS_ASSERT(params.sigma < params.r_list);

        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        const amrex::DistributionMapping dm(ba);

        const int num_neighbor_cells = 1;
        MDParticleContainer pc(geom, dm, ba, num_neighbor_cells);
        pc.initParticles(nppc, v_thermal);

        const CheckPair check_pair{static_cast<amrex::ParticleReal>(params.r_list*params.r_list)};

        // Copy particles near tile boundaries to their neighbors, build the
        // lists, and compute the initial forces.
        pc.fillNeighbors();
        pc.buildNeighborList(check_pair);
        pc.computeForces(params);

        const MDDiagnostics d0 = pc.diagnostics();
        const amrex::Real e0 = d0.kinetic + d0.potential;

        // At the initial maximum speed, a particle moves cfl*sigma per step.
        // Assuming speeds stay near that, the relative motion of a pair between
        // rebuilds must stay below the skin r_list - sigma.
        const amrex::Real dt = cfl * params.sigma / d0.max_speed;
        AMREX_ALWAYS_ASSERT(2*num_rebuild*cfl*params.sigma < params.r_list - params.sigma);
        amrex::Print() << "Number of particles: " << d0.np << ", dt = " << dt << "\n";

        // Velocity Verlet time stepping.
        for (int step = 1; step <= nsteps; ++step)
        {
            pc.kickDrift(0.5_rt*dt, dt);

            if (step % num_rebuild == 0) {
                // Particles moved less than one cell since the last rebuild.
                pc.Redistribute(0, -1, 0, 1);
                pc.fillNeighbors();
                pc.buildNeighborList(check_pair);
            } else {
                // Only refresh the data of the existing neighbor copies.
                pc.updateNeighbors();
            }

            pc.computeForces(params);
            pc.kickDrift(0.5_rt*dt, 0.0_rt);

            if (print_int > 0 && step % print_int == 0) {
                const MDDiagnostics d = pc.diagnostics();
                amrex::Print() << "Step " << step << ": kinetic " << d.kinetic
                               << ", potential " << d.potential
                               << ", total " << d.kinetic + d.potential << "\n";
            }
        }

        pc.Redistribute(0, -1, 0, 1);
        const MDDiagnostics d1 = pc.diagnostics();
        const amrex::Real e1 = d1.kinetic + d1.potential;

        // Pair forces are equal and opposite, so the total momentum is
        // conserved to round-off. Velocity Verlet keeps the energy error small.
        amrex::Real dp = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            dp = amrex::max(dp, std::abs(d1.momentum[d] - d0.momentum[d]));
        }
        const amrex::Real momentum_error = dp / d0.abs_momentum;
        const amrex::Real energy_error = std::abs(e1 - e0) / e0;
        amrex::Print() << "Relative momentum change: " << momentum_error << "\n"
                       << "Relative energy change:   " << energy_error << "\n";

        constexpr bool is_single = std::is_same_v<amrex::Real, float>
                                || std::is_same_v<amrex::ParticleReal, float>;
        constexpr amrex::Real momentum_tol = is_single ? 1.e-5_rt : 1.e-12_rt;
        const amrex::Real energy_tol = 1.e-2_rt;
        if (d1.np != d0.np || momentum_error > momentum_tol || energy_error > energy_tol) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
