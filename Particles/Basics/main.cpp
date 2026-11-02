#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <cmath>
#include <string>
#include <type_traits>

using namespace amrex::literals;

// Real components of each particle. In a pure struct-of-arrays (SoA)
// container, the first AMREX_SPACEDIM real components are the positions.
struct RealIdx
{
    enum {
        AMREX_D_DECL(x, y, z),
        AMREX_D_DECL(vx, vy, vz),
        w,      // weight, e.g., mass
        ncomps
    };
};

// All components, including position and id/cpu, are stored as separate arrays.
using PC = amrex::ParticleContainerPureSoA<RealIdx::ncomps, 0>;
using ParticleType = PC::ParticleType;

struct Stats
{
    amrex::Long np = 0;
    amrex::Real weight = 0;   // sum of w
    amrex::Real energy = 0;   // sum of w |v|^2 / 2
    amrex::Real path = 0;     // sum of the runtime component "dist"
    amrex::Real speed = 0;    // sum of |v|
};

// Global particle reductions.
Stats computeStats (PC const& pc)
{
    using ConstPTD = PC::ConstPTDType;
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                     amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_ops;
    auto r = amrex::ParticleReduce<amrex::ReduceData<amrex::Real, amrex::Real,
                                                     amrex::Real, amrex::Real>>(
        pc, [=] AMREX_GPU_DEVICE (const ConstPTD& ptd, int i) noexcept
            -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real, amrex::Real>
        {
            const amrex::Real w = ptd.rdata(RealIdx::w)[i];
            amrex::Real v2 = 0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                const amrex::Real v = ptd.rdata(RealIdx::vx+d)[i];
                v2 += v*v;
            }
            // Runtime components come after the compile-time ones.
            const amrex::Real dist = ptd.m_runtime_rdata[0][i];
            return {w, 0.5_rt*w*v2, dist, std::sqrt(v2)};
        }, reduce_ops);

    // ParticleReduce is local to this process; sum over MPI ranks.
    amrex::Real sums[4] = {amrex::get<0>(r), amrex::get<1>(r),
                           amrex::get<2>(r), amrex::get<3>(r)};
    amrex::ParallelDescriptor::ReduceRealSum(sums, 4);

    return Stats{pc.TotalNumberOfParticles(), sums[0], sums[1], sums[2], sums[3]};
}

// Add nppc randomly placed particles in every cell, directly on the device.
void initParticles (PC& pc, int nppc, amrex::Real vmax)
{
    const int lev = 0;
    const auto plo = pc.Geom(lev).ProbLoArray();
    const auto dx = pc.Geom(lev).CellSizeArray();
    const int myproc = amrex::ParallelDescriptor::MyProc();

    for (amrex::MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box bx = mfi.tilebox();
        const amrex::Long np = bx.numPts() * nppc;

        // Reserve a block of unique particle ids on this process.
        const amrex::Long id0 = ParticleType::NextID();
        ParticleType::NextID(id0 + np);

        // DefineAndReturnParticleTile also allocates the runtime components.
        auto& tile = pc.DefineAndReturnParticleTile(lev, mfi);
        tile.resize(np);
        auto ptd = tile.getParticleTileData();

        amrex::ParallelForRNG(np,
        [=] AMREX_GPU_DEVICE (amrex::Long ip, amrex::RandomEngine const& engine) noexcept
        {
            const amrex::IntVect iv = bx.atOffset(ip / nppc);
            ptd.idcpu(ip) = amrex::SetParticleIDandCPU(id0 + ip, myproc);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                ptd.pos(d, ip) = static_cast<amrex::ParticleReal>(
                    plo[d] + (iv[d] + amrex::Random(engine)) * dx[d]);
                ptd.rdata(RealIdx::vx+d)[ip] = static_cast<amrex::ParticleReal>(
                    vmax * (2.0_rt*amrex::Random(engine) - 1.0_rt));
            }
            ptd.rdata(RealIdx::w)[ip] = static_cast<amrex::ParticleReal>(
                0.5_rt + amrex::Random(engine));
            ptd.m_runtime_rdata[0][ip] = 0;
        });
    }

    // A particle may sit exactly on a cell face that belongs to another grid.
    // A global Redistribute moves every particle to its owning grid and rank.
    pc.Redistribute();
}

// Move particles in straight lines and record the distance traveled.
void pushParticles (PC& pc, amrex::Real dt, int dist_comp)
{
    const int lev = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (PC::ParIterType pti(pc, lev); pti.isValid(); ++pti)
    {
        const int np = pti.numParticles();
        auto ptd = pti.GetParticleTile().getParticleTileData();
        // Components can also be looked up as arrays, e.g., by index or by name.
        amrex::ParticleReal* dist = pti.GetStructOfArrays().GetRealData(dist_comp).data();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            amrex::ParticleReal v2 = 0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                const amrex::ParticleReal v = ptd.rdata(RealIdx::vx+d)[i];
                ptd.pos(d, i) += static_cast<amrex::ParticleReal>(v * dt);
                v2 += v*v;
            }
            dist[i] += static_cast<amrex::ParticleReal>(std::sqrt(v2) * dt);
        });
    }

    // Particles moved less than one cell, so only neighboring ranks need to
    // communicate (local = 1). Redistribute also applies periodic wrapping.
    pc.Redistribute(0, -1, 0, 1);
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int nppc = 2;
        int nsteps = 20;
        amrex::Real vmax = 1.0_rt;
        amrex::Real cfl = 0.5_rt;
        bool write_output = true;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nppc", nppc);
            pp.query("nsteps", nsteps);
            pp.query("vmax", vmax);
            pp.query("cfl", cfl);
            pp.query("write_output", write_output);
        }

        // Periodic unit domain.
        const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        const amrex::RealBox real_box({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                      {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
        const amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
        const amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        const amrex::DistributionMapping dm(ba);

        // The particle container uses the same grids as mesh data would.
        PC pc(geom, dm, ba);

        // Name the compile-time components and add a component at runtime.
        pc.SetSoACompileTimeNames({AMREX_D_DECL("x", "y", "z"),
                                   AMREX_D_DECL("vx", "vy", "vz"), "w"}, {});
        pc.AddRealComp("dist");
        const int dist_comp = pc.GetRealCompIndex("dist");

        amrex::Print() << "Real components:";
        for (auto const& name : pc.GetRealSoANames()) { amrex::Print() << " " << name; }
        amrex::Print() << "\n";

        initParticles(pc, nppc, vmax);
        const Stats s0 = computeStats(pc);
        amrex::Print() << "Initial: " << s0.np << " particles, total weight "
                       << s0.weight << ", kinetic energy " << s0.energy << "\n";

        // Each particle moves at most cfl cells per step in each direction.
        const amrex::Real dt = cfl * geom.CellSize(0) / vmax;
        for (int step = 1; step <= nsteps; ++step) {
            pushParticles(pc, dt, dist_comp);
        }

        const Stats s1 = computeStats(pc);
        amrex::Print() << "Final:   " << s1.np << " particles, total weight "
                       << s1.weight << ", kinetic energy " << s1.energy << "\n";

        constexpr bool is_single = std::is_same_v<amrex::Real, float>
                                || std::is_same_v<amrex::ParticleReal, float>;
        constexpr amrex::Real tol = is_single ? 1.e-4_rt : 1.e-12_rt;

        if (write_output) {
            // Mesh plotfile with the particle count per cell, plus the particles.
            amrex::MultiFab count(ba, dm, 1, 0);
            count.setVal(0.0);
            pc.Increment(count, 0);
            const std::string pltfile = amrex::Concatenate("plt", nsteps, 5);
            amrex::WriteSingleLevelPlotfile(pltfile, count, {"particle_count"},
                                            geom, nsteps*dt, nsteps);

            // Component names for I/O exclude the positions.
            amrex::Vector<std::string> real_names;
            auto const& names = pc.GetRealSoANames();
            for (int n = AMREX_SPACEDIM; n < static_cast<int>(names.size()); ++n) {
                real_names.push_back(names[n]);
            }
            pc.WritePlotFile(pltfile, "particles", real_names, {});

            // Write a checkpoint and read it back into a new container.
            const std::string chkfile = amrex::Concatenate("chk", nsteps, 5);
            pc.Checkpoint(chkfile, "particles", real_names, {});

            PC pc2(geom, dm, ba);
            pc2.AddRealComp("dist");
            pc2.Restart(chkfile, "particles");
            const Stats s2 = computeStats(pc2);
            amrex::Print() << "Restart: " << s2.np << " particles, total weight "
                           << s2.weight << "\n";
            if (s2.np != s1.np || std::abs(s2.weight - s1.weight) > tol*s1.weight
                || std::abs(s2.path - s1.path) > tol*s1.path) {
                amrex::Print() << "Restarted particles differ\n";
                status = 1;
            }
        }

        // No forces: count, weight and energy are conserved, and each particle
        // traveled nsteps*dt*|v|.
        const amrex::Real err_w = std::abs(s1.weight - s0.weight) / s0.weight;
        const amrex::Real err_e = std::abs(s1.energy - s0.energy) / s0.energy;
        const amrex::Real expected_path = nsteps * dt * s0.speed;
        const amrex::Real err_p = std::abs(s1.path - expected_path) / expected_path;
        amrex::Print() << "Relative errors: weight " << err_w << ", energy " << err_e
                       << ", path length " << err_p << "\n";

        if (s1.np != s0.np || err_w > tol || err_e > tol || err_p > tol) { status = 1; }
        if (status != 0) {
            amrex::Print() << "FAILED\n";
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
