#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleInterpolators.H>
#include <AMReX_ParticleMesh.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_AmrParticles.H>  // multi-level ParticleToMesh
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <cmath>
#include <type_traits>

using namespace amrex::literals;
using amrex::Real;

// Positions first, then the particle weight.
struct RealIdx
{
    enum { AMREX_D_DECL(x, y, z), w, ncomps };
};

using PC = amrex::ParticleContainerPureSoA<RealIdx::ncomps, 0>;
using PTD = PC::PTDType;
using ConstPTD = PC::ConstPTDType;

// Number of particles on level lev whose coarse cell lies in the region
// covered by level 1 (a Box in level-0 index space).
amrex::Long count_in_region (PC const& pc, int lev, amrex::Box const& region,
                             amrex::Geometry const& geom0)
{
    auto const plo = geom0.ProbLoArray();
    auto const dxi = geom0.InvCellSizeArray();
    amrex::Long n = amrex::ReduceSum(pc, lev,
        [=] AMREX_GPU_HOST_DEVICE (ConstPTD const& ptd, int i) noexcept -> amrex::Long
        {
            amrex::IntVect iv;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                iv[d] = int(std::floor((ptd.pos(d, i) - plo[d]) * dxi[d]));
            }
            return region.contains(iv) ? 1 : 0;
        });
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

// Every particle must be on the finest level that covers it: all level-1
// particles inside the fine region, no level-0 particle inside it.
bool check_levels (PC const& pc, amrex::Box const& region, amrex::Geometry const& geom0)
{
    amrex::Long const n0 = pc.NumberOfParticlesAtLevel(0);
    amrex::Long const n1 = pc.NumberOfParticlesAtLevel(1);
    amrex::Long const in0 = count_in_region(pc, 0, region, geom0);
    amrex::Long const in1 = count_in_region(pc, 1, region, geom0);
    amrex::Print() << "  level 0: " << n0 << " particles, " << in0 << " under level 1; "
                   << "level 1: " << n1 << " particles, " << in1 << " in its region\n";
    return in0 == 0 && in1 == n1;
}

// Deposit the particle weights on all levels with cloud-in-cell weights.
// ParticleToMesh handles the coarse/fine boundaries: contributions from
// fine particles to coarse cells and vice versa, and the average of the
// fine density replaces the coarse one under level 1.
void deposit (PC const& pc, amrex::Vector<amrex::MultiFab>& rho)
{
    amrex::ParticleToMesh(pc, amrex::GetVecOfPtrs(rho), 0, int(rho.size())-1,
        [=] AMREX_GPU_DEVICE (ConstPTD const& ptd, int i, amrex::Array4<Real> const& rho_arr,
                              amrex::GpuArray<Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<Real,AMREX_SPACEDIM> const& dxi) noexcept
        {
            auto const p = ptd[i];
            amrex::ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh(p, rho_arr, 0, 0, 1,
                [=] (PC::ParticleType::ConstType const&, int) noexcept
                {
                    return ptd.rdata(RealIdx::w)[i];
                });
        });
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        amrex::Long npart = 100000;
        int nsteps = 10;
        bool write_plotfile = true;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("npart", npart);
            pp.query("nsteps", nsteps);
            pp.query("write_plotfile", write_plotfile);
        }
        AMREX_ALWAYS_ASSERT(n_cell % 4 == 0 && max_grid_size % 2 == 0);

        // Two levels on the periodic unit domain: level 1 refines the middle
        // half of level 0 by 2.
        constexpr int nlevs = 2;
        amrex::IntVect const ratio(2);
        amrex::RealBox const rb({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
        amrex::Array<int,AMREX_SPACEDIM> const is_periodic{AMREX_D_DECL(1, 1, 1)};
        amrex::Box const domain0(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::Box region(amrex::IntVect(n_cell/4), amrex::IntVect(3*n_cell/4-1));

        amrex::Vector<amrex::Geometry> geom(nlevs);
        geom[0].define(domain0, rb, amrex::CoordSys::cartesian, is_periodic);
        geom[1].define(amrex::refine(domain0, ratio), rb, amrex::CoordSys::cartesian, is_periodic);
        amrex::Vector<amrex::BoxArray> grids(nlevs);
        grids[0].define(domain0);
        grids[1].define(amrex::refine(region, ratio));
        amrex::Vector<amrex::DistributionMapping> dmap(nlevs);
        for (int lev = 0; lev < nlevs; ++lev) {
            grids[lev].maxSize(max_grid_size);
            dmap[lev].define(grids[lev]);
        }

        // The container knows the whole hierarchy. A code built on AmrCore
        // constructs it from the AmrCore instead (AmrParticleContainer).
        PC pc(geom, dmap, grids, {ratio});
        pc.SetSoACompileTimeNames({AMREX_D_DECL("x", "y", "z"), "w"}, {});

        // Random particles of equal weight over the whole domain. Each one is
        // placed on the finest level that covers it, so the expected fraction
        // on level 1 is the volume fraction of the fine region.
        Real const weight = 1.0_rt / Real(npart);
        PC::ParticleInitData pdata = {{}, {}, {AMREX_D_DECL(0.0, 0.0, 0.0), weight}, {}};
        pc.InitRandom(npart, 451, pdata);
        amrex::Print() << "Initialized " << pc.TotalNumberOfParticles() << " particles on "
                       << pc.numLevels() << " levels\n";
        bool ok = check_levels(pc, region, geom[0]);
        Real const fine_fraction = Real(pc.NumberOfParticlesAtLevel(1)) / Real(npart);
        Real const expected_fraction = std::pow(0.5_rt, AMREX_SPACEDIM);
        amrex::Print() << "  fraction on level 1: " << fine_fraction
                       << " (expected about " << expected_fraction << ")\n";
        // Allow five standard deviations of the binomial count.
        Real const fraction_tol = 5.0_rt * std::sqrt(expected_fraction * (1.0_rt - expected_fraction)
                                                      / Real(npart));
        ok = ok && std::abs(fine_fraction - expected_fraction) < fraction_tol;

        // Density on both levels.
        amrex::Vector<amrex::MultiFab> rho(nlevs);
        for (int lev = 0; lev < nlevs; ++lev) {
            rho[lev].define(grids[lev], dmap[lev], 1, 0);
        }
        deposit(pc, rho);
        Real const vol0 = AMREX_D_TERM(geom[0].CellSize(0), *geom[0].CellSize(1), *geom[0].CellSize(2));
        Real const mass = rho[0].sum() * vol0;
        amrex::Print() << "Deposited mass on level 0 (fine data averaged down): " << mass << "\n";
        constexpr Real tol = std::is_same_v<Real,float> ? 1.e-5_rt : 1.e-12_rt;
        ok = ok && std::abs(mass - 1.0_rt) < tol;

        // Move the particles with a uniform velocity. After each step,
        // Redistribute moves them to their new grids and levels (and applies
        // the periodic boundaries).
        amrex::GpuArray<Real,AMREX_SPACEDIM> const vel{AMREX_D_DECL(0.7_rt, 0.4_rt, 0.2_rt)};
        Real const dt = 0.5_rt * geom[0].CellSize(0) / vel[0]; // half a coarse cell per step
        for (int step = 1; step <= nsteps; ++step) {
            for (int lev = 0; lev < nlevs; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
                for (PC::ParIterType pti(pc, lev); pti.isValid(); ++pti) {
                    PTD const ptd = pti.GetParticleTile().getParticleTileData();
                    int const np = pti.numParticles();
                    amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
                    {
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            ptd.pos(d, i) += vel[d] * dt;
                        }
                    });
                }
            }
            pc.Redistribute();
        }
        amrex::Print() << "After " << nsteps << " steps:\n";
        bool const levels_ok = check_levels(pc, region, geom[0]);
        ok = ok && levels_ok && (pc.TotalNumberOfParticles() == npart);

        // Regrid: level 1 moves to another region. The container gets the new
        // grids, and Redistribute moves the particles between the levels. A
        // container built from an AmrCore takes the new grids from it and
        // only needs the Redistribute.
        region.shift(0, n_cell/4);
        grids[1].define(amrex::refine(region, ratio));
        grids[1].maxSize(max_grid_size);
        dmap[1].define(grids[1]);
        pc.SetParticleBoxArray(1, grids[1]);
        pc.SetParticleDistributionMap(1, dmap[1]);
        pc.Redistribute();
        amrex::Print() << "After regridding level 1:\n";
        bool const regrid_ok = check_levels(pc, region, geom[0]);
        ok = ok && regrid_ok;

        rho[1].define(grids[1], dmap[1], 1, 0);
        deposit(pc, rho);
        Real const mass2 = rho[0].sum() * vol0;
        amrex::Print() << "Deposited mass on level 0: " << mass2 << "\n";
        ok = ok && std::abs(mass2 - 1.0_rt) < tol;

        if (write_plotfile) {
            // The density on both levels, and the particles of all levels.
            amrex::WriteMultiLevelPlotfile("plt_multilevel", nlevs, amrex::GetVecOfConstPtrs(rho),
                                           {"rho"}, geom, Real(nsteps)*dt, {nsteps, nsteps},
                                           {ratio});
            pc.WritePlotFile("plt_multilevel", "particles", {"w"}, {});
        }

        if (ok) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
