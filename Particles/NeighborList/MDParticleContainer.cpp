#include "MDParticleContainer.H"

#include <AMReX_ParticleReduce.H>
#include <AMReX_Random.H>

#include <cmath>

using namespace amrex::literals;

void
MDParticleContainer::initParticles (int nppc, amrex::Real v_thermal)
{
    const int lev = 0;
    const auto plo = Geom(lev).ProbLoArray();
    const auto dx = Geom(lev).CellSizeArray();
    const int myproc = amrex::ParallelDescriptor::MyProc();
    const int nppc_total = AMREX_D_TERM(nppc, *nppc, *nppc);
    // Index space of the particle lattice within one cell.
    const amrex::Box lattice(amrex::IntVect(0), amrex::IntVect(nppc-1));

    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box bx = mfi.tilebox();
        const amrex::Long np = bx.numPts() * nppc_total;
        const amrex::Long id0 = ParticleType::NextID();
        ParticleType::NextID(id0 + np);

        auto& tile = DefineAndReturnParticleTile(lev, mfi);
        tile.resize(np);
        auto ptd = tile.getParticleTileData();

        amrex::ParallelForRNG(np,
        [=] AMREX_GPU_DEVICE (amrex::Long ip, amrex::RandomEngine const& engine) noexcept
        {
            const amrex::IntVect iv = bx.atOffset(ip / nppc_total);
            const amrex::IntVect il = lattice.atOffset(ip % nppc_total);
            ptd.idcpu(ip) = amrex::SetParticleIDandCPU(id0 + ip, myproc);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                ptd.pos(d, ip) = static_cast<amrex::ParticleReal>(
                    plo[d] + (iv[d] + (il[d] + 0.5_rt) / nppc) * dx[d]);
                ptd.rdata(RealIdx::vx+d)[ip] = static_cast<amrex::ParticleReal>(
                    amrex::RandomNormal(0.0_rt, v_thermal, engine));
                ptd.rdata(RealIdx::ax+d)[ip] = 0;
            }
            ptd.rdata(RealIdx::pe)[ip] = 0;
        });
    }
}

void
MDParticleContainer::computeForces (MDParams const& params)
{
    const int lev = 0;
    const amrex::Real sigma = params.sigma;
    const amrex::Real k = params.stiffness;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        const auto index = std::make_pair(pti.index(), pti.LocalTileIndex());
        const int np = pti.numParticles(); // real particles, without neighbors
        auto ptd = pti.GetParticleTile().getParticleTileData();

        // For each real particle i, the list holds the indices of nearby
        // particles in this tile, which may be real particles or neighbor copies.
        const auto nbor_data = m_neighbor_list[lev].at(index).data();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            amrex::Real acc[AMREX_SPACEDIM] = {AMREX_D_DECL(0._rt, 0._rt, 0._rt)};
            amrex::Real pe = 0;
            for (const auto& p2 : nbor_data.getNeighbors(i))
            {
                amrex::Real dr[AMREX_SPACEDIM];
                amrex::Real r2 = 0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    dr[d] = ptd.pos(d, i) - p2.pos(d);
                    r2 += dr[d]*dr[d];
                }
                if (r2 >= sigma*sigma) { continue; }

                // Repulsive force k (sigma - r) along the unit vector dr/r.
                const amrex::Real r = amrex::max(std::sqrt(r2), 1.e-6_rt*sigma);
                const amrex::Real overlap = sigma - r;
                const amrex::Real coef = k * overlap / r;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    acc[d] += coef * dr[d];
                }
                // Each pair is visited from both sides; take half the energy.
                pe += 0.25_rt * k * overlap * overlap;
            }
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                ptd.rdata(RealIdx::ax+d)[i] = static_cast<amrex::ParticleReal>(acc[d]);
            }
            ptd.rdata(RealIdx::pe)[i] = static_cast<amrex::ParticleReal>(pe);
        });
    }
}

void
MDParticleContainer::kickDrift (amrex::Real dt_kick, amrex::Real dt_drift)
{
    const int lev = 0;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        const int np = pti.numParticles();
        auto ptd = pti.GetParticleTile().getParticleTileData();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                auto& v = ptd.rdata(RealIdx::vx+d)[i];
                v += static_cast<amrex::ParticleReal>(ptd.rdata(RealIdx::ax+d)[i] * dt_kick);
                ptd.pos(d, i) += static_cast<amrex::ParticleReal>(v * dt_drift);
            }
        });
    }
}

MDDiagnostics
MDParticleContainer::diagnostics () const
{
    using ConstPTD = ConstPTDType;
    // One pass over the particles per quantity for clarity; ParticleReduce
    // can do several in one pass.
    MDDiagnostics diag;
    diag.np = TotalNumberOfParticles();

    diag.kinetic = amrex::ReduceSum(*this,
        [=] AMREX_GPU_HOST_DEVICE (const ConstPTD& ptd, int i) noexcept -> amrex::Real
        {
            amrex::Real v2 = 0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                v2 += ptd.rdata(RealIdx::vx+d)[i] * ptd.rdata(RealIdx::vx+d)[i];
            }
            return 0.5_rt * v2;
        });

    diag.potential = amrex::ReduceSum(*this,
        [=] AMREX_GPU_HOST_DEVICE (const ConstPTD& ptd, int i) noexcept -> amrex::Real
        {
            return ptd.rdata(RealIdx::pe)[i];
        });

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        diag.momentum[d] = amrex::ReduceSum(*this,
            [=] AMREX_GPU_HOST_DEVICE (const ConstPTD& ptd, int i) noexcept -> amrex::Real
            {
                return ptd.rdata(RealIdx::vx+d)[i];
            });
    }

    diag.max_speed = amrex::ReduceMax(*this,
        [=] AMREX_GPU_HOST_DEVICE (const ConstPTD& ptd, int i) noexcept -> amrex::Real
        {
            amrex::Real v2 = 0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                v2 += ptd.rdata(RealIdx::vx+d)[i] * ptd.rdata(RealIdx::vx+d)[i];
            }
            return std::sqrt(v2);
        });

    // The reductions above are local to this process.
    amrex::ParallelDescriptor::ReduceRealSum(diag.kinetic);
    amrex::ParallelDescriptor::ReduceRealSum(diag.potential);
    amrex::ParallelDescriptor::ReduceRealSum(diag.momentum.data(), AMREX_SPACEDIM);
    amrex::ParallelDescriptor::ReduceRealMax(diag.max_speed);
    return diag;
}
