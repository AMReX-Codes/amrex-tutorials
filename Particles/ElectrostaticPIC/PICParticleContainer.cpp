#include "PICParticleContainer.H"

#include <AMReX_ParticleInterpolators.H>
#include <AMReX_ParticleMesh.H>
#include <AMReX_ParticleReduce.H>

#include <cmath>

using namespace amrex::literals;

PICParticleContainer::PICParticleContainer (const amrex::Geometry& geom,
                                            const amrex::DistributionMapping& dm,
                                            const amrex::BoxArray& ba,
                                            amrex::Real charge, amrex::Real mass)
    : amrex::ParticleContainerPureSoA<PIdx::ncomps, 0>(geom, dm, ba),
      m_charge(charge), m_mass(mass)
{
    SetSoACompileTimeNames({AMREX_D_DECL("x", "y", "z"), AMREX_D_DECL("vx", "vy", "vz"),
                            AMREX_D_DECL("Ex", "Ey", "Ez"), "w"}, {});
}

void
PICParticleContainer::initLangmuirWave (const amrex::IntVect& nppc, amrex::Real density,
                                        amrex::Real k, amrex::Real amplitude)
{
    const int lev = 0;
    const auto plo = Geom(lev).ProbLoArray();
    const auto dx = Geom(lev).CellSizeArray();
    const int myproc = amrex::ParallelDescriptor::MyProc();
    const amrex::Box lattice(amrex::IntVect(0), nppc - 1);
    const auto nppc_total = static_cast<int>(lattice.numPts());
    const amrex::Real cell_volume = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    const amrex::Real weight = density * cell_volume / nppc_total;

    for (amrex::MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box bx = mfi.tilebox();
        const amrex::Long np = bx.numPts() * nppc_total;
        const amrex::Long id0 = ParticleType::NextID();
        ParticleType::NextID(id0 + np);

        auto& tile = DefineAndReturnParticleTile(lev, mfi);
        tile.resize(np);
        auto ptd = tile.getParticleTileData();

        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (amrex::Long ip) noexcept
        {
            const amrex::IntVect iv = bx.atOffset(ip / nppc_total);
            const amrex::IntVect il = lattice.atOffset(ip % nppc_total);
            ptd.idcpu(ip) = amrex::SetParticleIDandCPU(id0 + ip, myproc);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                amrex::Real pos = plo[d] + (iv[d] + (il[d] + 0.5_rt) / nppc[d]) * dx[d];
                if (d == 0) { pos += amplitude * std::sin(k * (pos - plo[0])); }
                ptd.pos(d, ip) = static_cast<amrex::ParticleReal>(pos);
                ptd.rdata(PIdx::vx+d)[ip] = 0;
                ptd.rdata(PIdx::Ex+d)[ip] = 0;
            }
            ptd.rdata(PIdx::w)[ip] = static_cast<amrex::ParticleReal>(weight);
        });
    }

    // Displaced particles may have left their grid.
    Redistribute();
}

void
PICParticleContainer::depositCharge (amrex::MultiFab& rho) const
{
    const int lev = 0;
    const auto plo = Geom(lev).ProbLoArray();
    const auto dxi = Geom(lev).InvCellSizeArray();
    const amrex::Real q_over_vol = m_charge * AMREX_D_TERM(dxi[0], *dxi[1], *dxi[2]);
    using ConstParticle = ParticleType::ConstType;

    // ParticleToMesh zeros rho, then adds ghost cell contributions to the
    // valid cells they overlap, including across periodic boundaries.
    amrex::ParticleToMesh(*this, rho, lev,
        [=] AMREX_GPU_DEVICE (const ConstPTDType& ptd, int i,
                              amrex::Array4<amrex::Real> const& rho_arr) noexcept
        {
            const ConstParticle p = ptd[i];
            amrex::ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh(p, rho_arr, 0, 0, 1,
                [=] (ConstParticle const&, int) noexcept
                {
                    return q_over_vol * ptd.rdata(PIdx::w)[i];
                });
        });
}

void
PICParticleContainer::gatherField (const amrex::MultiFab& efield)
{
    const int lev = 0;
    const auto plo = Geom(lev).ProbLoArray();
    const auto dxi = Geom(lev).InvCellSizeArray();

    amrex::MeshToParticle(*this, efield, lev,
        [=] AMREX_GPU_DEVICE (const PTDType& ptd, int i,
                              amrex::Array4<const amrex::Real> const& e_arr) noexcept
        {
            auto p = ptd[i];
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                ptd.rdata(PIdx::Ex+d)[i] = 0;
            }
            amrex::ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.MeshToParticle(p, e_arr, 0, 0, AMREX_SPACEDIM,
                [=] (amrex::Array4<const amrex::Real> const& arr,
                     int ii, int jj, int kk, int comp) noexcept
                {
                    return arr(ii, jj, kk, comp);
                },
                [=] (ParticleType&, int comp, amrex::Real val) noexcept
                {
                    ptd.rdata(PIdx::Ex+comp)[i] += static_cast<amrex::ParticleReal>(val);
                });
        });
}

void
PICParticleContainer::pushVelocity (amrex::Real dt)
{
    const int lev = 0;
    const amrex::Real coef = m_charge / m_mass * dt;
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
                ptd.rdata(PIdx::vx+d)[i] += static_cast<amrex::ParticleReal>(
                    coef * ptd.rdata(PIdx::Ex+d)[i]);
            }
        });
    }
}

void
PICParticleContainer::pushPosition (amrex::Real dt)
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
                ptd.pos(d, i) += static_cast<amrex::ParticleReal>(ptd.rdata(PIdx::vx+d)[i] * dt);
            }
        });
    }

    // The time step keeps particles within one cell per step (see main.cpp).
    Redistribute(0, -1, 0, 1);
}

amrex::Real
PICParticleContainer::kineticEnergy () const
{
    amrex::Real ke = amrex::ReduceSum(*this,
        [=] AMREX_GPU_HOST_DEVICE (const ConstPTDType& ptd, int i) noexcept -> amrex::Real
        {
            amrex::Real v2 = 0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                v2 += ptd.rdata(PIdx::vx+d)[i] * ptd.rdata(PIdx::vx+d)[i];
            }
            return ptd.rdata(PIdx::w)[i] * v2;
        });
    amrex::ParallelDescriptor::ReduceRealSum(ke);
    return 0.5_rt * m_mass * ke;
}

amrex::Real
PICParticleContainer::totalCharge () const
{
    amrex::Real sum_w = amrex::ReduceSum(*this,
        [=] AMREX_GPU_HOST_DEVICE (const ConstPTDType& ptd, int i) noexcept -> amrex::Real
        {
            return ptd.rdata(PIdx::w)[i];
        });
    amrex::ParallelDescriptor::ReduceRealSum(sum_w);
    return m_charge * sum_w;
}
