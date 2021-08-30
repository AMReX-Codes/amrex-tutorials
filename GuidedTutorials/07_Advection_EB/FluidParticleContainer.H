#ifndef BL_FLUIDPARTICLES_H_
#define BL_FLUIDPARTICLES_H_

#include <AMReX_Particles.H>
#include <AMReX_ParticleMesh.H>
#include <AMReX_TracerParticle_mod_K.H>

#include "Indexing.H"

namespace Interpolation {
    enum {NGP=0, CIC};
}

namespace PIdx {
    enum {AMREX_D_DECL(BufferX=0, BufferY, BufferZ), Weight, NStructReal};
}

namespace amrex {

class FluidParticleContainer
    : public ParticleContainer<PIdx::NStructReal>
{
private:
    int m_number_particles_per_cell;

public:

    FluidParticleContainer (ParGDBBase* gdb)
        : ParticleContainer<PIdx::NStructReal>(gdb),
          m_number_particles_per_cell(0)
        {}

    FluidParticleContainer (const Geometry            & geom,
                         const DistributionMapping & dmap,
                         const BoxArray            & ba)
        : ParticleContainer<PIdx::NStructReal>(geom,dmap,ba),
          m_number_particles_per_cell(0)
        {}

    ~FluidParticleContainer () {}

    void InitParticles(const MultiFab& phi, const MultiFab& ebvol, Real density_cutoff, int nppc, int interpolation=Interpolation::CIC);

    int NumParticlesPerCell() { return m_number_particles_per_cell; }

    Real SumPhi();

    void AdvectWithUmac (MultiFab* umac, int lev, Real dt);

    void DepositToMesh (MultiFab& phi, int interpolation=Interpolation::CIC);

    void InterpolateFromMesh (const MultiFab& phi, int interpolation=Interpolation::CIC);

    void RemoveCoveredParticles (const MultiFab& ebvol, Real density_cutoff);
};

}

#endif

