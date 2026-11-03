#ifndef PIC_PARTICLE_CONTAINER_H_
#define PIC_PARTICLE_CONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_Particles.H>

// Particle components: position, velocity, electric field at the particle,
// and weight (number of physical particles per macro-particle).
struct PIdx
{
    enum {
        AMREX_D_DECL(x, y, z),
        AMREX_D_DECL(vx, vy, vz),
        AMREX_D_DECL(Ex, Ey, Ez),
        w,
        ncomps
    };
};

// One particle species with charge q and mass m.
class PICParticleContainer
    : public amrex::ParticleContainerPureSoA<PIdx::ncomps, 0>
{
public:
    PICParticleContainer (const amrex::Geometry& geom,
                          const amrex::DistributionMapping& dm,
                          const amrex::BoxArray& ba,
                          amrex::Real charge, amrex::Real mass);

    // A cold plasma of the given density, with nppc particles per cell on a
    // lattice, displaced by amplitude*sin(k x) in the x-direction.
    void initLangmuirWave (const amrex::IntVect& nppc, amrex::Real density,
                           amrex::Real k, amrex::Real amplitude);

    // Charge density with cloud-in-cell weights. rho needs one ghost cell.
    void depositCharge (amrex::MultiFab& rho) const;

    // Interpolate the cell-centered E (filled ghost cells) to the particles.
    void gatherField (const amrex::MultiFab& efield);

    // v += (q/m) E dt
    void pushVelocity (amrex::Real dt);

    // x += v dt, followed by Redistribute.
    void pushPosition (amrex::Real dt);

    [[nodiscard]] amrex::Real kineticEnergy () const;
    [[nodiscard]] amrex::Real totalCharge () const;

    [[nodiscard]] amrex::Real charge () const { return m_charge; }
    [[nodiscard]] amrex::Real mass () const { return m_mass; }

private:
    amrex::Real m_charge;
    amrex::Real m_mass;
};

#endif
