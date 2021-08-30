#ifndef BL_MYPARTICLES_H_
#define BL_MYPARTICLES_H_

#include <AMReX_TracerParticles.H>

namespace amrex {

struct PIdx
{
    enum {
        vx = 0, vy,
#if (AMREX_SPACEDIM == 3)
        vz,
#endif
        ncomps
    };
};

class MyParticleContainer
    : public amrex::ParticleContainer<PIdx::ncomps, 1>
{

public:

    MyParticleContainer (const amrex::Geometry            & a_geom,
                         const amrex::DistributionMapping & a_dmap,
                         const amrex::BoxArray            & a_ba)
        : ParticleContainer<PIdx::ncomps, 1>(a_geom, a_dmap, a_ba)
    {}

    ~MyParticleContainer () {};

    void InitPachinko (std::string initial_particle_file, Real zlen);

    void AdvectPachinko (Real dt, amrex::Vector<amrex::RealArray>& obstacle_center,
                         Real obstacle_radius, Real particle_radius);
};

}

#endif
