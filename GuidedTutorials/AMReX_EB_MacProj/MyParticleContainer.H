#ifndef BL_MYPARTICLES_H_
#define BL_MYPARTICLES_H_

#include <AMReX_TracerParticles.H>

namespace amrex {

class MyParticleContainer
    : public TracerParticleContainer
{
public:

    MyParticleContainer (ParGDBBase* gdb)
        : TracerParticleContainer(gdb)
        {}

    MyParticleContainer (const Geometry            & geom,
                         const DistributionMapping & dmap,
                         const BoxArray            & ba)
        : TracerParticleContainer(geom,dmap,ba)
        {}

    ~MyParticleContainer () {}

    Real FindWinner (int n);

    void InitParticles (std::string initial_particle_file, Real zlen);

};

}

#endif

