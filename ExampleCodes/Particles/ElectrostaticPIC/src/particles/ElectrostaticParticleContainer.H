#ifndef _ELECTROSTATICPARTICLECONTAINER_H_
#define _ELECTROSTATICPARTICLECONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>

#include "PhysConst.H"

struct PIdx
{
    enum { // Particle Attributes stored in amrex::ParticleContainer's struct-of-arrays
          w = 0,
          vx, vy,
          Ex, Ey,
          nattribs
    };
};

class MyParIter
    : public amrex::ParIter<0, 0, PIdx::nattribs>
{
public:
    using amrex::ParIter<0, 0, PIdx::nattribs>::ParIter;
    using RealVector = ParIter::ContainerType::RealVector;

    const std::array<RealVector, PIdx::nattribs>& GetAttribs () const {
        return GetStructOfArrays().GetRealData();
    }

    std::array<RealVector, PIdx::nattribs>& GetAttribs () {
        return GetStructOfArrays().GetRealData();
    }

    const RealVector& GetAttribs (int comp) const {
        return GetStructOfArrays().GetRealData(comp);
    }

    RealVector& GetAttribs (int comp) {
        return GetStructOfArrays().GetRealData(comp);
    }

};

class ElectrostaticParticleContainer
    : public amrex::ParticleContainer<0, 0, PIdx::nattribs>
{
public:

    ElectrostaticParticleContainer (const amrex::Vector<amrex::Geometry>            & geom,
                                    const amrex::Vector<amrex::DistributionMapping> & dmap,
                                    const amrex::Vector<amrex::BoxArray>            & ba,
                                    const amrex::Vector<int>                        & rr)
        : amrex::ParticleContainer<0, 0, PIdx::nattribs> (geom, dmap, ba, rr)
        {
            charge = -PhysConst::q_e;
            mass   =  PhysConst::m_e;
        }

    void InitParticles ();

    void pushX (const amrex::Real& dt);

    void DepositCharge (const amrex::Vector<amrex::MultiFab*>& rho);

    void FieldGather (const amrex::Vector<std::array<const amrex::MultiFab*, AMREX_SPACEDIM> >& E,
                      const amrex::Vector<const amrex::iMultiFab*>& masks);

    void Evolve (const amrex::Vector<std::array<const amrex::MultiFab*, AMREX_SPACEDIM> >& E,
                 const amrex::Vector<const amrex::MultiFab*>& rho, const amrex::Real& dt);

    void writeParticles (int n);

private:

    amrex::Real charge, mass;
};

#endif // _ELECTROSTATICPARTICLECONTAINER_H_
