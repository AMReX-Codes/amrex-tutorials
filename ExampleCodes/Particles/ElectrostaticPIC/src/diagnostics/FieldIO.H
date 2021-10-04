#ifndef FIELD_IO_H_
#define FIELD_IO_H_

#include "particles/ElectrostaticParticleContainer.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <array>

void WritePlotFile (const amrex::Vector<const amrex::MultiFab*>& rhs,
                    const amrex::Vector<const amrex::MultiFab*>& phi,
                    const amrex::Vector<std::array<const amrex::MultiFab*, AMREX_SPACEDIM> >& E,
                    const ElectrostaticParticleContainer& pc,
                    const amrex::Vector<amrex::Geometry>& geom,
                    int nstep);

#endif
