#ifndef ELECTROSTATIC_PIC_UTIL_H_
#define ELECTROSTATIC_PIC_UTIL_H_

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

void sumFineToCrseNodal (const amrex::MultiFab& fine, amrex::MultiFab& crse,
                         const amrex::Geometry& cgeom, const amrex::IntVect& ratio);

#endif
