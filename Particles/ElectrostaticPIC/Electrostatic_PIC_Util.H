#ifndef ELECTROSTATIC_PIC_UTIL_H_
#define ELECTROSTATIC_PIC_UTIL_H_

#include <AMReX_BaseFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

void sumFineToCrseNodal (const amrex::MultiFab& fine, amrex::MultiFab& crse,
                         const amrex::Geometry& cgeom, const amrex::IntVect& ratio);

void zeroOutBoundary (amrex::MultiFab& input_data,
                      amrex::MultiFab& bndry_data,
                      const amrex::FabArray<amrex::BaseFab<int> >& mask);

#endif
