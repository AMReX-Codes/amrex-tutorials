#ifndef ELECTROSTATIC_PIC_UTIL_H_
#define ELECTROSTATIC_PIC_UTIL_H_

#include <AMReX_BaseFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

void sumFineToCrseNodal (const amrex::MultiFab& fine, amrex::MultiFab& crse,
                         const amrex::Geometry& cgeom, const amrex::IntVect& ratio);

void zeroOutBoundary (amrex::MultiFab& input_data,
                      amrex::MultiFab& bndry_data,
                      const amrex::FabArray<amrex::BaseFab<int> >& mask);

void getLevelMasks (amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<int> > > >& masks,
                    const amrex::Vector<amrex::BoxArray>& grids,
                    const amrex::Vector<amrex::DistributionMapping>& dmap,
                    const amrex::Vector<amrex::Geometry>& geom,
                    const int ncells = 1);

#endif
