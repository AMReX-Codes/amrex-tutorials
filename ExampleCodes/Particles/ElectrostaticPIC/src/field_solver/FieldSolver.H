#ifndef FIELD_SOLVER_H_
#define FIELD_SOLVER_H_

#include "FieldSolver_K.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_IntVect.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <memory>

namespace FieldSolver {

    void fixRHSForSolve (const amrex::Vector<amrex::MultiFab*>& rhs,
                         const amrex::Vector<const amrex::iMultiFab*>& masks,
                         const amrex::Vector<amrex::Geometry>& geom,
                         const amrex::IntVect& ratio);

    void computePhi (const amrex::Vector<const amrex::MultiFab*>& rhs,
                     const amrex::Vector<amrex::MultiFab*>& phi,
                     amrex::Vector<amrex::BoxArray>& grids,
                     amrex::Vector<amrex::DistributionMapping>& dm,
                     amrex::Vector<amrex::Geometry>& geom,
                     const amrex::Vector<const amrex::iMultiFab*>& masks);

    void computeE (const amrex::Vector<std::array<amrex::MultiFab*, AMREX_SPACEDIM> >& E,
                   const amrex::Vector<const amrex::MultiFab*>& phi,
                   const amrex::Vector<amrex::Geometry>& geom);

    void sumFineToCrseNodal (const amrex::MultiFab& fine, amrex::MultiFab& crse,
                             const amrex::Geometry& cgeom, const amrex::IntVect& ratio);

    void zeroOutBoundary (amrex::MultiFab& input_data,
                          amrex::MultiFab& bndry_data,
                          const amrex::iMultiFab& mask);

    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > getLevelMasks
         (const amrex::Vector<amrex::BoxArray>& grids,
          const amrex::Vector<amrex::DistributionMapping>& dmap,
          const amrex::Vector<amrex::Geometry>& geom,
          const int ncells = 1);

};

#endif
