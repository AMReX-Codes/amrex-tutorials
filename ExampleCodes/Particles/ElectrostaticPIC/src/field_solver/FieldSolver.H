#ifndef FIELD_SOLVER_H_
#define FIELD_SOLVER_H_

#include "Electrostatic_PIC_Util.H"

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

    void fixRHSForSolve (amrex::Vector<amrex::MultiFab*>& rhs,
                         const amrex::Vector<const amrex::iMultiFab*>& masks,
                         const amrex::Vector<amrex::Geometry>& geom,
                         const amrex::IntVect& ratio);

    void computePhi (VectorMeshData& rhs, ScalarMeshData& phi,
                     amrex::Vector<amrex::BoxArray>& grids,
                     amrex::Vector<amrex::DistributionMapping>& dm,
                     amrex::Vector<amrex::Geometry>& geom,
                     amrex::Vector<const amrex::iMultiFab*>& masks);

    void computeE (VectorMeshData& E, const ScalarMeshData& phi,
                   const amrex::Vector<amrex::Geometry>& geom);

    void sumFineToCrseNodal (const amrex::MultiFab& fine, amrex::MultiFab& crse,
                         const amrex::Geometry& cgeom, const amrex::IntVect& ratio);

    void zeroOutBoundary (amrex::MultiFab& input_data,
                          amrex::MultiFab& bndry_data,
                          const amrex::FabArray<amrex::BaseFab<int> >& mask);

    void getLevelMasks (amrex::Vector<amrex::iMultiFab*>& masks,
                        const amrex::Vector<amrex::BoxArray>& grids,
                        const amrex::Vector<amrex::DistributionMapping>& dmap,
                        const amrex::Vector<amrex::Geometry>& geom,
                        const int ncells = 1);

};

#endif
