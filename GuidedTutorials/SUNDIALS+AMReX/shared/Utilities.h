/*--------------------------------------------------------------------
  Time Integration and Nonlinear Solvers
  Hands-on Lessons with SUNDIALS + AMReX
  2019 Argonne Training Program in Extreme-Scale Computing

  Authors (alphabetical):
    David Gardner (gardner48@llnl.gov)
    John Loffeld (loffeld1@llnl.gov)
    Daniel Reynolds (reynolds@smu.edu)
    Donald Willcox (dewillcox@lbl.gov)

  --------------------------------------------------------------------
  Header file for general utility routines shared between hands-on
  lessons.
  --------------------------------------------------------------------*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Array.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>

// user-data structure passed through SUNDIALS to RHS functions
struct ProblemData
{
   // Requested grid options
   int n_cell;
   int max_grid_size;

   // AMReX grid data structures
   amrex::Geometry* geom;
   amrex::BoxArray* grid;
   amrex::DistributionMapping* dmap;

   // AMReX MLMG data and parameters
   amrex::MultiFab* acoef;
   amrex::MultiFab* bcoef;

   int mg_agglomeration;
   int mg_consolidation;
   int mg_max_coarsening_level;
   int mg_linop_maxorder;
   int mg_max_iter;
   int mg_max_fmg_iter;
   int mg_verbose;
   int mg_bottom_verbose;
   int mg_use_hypre;
   int mg_hypre_interface;
   int mg_use_petsc;
   amrex::Real mg_tol_rel;

   // Problem data
   amrex::Real advCoeffx;
   amrex::Real advCoeffy;
   amrex::Real diffCoeffx;
   amrex::Real diffCoeffy;
   amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* flux;
};

// Set the ODE initial condition
void FillInitConds2D(amrex::MultiFab& sol, const amrex::Geometry& geom);

// Decompose the problem in space
void SetUpGeometry(amrex::BoxArray& ba, amrex::Geometry& geom, ProblemData& prob_data);

// SUNDIALS ODE RHS functions
int ComputeRhsAdv(amrex::Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data);
int ComputeRhsDiff(amrex::Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data);
int ComputeRhsAdvDiff(amrex::Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data);

// Advective portion of ODE RHS
void ComputeAdvectionUpwind(amrex::MultiFab& sol,
                            amrex::MultiFab& advection,
                            amrex::Geometry& geom,
                            int comp,
                            amrex::Real advCoeffx,
                            amrex::Real advCoeffy);

// Diffusive portion of ODE RHS
void ComputeDiffusion(amrex::MultiFab& sol,
                      amrex::MultiFab& diff_mf,
                      amrex::MultiFab& fx_mf,
                      amrex::MultiFab& fy_mf,
                      amrex::Geometry& geom,
                      int comp,
                      amrex::Real diffCoeffx,
                      amrex::Real diffCoeffy);

// Utility functions for computing diffusion
void ComputeDiffFlux(amrex::MultiFab& sol,
                     amrex::MultiFab& fx,
                     amrex::MultiFab& fy,
                     amrex::Geometry& geom,
                     int comp,
                     amrex::Real diffCoeffx,
                     amrex::Real diffCoeffy);

void ComputeDivergence(amrex::MultiFab& div,
                       amrex::MultiFab& fx,
                       amrex::MultiFab& fy,
                       amrex::Geometry& geom,
                       int comp);

// Preconditioner routines
int precondition_setup(realtype tn, N_Vector u, N_Vector fu,
                       booleantype jok, booleantype *jcurPtr,
                       realtype gamma, void *user_data);

int precondition_solve(realtype tn, N_Vector u, N_Vector fu,
                       N_Vector r, N_Vector z,
                       realtype gamma, realtype delta,
                       int lr, void *user_data);


#endif
