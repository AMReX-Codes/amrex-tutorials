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
  Header file for 'simplest' ARKode integration of 2D
  Advection-Diffusion example problem.  This allows for either
  fixed-step or temporal adaptivity, but requires explicit time
  integration.
  --------------------------------------------------------------------*/


#ifndef HANDSON1_H
#define HANDSON1_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include "Utilities.h"

// user-data structure for problem options
struct ProblemOpt
{
   int help;
   int plot_int;
   int arkode_order;
   amrex::Real rtol;
   amrex::Real atol;
   amrex::Real fixed_dt;
   amrex::Real tfinal;
   amrex::Real dtout;
   int max_steps;
   int write_diag;
};

// Run problem
void DoProblem();

// Parse the problem input file
void ParseInputs(ProblemOpt& prob_opt, ProblemData& prob_data);

// Advance the solution in time with ARKode ARKStep
void ComputeSolutionARK(N_Vector nv_sol, ProblemOpt* prob_opt,
                        ProblemData* prob_data);

#endif
