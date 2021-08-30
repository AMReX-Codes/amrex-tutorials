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
  Header file for 'advanced' ARKode integration of 2D
  Advection-Diffusion example problem.  This treats the diffusion
  portion of the problem implicitly, but allows advection to be either
  implicit (via DIRK methods) or explicit (via ARK-IMEX methods).
  The implicit portion may be solved using either Newton's method or
  an Anderson-accelerated fixed-point nonlinear solver.  When using
  Newton's method, the underlying linear systems are solved using a
  GMRES iterative linear solver, with optional preconditioning.
  This preconditioner treats only the diffusion portion (exact for
  ARK-IMEX time integration, or only approximate for DIRK), using the
  scalable MLMG linear solver from AMReX.  As with the other hands-on
  lessons, this program allows for either fixed time-step sizes or
  temporal adaptivity.
  --------------------------------------------------------------------*/


#ifndef HANDSON3_H
#define HANDSON3_H

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
   int nls_method;
   int nls_max_iter;
   int nls_fp_acc;
   int ls_max_iter;
   int rhs_adv;
   amrex::Real rtol;
   amrex::Real atol;
   amrex::Real fixed_dt;
   amrex::Real tfinal;
   amrex::Real dtout;
   int max_steps;
   int write_diag;
   int use_preconditioner;
};

// Run problem
void DoProblem();

// Parse the problem input file
void ParseInputs(ProblemOpt& prob_opt, ProblemData& prob_data);

// Advance the solution in time with ARKode ARKStep
void ComputeSolutionARK(N_Vector nv_sol, ProblemOpt* prob_opt,
                        ProblemData* prob_data);

#endif
