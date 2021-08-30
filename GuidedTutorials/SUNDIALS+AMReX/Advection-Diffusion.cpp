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
  Implementation file for 'general' SUNDIALS interface (CVODE and
  ARKode, with many configuration options) for 2D Advection-Diffusion
  example problem.
  --------------------------------------------------------------------*/

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cvode/cvode.h>
#include <arkode/arkode_arkstep.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "Advection-Diffusion.h"
#include "Utilities.h"

#include "NVector_Multifab.h"

using namespace amrex;


void ComputeSolutionCV(N_Vector nv_sol, ProblemOpt* prob_opt,
                       ProblemData* prob_data)
{
   // Extract problem data and options
   Geometry* geom         = prob_data->geom;
   int       plot_int     = prob_opt->plot_int;
   int       cvode_method = prob_opt->cvode_method;
   int       nls_method   = prob_opt->nls_method;
   int       nls_max_iter = prob_opt->nls_max_iter;
   int       nls_fp_acc   = prob_opt->nls_fp_acc;
   int       ls_max_iter  = prob_opt->ls_max_iter;
   int       rhs_adv      = prob_opt->rhs_adv;
   int       rhs_diff     = prob_opt->rhs_diff;
   Real      rtol         = prob_opt->rtol;
   Real      atol         = prob_opt->atol;
   Real      tfinal       = prob_opt->tfinal;
   Real      dtout        = prob_opt->dtout;
   int       max_steps    = prob_opt->max_steps;
   int use_preconditioner = prob_opt->use_preconditioner;

   // initial time, number of outputs, and error flag
   Real time = 0.0;
   int  nout = ceil(tfinal/dtout);
   int  ier  = 0;

   // Write a plotfile of the initial data
   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, time, 0);
   }

   // Create CVODE memory
   void* cvode_mem = NULL;
   if (cvode_method == 0)
      cvode_mem = CVodeCreate(CV_BDF);
   else
      cvode_mem = CVodeCreate(CV_ADAMS);

   // Initialize CVODE
   if (rhs_adv > 0 && rhs_diff > 0)
   {
      // implicit Advection and Diffusion
      CVodeInit(cvode_mem, ComputeRhsAdvDiff, time, nv_sol);
   }
   else if (rhs_adv > 0)
   {
      // implicit Advection
      CVodeInit(cvode_mem, ComputeRhsAdv, time, nv_sol);
   }
   else if (rhs_diff > 0)
   {
      // implicit Diffusion
      CVodeInit(cvode_mem, ComputeRhsDiff, time, nv_sol);
   }
   else
   {
      amrex::Print() << "Invalid RHS options for CVODE" << std::endl;
      return;
   }

   // Attach the user data structure to CVODE
   CVodeSetUserData(cvode_mem, prob_data);

   // Set integration tolerances
   CVodeSStolerances(cvode_mem, atol, rtol);

   // Set the max number of steps between outputs
   CVodeSetMaxNumSteps(cvode_mem, max_steps);

   // Attach nonlinear/linear solvers as needed
   if (nls_method == 0)
   {
      // Create and attach GMRES linear solver for Newton
      SUNLinearSolver LS;
      if (use_preconditioner)
         LS = SUNLinSol_SPGMR(nv_sol, PREC_LEFT, ls_max_iter);
      else
         LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, ls_max_iter);

      ier = CVodeSetLinearSolver(cvode_mem, LS, NULL);
      if (ier != CVLS_SUCCESS)
      {
         amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
         return;
      }

      if (use_preconditioner) {
         // Attach preconditioner setup/solve functions
         ier = CVodeSetPreconditioner(cvode_mem, precondition_setup, precondition_solve);
         if (ier != CVLS_SUCCESS)
         {
            amrex::Print() << "Attachment of preconditioner unsuccessful" << std::endl;
            return;
         }
      }
   }
   else
   {
      // Create and attach fixed point solver
      SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(nv_sol, nls_fp_acc);
      ier = CVodeSetNonlinearSolver(cvode_mem, NLS);
      if (ier != CV_SUCCESS)
      {
         amrex::Print() << "Creation of nonlinear solver unsuccessful" << std::endl;
         return;
      }
   }

   // Set max number of nonlinear iterations
   ier = CVodeSetMaxNonlinIters(cvode_mem, nls_max_iter);
   if (ier != CV_SUCCESS)
   {
      amrex::Print() << "Error setting max number of nonlinear iterations" << std::endl;
      return;
   }

   // Advance the solution in time
   Real tout = time + dtout; // first output time
   Real tret;                // return time
   for (int iout=0; iout < nout; iout++)
   {
      ier = CVode(cvode_mem, tout, nv_sol, &tret, CV_NORMAL);
      if (ier < 0)
      {
         amrex::Print() << "Error in CVODE" << std::endl;
         return;
      }

      // Get integration stats
      long nf_evals;
      CVodeGetNumRhsEvals(cvode_mem, &nf_evals);
      amrex::Print() << "t = " << std::setw(5) << tret
                     << "  rhs evals = " << std::setw(7) << nf_evals
                     << std::endl;

      // Write output
      if (plot_int > 0)
      {
         const std::string& pltfile = amrex::Concatenate("plt", iout+1, 5);
         MultiFab* sol = NV_MFAB(nv_sol);
         WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, tret, iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }
}


void ComputeSolutionARK(N_Vector nv_sol, ProblemOpt* prob_opt,
                        ProblemData* prob_data)
{
   // Extract problem data and options
   Geometry* geom         = prob_data->geom;
   int       plot_int     = prob_opt->plot_int;
   int       arkode_order = prob_opt->arkode_order;
   int       nls_method   = prob_opt->nls_method;
   int       nls_max_iter = prob_opt->nls_max_iter;
   int       nls_fp_acc   = prob_opt->nls_fp_acc;
   int       ls_max_iter  = prob_opt->ls_max_iter;
   int       rhs_adv      = prob_opt->rhs_adv;
   int       rhs_diff     = prob_opt->rhs_diff;
   Real      rtol         = prob_opt->rtol;
   Real      atol         = prob_opt->atol;
   Real      fixed_dt     = prob_opt->fixed_dt;
   Real      tfinal       = prob_opt->tfinal;
   Real      dtout        = prob_opt->dtout;
   int       max_steps    = prob_opt->max_steps;
   int       write_diag   = prob_opt->write_diag;
   int use_preconditioner = prob_opt->use_preconditioner;

   // initial time, number of outputs, and error flag
   Real time = 0.0;
   int  nout = ceil(tfinal/dtout);
   int  ier  = 0;

   // Write a plotfile of the initial data
   if (plot_int > 0)
   {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, time, 0);
   }

   // Create the ARK stepper
   void* arkode_mem = NULL;

   if (rhs_adv > 0 && rhs_diff > 0)
   {
      if (rhs_adv > 1 && rhs_diff > 1)
      {
         // explicit advection and diffusion
         arkode_mem = ARKStepCreate(ComputeRhsAdvDiff, NULL,
                                    time, nv_sol);
      }
      else if (rhs_adv > 1)
      {
         // explicit advection and implicit diffusion
         arkode_mem = ARKStepCreate(ComputeRhsAdv, ComputeRhsDiff,
                                    time, nv_sol);
      }
      else if (rhs_diff > 1)
      {
         // implicit advection and explicit diffusion
         arkode_mem = ARKStepCreate(ComputeRhsDiff, ComputeRhsAdv,
                                    time, nv_sol);
      }
      else
      {
         // implicit advection and diffusion
         arkode_mem = ARKStepCreate(NULL, ComputeRhsAdvDiff,
                                    time, nv_sol);
      }
   }
   else if (rhs_adv > 0)
   {
      if (rhs_adv > 1)
      {
         // explicit advection
         arkode_mem = ARKStepCreate(ComputeRhsAdv, NULL,
                                    time, nv_sol);
      }
      else
      {
         // implicit advection
         arkode_mem = ARKStepCreate(NULL, ComputeRhsAdv,
                                    time, nv_sol);
      }
   }
   else if (rhs_diff > 0)
   {
      if (rhs_diff > 1)
      {
         // explicit diffusion
         arkode_mem = ARKStepCreate(ComputeRhsDiff, NULL,
                                    time, nv_sol);
      }
      else
      {
         // implicit diffusion
         arkode_mem = ARKStepCreate(NULL, ComputeRhsDiff,
                                    time, nv_sol);
      }
   }
   else
   {
      amrex::Print() << "Invalid RHS options for ARKode" << std::endl;
      return;
   }

   // Attach the user data structure to ARKStep
   ARKStepSetUserData(arkode_mem, prob_data);

   // Set the method order
   ARKStepSetOrder(arkode_mem, arkode_order);

   // Set the time step size or integration tolerances
   if (fixed_dt > 0.0)
      ARKStepSetFixedStep(arkode_mem, fixed_dt);
   else
      ARKStepSStolerances(arkode_mem, atol, rtol);

   // Set the max number of steps between outputs
   ARKStepSetMaxNumSteps(arkode_mem, max_steps);

   // Set file for writing ARKStep diagnostics
   FILE* diagfp = NULL;
   if (write_diag) {
      diagfp = fopen("ARKStep_diagnostics.txt", "w");
      ARKStepSetDiagnostics(arkode_mem, diagfp);
   }

   // Attach nonlinear/linear solvers as needed
   if (rhs_adv == 1 || rhs_diff == 1)
   {
      if (nls_method == 0)
      {
         // Create and attach GMRES linear solver for Newton
         SUNLinearSolver LS;
         if (use_preconditioner)
            LS = SUNLinSol_SPGMR(nv_sol, PREC_LEFT, ls_max_iter);
         else
            LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, ls_max_iter);

         ier = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
         if (ier != ARKLS_SUCCESS)
         {
            amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
            return;
         }

         if (use_preconditioner) {
            // Attach preconditioner setup/solve functions
            ier = ARKStepSetPreconditioner(arkode_mem, precondition_setup, precondition_solve);
            if (ier != ARKLS_SUCCESS)
            {
               amrex::Print() << "Attachment of preconditioner unsuccessful" << std::endl;
               return;
            }
         }
      }
      else
      {
         // Create and attach GMRES linear solver (if implicit and using Newton)
         SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(nv_sol, nls_fp_acc);
         ier = ARKStepSetNonlinearSolver(arkode_mem, NLS);
         if (ier != ARK_SUCCESS)
         {
            amrex::Print() << "Creation of nonlinear solver unsuccessful" << std::endl;
            return;
         }
      }

      // Set max number of nonlinear iterations
      ier = ARKStepSetMaxNonlinIters(arkode_mem, nls_max_iter);
      if (ier != ARK_SUCCESS)
      {
         amrex::Print() << "Error setting max number of nonlinear iterations" << std::endl;
         return;
      }
   }

   // Advance the solution in time
   Real tout = time + dtout; // first output time
   Real tret;                // return time
   for (int iout=0; iout < nout; iout++)
   {
      ier = ARKStepEvolve(arkode_mem, tout, nv_sol, &tret, ARK_NORMAL);
      if (ier < 0)
      {
         amrex::Print() << "Error in ARKStepEvolve" << std::endl;
         return;
      }

      // Get integration stats
      long nfe_evals, nfi_evals;
      ARKStepGetNumRhsEvals(arkode_mem, &nfe_evals, &nfi_evals);
      amrex::Print() << "t = " << std::setw(5) << tret
                     << "  explicit evals = " << std::setw(7) << nfe_evals
                     << "  implicit evals = " << std::setw(7) << nfi_evals
                     << std::endl;

      // Write output
      if (plot_int > 0)
      {
         const std::string& pltfile = amrex::Concatenate("plt", iout+1, 5);
         MultiFab* sol = NV_MFAB(nv_sol);
         WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, tret, iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }

   // Output final solution statistics
   long int nst, nst_a, nfe, nfi, nsetups, nli, nJv, nlcf, nni, ncfn, netf, npe, nps;
   nst = nst_a = nfe = nfi = nsetups = nli = nJv = nlcf = nni = ncfn = netf = npe = nps = 0;
   ARKStepGetNumSteps(arkode_mem, &nst);
   ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
   ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
   ARKStepGetNumErrTestFails(arkode_mem, &netf);
   ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
   ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
   if (nls_method == 0) {
      ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
      ARKStepGetNumLinIters(arkode_mem, &nli);
      ARKStepGetNumJtimesEvals(arkode_mem, &nJv);
      ARKStepGetNumLinConvFails(arkode_mem, &nlcf);
      ARKStepGetNumPrecEvals(arkode_mem, &npe);
      ARKStepGetNumPrecSolves(arkode_mem, &nps);
   }
   amrex::Print() << "\nFinal Solver Statistics:\n"
                  << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n"
                  << "   Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << "\n"
                  << "   Total number of nonlinear iterations = " << nni << "\n"
                  << "   Total number of nonlinear solver convergence failures = " << ncfn << "\n"
                  << "   Total number of error test failures = " << netf << "\n";
   if (nls_method == 0) {
     amrex::Print() << "   Total linear solver setups = " << nsetups << "\n"
                    << "   Total linear iterations = " << nli << "\n"
                    << "   Total number of Jacobian-vector products = " << nJv << "\n"
                    << "   Total number of linear solver convergence failures = " << nlcf << "\n";
     if ((npe != 0) || (nps != 0))
       amrex::Print() << "   Total number of Preconditioner setups = " << npe << "\n"
                      << "   Total number of Preconditioner solves = " << nps << "\n";
   }


   // Close diagnostics file
   if (write_diag)
      fclose(diagfp);
}


void ParseInputs(ProblemOpt& prob_opt, ProblemData& prob_data)
{
   // ParmParse is way of reading inputs from the inputs file
   ParmParse pp;

   // --------------------------------------------------------------------------
   // Problem options
   // --------------------------------------------------------------------------

   // Enable (>0) or disable (<0) writing output files
   int plot_int = -1; // plots off
   pp.query("plot_int", plot_int);
   prob_opt.plot_int = plot_int;

   // Specify which integration method to use
   // 0 = CVODE
   // 1 = ARKStep
   int stepper = 1;
   pp.query("stepper", stepper);
   prob_opt.stepper = stepper;

   // Specify which CVODE method to use
   int cvode_method = 0; // BDF
   pp.query("cvode_method", cvode_method);
   prob_opt.cvode_method = cvode_method;

   // Specify the ARKode method order
   int arkode_order = 4; // 4th order
   pp.query("arkode_order", arkode_order);
   prob_opt.arkode_order = arkode_order;

   // Specify the nonlinear solver
   int nls_method = 0; // Newton
   pp.query("nls_method", nls_method);
   prob_opt.nls_method = nls_method;

   // Specify the max number of nonlinear iterations
   int nls_max_iter = 3;
   pp.query("nls_max_iter", nls_max_iter);
   prob_opt.nls_max_iter = nls_max_iter;

   // Specify the number of fixed point acceleration vectors
   int nls_fp_acc = 3; // no acceleration
   pp.query("nls_fp_acc", nls_fp_acc);
   prob_opt.nls_fp_acc = nls_fp_acc;

   // Specify the max number of linear iterations
   int ls_max_iter = 5;
   pp.query("ls_max_iter", ls_max_iter);
   prob_opt.ls_max_iter = ls_max_iter;

   // Specify RHS functions/splitting
   int rhs_adv  = 1; // implicit advection
   int rhs_diff = 1; // implicit diffusion
   pp.query("rhs_adv", rhs_adv);
   pp.query("rhs_diff", rhs_diff);
   prob_opt.rhs_adv  = rhs_adv;
   prob_opt.rhs_diff = rhs_diff;

   // Specify relative and absolute tolerances
   Real rtol = 1.0e-4;
   Real atol = 1.0e-9;
   pp.query("rtol", rtol);
   pp.query("atol", atol);
   prob_opt.rtol = rtol;
   prob_opt.atol = atol;

   // Specify a fixed time step size
   Real fixed_dt = -1.0; // diabled by default (use adaptive steps)
   pp.query("fixed_dt", fixed_dt);
   prob_opt.fixed_dt = fixed_dt;

   // Specify final time for integration
   Real tfinal = 1.0e4;
   pp.query("tfinal", tfinal);
   prob_opt.tfinal = tfinal;

   // Specify output frequency
   Real dtout = tfinal;
   pp.query("dtout", dtout);
   prob_opt.dtout = dtout;

   // Specify maximum number of steps between outputs
   int max_steps = 1000;
   pp.query("max_steps", max_steps);
   prob_opt.max_steps = max_steps;

   // Output integrator diagnostics to a file
   int write_diag = 1;
   pp.query("write_diag", write_diag);
   prob_opt.write_diag = write_diag;

   // Decide whether to use a preconditioner or not
   int use_preconditioner = 0;
   pp.query("use_preconditioner", use_preconditioner);
   if (nls_method==0)
     prob_opt.use_preconditioner = use_preconditioner;

   // --------------------------------------------------------------------------
   // Problem data
   // --------------------------------------------------------------------------

   // The number of cells on each side of a square domain.
   int n_cell = 256;
   pp.query("n_cell", n_cell);
   prob_data.n_cell = n_cell;

   // The domain is broken into boxes of size max_grid_size
   int max_grid_size = 64;
   pp.query("max_grid_size", max_grid_size);
   prob_data.max_grid_size = max_grid_size;

   // Advection coefficients
   Real advCoeffx = 5.0e-4;
   Real advCoeffy = 2.5e-4;
   pp.query("advCoeffx", advCoeffx);
   pp.query("advCoeffy", advCoeffy);
   prob_data.advCoeffx = advCoeffx;
   prob_data.advCoeffy = advCoeffy;

   // Diffusion coefficients
   Real diffCoeffx = 1.0e-6;
   Real diffCoeffy = 1.0e-6;
   pp.query("diffCoeffx", diffCoeffx);
   pp.query("diffCoeffy", diffCoeffy);
   prob_data.diffCoeffx = diffCoeffx;
   prob_data.diffCoeffy = diffCoeffy;

   // MLMG options
   ParmParse ppmg("mlmg");
   prob_data.mg_agglomeration = 1;
   ppmg.query("agglomeration", prob_data.mg_agglomeration);
   prob_data.mg_consolidation = 1;
   ppmg.query("consolidation", prob_data.mg_consolidation);
   prob_data.mg_max_coarsening_level = 1000;
   ppmg.query("max_coarsening_level", prob_data.mg_max_coarsening_level);
   prob_data.mg_linop_maxorder = 2;
   ppmg.query("linop_maxorder", prob_data.mg_linop_maxorder);
   prob_data.mg_max_iter = 1000;
   ppmg.query("max_iter", prob_data.mg_max_iter);
   prob_data.mg_max_fmg_iter = 1000;
   ppmg.query("max_fmg_iter", prob_data.mg_max_fmg_iter);
   prob_data.mg_verbose = 0;
   ppmg.query("verbose", prob_data.mg_verbose);
   prob_data.mg_bottom_verbose = 0;
   ppmg.query("bottom_verbose", prob_data.mg_bottom_verbose);
   prob_data.mg_use_hypre = 1;
   ppmg.query("use_hypre", prob_data.mg_use_hypre);
   prob_data.mg_hypre_interface = 3;
   ppmg.query("hypre_interface", prob_data.mg_hypre_interface);
   prob_data.mg_use_petsc = 0;
   ppmg.query("use_petsc", prob_data.mg_use_petsc);
   prob_data.mg_tol_rel = 1.0e-6;
   ppmg.query("tol_rel", prob_data.mg_tol_rel);

   // Ouput problem options and parameters
   amrex::Print()
      << "n_cell        = " << n_cell        << std::endl
      << "max_grid_size = " << max_grid_size << std::endl
      << "plot_int      = " << plot_int      << std::endl
      << "stepper       = " << stepper       << std::endl;
   if (stepper == 0)
      amrex::Print()
         << "cvode_method  = " << cvode_method  << std::endl;
   else
      amrex::Print()
         << "arkode_order  = " << arkode_order << std::endl;
   amrex::Print()
      << "rhs_adv       = " << rhs_adv  << std::endl
      << "rhs_diff      = " << rhs_diff << std::endl;
   if (fixed_dt > 0.0)
      amrex::Print()
         << "fixed_dt      = " << fixed_dt << std::endl;
   else
      amrex::Print()
         << "rtol          = " << rtol << std::endl
         << "atol          = " << atol << std::endl;
   amrex::Print()
      << "tfinal        = " << tfinal        << std::endl
      << "dtout         = " << dtout         << std::endl
      << "write_diag    = " << write_diag    << std::endl;
   if (rhs_adv > 0)
      amrex::Print()
         << "advCoeffx     = " << advCoeffx << std::endl
         << "advCoeffy     = " << advCoeffy << std::endl;
   if (rhs_diff > 0)
      amrex::Print()
         << "diffCoeffx    = " << diffCoeffx << std::endl
         << "diffCoeffy    = " << diffCoeffy << std::endl;
   if ((rhs_adv > 0) && (rhs_diff > 0) && (rhs_adv != rhs_diff))
     if (rhs_diff > 1) {
      amrex::Print() << "ImEx treatment: implicit advection and explicit diffusion" << std::endl;
     } else {
      amrex::Print() << "ImEx treatment: implicit diffusion and explicit advection" << std::endl;
     }
   if (use_preconditioner)
      amrex::Print()
        << "preconditioning enabled" << std::endl
        << "  mlmg.agglomeration        = " << prob_data.mg_agglomeration << std::endl
        << "  mlmg.consolidation        = " << prob_data.mg_consolidation << std::endl
        << "  mlmg.max_coarsening_level = " << prob_data.mg_max_coarsening_level << std::endl
        << "  mlmg.linop_maxorder       = " << prob_data.mg_linop_maxorder << std::endl
        << "  mlmg.max_iter             = " << prob_data.mg_max_iter << std::endl
        << "  mlmg.max_fmg_iter         = " << prob_data.mg_max_fmg_iter << std::endl
        << "  mlmg.verbose              = " << prob_data.mg_verbose << std::endl
        << "  mlmg.bottom_verbose       = " << prob_data.mg_bottom_verbose << std::endl
        << "  mlmg.use_hypre            = " << prob_data.mg_use_hypre << std::endl
        << "  mlmg.hypre_interface      = " << prob_data.mg_hypre_interface << std::endl
        << "  mlmg.use_petsc            = " << prob_data.mg_use_petsc << std::endl
        << "  mlmg.tol_rel              = " << prob_data.mg_tol_rel << std::endl;

}


int main(int argc, char* argv[])
{
   amrex::Initialize(argc,argv);

   DoProblem();

   amrex::Finalize();
   return 0;
}

void DoProblem()
{
   // What time is it now?  We'll use this to compute total run time.
   Real strt_time = amrex::second();

   // Set problem data and options
   ProblemData prob_data;
   ProblemOpt  prob_opt;
   ParseInputs(prob_opt, prob_data);

   // Make BoxArray and Geometry
   BoxArray ba;
   Geometry geom;
   SetUpGeometry(ba, geom, prob_data);

   // How Boxes are distrubuted among MPI processes
   DistributionMapping dm(ba);
   prob_data.dmap = &dm;

   // Allocate the solution MultiFab
   int nGhost = 1;  // number of ghost cells for each array
   int nComp  = 1;  // number of components for each array
   MultiFab sol(ba, dm, nComp, nGhost);

   // Allocate the linear solver coefficient MultiFabs
   MultiFab acoef(ba, dm, nComp, nGhost);
   MultiFab bcoef(ba, dm, nComp, nGhost);
   acoef = 1.0;
   bcoef = 1.0;
   prob_data.acoef = &acoef;
   prob_data.bcoef = &bcoef;

   // Build the flux MultiFabs
   Array<MultiFab, AMREX_SPACEDIM> flux;
   for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
   {
      // flux(dir) has one component, zero ghost cells, and is nodal in
      // direction dir
      BoxArray edge_ba = ba;
      edge_ba.surroundingNodes(dir);
      flux[dir].define(edge_ba, dm, 1, 0);
   }
   prob_data.flux = &flux;

   // Create an N_Vector wrapper for the solution MultiFab
   sunindextype length = nComp * prob_data.n_cell * prob_data.n_cell;
   N_Vector nv_sol     = N_VMake_Multifab(length, &sol);

   // Set the initial condition
   FillInitConds2D(sol, geom);

   // Integrate in time
   switch (prob_opt.stepper)
   {
   case 0:
      ComputeSolutionCV(nv_sol, &prob_opt, &prob_data);
      break;
   case 1:
      ComputeSolutionARK(nv_sol, &prob_opt, &prob_data);
      break;
   default:
      amrex::Print() << "Invalid stepper option" << std::endl;
      return;
   }

   // Call the timer again and compute the maximum difference between the start
   // time and stop time over all processors
   Real stop_time = amrex::second() - strt_time;
   const int IOProc = ParallelDescriptor::IOProcessorNumber();
   ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

   // Tell the I/O Processor to write out the "run time"
   amrex::Print() << "Run time = " << stop_time << std::endl;
}
