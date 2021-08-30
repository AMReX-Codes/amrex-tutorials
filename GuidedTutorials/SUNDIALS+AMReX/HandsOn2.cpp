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
  Implementation file for 'intermediate' ARKode integration of 2D
  Advection-Diffusion example problem.  This treats the diffusion
  portion of the problem implicitly, but allows advection to be either
  implicit (via DIRK methods) or explicit (via ARK-IMEX methods).
  The implicit portion may be solved using either Newton's method or
  an Anderson-accelerated fixed-point nonlinear solver.  When using
  Newton's method, the underlying linear systems are solved using a
  basic GMRES iterative linear solver.  This program allows for either
  fixed time-step sizes or temporal adaptivity.
  --------------------------------------------------------------------*/


#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <arkode/arkode_arkstep.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "HandsOn2.h"
#include "Utilities.h"

#include "NVector_Multifab.h"

using namespace amrex;


void ComputeSolutionARK(N_Vector nv_sol, ProblemOpt* prob_opt,
                        ProblemData* prob_data)
{
   // Exit immediately if "help" was requested
   if (prob_opt->help > 0)  return;

  // Extract problem data and options
   Geometry* geom         = prob_data->geom;
   int       plot_int     = prob_opt->plot_int;
   int       arkode_order = prob_opt->arkode_order;
   int       nls_method   = prob_opt->nls_method;
   int       nls_max_iter = prob_opt->nls_max_iter;
   int       nls_fp_acc   = prob_opt->nls_fp_acc;
   int       ls_max_iter  = prob_opt->ls_max_iter;
   int       rhs_adv      = prob_opt->rhs_adv;
   Real      rtol         = prob_opt->rtol;
   Real      atol         = prob_opt->atol;
   Real      fixed_dt     = prob_opt->fixed_dt;
   Real      tfinal       = prob_opt->tfinal;
   Real      dtout        = prob_opt->dtout;
   int       max_steps    = prob_opt->max_steps;
   int       write_diag   = prob_opt->write_diag;

   // initial time, number of outputs, and error flag
   Real time = 0.0;
   int  nout = ceil(tfinal/dtout);
   int  ier  = 0;

   // Write a plotfile of the initial data
   if (plot_int > 0) {
      const std::string& pltfile = amrex::Concatenate("plt", 0, 5);
      MultiFab* sol = NV_MFAB(nv_sol);
      WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, time, 0);
   }

   // Create the ARK stepper     ***** UPDATED FROM HandsOn1 *****
   void* arkode_mem = NULL;
   if (rhs_adv > 1) {
      // explicit advection and implicit diffusion
      arkode_mem = ARKStepCreate(ComputeRhsAdv, ComputeRhsDiff, time, nv_sol);
   } else {
      // implicit advection and diffusion
      arkode_mem = ARKStepCreate(NULL, ComputeRhsAdvDiff, time, nv_sol);
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
      diagfp = fopen("HandsOn2_diagnostics.txt", "w");
      ARKStepSetDiagnostics(arkode_mem, diagfp);
   }

   // Attach nonlinear/linear solvers as needed     ***** UPDATED FROM HandsOn1 *****
   if (nls_method == 0) {
      // Create and attach GMRES linear solver for Newton
      SUNLinearSolver LS = SUNLinSol_SPGMR(nv_sol, PREC_NONE, ls_max_iter);
      ier = ARKStepSetLinearSolver(arkode_mem, LS, NULL);
      if (ier != ARKLS_SUCCESS) {
         amrex::Print() << "Creation of linear solver unsuccessful" << std::endl;
         return;
      }
   } else {
      // Create and attach fixed-point nonlinear solver
      SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(nv_sol, nls_fp_acc);
      ier = ARKStepSetNonlinearSolver(arkode_mem, NLS);
      if (ier != ARK_SUCCESS) {
         amrex::Print() << "Creation of nonlinear solver unsuccessful" << std::endl;
         return;
      }
   }

   // Set max number of nonlinear iterations     ***** UPDATED FROM HandsOn1 *****
   ier = ARKStepSetMaxNonlinIters(arkode_mem, nls_max_iter);
   if (ier != ARK_SUCCESS) {
      amrex::Print() << "Error setting max number of nonlinear iterations" << std::endl;
      return;
   }

   // Advance the solution in time
   Real tout = time + dtout; // first output time
   Real tret;                // return time
   for (int iout=0; iout < nout; iout++) {

      ier = ARKStepEvolve(arkode_mem, tout, nv_sol, &tret, ARK_NORMAL);
      if (ier < 0) {
         amrex::Print() << "Error in ARKStepEvolve" << std::endl;
         return;
      }

      // Get integration stats
      long nfe_evals, nfi_evals;
      ARKStepGetNumRhsEvals(arkode_mem, &nfe_evals, &nfi_evals);
      if (nfe_evals > 0)     /***** UPDATED FROM HandsOn1 *****/
         amrex::Print() << "t = " << std::setw(5) << tret
                        << "  explicit evals = " << std::setw(7) << nfe_evals
                        << "  implicit evals = " << std::setw(7) << nfi_evals
                        << std::endl;
      else
         amrex::Print() << "t = " << std::setw(5) << tret
                        << "  RHS evals = " << std::setw(7) << nfi_evals
                        << std::endl;

      // Write output
      if (plot_int > 0) {
         const std::string& pltfile = amrex::Concatenate("plt", iout+1, 5);
         MultiFab* sol = NV_MFAB(nv_sol);
         WriteSingleLevelPlotfile(pltfile, *sol, {"u"}, *geom, tret, iout+1);
      }

      // Update output time
      tout += dtout;
      if (tout > tfinal) tout = tfinal;
   }

   // Output final solution statistics     ***** UPDATED FROM HandsOn1 *****
   long int nst, nst_a, nfe, nfi, nsetups, nli, nJv, nlcf, nni, ncfn, netf;
   nst = nst_a = nfe = nfi = nsetups = nli = nJv = nlcf = nni = ncfn = netf = 0;
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
   }
   amrex::Print() << "\nFinal Solver Statistics:\n"
                  << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n";
   if (nfe > 0)
      amrex::Print() << "   Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << "\n";
   else
      amrex::Print() << "   Total RHS evals = " << nfi << "\n";
   amrex::Print() << "   Total number of nonlinear iterations = " << nni << "\n"
                  << "   Total number of nonlinear solver convergence failures = " << ncfn << "\n"
                  << "   Total number of error test failures = " << netf << "\n";
   if (nls_method == 0) {
     amrex::Print() << "   Total linear solver setups = " << nsetups << "\n"
                    << "   Total linear iterations = " << nli << "\n"
                    << "   Total number of Jacobian-vector products = " << nJv << "\n"
                    << "   Total number of linear solver convergence failures = " << nlcf << "\n";
   }

   // Close diagnostics file
   if (write_diag)  fclose(diagfp);
}


void ParseInputs(ProblemOpt& prob_opt, ProblemData& prob_data) {

   // ParmParse is way of reading inputs from the inputs file
   ParmParse pp;

   // --------------------------------------------------------------------------
   // Problem options     ***** UPDATED FROM HandsOn1 *****
   // --------------------------------------------------------------------------

   // Store 'help' request
   int help = 0;  // do not output help message
   pp.query("help", help);
   prob_opt.help = help;

   // Enable (>0) or disable (<0) writing output files
   int plot_int = -1; // plots off
   pp.query("plot_int", plot_int);
   prob_opt.plot_int = plot_int;

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

   // Specify RHS splitting (diffusion is always implicit)
   int rhs_adv  = 1;        // implicit advection
   pp.query("rhs_adv", rhs_adv);
   if (rhs_adv > 0)
      prob_opt.rhs_adv = rhs_adv;

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
   int max_steps = 10000;
   pp.query("max_steps", max_steps);
   prob_opt.max_steps = max_steps;

   // Output integrator diagnostics to a file
   int write_diag = 1;
   pp.query("write_diag", write_diag);
   prob_opt.write_diag = write_diag;

   // --------------------------------------------------------------------------
   // Problem data
   // --------------------------------------------------------------------------

   // The number of cells on each side of a square domain.
   int n_cell = 128;
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

   // Handle 'help' request
   if (help > 0) {
     amrex::Print() << std:: endl
       << "Usage: HandsOn2.exe [fname] [options]" << std::endl
       << "Options:" << std::endl
       << "  help=1" << std::endl
       << "    Print this help message and exit." << std::endl
       << "  plot_int=<int>" << std::endl
       << "    enable (1) or disable (0) plots [default=0]." << std::endl
       << "  arkode_order=<int>" << std::endl
       << "    ARKStep method order [default=4]." << std::endl
       << "  nls_method=<int>" << std::endl
       << "    use Newton (0) or fixed-point (1) solver [default=0]." << std::endl
       << "  nls_max_iter=<int>" << std::endl
       << "    maximum number of nonlinear iterations [default=3]." << std::endl
       << "  nls_fp_acc=<int>" << std::endl
       << "    number of fixed-point acceleration vectors [default=3]." << std::endl
       << "  ls_max_iter=<int>" << std::endl
       << "    maximum number of linear iterations [default=5]." << std::endl
       << "  rhs_adv=<int>" << std::endl
       << "    treat advection implicitly (1) or explicitly (2) [default=1]." << std::endl
       << "  fixed_dt=<float>" << std::endl
       << "    use a fixed time step size (if value > 0.0) [default=-1.0]." << std::endl
       << "  rtol=<float>" << std::endl
       << "    relative tolerance for time step adaptivity [default=1e-4]." << std::endl
       << "  atol=<float>" << std::endl
       << "    absolute tolerance for time step adaptivity [default=1e-9]." << std::endl
       << "  tfinal=<float>" << std::endl
       << "    final integration time [default=1e4]." << std::endl
       << "  dtout=<float>" << std::endl
       << "    time between outputs [default=tfinal]." << std::endl
       << "  max_steps=<int>" << std::endl
       << "    maximum number of internal steps between outputs [default=10000]." << std::endl
       << "  write_diag=<int>" << std::endl
       << "    output ARKStep time step adaptivity diagnostics to a file [default=1]." << std::endl
       << "  n_cell=<int>" << std::endl
       << "    number of cells on each side of the square domain [default=128]." << std::endl
       << "  max_grid_size=<int>" << std::endl
       << "    max size of boxes in box array [default=64]." << std::endl
       << "  advCoeffx=<float>" << std::endl
       << "    advection speed in the x-direction [default=5e-4]." << std::endl
       << "  advCoeffy=<float>" << std::endl
       << "    advection speed in the y-direction [default=2.5e-4]." << std::endl
       << "  diffCoeffx=<float>" << std::endl
       << "    diffusion coefficient in the x-direction [default=1e-6]." << std::endl
       << "  diffCoeffy=<float>" << std::endl
       << "    diffusion coefficient in the y-direction [default=1e-6]." << std::endl << std::endl
       << "If a file name 'fname' is provided, it will be parsed for each of the above" << std::endl
       << "options.  If an option is specified in both the input file and on the" << std::endl
       << "command line, then the command line option takes precedence." << std::endl << std::endl;
     return;
   }

   // Ouput problem options and parameters
   amrex::Print()
      << "n_cell        = " << n_cell        << std::endl
      << "max_grid_size = " << max_grid_size << std::endl
      << "plot_int      = " << plot_int      << std::endl
      << "arkode_order  = " << arkode_order  << std::endl;
   if (rhs_adv != 1)
      amrex::Print()
         << "ImEx treatment (implicit diffusion, explicit advection)" << std::endl;
   else
      amrex::Print()
         << "fully implicit treatment" << std::endl;
   if (fixed_dt > 0.0)
      amrex::Print()
         << "fixed_dt      = " << fixed_dt << std::endl;
   else
      amrex::Print()
         << "rtol          = " << rtol << std::endl
         << "atol          = " << atol << std::endl;
   amrex::Print()
      << "tfinal        = " << tfinal     << std::endl
      << "dtout         = " << dtout      << std::endl
      << "write_diag    = " << write_diag << std::endl
      << "advCoeffx     = " << advCoeffx  << std::endl
      << "advCoeffy     = " << advCoeffy  << std::endl
      << "diffCoeffx    = " << diffCoeffx << std::endl
      << "diffCoeffy    = " << diffCoeffy << std::endl;
   if (nls_method == 0)
      amrex::Print()
         << "Newton nonlinear solver:" << std::endl
         << "  max_iter    = " << nls_max_iter << std::endl
         << "  ls_max_iter = " << ls_max_iter << std::endl;
   else
      amrex::Print()
         << "Accelerated fixed-point nonlinear solver:" << std::endl
         << "  max_iter   = " << nls_max_iter << std::endl
         << "  nls_fp_acc = " << nls_fp_acc << std::endl;
}


int main(int argc, char* argv[]) {

   amrex::Initialize(argc,argv);

   DoProblem();

   amrex::Finalize();
   return 0;
}


void DoProblem() {

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
   for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
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
   ComputeSolutionARK(nv_sol, &prob_opt, &prob_data);

   // Call the timer again and compute the maximum difference between the start
   // time and stop time over all processors
   Real stop_time = amrex::second() - strt_time;
   const int IOProc = ParallelDescriptor::IOProcessorNumber();
   ParallelDescriptor::ReduceRealMax(stop_time, IOProc);

   // Tell the I/O Processor to write out the "run time"
   amrex::Print() << "Run time = " << stop_time << std::endl;
}
