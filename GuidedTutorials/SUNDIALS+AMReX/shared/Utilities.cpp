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
  Implementation file for general utility routines shared between
  hands-on lessons.
  --------------------------------------------------------------------*/

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cvode/cvode.h>
#include <arkode/arkode_arkstep.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#include "Utilities.h"

#include "NVector_Multifab.h"



using namespace amrex;

void FillInitConds2D(MultiFab& sol, const Geometry& geom)
{
   const auto dx = geom.CellSize();
   const auto prob_lo = geom.ProbLo();
   const auto prob_hi = geom.ProbHi();

   Real sigma = 0.1;
   Real a = 1.0/(sigma*sqrt(2*M_PI));
   Real b = -0.5/(sigma*sigma);

   for (MFIter mfi(sol); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& fab = sol.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);
      for (int j = lo.y; j <= hi.y; ++j) {
         Real y = prob_lo[1] + (((Real) j) + 0.5) * dx[1];

         for (int i = lo.x; i <= hi.x; ++i) {
            Real x = prob_lo[0] + (((Real) i) + 0.5) * dx[0];

            Real r = x*x + y*y;
            fab(i,j,0,0) = a * exp(b*r);
         }
      }
   }
}

void SetUpGeometry(BoxArray& ba, Geometry& geom, ProblemData& prob_data)
{
   // Extract problem options
   int n_cell = prob_data.n_cell;
   int max_grid_size = prob_data.max_grid_size;

   IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
   IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
   Box domain(dom_lo, dom_hi); // cell-centered

   // Initialize the boxarray "ba" from the single box "domain"
   ba.define(domain);

   // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a
   // direction
   ba.maxSize(max_grid_size);

   // This defines the physical box, [-1,1] in each direction.
   RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                    {AMREX_D_DECL(1.0, 1.0, 1.0)});

   // This defines a Geometry object
   Vector<int> is_periodic(AMREX_SPACEDIM, 1);  // periodic in all direction
   geom.define(domain, &real_box, CoordSys::cartesian, is_periodic.data());

   prob_data.geom = &geom;
   prob_data.grid = &ba;
}


/* ---------------------------------------------------------------------------
 * SUNDIALS RHS functions
 * ---------------------------------------------------------------------------*/

int ComputeRhsAdv(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
   // extract MultiFabs
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);

   // extract problem data
   ProblemData *prob_data = (ProblemData*) data;
   Geometry* geom = prob_data->geom;
   Real advCoeffx = prob_data->advCoeffx;
   Real advCoeffy = prob_data->advCoeffy;

   // clear the RHS
   *rhs = 0.0;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // compute advection
   ComputeAdvectionUpwind(*sol, *rhs, *geom, 0, advCoeffx, advCoeffy);

   return 0;
}

int ComputeRhsDiff(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
   // extract MultiFabs
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);

   // extract problem data
   ProblemData *prob_data = (ProblemData*) data;
   Geometry* geom = prob_data->geom;
   Array<MultiFab, AMREX_SPACEDIM>& flux = *(prob_data->flux);
   Real diffCoeffx = prob_data->diffCoeffx;
   Real diffCoeffy = prob_data->diffCoeffy;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // clear the RHS
   *rhs = 0.0;

   // compute diffusion
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 0,
                    diffCoeffx, diffCoeffy);

   return 0;
}

int ComputeRhsAdvDiff(Real t, N_Vector nv_sol, N_Vector nv_rhs, void* data)
{
   // extract MultiFabs
   MultiFab* sol = NV_MFAB(nv_sol);
   MultiFab* rhs = NV_MFAB(nv_rhs);

   // extract problem data
   ProblemData *prob_data = (ProblemData*) data;
   Geometry* geom = prob_data->geom;
   Array<MultiFab, AMREX_SPACEDIM>& flux = *(prob_data->flux);
   Real advCoeffx = prob_data->advCoeffx;
   Real advCoeffy = prob_data->advCoeffy;
   Real diffCoeffx = prob_data->diffCoeffx;
   Real diffCoeffy = prob_data->diffCoeffy;

   // clear the RHS
   *rhs = 0.0;

   // fill ghost cells
   sol->FillBoundary(geom->periodicity());

   // compute advection
   ComputeAdvectionUpwind(*sol, *rhs, *geom, 0, advCoeffx, advCoeffy);

   // compute diffusion
   ComputeDiffusion(*sol, *rhs, flux[0], flux[1], *geom, 0,
                    diffCoeffx, diffCoeffy);

   return 0;
}

/* ---------------------------------------------------------------------------
 * Advection RHS functions
 * ---------------------------------------------------------------------------*/

// Assumes ghost cells already filled
// Adds result to adv_mf MultiFab
void ComputeAdvectionUpwind(MultiFab& sol_mf, MultiFab& adv_mf, Geometry& geom,
                            int comp, Real advCoeffx, Real advCoeffy)
{
   const auto dx = geom.CellSize();
   Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
   Real dyInv = 1.0 / dx[1]; // assume same over entire mesh
   Real sideCoeffx = advCoeffx * dxInv;
   Real sideCoeffy = advCoeffy * dyInv;

   int c = comp;  // for brevity
   for (MFIter mfi(sol_mf); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& sol_fab = sol_mf.array(mfi);
      Array4<Real> const& adv_fab = adv_mf.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      // x-direction
      if (advCoeffx > 0) {
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               adv_fab(i,j,0,c) -= sideCoeffx * (sol_fab(i,j,0,c) - sol_fab(i-1,j,0,c));
            }
         }
      } else {
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               adv_fab(i,j,0,c) -= sideCoeffx * (sol_fab(i+1,j,0,c) - sol_fab(i,j,0,c));
            }
         }
      }

      // y-direction
      if (advCoeffy > 0) {
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               adv_fab(i,j,0,c) -= sideCoeffy * (sol_fab(i,j,0,c) - sol_fab(i,j-1,0,c));
            }
         }
      } else {
         for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
               adv_fab(i,j,0,c) -= sideCoeffy * (sol_fab(i,j+1,0,c) - sol_fab(i,j,0,c));
            }
         }
      }
   }
}


/* ---------------------------------------------------------------------------
 * Diffusion RHS functions
 * ---------------------------------------------------------------------------*/

// Assumes ghots cells are already filled
// Adds result to diff_mf
void ComputeDiffusion(MultiFab& sol, MultiFab& diff_mf, MultiFab& fx_mf,
                      MultiFab& fy_mf, Geometry& geom, int comp,
                      Real diffCoeffx, Real diffCoeffy)
{
   ComputeDiffFlux(sol, fx_mf, fy_mf, geom, comp, diffCoeffx, diffCoeffy);
   ComputeDivergence(diff_mf, fx_mf, fy_mf, geom, comp);
}

// Assumes ghost cells already filled
// Overwrites fx_mf and fy_mf MultiFabs
void ComputeDiffFlux(MultiFab& sol_mf, MultiFab& fx_mf, MultiFab& fy_mf,
                     Geometry& geom, int comp, Real diffCoeffx, Real diffCoeffy)
{
   const auto dx = geom.CellSize();
   Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
   Real dyInv = 1.0 / dx[1]; // assume same over entire mesh
   Real coeffX = diffCoeffx * dxInv;
   Real coeffY = diffCoeffy * dyInv;

   int c = comp;  // for brevity
   for (MFIter mfi(sol_mf); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& sol = sol_mf.array(mfi);
      Array4<Real> const& fx = fx_mf.array(mfi);
      Array4<Real> const& fy = fy_mf.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      // x-flux
      for (int j = lo.y; j <= hi.y; ++j) {
         for (int i = lo.x; i <= hi.x+1; ++i) {
            // always use zero component for flux
            fx(i,j,0,0) = coeffX * (sol(i,j,0,c) - sol(i-1,j,0,c));
         }
      }

      // y-flux
      for (int j = lo.y; j <= hi.y+1; ++j) {
         for (int i = lo.x; i <= hi.x; ++i) {
            // always use zero component for flux
            fy(i,j,0,0) = coeffY * (sol(i,j,0,c) - sol(i,j-1,0,c));
         }
      }
   }
}

// Assumes ghost cells already filled
// Adds result to div_mf MultiFab
void ComputeDivergence(MultiFab& div_mf, MultiFab& fx_mf,
                       MultiFab& fy_mf, Geometry& geom, int comp)
{
   const auto dx = geom.CellSize();
   Real dxInv = 1.0 / dx[0]; // assume same over entire mesh
   Real dyInv = 1.0 / dx[1]; // assume same over entire mesh

   int c = comp;  // for brevity
   for (MFIter mfi(div_mf); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.validbox();
      Array4<Real> const& div = div_mf.array(mfi);
      Array4<Real> const& fx = fx_mf.array(mfi);
      Array4<Real> const& fy = fy_mf.array(mfi);
      const auto lo = lbound(bx);
      const auto hi = ubound(bx);

      for (int j = lo.y; j <= hi.y; ++j) {
         for (int i = lo.x; i <= hi.x; ++i) {
            // always use zero component for flux
            div(i,j,0,c) += dxInv * (fx(i+1,j,0,0) - fx(i,j,0,0)) +
                            dyInv * (fy(i,j+1,0,0) - fy(i,j,0,0));
         }
      }
   }
}


/* ---------------------------------------------------------------------------
 * Preconditioning routines
 * ---------------------------------------------------------------------------*/

int precondition_setup(realtype tn, N_Vector u, N_Vector fu,
                       booleantype jok, booleantype *jcurPtr,
                       realtype gamma, void *user_data)
{
  return(0);
}

int precondition_solve(realtype tn, N_Vector u, N_Vector fu,
                       N_Vector r, N_Vector z,
                       realtype gamma, realtype delta,
                       int lr, void *user_data)
{
  ProblemData *prob_data = (ProblemData*) user_data;

  auto geom = *(prob_data->geom);
  auto grid = *(prob_data->grid);
  auto dmap = *(prob_data->dmap);
  auto& acoef = *(prob_data->acoef);
  auto& bcoef = *(prob_data->acoef);

  MultiFab* solution = NV_MFAB(z);
  MultiFab* rhs = NV_MFAB(r);

  LPInfo info;
  info.setAgglomeration(prob_data->mg_agglomeration);
  info.setConsolidation(prob_data->mg_consolidation);
  info.setMaxCoarseningLevel(prob_data->mg_max_coarsening_level);

  const Real tol_rel = 1.e-10;
  const Real tol_abs = 0.0;

  const int nlevels = 1;

  const Real ascalar = 1.0;
  const Real bscalar = gamma;

  MLABecLaplacian mlabec({geom}, {grid}, {dmap}, info);

  mlabec.setMaxOrder(prob_data->mg_linop_maxorder);

  // Set periodic BC
  mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                   LinOpBCType::Periodic,
                                   LinOpBCType::Periodic)},
    {AMREX_D_DECL(LinOpBCType::Periodic,
                  LinOpBCType::Periodic,
                  LinOpBCType::Periodic)});

  mlabec.setLevelBC(0, nullptr);

  mlabec.setScalars(ascalar, bscalar);

  mlabec.setACoeffs(0, acoef);

  Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      const BoxArray& ba = amrex::convert(bcoef.boxArray(),
                                          IntVect::TheDimensionVector(idim));
      face_bcoef[idim].define(ba, bcoef.DistributionMap(), 1, 0);

      switch (idim)
         {
            case 0:
               face_bcoef[idim] = prob_data->diffCoeffx;
            case 1:
               face_bcoef[idim] = prob_data->diffCoeffy;
         }
    }

  mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

  MLMG mlmg(mlabec);
  mlmg.setMaxIter(prob_data->mg_max_iter);
  mlmg.setMaxFmgIter(prob_data->mg_max_fmg_iter);
  mlmg.setVerbose(prob_data->mg_verbose);
  mlmg.setBottomVerbose(prob_data->mg_bottom_verbose);
#ifdef AMREX_USE_HYPRE
  if (prob_data->mg_use_hypre) {
    mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
    if (prob_data->mg_hypre_interface == 1)
       mlmg.setHypreInterface(amrex::Hypre::Interface::structed);
    else if (prob_data->mg_hypre_interface == 2)
       mlmg.setHypreInterface(amrex::Hypre::Interface::semi_structed);
    else
       mlmg.setHypreInterface(amrex::Hypre::Interface::ij);
  }
#endif
#ifdef AMREX_USE_PETSC
  if (prob_data->mg_use_petsc) {
    mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
  }
#endif

  mlmg.solve({solution}, {rhs}, prob_data->mg_tol_rel, tol_abs);

  return(0);
}
