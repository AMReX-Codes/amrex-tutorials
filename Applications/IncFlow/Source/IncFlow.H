#ifndef INCFLOW_H_
#define INCFLOW_H_

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>
#include <AMReX_ErrorList.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MultiFab.H>

#include <string>

// Incompressible Navier-Stokes equations with constant density (rho = 1),
//
//     du/dt + div(u u) = -grad(p) + nu lap(u),   div(u) = 0,
//
// on an adaptive hierarchy. All levels advance together with the same dt.
// The domain is periodic or has no-slip walls with a given wall velocity.
class IncFlow
    : public amrex::AmrCore
{
public:

    IncFlow ();

    // Build the initial hierarchy, or read it from a checkpoint.
    void InitData ();

    // Advance to stop_time or max_step.
    void Evolve ();

    [[nodiscard]] amrex::MultiFab const& Velocity (int lev) const { return vel_new[lev]; }
    [[nodiscard]] amrex::Real Time () const { return t_new; }
    [[nodiscard]] amrex::Real Viscosity () const { return nu; }
    [[nodiscard]] amrex::Real KineticEnergy () const;

    // The functions below are public because CUDA does not allow device
    // lambdas in private or protected member functions.

    // AmrCore interface
    void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                  const amrex::DistributionMapping& dm) override;
    void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                 const amrex::DistributionMapping& dm) override;
    void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                      const amrex::DistributionMapping& dm) override;
    void ClearLevel (int lev) override;
    void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    void ReadParameters ();
    void DefineLevel (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    // Time stepping (Advance.cpp)
    void ComputeDt ();
    void Advance ();
    // rhs = -div(u u) + nu lap(u), without the pressure gradient.
    void ComputeRHS (amrex::Vector<amrex::MultiFab>& rhs, amrex::Real time);
    void ComputeAdvection (amrex::Vector<amrex::MultiFab>& conv,
                           amrex::Vector<amrex::MultiFab> const& vel);
    void ComputeViscousTerm (amrex::Vector<amrex::MultiFab>& visc,
                             amrex::Vector<amrex::MultiFab>& vel);

    // Projections (Projection.cpp)
    void MacProjection (amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>>& umac);
    void NodalProjection (amrex::Vector<amrex::MultiFab>& vel, amrex::Real dt_proj);

    // Fill a MultiFab including ghost cells from this level and the level below.
    void FillPatch (int lev, amrex::Real time, amrex::MultiFab& mf);
    // Fill a MultiFab entirely from the level below (for a new level).
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf);
    // vel_old and/or vel_new, whichever brackets time.
    void GetData (int lev, amrex::Real time, amrex::Vector<amrex::MultiFab*>& data,
                  amrex::Vector<amrex::Real>& datatime);

    // Vorticity (2D) or its magnitude (3D) at time
    void ComputeVorticity (int lev, amrex::Real time, amrex::MultiFab& vort);

    // Geometry, grids and distribution of levels 0 to finest_level, for the
    // linear solvers
    [[nodiscard]] amrex::Vector<amrex::Geometry> LevelGeom () const {
        return {geom.begin(), geom.begin()+finest_level+1};
    }
    [[nodiscard]] amrex::Vector<amrex::BoxArray> LevelGrids () const {
        return {grids.begin(), grids.begin()+finest_level+1};
    }
    [[nodiscard]] amrex::Vector<amrex::DistributionMapping> LevelDmap () const {
        return {dmap.begin(), dmap.begin()+finest_level+1};
    }

    // Domain boundary types for the linear solvers: Periodic, or wall
    // (the same on both sides) in non-periodic directions
    [[nodiscard]] amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM> LinOpBC (amrex::LinOpBCType wall) const;

    // I/O (IO.cpp)
    void WritePlotFile ();
    void WriteCheckpointFile () const;
    void ReadCheckpointFile ();

private:

    int istep = 0;
    int last_regrid_step = 0;
    amrex::Real t_new = 0.0;
    amrex::Real t_old = 0.0;
    amrex::Real dt = 0.0;

    // Cell-centered velocity with one ghost cell, needed by the nodal
    // projection. Pressure on the nodes, for output only.
    amrex::Vector<amrex::MultiFab> vel_new;
    amrex::Vector<amrex::MultiFab> vel_old;
    amrex::Vector<amrex::MultiFab> pressure;

    // int_dir in periodic directions, ext_dir (no-slip wall) otherwise
    amrex::Vector<amrex::BCRec> bcs;
    // Wall velocity on the low and high side of each direction
    amrex::GpuArray<amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> wall_vel_lo{};
    amrex::GpuArray<amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> wall_vel_hi{};

    // Refine where |vorticity| exceeds a threshold
    amrex::Vector<amrex::AMRErrorTag> error_tags;

    // Runtime parameters
    int max_step = 20;
    amrex::Real stop_time = 1.0;
    amrex::Real nu = 0.01;
    amrex::Real cfl = 0.5;
    int regrid_int = 2;
    int mg_verbose = 0;
    std::string restart_chkfile;
    std::string plot_file {"plt"};
    int plot_int = -1;
    std::string chk_file {"chk"};
    int chk_int = -1;
};

#endif
