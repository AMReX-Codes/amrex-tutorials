#ifndef AMR_CORE_ADV_H_
#define AMR_CORE_ADV_H_

#include <AMReX_AmrCore.H>
#include <AMReX_BCRec.H>
#include <AMReX_ErrorList.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_MultiFab.H>

#include <limits>
#include <memory>
#include <string>

// Advection of a scalar phi on an adaptive hierarchy. AmrCore owns the grids
// and geometry of all levels; this class owns the data and implements the
// virtual functions AmrCore calls when it creates, remakes or deletes levels.
class AmrCoreAdv
    : public amrex::AmrCore
{
public:

    AmrCoreAdv ();

    // Build the initial hierarchy, or read it from a checkpoint.
    void InitData ();

    // Advance to stop_time or max_step.
    void Evolve ();

    // Integral of phi over the domain.
    [[nodiscard]] amrex::Real TotalMass () const;

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

    // Allocate phi and the flux register of a level.
    void DefineLevel (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    // Time stepping
    void ComputeDt ();
    [[nodiscard]] amrex::Real EstTimeStep (int lev, amrex::Real time) const;
    void TimeStepWithSubcycling (int lev, amrex::Real time);
    void TimeStepNoSubcycling (amrex::Real time);
    void AdvanceLevel (int lev, amrex::Real time, amrex::Real dt_lev);
    void ComputeFluxes (int lev, amrex::MultiFab const& state, amrex::Real time,
                        amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& flux, bool average) const;
    void AverageDownTo (int crse_lev);

    // Fill a MultiFab including ghost cells from this level and the level below.
    void FillPatch (int lev, amrex::Real time, amrex::MultiFab& mf);
    // Fill a MultiFab entirely from the level below (for a new level).
    void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf);
    // phi_old and/or phi_new, whichever brackets time.
    void GetData (int lev, amrex::Real time, amrex::Vector<amrex::MultiFab*>& data,
                  amrex::Vector<amrex::Real>& datatime);

    // I/O
    void WritePlotFile () const;
    void WriteCheckpointFile () const;
    void ReadCheckpointFile ();

private:

    amrex::Vector<int> istep;      // steps taken on each level
    amrex::Vector<int> nsubsteps;  // substeps per step of the next coarser level
    amrex::Vector<int> last_regrid_step;

    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    amrex::Vector<amrex::MultiFab> phi_new;
    amrex::Vector<amrex::MultiFab> phi_old;

    // Physical boundary types of phi (periodic or first-order extrapolation)
    amrex::Vector<amrex::BCRec> bcs;

    // flux_reg[lev] holds the flux mismatch on the boundary between
    // levels lev and lev-1. flux_reg[0] is never used.
    amrex::Vector<std::unique_ptr<amrex::FluxRegister>> flux_reg;

    // Refinement criteria
    amrex::Vector<amrex::AMRErrorTag> error_tags;

    // Runtime parameters
    int max_step = 10;
    amrex::Real stop_time = 2.0;
    amrex::Real cfl = 0.5;
    int regrid_int = 2;
    int do_reflux = 1;
    int do_subcycle = 1;
    std::string restart_chkfile;
    std::string plot_file {"plt"};
    int plot_int = -1;
    std::string chk_file {"chk"};
    int chk_int = -1;
};

#endif
