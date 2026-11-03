#ifndef CNS_H_
#define CNS_H_

#include "CNS_index.H"
#include "CNS_parm.H"
#include "cns_prob.H"

#include <AMReX_AmrLevel.H>
#include <AMReX_ErrorList.H>
#include <AMReX_FluxRegister.H>

#include <array>
#include <memory>

// One level of the AMR hierarchy.  Amr owns a CNS object per level and calls
// the virtual functions below to initialize, advance, regrid and write data.
class CNS
    :
    public amrex::AmrLevel
{
public:

    CNS () = default;
    CNS (amrex::Amr&            papa,
         int                    lev,
         const amrex::Geometry& level_geom,
         const amrex::BoxArray& bl,
         const amrex::DistributionMapping& dm,
         amrex::Real            time);
    ~CNS () override = default;

    CNS (const CNS& rhs) = delete;
    CNS (CNS&& rhs) = delete;
    CNS& operator= (const CNS& rhs) = delete;
    CNS& operator= (CNS&& rhs) = delete;

    // Restart from a checkpoint file.
    void restart (amrex::Amr& papa, std::istream& is, bool bReadSpecial = false) override;

    // Initialize data on this level from another CNS (during regrid).
    void init (amrex::AmrLevel& old) override;

    // Initialize data on a new level from the coarser level (during regrid).
    void init () override;

    // Initialize grid data at problem start-up.
    void initData () override;

    // Advance this level by one time step.
    amrex::Real advance (amrex::Real time, amrex::Real dt, int iteration, int ncycle) override;

    void computeInitialDt (int                                  finest_level,
                           int                                  sub_cycle,
                           amrex::Vector<int>&                  n_cycle,
                           const amrex::Vector<amrex::IntVect>& ref_ratio,
                           amrex::Vector<amrex::Real>&          dt_level,
                           amrex::Real                          stop_time) override;

    void computeNewDt (int                                  finest_level,
                       int                                  sub_cycle,
                       amrex::Vector<int>&                  n_cycle,
                       const amrex::Vector<amrex::IntVect>& ref_ratio,
                       amrex::Vector<amrex::Real>&          dt_min,
                       amrex::Vector<amrex::Real>&          dt_level,
                       amrex::Real                          stop_time,
                       int                                  post_regrid_flag) override;

    // Reflux and average down after this level has caught up with the next finer level.
    void post_timestep (int iteration) override;

    // Called after regridding.  Nothing to do here.
    void post_regrid (int /*lbase*/, int /*new_finest*/) override {}

    // Called after each coarse time step.
    void postCoarseTimeStep (amrex::Real time) override;

    // Called after the initial hierarchy is built.
    void post_init (amrex::Real stop_time) override;

    // Tag cells for refinement.
    void errorEst (amrex::TagBoxArray& tb, int clearval, int tagval, amrex::Real time,
                   int n_error_buf = 0, int ngrow = 0) override;

    // Define the state data, boundary conditions and derived quantities.
    static void variableSetUp ();
    static void variableCleanUp ();

    // Volume-integrated rho, momentum and rho*E on level 0, which holds the
    // average of the finer levels.
    [[nodiscard]] std::array<amrex::Real,NEQNS> conservedTotals () const;

    // True if no mass can enter or leave the domain.
    static bool isClosedDomain ();

    // The functions below launch GPU kernels.  CUDA does not allow device
    // lambdas in private or protected member functions, so they are public.

    void avgDown ();

    amrex::Real estTimeStep ();

    // Recompute rho*e and temperature from the evolved variables.
    void computeTemp (amrex::MultiFab& State, int ng);

    void compute_dSdt (const amrex::MultiFab& S, amrex::MultiFab& dSdt, amrex::Real dt,
                       amrex::FluxRegister* fr_as_crse, amrex::FluxRegister* fr_as_fine);

    void printTotal () const;

    CNS& getLevel (int lev) { return dynamic_cast<CNS&>(parent->getLevel(lev)); }

    static constexpr int NUM_GROW = 2;  // ghost cells needed by the hydro scheme

    enum StateDataType {
        State_Type = 0,
        NUM_STATE_DATA_TYPE
    };

    static Parm parm;
    static ProbParm prob_parm;

private:

    static void read_params ();

    // Holds coarse and fine fluxes at the boundary between this level and the
    // next coarser level.
    std::unique_ptr<amrex::FluxRegister> flux_reg;

    // Runtime parameters
    static amrex::BCRec phys_bc;
    static int verbose;
    static amrex::Real cfl;
    static int rk_order;
    static int do_reflux;
    static amrex::Vector<amrex::AMRErrorTag> error_tags;
};

// Implemented by each problem in Exec/<Problem>/cns_prob.cpp

// Set default runtime parameters (domain, boundary conditions, ...).  It runs
// inside amrex::Initialize, so the inputs file and command line override them.
void cns_prob_defaults ();

// Read the problem parameters (prob.*).
void cns_prob_read_params (ProbParm& prob_parm);

#endif
