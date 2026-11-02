#ifndef AMRLEVEL_WAVE_H_
#define AMRLEVEL_WAVE_H_

#include <AMReX_AmrLevel.H>

#include <cmath>

enum StateType { State_Type = 0, NUM_STATE_TYPE };

// Initial velocity v0(x) = exp(-16 r^2) cos^6(pi r) with r = x - 1/2. It
// vanishes at x = 0 and 1 with all derivatives up to the fifth.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wave_v0 (amrex::Real x)
{
    amrex::Real const r = x - amrex::Real(0.5);
    return std::exp(amrex::Real(-16)*r*r)
        * amrex::Math::powi<6>(std::cos(amrex::Math::pi<amrex::Real>()*r));
}

// Exact velocity at time t. The solution does not depend on y and z, and
// with u = 0 at x = 0 and 1, d'Alembert's formula gives
// v(x,t) = (V(x+t) + V(x-t))/2, where V is v0 extended to an odd function
// about x = 0 and x = 1 (period 2).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wave_v_exact (amrex::Real x, amrex::Real t)
{
    auto odd_ext = [] (amrex::Real s) {
        s -= amrex::Real(2)*std::floor(s*amrex::Real(0.5)); // s in [0,2)
        return (s <= amrex::Real(1)) ? wave_v0(s) : -wave_v0(amrex::Real(2)-s);
    };
    return amrex::Real(0.5)*(odd_ext(x+t) + odd_ext(x-t));
}

// One AMR level of the wave equation solver. Amr creates one object of
// this class per level and calls its virtual functions to initialize,
// advance, regrid and write the level.
class AmrLevelWave
    : public amrex::AmrLevel
{
public:

    AmrLevelWave () = default;

    AmrLevelWave (amrex::Amr& parent, int lev, const amrex::Geometry& gm,
                  const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                  amrex::Real time);

    //! Define the state data and its boundary conditions.
    static void variableSetUp ();

    //! Clean up the state descriptors at the end of the run.
    static void variableCleanUp ();

    //! Initialize data at the start of the run.
    void initData () override;

    //! Initialize data on this level from the old level during regridding.
    void init (amrex::AmrLevel& old) override;

    //! Initialize data on a new level during regridding.
    void init () override;

    /**
     * \brief Advance this level by one step.
     *
     * \param time      time at the beginning of the step
     * \param dt        time step
     * \param iteration substep number on this level during a coarser
     *                  level step: 1 or 2 with subcycling and a refinement
     *                  ratio of 2; always 1 on level 0.
     * \param ncycle    number of substeps: 1 without subcycling or on
     *                  level 0, usually 2 or 4 on fine levels.
     */
    amrex::Real advance (amrex::Real time, amrex::Real dt,
                         int iteration, int ncycle) override;

    //! Compute the initial time steps of all levels (dt_level).
    void computeInitialDt (int finest_level, int sub_cycle,
                           amrex::Vector<int>& n_cycle,
                           const amrex::Vector<amrex::IntVect>& ref_ratio,
                           amrex::Vector<amrex::Real>& dt_level,
                           amrex::Real stop_time) override;

    //! Compute the new time steps of all levels (dt_level).
    void computeNewDt (int finest_level, int sub_cycle,
                       amrex::Vector<int>& n_cycle,
                       const amrex::Vector<amrex::IntVect>& ref_ratio,
                       amrex::Vector<amrex::Real>& dt_min,
                       amrex::Vector<amrex::Real>& dt_level,
                       amrex::Real stop_time, int post_regrid_flag) override;

    //! Synchronize with the finer level after each step.
    void post_timestep (int iteration) override;

    void post_regrid (int /*lbase*/, int /*new_finest*/) override {}

    void post_init (amrex::Real /*stop_time*/) override {}

    //! Tag cells for refinement.
    void errorEst (amrex::TagBoxArray& tags, int clearval, int tagval,
                   amrex::Real time, int n_error_buf = 0, int ngrow = 0) override;

    //! Right-hand side of the system; public for CUDA lambdas.
    void computeRHS (amrex::MultiFab& dSdt, amrex::MultiFab const& S);

    //! Max error of v on this level against the exact solution at time t.
    amrex::Real maxError (amrex::Real t);

private:

    static constexpr int ncomp  = 2; // u and v
    static constexpr int nghost = 2; // for the fourth-order Laplacian

    static int rk_order;
    static amrex::Real cfl;

    static void read_params ();

    AmrLevelWave& getLevel (int lev) {
        return static_cast<AmrLevelWave&>(parent->getLevel(lev));
    }
};

#endif
