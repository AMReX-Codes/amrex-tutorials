#include "AmrLevelWave.H"
#include <AMReX_ParmParse.H>
#include <numeric>

using namespace amrex;

constexpr int AmrLevelWave::ncomp;
constexpr int AmrLevelWave::nghost;
int  AmrLevelWave::verbose = 0;
int  AmrLevelWave::rk_order = 4;
Real AmrLevelWave::cfl = 0.2;

namespace {
    struct WaveBCFill {
        AMREX_GPU_DEVICE
        void operator() (const IntVect& iv, Array4<Real> const& dest,
                         const int /*dcomp*/, const int /*numcomp*/,
                         GeometryData const& geom, const Real /*time*/,
                         const BCRec* /*bcr*/, const int /*bcomp*/,
                         const int /*orig_comp*/) const
        {
            // In this test, we only need to fill the x-direction bounary,
            // because it's periodic in other directions.  We also could
            // have used BCType::reflect_odd, and then we would not need to
            // do anything here.  However, this is an example of how to fill
            // external Dirichlet BC.
            const int ilo = geom.Domain().smallEnd(0);
            const int ihi = geom.Domain().bigEnd(0);
            int i = iv[0];
            int j = iv[1];
#if (AMREX_SPACEDIM == 3)
            int k = iv[2];
#else
            int k = 0;
#endif
            if (i < ilo) {
                dest(i,j,k,0) = -dest(2*ilo-i-1,j,k,0);
                dest(i,j,k,1) = -dest(2*ilo-i-1,j,k,1);
            }
            if (i > ihi) {
                dest(i,j,k,0) = -dest(2*ihi-i+1,j,k,0);
                dest(i,j,k,1) = -dest(2*ihi-i+1,j,k,1);
            }
        }
    };

    void wave_bcfill (Box const& bx, FArrayBox& data,
                      const int dcomp, const int numcomp,
                      Geometry const& geom, const Real time,
                      const Vector<BCRec>& bcr, const int bcomp,
                      const int scomp)
    {
        GpuBndryFuncFab<WaveBCFill> gpu_bndry_func(WaveBCFill{});
        gpu_bndry_func(bx,data,dcomp,numcomp,geom,time,bcr,bcomp,scomp);
    }
}

AmrLevelWave::AmrLevelWave (Amr& amr, int lev, const Geometry& gm,
                            const BoxArray& ba, const DistributionMapping& dm,
                            Real time)
    : AmrLevel(amr,lev,gm,ba,dm,time)
{}

AmrLevelWave::~AmrLevelWave () {}

void
AmrLevelWave::variableSetUp ()
{
    read_params();

    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, nghost, ncomp,
                           &cell_quartic_interp);

    int lo_bc[BL_SPACEDIM] = {AMREX_D_DECL(BCType::ext_dir,    // external Dirichlet
                                           BCType::int_dir,    // periodic
                                           BCType::int_dir) }; // periodic
    int hi_bc[BL_SPACEDIM] = {AMREX_D_DECL(BCType::ext_dir,
                                           BCType::int_dir,
                                           BCType::int_dir) };
    Vector<BCRec> bcs(ncomp, BCRec(lo_bc, hi_bc));

    StateDescriptor::BndryFunc bndryfunc(wave_bcfill);
    bndryfunc.setRunOnGPU(true);

    desc_lst.setComponent(State_Type, 0, {"u", "v"}, bcs, bndryfunc);
}

void
AmrLevelWave::variableCleanUp ()
{
    desc_lst.clear();
}


void
AmrLevelWave::init (AmrLevel &old)
{
    Real dt_new    = parent->dtLevel(Level());
    Real cur_time  = old.get_state_data(State_Type).curTime();
    Real prev_time = old.get_state_data(State_Type).prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    for (int k = 0; k < NUM_STATE_TYPE; ++k) {
        MultiFab& S_new = get_new_data(k);
        FillPatch(old, S_new, 0, cur_time, k, 0, ncomp);
    }
}

void
AmrLevelWave::init ()
{
    Real dt        = parent->dtLevel(Level());
    Real cur_time  = getLevel(Level()-1).state[State_Type].curTime();
    Real prev_time = getLevel(Level()-1).state[State_Type].prevTime();
    Real dt_old = (cur_time - prev_time)/(Real)parent->MaxRefRatio(Level()-1);
    setTimeLevel(cur_time,dt_old,dt);

    for (int k = 0; k < NUM_STATE_TYPE; ++k) {
        MultiFab& S_new = get_new_data(k);
        FillCoarsePatch(S_new, 0, cur_time, k, 0, ncomp);
    }
}

void
AmrLevelWave::computeInitialDt (int finest_level, int /*sub_cycle*/,
                                Vector<int>& n_cycle,
                                const Vector<IntVect>& /*ref_ratio*/,
                                Vector<Real>& dt_level, Real stop_time)
{
    if (Level() > 0) { return; } // Level 0 does this for every level.

    Vector<int> nsteps(n_cycle.size()); // Total number of steps in one level 0 step
    std::partial_sum(n_cycle.begin(), n_cycle.end(), nsteps.begin(),
                     std::multiplies<int>());

    Real dt_0 = std::numeric_limits<Real>::max();
    for (int ilev = 0; ilev <= finest_level; ++ilev) {
        const auto dx = parent->Geom(ilev).CellSizeArray();
        Real dtlev = cfl * std::min({AMREX_D_DECL(dx[0],dx[1],dx[2])});
        dt_0 = std::min(dt_0, nsteps[ilev] * dtlev);
    }
    // dt_0 will be the time step on level 0 (unless limited by stop_time).

    if (stop_time > 0) {
        // Limit dt's by the value of stop_time.
        const Real eps = 0.001 * dt_0;
        const Real cur_time = get_state_data(State_Type).curTime();
        if ((cur_time + dt_0) > (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
        }
    }

    for (int ilev = 0; ilev <= finest_level; ++ilev) {
        dt_level[ilev] = dt_0 / Real(nsteps[ilev]);
    }
}

void
AmrLevelWave::computeNewDt (int finest_level, int sub_cycle,
                            Vector<int>& n_cycle,
                            const Vector<IntVect>& ref_ratio,
                            Vector<Real>& /*dt_min*/, Vector<Real>& dt_level,
                            Real stop_time, int /*post_regrid_flag*/)
{
    // For this code we can just call computeInitialDt.
    computeInitialDt(finest_level, sub_cycle, n_cycle, ref_ratio, dt_level,
                     stop_time);
}

void
AmrLevelWave::post_timestep (int iteration)
{
    if (Level() < parent->finestLevel()) {
        auto& fine_level = getLevel(Level()+1);
        MultiFab& S_fine = fine_level.get_new_data(State_Type);
        MultiFab& S_crse =      this->get_new_data(State_Type);
        Real t = get_state_data(State_Type).curTime();

        IntVect ratio = parent->refRatio(Level());
        AMREX_ASSERT(ratio == 2 || ratio == 4);
        if (ratio == 2) {
            // Need to fill one ghost cell for the high-order interpolation below
            FillPatch(fine_level, S_fine, 1, t, State_Type, 0, ncomp);
        }

        FourthOrderInterpFromFineToCoarse(S_crse, 0, 2, S_fine, ratio);
    }

    AmrLevel::post_timestep(iteration);
}

void
AmrLevelWave::errorEst (TagBoxArray& tags, int /*clearval*/, int /*tagval*/,
                        Real /*time*/, int /*n_error_buf*/, int /*ngrow*/)
{
    auto const& S_new = get_new_data(State_Type);

    const char tagval = TagBox::SET;
    auto const& a = tags.arrays();
    auto const& s = S_new.const_arrays();
    amrex::ParallelFor(tags, [=] AMREX_GPU_DEVICE
        (int bi, int i, int j, int k)
    {
        // Just an example, not necessarily good choice.
        if (amrex::Math::abs(s[bi](i,j,k,1)) > 0.75) {
            a[bi](i,j,k) = tagval;
        }
    });
}

void
AmrLevelWave::read_params ()
{
    ParmParse pp("wave");
    pp.query("v", verbose); // Could use this to control verbosity during the run
    pp.query("rk_order", rk_order);
    pp.query("cfl", cfl);
}
