#include "CNS.H"
#include "CNS_K.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <limits>

using namespace amrex;

Parm      CNS::parm;
ProbParm  CNS::prob_parm;
BCRec     CNS::phys_bc;
int       CNS::verbose   = 0;
Real      CNS::cfl       = Real(0.3);
int       CNS::rk_order  = 2;
int       CNS::do_reflux = 1;
Vector<AMRErrorTag> CNS::error_tags;

CNS::CNS (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          const DistributionMapping& dm,
          Real            time)
    : AmrLevel(papa,lev,level_geom,bl,dm,time)
{
    if (do_reflux && level > 0) {
        flux_reg = std::make_unique<FluxRegister>(grids,dmap,crse_ratio,level,NCONS);
    }
}

void
CNS::restart (Amr& papa, std::istream& is, bool bReadSpecial)
{
    AmrLevel::restart(papa,is,bReadSpecial);

    if (do_reflux && level > 0) {
        flux_reg = std::make_unique<FluxRegister>(grids,dmap,crse_ratio,level,NCONS);
    }
}

void
CNS::read_params ()
{
    ParmParse pp("cns");

    pp.query("v", verbose);
    pp.query("cfl", cfl);
    pp.query("rk_order", rk_order);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(rk_order == 2 || rk_order == 3,
                                     "cns.rk_order must be 2 or 3");
    pp.query("do_reflux", do_reflux);

    // Physical boundary types: 0 Interior (periodic), 1 Inflow, 2 Outflow,
    // 3 Symmetry, 4 SlipWall, 5 NoSlipWall
    Vector<int> lo_bc(AMREX_SPACEDIM, 0), hi_bc(AMREX_SPACEDIM, 0);
    pp.queryarr("lo_bc", lo_bc, 0, AMREX_SPACEDIM);
    pp.queryarr("hi_bc", hi_bc, 0, AMREX_SPACEDIM);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int bc : {lo_bc[dir], hi_bc[dir]}) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bc >= 0 && bc <= 5, "cns.lo_bc/hi_bc must be in [0,5]");
            // Inflow needs problem-specific boundary values. See CnsFillExtDir.
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(bc != 1, "cns: Inflow is not implemented");
        }
        phys_bc.setLo(dir, lo_bc[dir]);
        phys_bc.setHi(dir, hi_bc[dir]);
    }

    // Refine where the density jump between neighboring cells exceeds
    // refine_dengrad, on levels below refine_max_dengrad_lev.
    int refine_max_dengrad_lev = -1;
    Real refine_dengrad = Real(1.e10);
    pp.query("refine_max_dengrad_lev", refine_max_dengrad_lev);
    pp.query("refine_dengrad", refine_dengrad);
    error_tags.clear();
    if (refine_max_dengrad_lev > 0) {
        error_tags.emplace_back(refine_dengrad, AMRErrorTag::GRAD, "density",
                                AMRErrorTagInfo().SetMaxLevel(refine_max_dengrad_lev));
    }

    pp.query("gravity", parm.gravity);
    pp.query("eos_gamma", parm.eos_gamma);
    pp.query("eos_mu", parm.eos_mu);
    parm.Initialize();
}

bool
CNS::isClosedDomain ()
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        // 1: Inflow, 2: Outflow
        for (int bc : {phys_bc.lo(dir), phys_bc.hi(dir)}) {
            if (bc == 1 || bc == 2) { return false; }
        }
    }
    return true;
}

void
CNS::init (AmrLevel& old)
{
    auto& oldlev = dynamic_cast<CNS&>(old);

    Real dt_new    = parent->dtLevel(level);
    Real cur_time  = oldlev.state[State_Type].curTime();
    Real prev_time = oldlev.state[State_Type].prevTime();
    Real dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time,dt_old,dt_new);

    MultiFab& S_new = get_new_data(State_Type);
    FillPatch(old,S_new,0,cur_time,State_Type,0,NCONS);
}

void
CNS::init ()
{
    Real dt        = parent->dtLevel(level);
    Real cur_time  = getLevel(level-1).state[State_Type].curTime();
    Real prev_time = getLevel(level-1).state[State_Type].prevTime();
    Real dt_old = (cur_time - prev_time)/static_cast<Real>(parent->MaxRefRatio(level-1));
    setTimeLevel(cur_time,dt_old,dt);

    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NCONS);
}

void
CNS::initData ()
{
    BL_PROFILE("CNS::initData()");

    MultiFab& S_new = get_new_data(State_Type);
    auto const& sa = S_new.arrays();
    auto const geomdata = geom.data();
    auto const lparm = parm;
    auto const lprobparm = prob_parm;

    ParallelFor(S_new, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        cns_initdata(i, j, k, sa[b], geomdata, lparm, lprobparm);
    });

    computeTemp(S_new, 0);
}

void
CNS::computeInitialDt (int                    finest_level,
                       int                    /*sub_cycle*/,
                       Vector<int>&           n_cycle,
                       const Vector<IntVect>& /*ref_ratio*/,
                       Vector<Real>&          dt_level,
                       Real                   stop_time)
{
    // Level 0 computes dt for all levels.
    if (level > 0) { return; }

    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        dt_level[i] = getLevel(i).estTimeStep();
        n_factor   *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_level[i]);
    }

    // Do not step past stop_time.
    const Real eps = Real(0.001)*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= Real(0.0)) {
        if ((cur_time + dt_0) > (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
        }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
CNS::computeNewDt (int                    finest_level,
                   int                    /*sub_cycle*/,
                   Vector<int>&           n_cycle,
                   const Vector<IntVect>& /*ref_ratio*/,
                   Vector<Real>&          dt_min,
                   Vector<Real>&          dt_level,
                   Real                   stop_time,
                   int                    post_regrid_flag)
{
    // Level 0 computes dt for all levels at the end of a coarse step.
    if (level > 0) { return; }

    for (int i = 0; i <= finest_level; i++)
    {
        dt_min[i] = getLevel(i).estTimeStep();
    }

    if (post_regrid_flag == 1)
    {
        // Limit dt's by pre-regrid dt
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],dt_level[i]);
        }
    }
    else
    {
        // Limit the growth of dt
        constexpr Real change_max = Real(1.1);
        for (int i = 0; i <= finest_level; i++)
        {
            dt_min[i] = std::min(dt_min[i],change_max*dt_level[i]);
        }
    }

    // Find the minimum over all levels
    Real dt_0 = std::numeric_limits<Real>::max();
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_0 = std::min(dt_0,n_factor*dt_min[i]);
    }

    // Do not step past stop_time.
    const Real eps = Real(0.001)*dt_0;
    Real cur_time  = state[State_Type].curTime();
    if (stop_time >= Real(0.0)) {
        if ((cur_time + dt_0) > (stop_time - eps)) {
            dt_0 = stop_time - cur_time;
        }
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0/n_factor;
    }
}

void
CNS::post_timestep (int iteration)
{
    BL_PROFILE("CNS::post_timestep()");

    if (level < parent->finestLevel()) {
        // Replace the coarse fluxes at the coarse/fine interface with the
        // sum of the fine fluxes, which keeps the scheme conservative.
        if (do_reflux) {
            MultiFab& S = get_new_data(State_Type);
            getLevel(level+1).flux_reg->Reflux(S, Real(1.0), 0, 0, NCONS, geom);
        }
        avgDown();
    }

    // Resets the FillPatcher of the next finer level used by RK.
    AmrLevel::post_timestep(iteration);
}

void
CNS::postCoarseTimeStep (Real /*time*/)
{
    if (verbose >= 2) {
        printTotal();
    }
}

void
CNS::post_init (Real /*stop_time*/)
{
    if (level > 0) { return; }
    for (int k = parent->finestLevel()-1; k >= 0; --k) {
        getLevel(k).avgDown();
    }

    if (verbose >= 2) {
        printTotal();
    }
}

std::array<Real,NEQNS>
CNS::conservedTotals () const
{
    const MultiFab& S_new = get_new_data(State_Type);
    const Real dv = AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
    std::array<Real,NEQNS> tot{};
    for (int comp = 0; comp < NEQNS; ++comp) {
        tot[comp] = S_new.sum(comp,true) * dv;
    }
    ParallelDescriptor::ReduceRealSum(tot.data(), NEQNS);
    return tot;
}

void
CNS::printTotal () const
{
    auto tot = conservedTotals();
    amrex::Print().SetPrecision(17) << "\n[CNS] Total mass       is " << tot[0] << "\n"
                                    <<   "      Total x-momentum is " << tot[1] << "\n"
                                    <<   "      Total y-momentum is " << tot[2] << "\n"
                                    <<   "      Total z-momentum is " << tot[3] << "\n"
                                    <<   "      Total energy     is " << tot[4] << "\n";
}

void
CNS::errorEst (TagBoxArray& tags, int clearval, int tagval, Real time, int /*n_error_buf*/,
               int /*ngrow*/)
{
    BL_PROFILE("CNS::errorEst()");

    if (error_tags.empty()) { return; }

    // The gradient test needs one ghost cell.
    MultiFab rho(grids, dmap, 1, 1);
    FillPatch(*this, rho, rho.nGrow(), time, State_Type, URHO, 1, 0);

    for (auto const& tag : error_tags) {
        tag(tags, &rho, static_cast<char>(clearval), static_cast<char>(tagval),
            time, level, geom);
    }
}

void
CNS::avgDown ()
{
    BL_PROFILE("CNS::avgDown()");

    if (level == parent->finestLevel()) { return; }

    auto& fine_lev = getLevel(level+1);

    MultiFab& S_crse =          get_new_data(State_Type);
    MultiFab& S_fine = fine_lev.get_new_data(State_Type);

    amrex::average_down(S_fine, S_crse, fine_lev.geom, geom,
                        0, S_fine.nComp(), parent->refRatio(level));

    computeTemp(S_crse, 0);
}

Real
CNS::estTimeStep ()
{
    BL_PROFILE("CNS::estTimeStep()");

    const auto dx = geom.CellSizeArray();
    const MultiFab& S = get_new_data(State_Type);
    auto const& sa = S.const_arrays();
    auto const lparm = parm;

    Real estdt = ParReduce(TypeList<ReduceOpMin>{}, TypeList<Real>{}, S, IntVect(0),
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Real>
    {
        return cns_estdt(i, j, k, sa[b], dx, lparm);
    });

    ParallelDescriptor::ReduceRealMin(estdt);
    return cfl * estdt;
}

void
CNS::computeTemp (MultiFab& State, int ng)
{
    BL_PROFILE("CNS::computeTemp()");

    auto const& sa = State.arrays();
    auto const lparm = parm;
    ParallelFor(State, IntVect(ng), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        cns_compute_temperature(i, j, k, sa[b], lparm);
    });
}
