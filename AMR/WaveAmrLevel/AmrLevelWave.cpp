#include "AmrLevelWave.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>

using namespace amrex;
using namespace amrex::literals;

int  AmrLevelWave::rk_order = 4;
Real AmrLevelWave::cfl = 0.2_rt;

namespace {

// Fills ghost cells outside the domain in x, where u = 0 on the boundary.
// The other directions are periodic. With u = 0 on the face, the ghost
// values are the negated mirror images (BCType::reflect_odd would do the
// same); this shows how to write a boundary function for external
// Dirichlet data.
struct WaveBCFill
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int /*dcomp*/, const int /*numcomp*/,
                     GeometryData const& geom, const Real /*time*/,
                     const BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
    {
        const int ilo = geom.Domain().smallEnd(0);
        const int ihi = geom.Domain().bigEnd(0);
        const auto [i,j,k] = iv.dim3();
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

void wave_bcfill (Box const& bx, FArrayBox& data, const int dcomp, const int numcomp,
                  Geometry const& geom, const Real time, const Vector<BCRec>& bcr,
                  const int bcomp, const int scomp)
{
    GpuBndryFuncFab<WaveBCFill> gpu_bndry_func(WaveBCFill{});
    gpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}

}

AmrLevelWave::AmrLevelWave (Amr& amr, int lev, const Geometry& gm,
                            const BoxArray& ba, const DistributionMapping& dm,
                            Real time)
    : AmrLevel(amr, lev, gm, ba, dm, time)
{}

void
AmrLevelWave::variableSetUp ()
{
    read_params();

    // One cell-centered StateData with two components and the ghost cells
    // of the stencil. Fine-level ghost cells at the coarse/fine boundary
    // are interpolated from the coarse level with a fourth-order
    // interpolater, to match the fourth-order spatial discretization.
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, nghost, ncomp,
                           &cell_quartic_interp);

    int lo_bc[AMREX_SPACEDIM] = {AMREX_D_DECL(BCType::ext_dir,    // external Dirichlet
                                              BCType::int_dir,    // periodic
                                              BCType::int_dir)};  // periodic
    int hi_bc[AMREX_SPACEDIM] = {AMREX_D_DECL(BCType::ext_dir,
                                              BCType::int_dir,
                                              BCType::int_dir)};
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
AmrLevelWave::initData ()
{
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    MultiFab& S_new = get_new_data(State_Type);
    auto const& snew = S_new.arrays();
    ParallelFor(S_new, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k)
    {
        Real const x = problo[0] + (Real(i)+0.5_rt)*dx[0];
        snew[bi](i,j,k,0) = 0.0_rt;      // u
        snew[bi](i,j,k,1) = wave_v0(x);  // v = du/dt
    });
}

void
AmrLevelWave::init (AmrLevel& old)
{
    Real const dt_new    = parent->dtLevel(Level());
    Real const cur_time  = old.get_state_data(State_Type).curTime();
    Real const prev_time = old.get_state_data(State_Type).prevTime();
    Real const dt_old    = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    MultiFab& S_new = get_new_data(State_Type);
    FillPatch(old, S_new, 0, cur_time, State_Type, 0, ncomp);
}

void
AmrLevelWave::init ()
{
    Real const dt        = parent->dtLevel(Level());
    Real const cur_time  = getLevel(Level()-1).state[State_Type].curTime();
    Real const prev_time = getLevel(Level()-1).state[State_Type].prevTime();
    Real const dt_old    = (cur_time - prev_time)/Real(parent->MaxRefRatio(Level()-1));
    setTimeLevel(cur_time, dt_old, dt);

    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, ncomp);
}

Real
AmrLevelWave::advance (Real time, Real dt, int iteration, int ncycle)
{
    // The new data of the previous step become the old data of this step.
    for (int k = 0; k < NUM_STATE_TYPE; ++k) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    // AmrLevel::RK advances the StateData from old to new with RK3 or RK4.
    // It fills the ghost cells of every stage, including those at the
    // coarse/fine boundary, which on fine levels need coarse data at the
    // intermediate times of the RK stages. We only provide the right-hand
    // side.
    AmrLevel::RK(rk_order, State_Type, time, dt, iteration, ncycle,
                 [&] (int /*stage*/, MultiFab& dSdt, MultiFab const& S,
                      Real /*t*/, Real /*dtsub*/)
    {
        computeRHS(dSdt, S);
    });

    return dt;
}

// du/dt = v, dv/dt = Laplacian(u) with a fourth-order stencil
void
AmrLevelWave::computeRHS (MultiFab& dSdt, MultiFab const& S)
{
    const auto dxinv = Geom().InvCellSizeArray();
    AMREX_D_TERM(Real const dx2inv = dxinv[0]*dxinv[0];,
                 Real const dy2inv = dxinv[1]*dxinv[1];,
                 Real const dz2inv = dxinv[2]*dxinv[2]);
    constexpr Real c0 = -2.5_rt, c1 = 4._rt/3._rt, c2 = -1._rt/12._rt;
    auto const& sa = S.const_arrays();
    auto const& sdot = dSdt.arrays();
    ParallelFor(S, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k)
    {
        auto const& s = sa[bi];
        auto const& f = sdot[bi];
        f(i,j,k,0) = s(i,j,k,1);
        AMREX_D_TERM(Real const lapx = dx2inv*(c0*s(i,j,k,0) + c1*(s(i-1,j,k,0)+s(i+1,j,k,0))
                                               + c2*(s(i-2,j,k,0)+s(i+2,j,k,0)));,
                     Real const lapy = dy2inv*(c0*s(i,j,k,0) + c1*(s(i,j-1,k,0)+s(i,j+1,k,0))
                                               + c2*(s(i,j-2,k,0)+s(i,j+2,k,0)));,
                     Real const lapz = dz2inv*(c0*s(i,j,k,0) + c1*(s(i,j,k-1,0)+s(i,j,k+1,0))
                                               + c2*(s(i,j,k-2,0)+s(i,j,k+2,0))));
        f(i,j,k,1) = AMREX_D_TERM(lapx, +lapy, +lapz);
    });
}

void
AmrLevelWave::computeInitialDt (int finest_level, int /*sub_cycle*/,
                                Vector<int>& n_cycle,
                                const Vector<IntVect>& /*ref_ratio*/,
                                Vector<Real>& dt_level, Real stop_time)
{
    if (Level() > 0) { return; } // level 0 sets dt for all levels

    // Number of steps of each level in one level-0 step
    Vector<int> nsteps(n_cycle.size());
    std::partial_sum(n_cycle.begin(), n_cycle.end(), nsteps.begin(),
                     std::multiplies<>());

    // The wave speed is 1, so dt = cfl * dx on each level.
    Real dt_0 = std::numeric_limits<Real>::max();
    for (int ilev = 0; ilev <= finest_level; ++ilev) {
        const auto dx = parent->Geom(ilev).CellSizeArray();
        Real const dtlev = cfl * std::min({AMREX_D_DECL(dx[0],dx[1],dx[2])});
        dt_0 = std::min(dt_0, Real(nsteps[ilev]) * dtlev);
    }

    // Do not step past stop_time.
    if (stop_time > 0) {
        Real const eps = 0.001_rt * dt_0;
        Real const cur_time = get_state_data(State_Type).curTime();
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
    computeInitialDt(finest_level, sub_cycle, n_cycle, ref_ratio, dt_level, stop_time);
}

void
AmrLevelWave::post_timestep (int iteration)
{
    // Replace the coarse data under the fine level with a fourth-order
    // average of the fine data. (This code does not reflux, so the
    // solution is not conservative across levels.)
    if (Level() < parent->finestLevel()) {
        auto& fine_level = getLevel(Level()+1);
        MultiFab& S_fine = fine_level.get_new_data(State_Type);
        MultiFab& S_crse = get_new_data(State_Type);
        Real const t = get_state_data(State_Type).curTime();

        IntVect const ratio = parent->refRatio(Level());
        AMREX_ALWAYS_ASSERT(ratio == 2 || ratio == 4);
        if (ratio == 2) {
            // The fourth-order average needs one fine ghost cell.
            FillPatch(fine_level, S_fine, 1, t, State_Type, 0, ncomp);
        }
        FourthOrderInterpFromFineToCoarse(S_crse, 0, ncomp, S_fine, ratio);
    }

    AmrLevel::post_timestep(iteration);
}

void
AmrLevelWave::errorEst (TagBoxArray& tags, int /*clearval*/, int /*tagval*/,
                        Real /*time*/, int /*n_error_buf*/, int /*ngrow*/)
{
    // Refine where |v| is large, i.e., around the pulses.
    auto const& S_new = get_new_data(State_Type);
    auto const& a = tags.arrays();
    auto const& s = S_new.const_arrays();
    ParallelFor(tags, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k)
    {
        if (std::abs(s[bi](i,j,k,1)) > 0.75_rt) {
            a[bi](i,j,k) = TagBox::SET;
        }
    });
}

Real
AmrLevelWave::maxError (Real t)
{
    MultiFab const& S = get_new_data(State_Type);
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    MultiFab err(S.boxArray(), S.DistributionMap(), 1, 0);
    auto const& ea = err.arrays();
    auto const& sa = S.const_arrays();
    ParallelFor(err, [=] AMREX_GPU_DEVICE (int bi, int i, int j, int k)
    {
        Real const x = problo[0] + (Real(i)+0.5_rt)*dx[0];
        ea[bi](i,j,k) = sa[bi](i,j,k,1) - wave_v_exact(x, t);
    });
    return err.norm0();
}

void
AmrLevelWave::read_params ()
{
    ParmParse pp("wave");
    pp.query("rk_order", rk_order); // 3 or 4
    pp.query("cfl", cfl);
    AMREX_ALWAYS_ASSERT(rk_order == 3 || rk_order == 4);
}
