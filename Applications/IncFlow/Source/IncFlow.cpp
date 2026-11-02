#include "IncFlow.H"
#include "prob.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>

#include <algorithm>

using namespace amrex;

namespace {
// Sets the ghost cells outside no-slip walls to the wall velocity. With
// ext_dir, a ghost cell next to the domain holds the value on the boundary.
struct WallVelocityFill
{
    GpuArray<GpuArray<Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> lo;
    GpuArray<GpuArray<Real,AMREX_SPACEDIM>,AMREX_SPACEDIM> hi;

    AMREX_GPU_DEVICE
    void operator() (IntVect const& iv, Array4<Real> const& dest, int dcomp, int numcomp,
                     GeometryData const& geom, Real /*time*/, BCRec const* bcr, int bcomp,
                     int orig_comp) const
    {
        Box const& domain = geom.Domain();
        for (int n = 0; n < numcomp; ++n) {
            BCRec const& bc = bcr[bcomp+n];
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (bc.lo(d) == BCType::ext_dir && iv[d] < domain.smallEnd(d)) {
                    dest(iv, dcomp+n) = lo[d][orig_comp+n];
                } else if (bc.hi(d) == BCType::ext_dir && iv[d] > domain.bigEnd(d)) {
                    dest(iv, dcomp+n) = hi[d][orig_comp+n];
                }
            }
        }
    }
};
using BndryFunc = GpuBndryFuncFab<WallVelocityFill>;
}

IncFlow::IncFlow ()
{
    // The AmrCore constructor has read amr.* and geometry.* and defined the
    // Geometry of every level. No grids exist yet.
    ReadParameters();

    int const nlevs_max = max_level + 1;
    vel_new.resize(nlevs_max);
    vel_old.resize(nlevs_max);
    pressure.resize(nlevs_max);

    bcs.resize(AMREX_SPACEDIM);
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            int const bc = Geom(0).isPeriodic(d) ? BCType::int_dir : BCType::ext_dir;
            bcs[n].setLo(d, bc);
            bcs[n].setHi(d, bc);
        }
    }
}

void
IncFlow::ReadParameters ()
{
    {
        ParmParse pp;
        pp.query("max_step", max_step);
        pp.query("stop_time", stop_time);
    }
    {
        ParmParse pp("amr");
        pp.query("regrid_int", regrid_int);
        pp.query("plot_file", plot_file);
        pp.query("plot_int", plot_int);
        pp.query("chk_file", chk_file);
        pp.query("chk_int", chk_int);
        pp.query("restart", restart_chkfile);
    }
    {
        ParmParse pp("incflow");
        pp.query("nu", nu);
        pp.query("cfl", cfl);
        pp.query("mg_verbose", mg_verbose);

        // Refine where |vorticity| > vort_tag[lev]. The last value is used
        // for finer levels.
        Vector<Real> vort_tag{10.0_rt};
        pp.queryarr("vort_tag", vort_tag);
        error_tags.emplace_back(vort_tag, AMRErrorTag::GREATER, "vort");
    }

    // Wall velocity, e.g., yhi.velocity = 1 0 for a lid moving in x
    char const* names[3][2] = {{"xlo","xhi"}, {"ylo","yhi"}, {"zlo","zhi"}};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Vector<Real> v;
        if (ParmParse(names[d][0]).queryarr("velocity", v)) {
            for (int n = 0; n < std::min(int(v.size()), AMREX_SPACEDIM); ++n) {
                wall_vel_lo[d][n] = v[n];
            }
        }
        v.clear();
        if (ParmParse(names[d][1]).queryarr("velocity", v)) {
            for (int n = 0; n < std::min(int(v.size()), AMREX_SPACEDIM); ++n) {
                wall_vel_hi[d][n] = v[n];
            }
        }
        // Walls do not move in the normal direction.
        wall_vel_lo[d][d] = wall_vel_hi[d][d] = 0.0_rt;
    }
}

void
IncFlow::InitData ()
{
    if (restart_chkfile.empty()) {
        // Calls ErrorEst and MakeNewLevelFromScratch level by level.
        InitFromScratch(0.0_rt);
        for (int lev = finest_level-1; lev >= 0; --lev) {
            amrex::average_down(vel_new[lev+1], vel_new[lev], geom[lev+1], geom[lev],
                                0, AMREX_SPACEDIM, refRatio(lev));
        }
        // Make the initial velocity divergence free. The pressure from this
        // projection is not meaningful.
        NodalProjection(vel_new, 0.0_rt);
        if (chk_int > 0) {
            WriteCheckpointFile();
        }
    } else {
        ReadCheckpointFile();
    }

    if (plot_int > 0) {
        WritePlotFile();
    }
}

void
IncFlow::Evolve ()
{
    int last_plot_file_step = istep;

    while (istep < max_step && t_new < stop_time * (1.0_rt - 1.e-6_rt))
    {
        if (regrid_int > 0 && max_level > 0 && istep > last_regrid_step
            && istep % regrid_int == 0)
        {
            regrid(0, t_new);
            last_regrid_step = istep;
        }

        ComputeDt();
        Advance();
        ++istep;

        amrex::Print() << "Step " << istep << " time = " << t_new << " dt = " << dt
                       << " levels = " << finest_level+1
                       << " kinetic energy = " << KineticEnergy() << "\n";

        if (plot_int > 0 && istep % plot_int == 0) {
            last_plot_file_step = istep;
            WritePlotFile();
        }
        if (chk_int > 0 && istep % chk_int == 0) {
            WriteCheckpointFile();
        }
    }

    if (plot_int > 0 && istep > last_plot_file_step) {
        WritePlotFile();
    }
}

Real
IncFlow::KineticEnergy () const
{
    // Level 0 covers the domain and holds the average of finer levels.
    Real const vol = AMREX_D_TERM(Geom(0).CellSize(0), *Geom(0).CellSize(1), *Geom(0).CellSize(2));
    Real ke = 0.0_rt;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        ke += MultiFab::Dot(vel_new[0], n, vel_new[0], n, 1, 0);
    }
    return 0.5_rt * ke * vol;
}

void
IncFlow::DefineLevel (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    vel_new[lev].define(ba, dm, AMREX_SPACEDIM, 1);
    vel_old[lev].define(ba, dm, AMREX_SPACEDIM, 1);
    pressure[lev].define(amrex::convert(ba, IntVect::TheNodeVector()), dm, 1, 0);
    vel_new[lev].setVal(0.0_rt);
    vel_old[lev].setVal(0.0_rt);
    pressure[lev].setVal(0.0_rt);
}

void
IncFlow::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                  const DistributionMapping& dm)
{
    DefineLevel(lev, ba, dm);
    t_new = t_old = time;

    auto const& vel = vel_new[lev].arrays();
    auto const geomdata = Geom(lev).data();
    ParallelFor(vel_new[lev], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        prob_initdata(i, j, k, vel[b], geomdata);
    });
}

void
IncFlow::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                 const DistributionMapping& dm)
{
    DefineLevel(lev, ba, dm);
    FillCoarsePatch(lev, time, vel_new[lev]);
}

void
IncFlow::RemakeLevel (int lev, Real time, const BoxArray& ba,
                      const DistributionMapping& dm)
{
    MultiFab new_vel(ba, dm, AMREX_SPACEDIM, 1);
    FillPatch(lev, time, new_vel);

    DefineLevel(lev, ba, dm);
    std::swap(new_vel, vel_new[lev]);
}

void
IncFlow::ClearLevel (int lev)
{
    vel_new[lev].clear();
    vel_old[lev].clear();
    pressure[lev].clear();
}

void
IncFlow::ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    MultiFab vort(grids[lev], dmap[lev], 1, 0);
    ComputeVorticity(lev, time, vort);
    vort.abs(0, 1); // AMRErrorTag::GREATER tests the value, so use the magnitude.
    for (auto const& tag : error_tags) {
        tag(tags, &vort, TagBox::CLEAR, TagBox::SET, time, lev, Geom(lev));
    }
}

void
IncFlow::ComputeVorticity (int lev, Real time, MultiFab& vort)
{
    MultiFab vel(vort.boxArray(), vort.DistributionMap(), AMREX_SPACEDIM, 1);
    FillPatch(lev, time, vel);

    auto const dxinv = Geom(lev).InvCellSizeArray();
    auto const& u = vel.const_arrays();
    auto const& w = vort.arrays();
    ParallelFor(vort, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        auto const& v = u[b];
        // Central differences of velocity component n in direction d
        auto ddx = [&] (int n, int d) {
            auto const e = IntVect::TheDimensionVector(d);
            IntVect const iv(AMREX_D_DECL(i,j,k));
            return 0.5_rt * (v(iv+e,n) - v(iv-e,n)) * dxinv[d];
        };
#if (AMREX_SPACEDIM == 2)
        w[b](i,j,k) = ddx(1,0) - ddx(0,1);
#else
        Real const wx = ddx(2,1) - ddx(1,2);
        Real const wy = ddx(0,2) - ddx(2,0);
        Real const wz = ddx(1,0) - ddx(0,1);
        w[b](i,j,k) = std::sqrt(wx*wx + wy*wy + wz*wz);
#endif
    });
}

Array<LinOpBCType,AMREX_SPACEDIM>
IncFlow::LinOpBC (LinOpBCType wall) const
{
    Array<LinOpBCType,AMREX_SPACEDIM> r;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        r[d] = Geom(0).isPeriodic(d) ? LinOpBCType::Periodic : wall;
    }
    return r;
}

void
IncFlow::GetData (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    Real const teps = (t_new - t_old) * 1.e-3_rt;
    if (t_new == t_old || std::abs(time - t_new) <= teps) {
        data.push_back(&vel_new[lev]);
        datatime.push_back(t_new);
    } else if (std::abs(time - t_old) <= teps) {
        data.push_back(&vel_old[lev]);
        datatime.push_back(t_old);
    } else {
        // FillPatch interpolates linearly in time.
        data.push_back(&vel_old[lev]);
        data.push_back(&vel_new[lev]);
        datatime.push_back(t_old);
        datatime.push_back(t_new);
    }
}

void
IncFlow::FillPatch (int lev, Real time, MultiFab& mf)
{
    BndryFunc const bndry_func(WallVelocityFill{wall_vel_lo, wall_vel_hi});
    if (lev == 0) {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        PhysBCFunct<BndryFunc> physbc(geom[0], bcs, bndry_func);
        FillPatchSingleLevel(mf, time, smf, stime, 0, 0, AMREX_SPACEDIM, geom[0], physbc, 0);
    } else {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime);
        GetData(lev  , time, fmf, ftime);

        PhysBCFunct<BndryFunc> cphysbc(geom[lev-1], bcs, bndry_func);
        PhysBCFunct<BndryFunc> fphysbc(geom[lev  ], bcs, bndry_func);
        FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime, 0, 0, AMREX_SPACEDIM,
                           geom[lev-1], geom[lev], cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), &cell_cons_interp, bcs, 0);
    }
}

void
IncFlow::FillCoarsePatch (int lev, Real time, MultiFab& mf)
{
    AMREX_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    AMREX_ALWAYS_ASSERT(cmf.size() == 1);

    BndryFunc const bndry_func(WallVelocityFill{wall_vel_lo, wall_vel_hi});
    PhysBCFunct<BndryFunc> cphysbc(geom[lev-1], bcs, bndry_func);
    PhysBCFunct<BndryFunc> fphysbc(geom[lev  ], bcs, bndry_func);
    InterpFromCoarseLevel(mf, time, *cmf[0], 0, 0, AMREX_SPACEDIM, geom[lev-1], geom[lev],
                          cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                          &cell_cons_interp, bcs, 0);
}
