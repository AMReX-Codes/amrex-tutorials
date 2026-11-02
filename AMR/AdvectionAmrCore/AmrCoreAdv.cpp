#include "AmrCoreAdv.H"
#include "Kernels.H"

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

namespace {
// Fills physical-boundary ghost cells according to the BCRec types. The
// FabFillNoOp hook would fill external Dirichlet (BCType::ext_dir) values.
using BndryFunc = GpuBndryFuncFab<FabFillNoOp>;
}

AmrCoreAdv::AmrCoreAdv ()
{
    // The AmrCore constructor has read amr.* and geometry.* and defined the
    // Geometry of every level. No grids exist yet.
    ReadParameters();

    int const nlevs_max = max_level + 1;

    istep.resize(nlevs_max, 0);
    last_regrid_step.resize(nlevs_max, 0);
    nsubsteps.resize(nlevs_max, 1);
    if (do_subcycle) {
        for (int lev = 1; lev <= max_level; ++lev) {
            nsubsteps[lev] = MaxRefRatio(lev-1);
        }
    }

    t_new.resize(nlevs_max, 0.0_rt);
    t_old.resize(nlevs_max, 0.0_rt);
    dt.resize(nlevs_max, 0.0_rt);

    phi_new.resize(nlevs_max);
    phi_old.resize(nlevs_max);
    flux_reg.resize(nlevs_max);

    // Periodic directions use int_dir. Otherwise, walls with zero gradient.
    bcs.resize(1);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        int const bc = Geom(0).isPeriodic(idim) ? BCType::int_dir : BCType::foextrap;
        bcs[0].setLo(idim, bc);
        bcs[0].setHi(idim, bc);
    }
}

void
AmrCoreAdv::ReadParameters ()
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
        ParmParse pp("adv");
        pp.query("cfl", cfl);
        pp.query("do_reflux", do_reflux);
        pp.query("do_subcycle", do_subcycle);

        // Refine where phi > phierr[lev]. The last value is used for finer levels.
        Vector<Real> phierr{0.01_rt, 0.1_rt, 0.5_rt};
        pp.queryarr("phierr", phierr);
        error_tags.emplace_back(phierr, AMRErrorTag::GREATER, "phi");
    }
}

void
AmrCoreAdv::InitData ()
{
    if (restart_chkfile.empty()) {
        // Calls ErrorEst and MakeNewLevelFromScratch level by level.
        InitFromScratch(0.0_rt);
        for (int lev = finest_level-1; lev >= 0; --lev) {
            AverageDownTo(lev);
        }
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
AmrCoreAdv::Evolve ()
{
    Real cur_time = t_new[0];
    int last_plot_file_step = istep[0];

    for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
        ComputeDt();

        if (do_subcycle) {
            TimeStepWithSubcycling(0, cur_time);
        } else {
            TimeStepNoSubcycling(cur_time);
        }

        cur_time += dt[0];
        for (int lev = 0; lev <= finest_level; ++lev) {
            t_new[lev] = cur_time;  // remove roundoff differences between levels
        }

        amrex::Print() << "Coarse step " << step+1 << " time = " << cur_time
                       << " dt = " << dt[0] << " levels = " << finest_level+1
                       << " mass = " << TotalMass() << "\n";

        if (plot_int > 0 && (step+1) % plot_int == 0) {
            last_plot_file_step = step+1;
            WritePlotFile();
        }

        if (chk_int > 0 && (step+1) % chk_int == 0) {
            WriteCheckpointFile();
        }

        if (cur_time >= stop_time - 1.e-6_rt*dt[0]) { break; }
    }

    if (plot_int > 0 && istep[0] > last_plot_file_step) {
        WritePlotFile();
    }
}

Real
AmrCoreAdv::TotalMass () const
{
    // Level 0 covers the domain and holds the average of finer levels.
    Real const vol = AMREX_D_TERM(Geom(0).CellSize(0), *Geom(0).CellSize(1), *Geom(0).CellSize(2));
    return phi_new[0].sum(0) * vol;
}

void
AmrCoreAdv::DefineLevel (int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    phi_new[lev].define(ba, dm, 1, 0);
    phi_old[lev].define(ba, dm, 1, 0);
    if (lev > 0 && do_reflux) {
        flux_reg[lev] = std::make_unique<FluxRegister>(ba, dm, refRatio(lev-1), lev, 1);
    }
}

void
AmrCoreAdv::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                     const DistributionMapping& dm)
{
    DefineLevel(lev, ba, dm);
    t_new[lev] = t_old[lev] = time;

    auto const& phi = phi_new[lev].arrays();
    auto const problo = Geom(lev).ProbLoArray();
    auto const dx = Geom(lev).CellSizeArray();
    ParallelFor(phi_new[lev], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        Real const x = problo[0] + (i+0.5_rt)*dx[0];
        Real const y = problo[1] + (j+0.5_rt)*dx[1];
        phi[b](i,j,k) = adv::initial_phi(x, y);
    });
}

void
AmrCoreAdv::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                    const DistributionMapping& dm)
{
    DefineLevel(lev, ba, dm);
    t_new[lev] = t_old[lev] = time;
    FillCoarsePatch(lev, time, phi_new[lev]);
}

void
AmrCoreAdv::RemakeLevel (int lev, Real time, const BoxArray& ba,
                         const DistributionMapping& dm)
{
    MultiFab new_state(ba, dm, 1, 0);
    FillPatch(lev, time, new_state);

    DefineLevel(lev, ba, dm);
    std::swap(new_state, phi_new[lev]);
    t_new[lev] = t_old[lev] = time;
}

void
AmrCoreAdv::ClearLevel (int lev)
{
    phi_new[lev].clear();
    phi_old[lev].clear();
    flux_reg[lev].reset();
}

void
AmrCoreAdv::ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/)
{
    for (auto const& tag : error_tags) {
        tag(tags, &phi_new[lev], TagBox::CLEAR, TagBox::SET, time, lev, Geom(lev));
    }
}

void
AmrCoreAdv::GetData (int lev, Real time, Vector<MultiFab*>& data, Vector<Real>& datatime)
{
    data.clear();
    datatime.clear();

    Real const teps = (t_new[lev] - t_old[lev]) * 1.e-3_rt;
    if (t_new[lev] == t_old[lev] || std::abs(time - t_new[lev]) <= teps) {
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_new[lev]);
    } else if (std::abs(time - t_old[lev]) <= teps) {
        data.push_back(&phi_old[lev]);
        datatime.push_back(t_old[lev]);
    } else {
        // FillPatch interpolates linearly in time.
        data.push_back(&phi_old[lev]);
        data.push_back(&phi_new[lev]);
        datatime.push_back(t_old[lev]);
        datatime.push_back(t_new[lev]);
    }
}

void
AmrCoreAdv::FillPatch (int lev, Real time, MultiFab& mf)
{
    if (lev == 0) {
        Vector<MultiFab*> smf;
        Vector<Real> stime;
        GetData(0, time, smf, stime);

        PhysBCFunct<BndryFunc> physbc(geom[0], bcs, BndryFunc(FabFillNoOp{}));
        FillPatchSingleLevel(mf, time, smf, stime, 0, 0, mf.nComp(), geom[0], physbc, 0);
    } else {
        Vector<MultiFab*> cmf, fmf;
        Vector<Real> ctime, ftime;
        GetData(lev-1, time, cmf, ctime);
        GetData(lev  , time, fmf, ftime);

        PhysBCFunct<BndryFunc> cphysbc(geom[lev-1], bcs, BndryFunc(FabFillNoOp{}));
        PhysBCFunct<BndryFunc> fphysbc(geom[lev  ], bcs, BndryFunc(FabFillNoOp{}));
        // Conservative linear interpolation from the coarse level where
        // there is no fine data.
        FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime, 0, 0, mf.nComp(),
                           geom[lev-1], geom[lev], cphysbc, 0, fphysbc, 0,
                           refRatio(lev-1), &cell_cons_interp, bcs, 0);
    }
}

void
AmrCoreAdv::FillCoarsePatch (int lev, Real time, MultiFab& mf)
{
    AMREX_ASSERT(lev > 0);

    Vector<MultiFab*> cmf;
    Vector<Real> ctime;
    GetData(lev-1, time, cmf, ctime);
    AMREX_ALWAYS_ASSERT(cmf.size() == 1);

    PhysBCFunct<BndryFunc> cphysbc(geom[lev-1], bcs, BndryFunc(FabFillNoOp{}));
    PhysBCFunct<BndryFunc> fphysbc(geom[lev  ], bcs, BndryFunc(FabFillNoOp{}));
    InterpFromCoarseLevel(mf, time, *cmf[0], 0, 0, mf.nComp(), geom[lev-1], geom[lev],
                          cphysbc, 0, fphysbc, 0, refRatio(lev-1),
                          &cell_cons_interp, bcs, 0);
}

void
AmrCoreAdv::AverageDownTo (int crse_lev)
{
    amrex::average_down(phi_new[crse_lev+1], phi_new[crse_lev],
                        geom[crse_lev+1], geom[crse_lev],
                        0, phi_new[crse_lev].nComp(), refRatio(crse_lev));
}

// Advance level lev by one step, then its finer levels by nsubsteps steps
// each (recursively), then synchronize lev with lev+1.
void
AmrCoreAdv::TimeStepWithSubcycling (int lev, Real time)
{
    // Regrid levels finer than lev every regrid_int steps of lev, unless a
    // coarser level has just done it.
    if (regrid_int > 0 && lev < max_level && istep[lev] > last_regrid_step[lev]
        && istep[lev] % regrid_int == 0)
    {
        int const old_finest = finest_level;
        regrid(lev, time);
        for (int k = lev; k <= finest_level; ++k) {
            last_regrid_step[k] = istep[k];
        }
        for (int k = old_finest+1; k <= finest_level; ++k) {
            dt[k] = dt[k-1] / nsubsteps[k];
        }
    }

    AdvanceLevel(lev, time, dt[lev]);
    ++istep[lev];

    if (Verbose()) {
        amrex::Print() << "  [Level " << lev << " step " << istep[lev] << "] advanced "
                       << CountCells(lev) << " cells with dt = " << dt[lev] << "\n";
    }

    if (lev < finest_level)
    {
        for (int i = 1; i <= nsubsteps[lev+1]; ++i) {
            TimeStepWithSubcycling(lev+1, time + (i-1)*dt[lev+1]);
        }

        if (do_reflux) {
            // Replace the coarse flux with the time and area averaged fine flux
            // on the coarse/fine boundary.
            flux_reg[lev+1]->Reflux(phi_new[lev], 1.0_rt, 0, 0, 1, geom[lev]);
        }

        AverageDownTo(lev);
    }
}

// Advance all levels with the same dt.
void
AmrCoreAdv::TimeStepNoSubcycling (Real time)
{
    if (regrid_int > 0 && max_level > 0 && istep[0] > last_regrid_step[0]
        && istep[0] % regrid_int == 0)
    {
        regrid(0, time);
        last_regrid_step[0] = istep[0];
        for (int lev = 1; lev <= finest_level; ++lev) {
            dt[lev] = dt[0];
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev) {
        AdvanceLevel(lev, time, dt[0]);
        ++istep[lev];
    }

    if (Verbose()) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            amrex::Print() << "  [Level " << lev << " step " << istep[lev] << "] advanced "
                           << CountCells(lev) << " cells with dt = " << dt[0] << "\n";
        }
    }

    for (int lev = finest_level-1; lev >= 0; --lev) {
        if (do_reflux) {
            flux_reg[lev+1]->Reflux(phi_new[lev], 1.0_rt, 0, 0, 1, geom[lev]);
        }
        AverageDownTo(lev);
    }
}

void
AmrCoreAdv::ComputeDt ()
{
    Vector<Real> dt_tmp(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        dt_tmp[lev] = EstTimeStep(lev, t_new[lev]);
    }

    // Let dt grow by at most 10% per step. The coarse dt must be stable on
    // every level after dividing by the number of substeps.
    constexpr Real change_max = 1.1_rt;
    Real dt_0 = dt_tmp[0];
    int n_factor = 1;
    for (int lev = 0; lev <= finest_level; ++lev) {
        if (dt[lev] > 0.0_rt) {
            dt_tmp[lev] = std::min(dt_tmp[lev], change_max*dt[lev]);
        }
        n_factor *= nsubsteps[lev];
        dt_0 = std::min(dt_0, n_factor*dt_tmp[lev]);
    }

    // Do not step past stop_time.
    Real const eps = 1.e-3_rt*dt_0;
    if (t_new[0] + dt_0 > stop_time - eps) {
        dt_0 = stop_time - t_new[0];
    }

    dt[0] = dt_0;
    for (int lev = 1; lev <= finest_level; ++lev) {
        dt[lev] = dt[lev-1] / nsubsteps[lev];
    }
}
