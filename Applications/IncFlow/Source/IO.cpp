#include "IncFlow.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <fstream>
#include <limits>
#include <sstream>

using namespace amrex;

void
IncFlow::WritePlotFile ()
{
    std::string const plotfilename = amrex::Concatenate(plot_file, istep, 5);
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";

    int const nlevs = finest_level + 1;
    Vector<std::string> names{AMREX_D_DECL("velx", "vely", "velz")};
    names.emplace_back("pressure");
    names.emplace_back("vorticity");
    int const ncomp = int(names.size());

    Vector<MultiFab> plt(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        plt[lev].define(grids[lev], dmap[lev], ncomp, 0);
        MultiFab::Copy(plt[lev], vel_new[lev], 0, 0, AMREX_SPACEDIM, 0);
        amrex::average_node_to_cellcenter(plt[lev], AMREX_SPACEDIM, pressure[lev], 0, 1);
        MultiFab vort(plt[lev], amrex::make_alias, AMREX_SPACEDIM+1, 1);
        ComputeVorticity(lev, t_new, vort);
    }

    amrex::WriteMultiLevelPlotfile(plotfilename, nlevs, GetVecOfConstPtrs(plt), names,
                                   Geom(), t_new, Vector<int>(nlevs, istep), refRatio());
}

// A checkpoint is a directory, e.g., chk00010, containing a Header with the
// time stepping state and the BoxArrays, and Level_<lev>/vel MultiFab files.
// The pressure is not needed to restart.
void
IncFlow::WriteCheckpointFile () const
{
    std::string const checkpointname = amrex::Concatenate(chk_file, istep, 5);
    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    int const nlevs = finest_level + 1;

    // The I/O rank creates the directories, followed by a barrier.
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevs, true);

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream header(checkpointname + "/Header");
        if (!header.good()) {
            amrex::FileOpenFailed(checkpointname + "/Header");
        }
        header.precision(std::numeric_limits<Real>::max_digits10);

        header << "Checkpoint file for IncFlow\n"
               << finest_level << "\n"
               << istep << " " << last_regrid_step << "\n"
               << t_new << " " << dt << "\n";
        for (int lev = 0; lev < nlevs; ++lev) {
            boxArray(lev).writeOn(header);
            header << "\n";
        }
    }

    for (int lev = 0; lev < nlevs; ++lev) {
        VisMF::Write(vel_new[lev],
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "vel"));
    }
}

void
IncFlow::ReadCheckpointFile ()
{
    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // The I/O rank reads the Header and broadcasts it.
    Vector<char> file_chars;
    ParallelDescriptor::ReadAndBcastFile(restart_chkfile + "/Header", file_chars);
    std::istringstream is(std::string(file_chars.dataPtr()));

    std::string line;
    std::getline(is, line);  // title

    is >> finest_level >> istep >> last_regrid_step >> t_new >> dt;
    t_old = t_new;

    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxArray ba;
        ba.readFrom(is);
        DistributionMapping dm(ba, ParallelDescriptor::NProcs());

        // Tell AmrCore about the grids, then allocate the data.
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);
        DefineLevel(lev, ba, dm);

        VisMF::Read(vel_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "vel"));
    }
}
