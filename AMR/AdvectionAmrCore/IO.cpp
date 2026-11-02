#include "AmrCoreAdv.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <fstream>
#include <sstream>

using namespace amrex;

void
AmrCoreAdv::WritePlotFile () const
{
    std::string const plotfilename = amrex::Concatenate(plot_file, istep[0], 5);
    amrex::Print() << "Writing plotfile " << plotfilename << "\n";
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1,
                                   GetVecOfConstPtrs(phi_new), {"phi"},
                                   Geom(), t_new[0], istep, refRatio());
}

// A checkpoint is a directory, e.g., chk00010, containing a Header with the
// time stepping state and the BoxArrays, and Level_<lev>/phi MultiFab files.
void
AmrCoreAdv::WriteCheckpointFile () const
{
    std::string const checkpointname = amrex::Concatenate(chk_file, istep[0], 5);
    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    int const nlevels = finest_level+1;

    // The I/O rank creates the directories, followed by a barrier.
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    if (ParallelDescriptor::IOProcessor()) {
        std::ofstream header(checkpointname + "/Header");
        if (!header.good()) {
            amrex::FileOpenFailed(checkpointname + "/Header");
        }
        header.precision(std::numeric_limits<Real>::max_digits10);

        header << "Checkpoint file for AmrCoreAdv\n";
        header << finest_level << "\n";
        for (int lev = 0; lev < nlevels; ++lev) { header << istep[lev] << " "; }
        header << "\n";
        for (int lev = 0; lev < nlevels; ++lev) { header << dt[lev] << " "; }
        header << "\n";
        for (int lev = 0; lev < nlevels; ++lev) { header << t_new[lev] << " "; }
        header << "\n";
        for (int lev = 0; lev < nlevels; ++lev) {
            boxArray(lev).writeOn(header);
            header << "\n";
        }
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        VisMF::Write(phi_new[lev],
                     amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "phi"));
    }
}

void
AmrCoreAdv::ReadCheckpointFile ()
{
    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // The I/O rank reads the Header and broadcasts it.
    Vector<char> file_chars;
    ParallelDescriptor::ReadAndBcastFile(restart_chkfile + "/Header", file_chars);
    std::istringstream is(std::string(file_chars.dataPtr()));

    std::string line;
    std::getline(is, line);  // title

    is >> finest_level;
    for (int lev = 0; lev <= finest_level; ++lev) { is >> istep[lev]; }
    for (int lev = 0; lev <= finest_level; ++lev) { is >> dt[lev]; }
    for (int lev = 0; lev <= finest_level; ++lev) { is >> t_new[lev]; }

    for (int lev = 0; lev <= finest_level; ++lev) {
        BoxArray ba;
        ba.readFrom(is);
        DistributionMapping dm(ba, ParallelDescriptor::NProcs());

        // Tell AmrCore about the grids, then allocate the data.
        SetBoxArray(lev, ba);
        SetDistributionMap(lev, dm);
        DefineLevel(lev, ba, dm);
        t_old[lev] = t_new[lev];

        VisMF::Read(phi_new[lev],
                    amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "phi"));
    }
}
