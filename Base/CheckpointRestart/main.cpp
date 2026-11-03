#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>

using amrex::Real;
using namespace amrex::literals;

namespace {

// Write a checkpoint directory, e.g., chk00010/ with
//   Header           step, time and BoxArray
//   Level_0/phi_*    MultiFab data written by VisMF
void WriteCheckpoint (std::string const& name, int step, Real time, amrex::MultiFab const& phi)
{
    // Create name/ and name/Level_0/ (on the I/O rank, then barrier).
    amrex::PreBuildDirectorHierarchy(name, "Level_", 1, true);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream header(name + "/Header");
        if (!header.good()) { amrex::FileOpenFailed(name + "/Header"); }
        header.precision(17);
        header << "HeatEquation checkpoint\n" << step << "\n" << time << "\n";
        phi.boxArray().writeOn(header);
        header << "\n";
    }

    // Every rank writes the data it owns.
    amrex::VisMF::Write(phi, amrex::MultiFabFileFullPrefix(0, name, "Level_", "phi"));
}

// Read a checkpoint written by WriteCheckpoint and define phi from it.
void ReadCheckpoint (std::string const& name, int& step, Real& time, amrex::MultiFab& phi,
                     int nghost)
{
    // Rank 0 reads the header and broadcasts it.
    amrex::Vector<char> buffer;
    amrex::ParallelDescriptor::ReadAndBcastFile(name + "/Header", buffer);
    std::istringstream is(buffer.dataPtr());

    std::string title;
    std::getline(is, title);
    is >> step >> time;
    is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    amrex::BoxArray ba;
    ba.readFrom(is);

    // A new DistributionMapping; the run may use a different number of ranks.
    amrex::DistributionMapping dm(ba);
    phi.define(ba, dm, 1, nghost);
    amrex::VisMF::Read(phi, amrex::MultiFabFileFullPrefix(0, name, "Level_", "phi"));
}

// Advance the periodic heat equation from step to nsteps.
void Advance (amrex::MultiFab& phi, amrex::Geometry const& geom, int& step, Real& time,
              int nsteps, Real dt, int chk_int, std::string const& chk_file)
{
    auto const dxinv = geom.InvCellSizeArray();
    amrex::MultiFab phi_new(phi.boxArray(), phi.DistributionMap(), 1, phi.nGrowVect());

    while (step < nsteps) {
        phi.FillBoundary(geom.periodicity());
        auto const& po = phi.const_arrays();
        auto const& pn = phi_new.arrays();
        amrex::ParallelFor(phi_new, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            auto const& p = po[b];
            pn[b](i,j,k) = p(i,j,k) + dt * (
                  AMREX_D_TERM((p(i+1,j,k) - 2.0_rt*p(i,j,k) + p(i-1,j,k)) * dxinv[0]*dxinv[0],
                             + (p(i,j+1,k) - 2.0_rt*p(i,j,k) + p(i,j-1,k)) * dxinv[1]*dxinv[1],
                             + (p(i,j,k+1) - 2.0_rt*p(i,j,k) + p(i,j,k-1)) * dxinv[2]*dxinv[2]));
        });
        std::swap(phi, phi_new);
        ++step;
        time += dt;

        if (chk_int > 0 && step % chk_int == 0) {
            std::string const name = amrex::Concatenate(chk_file, step, 5);
            amrex::Print() << "Writing checkpoint " << name << " at step " << step << "\n";
            WriteCheckpoint(name, step, time, phi);
        }
    }
}

Real max_diff (amrex::MultiFab const& a, amrex::MultiFab const& b)
{
    // b may have a different DistributionMapping, so ParallelCopy it into a's layout.
    amrex::MultiFab diff(a.boxArray(), a.DistributionMap(), 1, 0);
    diff.ParallelCopy(b, 0, 0, 1);
    amrex::MultiFab::Subtract(diff, a, 0, 0, 1, 0);
    return diff.norm0();
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int nsteps = 20;
        int chk_int = 10;            // checkpoint interval; < 0 for none
        std::string chk_file = "chk";
        std::string restart;         // checkpoint to restart from, e.g., chk00010
        bool test_restart = false;   // check that a restarted run is bitwise identical
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nsteps", nsteps);
            pp.query("chk_int", chk_int);
            pp.query("chk_file", chk_file);
            pp.query("restart", restart);
            pp.query("test_restart", test_restart);
        }

        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::Geometry geom(domain,
                             amrex::RealBox({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                            {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}),
                             amrex::CoordSys::cartesian, {AMREX_D_DECL(1, 1, 1)});
        auto const dx = geom.CellSizeArray();
        Real const dt = 0.9_rt * dx[0] * dx[0] / (2.0_rt * AMREX_SPACEDIM);
        int const nghost = 1;

        amrex::MultiFab phi;
        int step = 0;
        Real time = 0.0_rt;

        if (restart.empty()) {
            amrex::BoxArray ba(domain);
            ba.maxSize(max_grid_size);
            phi.define(ba, amrex::DistributionMapping(ba), 1, nghost);
            auto const& pa = phi.arrays();
            auto const prob_lo = geom.ProbLoArray();
            amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real r2 = 0.0_rt;
                amrex::IntVect const iv(AMREX_D_DECL(i,j,k));
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    Real x = prob_lo[d] + (iv[d] + 0.5_rt) * dx[d] - 0.5_rt;
                    r2 += x * x;
                }
                pa[b](i,j,k) = 1.0_rt + std::exp(-r2 / 0.01_rt);
            });
        } else {
            ReadCheckpoint(restart, step, time, phi, nghost);
            amrex::Print() << "Restarted from " << restart << " at step " << step
                           << ", time " << time << "\n";
        }

        Advance(phi, geom, step, time, nsteps, dt, chk_int, chk_file);

        std::string const plotfile = amrex::Concatenate("plt", step, 5);
        amrex::WriteSingleLevelPlotfile(plotfile, phi, {"phi"}, geom, time, step);
        amrex::Print() << "Finished step " << step << ", time " << time
                       << ", wrote " << plotfile << "\n";

        if (test_restart && chk_int > 0 && chk_int < nsteps) {
            // Restart from the first checkpoint and redo the remaining steps.
            int step2 = 0;
            Real time2 = 0.0_rt;
            amrex::MultiFab phi2;
            ReadCheckpoint(amrex::Concatenate(chk_file, chk_int, 5), step2, time2, phi2, nghost);
            Advance(phi2, geom, step2, time2, nsteps, dt, -1, chk_file);
            Real const restart_diff = max_diff(phi, phi2);
            amrex::Print() << "Restarted run vs. uninterrupted run: max diff = "
                           << restart_diff << ", time diff = " << time2 - time << "\n";

            // Read the plotfile back. Its DistributionMapping may differ from
            // phi's, so copy it into phi's layout before comparing.
            amrex::PlotFileData pf(plotfile);
            amrex::MultiFab phi_pf(phi.boxArray(), phi.DistributionMap(), 1, 0);
            phi_pf.ParallelCopy(pf.get(0)); // all variables of level 0, here only "phi"
            Real const plotfile_diff = max_diff(phi, phi_pf);
            amrex::Print() << "Plotfile " << plotfile << " (step " << pf.levelStep(0)
                           << ", time " << pf.time() << ") vs. memory: max diff = "
                           << plotfile_diff << "\n";

            if (restart_diff != 0.0_rt || time2 != time || plotfile_diff != 0.0_rt ||
                pf.levelStep(0) != step) {
                amrex::Print() << "FAILED\n";
                status = 1;
            } else {
                amrex::Print() << "PASSED\n";
            }
        }
    }
    amrex::Finalize();
    return status;
}
