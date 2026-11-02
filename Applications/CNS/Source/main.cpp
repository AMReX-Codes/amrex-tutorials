#include "CNS.H"

#include <AMReX.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PROB_AMR_F.H>

#include <cmath>
#include <functional>
#include <type_traits>

using namespace amrex;

LevelBld* getLevelBld ();

int main (int argc, char* argv[])
{
    // cns_prob_defaults adds problem defaults to ParmParse after it has read
    // the inputs file and command line.
    amrex::Initialize(argc, argv, std::function<void()>(cns_prob_defaults));
    int status = 0;
    {
        BL_PROFILE("main()");

        const auto strt_total = amrex::second();

        int max_step = -1;
        Real strt_time = Real(0.0);
        Real stop_time = Real(-1.0);
        {
            ParmParse pp;
            pp.query("max_step", max_step);
            pp.query("strt_time", strt_time);
            pp.query("stop_time", stop_time);
        }
        if (max_step < 0 && stop_time < Real(0.0)) {
            amrex::Abort("Neither max_step nor stop_time is non-negative.");
        }

        // Amr manages the hierarchy of CNS levels: it builds the initial
        // grids, regrids, subcycles in time, and writes plotfiles and
        // checkpoints according to the amr.* parameters.
        Amr amr(getLevelBld());
        amr.init(strt_time, stop_time);

        const Real mass0 = dynamic_cast<CNS&>(amr.getLevel(0)).conservedTotals()[0];

        while ( amr.okToContinue() &&
                (amr.levelSteps(0) < max_step || max_step < 0) &&
                (amr.cumTime() < stop_time || stop_time < Real(0.0)) )
        {
            // One coarse time step, which recursively advances finer levels.
            amr.coarseTimeStep(stop_time);
        }

        // Write final checkpoint and plotfile
        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
            amr.checkPoint();
        }
        if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
            amr.writePlotFile();
        }

        // With reflux and average down, total mass is conserved to round-off
        // unless mass can flow through the domain boundary.
        if (CNS::isClosedDomain()) {
            const Real mass1 = dynamic_cast<CNS&>(amr.getLevel(0)).conservedTotals()[0];
            const Real rel_change = std::abs(mass1-mass0) / mass0;
            const Real tol = std::is_same_v<Real,float> ? Real(1.e-4) : Real(1.e-11);
            amrex::Print() << "Relative change in total mass: " << rel_change << "\n";
            if (rel_change > tol) {
                amrex::Print() << "FAILED: mass is not conserved\n";
                status = 1;
            }
        }

        auto run_time = amrex::second() - strt_total;
        ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "Run time = " << run_time << "\n";
    }
    amrex::Finalize();
    return status;
}

#ifndef AMREX_NO_PROBINIT
// Amr calls this hook when AMReX is built with probin support.  CNS reads its
// problem parameters in cns_prob_read_params instead.
extern "C" void amrex_probinit (const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                                const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{}
#endif
