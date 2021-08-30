#include <new>
#include <iostream>
#include <iomanip>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include "AmrLevelAdv.H"

#if defined(AMREX_NO_SENSEI_AMR_INST)
#error Incompatible AMReX library configuration! AMREX_NO_SENSEI_AMR_INST must not be defined
#endif

/**
 * This tutorial illustrates in situ processing in a simulation that computes
 * mesh based data. In this case the mesh based data is provided by an instance
 * of amrex::Amr. In situ processing is managed by an instrumentation built into
 * amrex::Amr
 */

amrex::LevelBld* getLevelBld();

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);


    Real dRunTime1 = amrex::second();

    int  max_step;
    Real strt_time;
    Real stop_time;

    {
        ParmParse pp;

        max_step  = -1;
        strt_time =  0.0;
        stop_time = -1.0;

        pp.query("max_step",max_step);
        pp.query("strt_time",strt_time);
        pp.query("stop_time",stop_time);
    }

    if (strt_time < 0.0) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0) {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    {
        Amr amr(getLevelBld());

        amr.init(strt_time,stop_time);

        while ( amr.okToContinue() && (amr.levelSteps(0) < max_step || max_step < 0) &&
                (amr.cumTime() < stop_time || stop_time < 0.0) )
        {
            //
            // Do a coarse timestep.  Recursively calls timeStep()
            //
            amr.coarseTimeStep(stop_time);
        }

        //
        // Write final checkpoint and plotfile
        //
        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
            amr.checkPoint();
        }

        if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
            amr.writePlotFile();
        }

    }

    Real dRunTime2 = amrex::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;

    amrex::Finalize();

    return 0;
}
