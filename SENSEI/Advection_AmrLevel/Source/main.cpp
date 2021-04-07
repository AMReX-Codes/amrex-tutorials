
#include <new>
#include <iostream>
#include <iomanip>

#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include "AmrLevelAdv.H"

#ifdef BL_USE_SENSEI_INSITU
#ifndef AMREX_USE_SENSEI_AUTO
#include <AMReX_AmrInSituBridge.H>
#ifdef AMREX_PARTICLES
#include <AMReX_AmrParticleMeshInSituBridge.H>
#endif
#endif
#endif

amrex::LevelBld* getLevelBld();

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

#ifdef BL_USE_SENSEI_INSITU
#ifndef AMREX_USE_SENSEI_AUTO
#ifdef AMREX_PARTICLES
    AmrParticleMeshInSituBridge<3,0,0,0>* insitu_bridge = nullptr;
#else
    AmrInSituBridge* insitu_bridge = nullptr;
#endif
#endif
#endif

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
#ifdef BL_USE_SENSEI_INSITU
#ifndef AMREX_USE_SENSEI_AUTO
#ifdef AMREX_PARTICLES
        insitu_bridge = new AmrParticleMeshInSituBridge<3,0,0,0>;

        // define specifications for fields on particles
        std::map<std::string, std::vector<int>> rStructs = {{"u", {0,1,2}}};
        std::map<std::string, int> iStructs;
        std::map<std::string, std::vector<int>> rArrays;
        std::map<std::string, int> iArrays;
#else
        insitu_bridge = new AmrInSituBridge;
#endif
        if (insitu_bridge->initialize())
        {
            amrex::ErrorStream() << "AmrInSituBridge::initialize : Failed to initialize." << std::endl;
            amrex::Abort();
        }
#endif
#endif

        Amr amr(getLevelBld());

        amr.init(strt_time,stop_time);

        while ( amr.okToContinue() &&
                 (amr.levelSteps(0) < max_step || max_step < 0) &&
               (amr.cumTime() < stop_time || stop_time < 0.0) )

        {
            //
            // Do a coarse timestep.  Recursively calls timeStep()
            //
            amr.coarseTimeStep(stop_time);
#ifdef BL_USE_SENSEI_INSITU
#ifndef AMREX_USE_SENSEI_AUTO
#ifdef AMREX_PARTICLES
            auto theTracer = AmrLevelAdv::theTracerPC();
            auto tracers = static_cast<amrex::ParticleContainer<AMREX_SPACEDIM,0,0,0> *>(theTracer);
            if (insitu_bridge && insitu_bridge->update(&amr, tracers, rStructs))
            {
                amrex::ErrorStream() << "AmrParticleMeshInSituBridge::update : Failed to update." << std::endl;
                amrex::Abort();
            }
#else
            if (insitu_bridge && insitu_bridge->update(&amr))
            {
                amrex::ErrorStream() << "AmrInSituBridge::update : Failed to update." << std::endl;
                amrex::Abort();
            }
#endif
#endif
#endif
        }

        // Write final checkpoint and plotfile
        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
            amr.checkPoint();
        }

        if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
            amr.writePlotFile();
        }
#ifdef BL_USE_SENSEI_INSITU
#ifndef AMREX_USE_SENSEI_AUTO
        if (insitu_bridge)
        {
            if (insitu_bridge->finalize())
                amrex::ErrorStream() << "AmrInSituBridge::finalizeInSitu : Failed to finalize." << std::endl;

            delete insitu_bridge;
            insitu_bridge = nullptr;
        }
#endif
#endif
    }

    Real dRunTime2 = amrex::second() - dRunTime1;

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Run time = " << dRunTime2 << std::endl;



    amrex::Finalize();

    return 0;
}
