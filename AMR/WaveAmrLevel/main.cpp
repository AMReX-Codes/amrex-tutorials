#include "AmrLevelWave.H"

#include <AMReX.H>
#include <AMReX_Amr.H>
#include <AMReX_LevelBld.H>
#include <AMReX_ParmParse.H>

#include <memory>

using namespace amrex;

// Amr uses a LevelBld to set up the state descriptors and create the
// AmrLevelWave object of each level.
class LevelBldWave
    : public LevelBld
{
    void variableSetUp () override { AmrLevelWave::variableSetUp(); }
    void variableCleanUp () override { AmrLevelWave::variableCleanUp(); }
    AmrLevel* operator() () override { return new AmrLevelWave; }
    AmrLevel* operator() (Amr& amr, int lev, const Geometry& gm,
                          const BoxArray& ba, const DistributionMapping& dm,
                          Real time) override
    {
        return new AmrLevelWave(amr, lev, gm, ba, dm, time);
    }
};

#ifndef AMREX_NO_PROBINIT
// Amr calls this hook when AMReX is built with probin support. This code
// reads its parameters with ParmParse instead.
extern "C" void amrex_probinit (const int* /*init*/, const int* /*name*/, const int* /*namelen*/,
                                const amrex_real* /*problo*/, const amrex_real* /*probhi*/)
{}
#endif

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int max_step = 100000;
        Real stop_time = 1.0;
        {
            ParmParse pp;
            pp.query("max_step", max_step);
            pp.query("stop_time", stop_time);
        }

        LevelBldWave level_bld;
        auto amr = std::make_unique<Amr>(&level_bld);
        amr->init(0.0, stop_time);

        while (amr->okToContinue() && amr->levelSteps(0) < max_step
               && amr->cumTime() < stop_time)
        {
            amr->coarseTimeStep(stop_time);
        }

        // Write a final checkpoint and plotfile if they are enabled.
        if (amr->stepOfLastCheckPoint() < amr->levelSteps(0)) {
            amr->checkPoint();
        }
        if (amr->stepOfLastPlotFile() < amr->levelSteps(0)) {
            amr->writePlotFile();
        }

        // Level 0 holds the averaged fine data where it is covered.
        Real const time = amr->cumTime();
        Real const err = static_cast<AmrLevelWave&>(amr->getLevel(0)).maxError(time);
        amrex::Print() << "Finished " << amr->levelSteps(0) << " level-0 steps at time "
                       << time << " with " << amr->finestLevel()+1 << " levels\n"
                       << "Max error of v: " << err << "\n";

        if (err > Real(1.e-4)) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
