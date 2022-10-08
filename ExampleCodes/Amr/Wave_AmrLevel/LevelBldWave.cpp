#include "AmrLevelWave.H"
#include <AMReX_LevelBld.H>

using namespace amrex;

class LevelBldWave
    : public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (Amr& amr, int lev, const Geometry& gm,
                                  const BoxArray& ba, const DistributionMapping& dm,
                                  Real time) override;
};

LevelBldWave Wave_bld;

LevelBld*
getLevelBld ()
{
    return &Wave_bld;
}

void
LevelBldWave::variableSetUp ()
{
    AmrLevelWave::variableSetUp();
}

void
LevelBldWave::variableCleanUp ()
{
    AmrLevelWave::variableCleanUp();
}

AmrLevel*
LevelBldWave::operator() ()
{
    return new AmrLevelWave;
}

AmrLevel*
LevelBldWave::operator() (Amr& amr, int lev, const Geometry& gm,
                          const BoxArray& ba, const DistributionMapping& dm,
                          Real time)
{
    return new AmrLevelWave(amr, lev, gm, ba, dm, time);
}
