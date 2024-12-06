#include "AMReX_PlotFileUtil.H"

#ifndef MYFUNC_H_
#define MYFUNC_H_

using namespace amrex;

void main_main ();

void WriteCheckpoint(const int& step,
                     const Real& time,
                     const MultiFab& phi);

void ReadCheckpoint(const int& restart,
                    Real& time,
                    MultiFab& phi,
                    BoxArray& ba,
                    DistributionMapping& dm);

#endif
