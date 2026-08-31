#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_MultiFab.H>

using namespace amrex;

void main_main ();

void InitializeData(MultiFab& phi,
                    const GpuArray<Real,AMREX_SPACEDIM> dx,
                    const GpuArray<Real,AMREX_SPACEDIM> prob_lo,
                    const GpuArray<Real,AMREX_SPACEDIM> prob_hi,
                    const Real& time,
                    const Real& Ax,
                    const Real& Ay);

void ComputeDiffusion(MultiFab& S_rhs,
                      MultiFab& S_data,
                      const Real& Dx,
                      const Real& Dy,
                      const GpuArray<Real,AMREX_SPACEDIM> dx);

void ComputeAdvection(MultiFab& S_rhs,
                      MultiFab& S_data,
                      const Real& Ax,
                      const Real& Ay,
                      const GpuArray<Real,AMREX_SPACEDIM> dx);

#endif
