#ifndef MYFUNC_H_
#define MYFUNC_H_

#include <AMReX_MultiFabUtil.H>

#ifdef AMREX_USE_CUDA
#include <cufft.h>
#endif

using namespace amrex;

#ifdef AMREX_USE_CUDA
std::string cufftErrorToString (const cufftResult& err);
#endif

#endif
