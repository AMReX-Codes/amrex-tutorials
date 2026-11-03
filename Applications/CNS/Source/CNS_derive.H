#ifndef CNS_DERIVE_H_
#define CNS_DERIVE_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

// Derived quantities.  data holds the state components registered with
// derive_lst.addComponent, and der receives the result starting at dcomp.

void cns_derpres (amrex::MultiFab& der, int dcomp, int ncomp,
                  const amrex::MultiFab& data, const amrex::Geometry& geom,
                  amrex::Real time, const int* bcrec, int level);

void cns_dervel (amrex::MultiFab& der, int dcomp, int ncomp,
                 const amrex::MultiFab& data, const amrex::Geometry& geom,
                 amrex::Real time, const int* bcrec, int level);

#endif
