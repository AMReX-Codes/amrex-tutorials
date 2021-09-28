#ifndef AMREX_COMPUTE_FORCE_F_H_
#define AMREX_COMPUTE_FORCE_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void interpolate_cic_two_levels(const amrex::ParticleReal* particles, int ns, int np,
                                    amrex::ParticleReal* Ex_p, amrex::ParticleReal* Ey_p,
#if AMREX_SPACEDIM == 3
                                    amrex::ParticleReal* Ez_p,
#endif
                                    const amrex::Real* Ex, const amrex::Real* Ey,
#if AMREX_SPACEDIM == 3
                                    const amrex::Real* Ez,
#endif
                                    const int* lo, const int* hi, const amrex::Real* dx,
                                    const amrex::Real* cEx, const amrex::Real* cEy,
#if AMREX_SPACEDIM == 3
                                    const amrex::Real* cEz,
#endif
                                    const int* mask,
                                    const int* clo, const int* chi, const amrex::Real* cdx,
                                    const amrex::Real* plo, const int* ng, const int* lev);
#ifdef __cplusplus
}
#endif

#endif /* _COMPUTE_FORCE_F_H_ */
