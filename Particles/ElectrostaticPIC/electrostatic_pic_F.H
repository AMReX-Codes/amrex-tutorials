#ifndef AMREX_COMPUTE_FORCE_F_H_
#define AMREX_COMPUTE_FORCE_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void compute_E_nodal (const int* lo, const int* hi,
                          const amrex::Real* phi,
                          amrex::Real* Ex, amrex::Real* Ey,
#if AMREX_SPACEDIM == 3
                          amrex::Real* Ez,
#endif
                          const amrex::Real* dx);

    void deposit_cic(const amrex::ParticleReal* particles, int ns, int np,
                     const amrex::ParticleReal* weights,
                     const amrex::Real* charge,
                     amrex::Real* rho, const int* lo, const int* hi,
                     const amrex::Real* plo, const amrex::Real* dx,
                     const int* ng);

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

    void interpolate_cic(const amrex::ParticleReal* particles, int ns, int np,
                         amrex::ParticleReal* Ex_p, amrex::ParticleReal* Ey_p,
#if AMREX_SPACEDIM == 3
                         amrex::ParticleReal* Ez_p,
#endif
                         const amrex::Real* Ex, const amrex::Real* Ey,
#if AMREX_SPACEDIM == 3
                         const amrex::Real* Ez,
#endif
                         const int* lo, const int* hi,
                         const amrex::Real* plo, const amrex::Real* dx,
                         const int* ng);

#ifdef __cplusplus
}
#endif

#endif /* _COMPUTE_FORCE_F_H_ */
