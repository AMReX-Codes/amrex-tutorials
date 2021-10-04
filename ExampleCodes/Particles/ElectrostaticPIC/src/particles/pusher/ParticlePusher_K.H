#ifndef PARTICLE_PUSHER_K_H_
#define PARTICLE_PUSHER_K_H_

#include <AMReX_REAL.H>

/**
 *
 * \brief This routine advances the particle positions using the given
 * velocity. This is used by the "push_leapfrog" routine below, and also directly to
 * desynchronize the particle positions from the velocities after particle initialization.
 *
 * Arguments:
 *     x         : the particle x position
 *     y         : the particle y position
 *     vx        : the particle x-velocity
 *     vy        : the particle y-velocity
 *     dt        : the time step
 *     plo       : the LH corner of the problem domain in physical coordinates
 *     phi       : the RH corner of the problem domain in physical coordinates
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void push_leapfrog_positions (amrex::ParticleReal& x,  amrex::ParticleReal& y,
                              amrex::ParticleReal& vx, amrex::ParticleReal& vy,
                              const amrex::Real dt,
                              const amrex::GpuArray<amrex::Real, 2>& plo,
                              const amrex::GpuArray<amrex::Real, 2>& phi) noexcept
{
    x += dt * vx;
    y += dt * vy;

    // bounce off the walls in the x...
    while (x < plo[0] or x >= phi[0])
    {
        if (x < plo[0]) {
            x = 2.0*plo[0] - x;
        }
        else {
            x = 2.0*phi[0] - x;
        }
        vx = -vx;
    }

    // .. and y directions.
    while (y < plo[1] or y >= phi[1])
    {
        if (y < plo[1]) {
            y = 2.0*plo[1] - y;
        }
        else {
            y = 2.0*phi[1] - y;
        }
        vy = -vy;
    }
}

/**
*
* \brief This routine updates the particle positions and velocities using the
* leapfrog time integration algorithm, given the electric fields at the
* particle positions. It also enforces specular reflection off the domain
* walls.
*
* Arguments:
*     x         : the particle x position
*     y         : the particle y position
*     vx        : the particle x-velocity
*     vy        : the particle y-velocity
*     Ex        : the electric field in the x-direction at the particle position
*     Ey        : the electric field in the y-direction at the particle position
*     charge    : the charge of this particle species
*     mass      : the mass of this particle species
*     dt        : the time step
*     plo       : the LH corner of the problem domain in physical coordinates
*     phi       : the RH corner of the problem domain in physical coordinates
*/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void push_leapfrog (amrex::ParticleReal& x,  amrex::ParticleReal& y,
                    amrex::ParticleReal& vx, amrex::ParticleReal& vy,
                    const amrex::ParticleReal Ex, const amrex::ParticleReal Ey,
                    const amrex::Real charge, const amrex::Real mass, const amrex::Real dt,
                    const amrex::GpuArray<amrex::Real, 2>& plo,
                    const amrex::GpuArray<amrex::Real, 2>& phi) noexcept
{
    const amrex::Real fac = charge * dt / mass;

    vx += fac * Ex;
    vy += fac * Ey;

    push_leapfrog_positions(x, y, vx, vy, dt, plo, phi);
}

#endif
