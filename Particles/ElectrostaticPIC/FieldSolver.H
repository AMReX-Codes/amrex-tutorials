#ifndef FIELD_SOLVER_H_
#define FIELD_SOLVER_H_

#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFab.H>

// Solves lap(phi) = -rho/epsilon0 on a periodic domain with MLMG, and
// computes E = -grad(phi) at cell centers.
class FieldSolver
{
public:
    FieldSolver (const amrex::Geometry& geom, const amrex::BoxArray& ba,
                 const amrex::DistributionMapping& dm, amrex::Real epsilon0);

    // phi needs one ghost cell; it is filled on return.
    void solvePhi (amrex::MultiFab& phi, const amrex::MultiFab& rho);

    // efield has AMREX_SPACEDIM components and one ghost cell, filled on return.
    void computeE (amrex::MultiFab& efield, const amrex::MultiFab& phi) const;

    // epsilon0/2 * integral of |E|^2
    [[nodiscard]] amrex::Real fieldEnergy (const amrex::MultiFab& efield) const;

private:
    amrex::Geometry m_geom;
    amrex::Real m_epsilon0;
    // The operator and the solver are set up once and reused every step.
    amrex::MLPoisson m_linop;
    amrex::MLMG m_mlmg;
};

#endif
