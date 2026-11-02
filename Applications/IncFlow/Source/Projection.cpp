#include "IncFlow.H"

#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFabUtil.H>

#include <type_traits>

using namespace amrex;

namespace {
Real solver_tolerance ()
{
    return std::is_same_v<Real,float> ? 1.e-4_rt : 1.e-10_rt;
}
}

// MAC projection of the face velocity: solve lap(phi) = div(umac), then
// umac -= grad(phi). Walls have umac.n = 0, i.e., Neumann for phi.
void
IncFlow::MacProjection (Vector<Array<MultiFab,AMREX_SPACEDIM>>& umac)
{
    int const nlevs = finest_level + 1;

    Vector<MultiFab> rhs(nlevs);
    Vector<MultiFab> phi(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        rhs[lev].define(grids[lev], dmap[lev], 1, 0);
        phi[lev].define(grids[lev], dmap[lev], 1, 1);
        phi[lev].setVal(0.0_rt);
        amrex::computeDivergence(rhs[lev], GetArrOfConstPtrs(umac[lev]), Geom(lev));
    }

    // All levels are solved together (composite solve).
    MLPoisson op(LevelGeom(), LevelGrids(), LevelDmap());
    op.setDomainBC(LinOpBC(LinOpBCType::Neumann), LinOpBC(LinOpBCType::Neumann));
    for (int lev = 0; lev < nlevs; ++lev) {
        op.setLevelBC(lev, nullptr);
    }

    MLMG mlmg(op);
    mlmg.setVerbose(mg_verbose);
    mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), solver_tolerance(), 0.0_rt);

    // The fluxes are -grad(phi) on the faces.
    Vector<Array<MultiFab,AMREX_SPACEDIM>> gradphi(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            gradphi[lev][d].define(umac[lev][d].boxArray(), dmap[lev], 1, 0);
        }
    }
    mlmg.getFluxes(GetVecOfArrOfPtrs(gradphi));
    for (int lev = 0; lev < nlevs; ++lev) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            MultiFab::Add(umac[lev][d], gradphi[lev][d], 0, 0, 1, 0);
        }
    }

    for (int lev = nlevs-1; lev > 0; --lev) {
        amrex::average_down_faces(GetArrOfConstPtrs(umac[lev]), GetArrOfPtrs(umac[lev-1]),
                                  refRatio(lev-1), Geom(lev-1));
    }
}

// Nodal projection of the cell-centered velocity: solve
// div(grad(phi)) = div(vel) with phi on the nodes, then vel -= grad(phi).
// The update has no pressure gradient, so vel = u + dt grad(p) with u
// divergence free, and p = phi/dt_proj if dt_proj > 0.
void
IncFlow::NodalProjection (Vector<MultiFab>& vel, Real dt_proj)
{
    int const nlevs = finest_level + 1;

    Vector<MultiFab> rhs(nlevs);
    Vector<MultiFab> phi(nlevs);
    Vector<MultiFab> gradphi(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        BoxArray const nba = amrex::convert(grids[lev], IntVect::TheNodeVector());
        rhs[lev].define(nba, dmap[lev], 1, 0);
        phi[lev].define(nba, dmap[lev], 1, 1);
        phi[lev].setVal(0.0_rt);
        gradphi[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
    }

    // Constant sigma = 1. The empty list of factories means no embedded
    // boundaries.
    Vector<FabFactory<FArrayBox> const*> const no_eb{};
    MLNodeLaplacian op(LevelGeom(), LevelGrids(), LevelDmap(), LPInfo(), no_eb, 1.0_rt);
    op.setDomainBC(LinOpBC(LinOpBCType::Neumann), LinOpBC(LinOpBCType::Neumann));

    // rhs = div(vel) on the nodes, including the coarse/fine synchronization.
    // This also sets the velocity outside walls to zero.
    op.compRHS(GetVecOfPtrs(rhs), GetVecOfPtrs(vel), {}, {});

    MLMG mlmg(op);
    mlmg.setVerbose(mg_verbose);
    mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), solver_tolerance(), 0.0_rt);

    // The fluxes are -grad(phi) at the cell centers.
    mlmg.getFluxes(GetVecOfPtrs(gradphi));
    for (int lev = 0; lev < nlevs; ++lev) {
        MultiFab::Add(vel[lev], gradphi[lev], 0, 0, AMREX_SPACEDIM, 0);
    }
    for (int lev = nlevs-2; lev >= 0; --lev) {
        amrex::average_down(vel[lev+1], vel[lev], geom[lev+1], geom[lev],
                            0, AMREX_SPACEDIM, refRatio(lev));
    }

    if (dt_proj > 0.0_rt) {
        for (int lev = 0; lev < nlevs; ++lev) {
            MultiFab::Copy(pressure[lev], phi[lev], 0, 0, 1, 0);
            pressure[lev].mult(1.0_rt/dt_proj);
        }
    }
}
