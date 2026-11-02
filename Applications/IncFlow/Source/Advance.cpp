#include "IncFlow.H"
#include "Advection_K.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>

#include <algorithm>
#include <limits>

using namespace amrex;

void
IncFlow::ComputeDt ()
{
    // Walls also move the fluid, e.g., at the start of a lid-driven cavity.
    GpuArray<Real,AMREX_SPACEDIM> umax{};
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            umax[n] = std::max({umax[n], std::abs(wall_vel_lo[d][n]), std::abs(wall_vel_hi[d][n])});
        }
    }

    Real dt_new = std::numeric_limits<Real>::max();
    for (int lev = 0; lev <= finest_level; ++lev) {
        auto const dx = Geom(lev).CellSizeArray();
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            Real const u = std::max(umax[n], vel_new[lev].norminf(n, 0));
            if (u > 0.0_rt) {
                dt_new = std::min(dt_new, cfl * dx[n] / u);
            }
        }
        if (nu > 0.0_rt) {
            Real const h = *std::min_element(dx.begin(), dx.end());
            dt_new = std::min(dt_new, cfl * h * h / (2.0_rt * AMREX_SPACEDIM * nu));
        }
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt_new < std::numeric_limits<Real>::max(),
                                     "No velocity and no viscosity to set dt");

    // Let dt grow by at most 10% per step, and do not step past stop_time.
    if (dt > 0.0_rt) {
        dt_new = std::min(dt_new, 1.1_rt * dt);
    }
    if (t_new + dt_new > stop_time - 1.e-3_rt * dt_new) {
        dt_new = stop_time - t_new;
    }
    dt = dt_new;
}

// MOL predictor-corrector. The update has no pressure gradient; the nodal
// projection removes the gradient part of the velocity, and its potential
// is the pressure.
void
IncFlow::Advance ()
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        std::swap(vel_old[lev], vel_new[lev]);
    }
    t_old = t_new;
    t_new = t_old + dt;

    int const nlevs = finest_level + 1;
    Vector<MultiFab> rhs_old(nlevs);
    Vector<MultiFab> rhs_new(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        rhs_old[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
        rhs_new[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
    }

    // Predictor: u* = P(u^n + dt R(u^n))
    ComputeRHS(rhs_old, t_old);
    for (int lev = 0; lev < nlevs; ++lev) {
        MultiFab::LinComb(vel_new[lev], 1.0_rt, vel_old[lev], 0, dt, rhs_old[lev], 0,
                          0, AMREX_SPACEDIM, 0);
    }
    NodalProjection(vel_new, dt);

    // Corrector: u^{n+1} = P(u^n + dt/2 (R(u^n) + R(u*)))
    ComputeRHS(rhs_new, t_new);
    for (int lev = 0; lev < nlevs; ++lev) {
        MultiFab::LinComb(vel_new[lev], 1.0_rt, vel_old[lev], 0, 0.5_rt*dt, rhs_old[lev], 0,
                          0, AMREX_SPACEDIM, 0);
        MultiFab::Saxpy(vel_new[lev], 0.5_rt*dt, rhs_new[lev], 0, 0, AMREX_SPACEDIM, 0);
    }
    NodalProjection(vel_new, dt);
}

void
IncFlow::ComputeRHS (Vector<MultiFab>& rhs, Real time)
{
    int const nlevs = finest_level + 1;

    // Velocity at time with two ghost cells for the slopes, filled from the
    // same level, the level below, or the wall values
    Vector<MultiFab> vel(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        vel[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 2);
        FillPatch(lev, time, vel[lev]);
    }

    ComputeAdvection(rhs, vel);

    if (nu > 0.0_rt) {
        Vector<MultiFab> visc(nlevs);
        for (int lev = 0; lev < nlevs; ++lev) {
            visc[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
        }
        ComputeViscousTerm(visc, vel);
        for (int lev = 0; lev < nlevs; ++lev) {
            MultiFab::Add(rhs[lev], visc[lev], 0, 0, AMREX_SPACEDIM, 0);
        }
    }
}

// conv = -div(umac u), with u on faces upwinded by the divergence-free
// face velocity umac
void
IncFlow::ComputeAdvection (Vector<MultiFab>& conv, Vector<MultiFab> const& vel)
{
    int const nlevs = finest_level + 1;

    Vector<Array<MultiFab,AMREX_SPACEDIM>> umac(nlevs);
    Vector<Array<MultiFab,AMREX_SPACEDIM>> flux(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            BoxArray const ba = amrex::convert(grids[lev], IntVect::TheDimensionVector(d));
            umac[lev][d].define(ba, dmap[lev], 1, 0);
            flux[lev][d].define(ba, dmap[lev], AMREX_SPACEDIM, 0);
        }
    }

    // Walls and domain extent of each level
    Vector<GpuArray<incflow::DirInfo,AMREX_SPACEDIM>> dirinfo(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        Box const& domain = Geom(lev).Domain();
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            bool const wall = !Geom(lev).isPeriodic(d);
            dirinfo[lev][d] = {wall, wall, domain.smallEnd(d), domain.bigEnd(d)};
        }
    }

    // Normal velocity on the faces, extrapolated from both sides
    for (int lev = 0; lev < nlevs; ++lev) {
        auto const& q = vel[lev].const_arrays();
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            auto const di = dirinfo[lev][d];
            auto const& um = umac[lev][d].arrays();
            ParallelFor(umac[lev][d], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real ul, ur;
                incflow::face_values(IntVect(AMREX_D_DECL(i,j,k)), d, d, q[b], di, ul, ur);
                um[b](i,j,k) = incflow::normal_velocity(ul, ur);
            });
        }
    }

    // Make it divergence free.
    MacProjection(umac);

    // Fluxes of all velocity components
    for (int lev = 0; lev < nlevs; ++lev) {
        auto const& q = vel[lev].const_arrays();
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            auto const di = dirinfo[lev][d];
            auto const& um = umac[lev][d].const_arrays();
            auto const& fx = flux[lev][d].arrays();
            ParallelFor(flux[lev][d], IntVect(0), AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
            {
                Real ql, qr;
                incflow::face_values(IntVect(AMREX_D_DECL(i,j,k)), n, d, q[b], di, ql, qr);
                Real const u = um[b](i,j,k);
                fx[b](i,j,k,n) = u * incflow::upwind(ql, qr, u);
            });
        }
    }

    // On faces covered by a finer level, use the average of the fine fluxes,
    // so that the coarse/fine update is conservative.
    for (int lev = nlevs-1; lev > 0; --lev) {
        amrex::average_down_faces(GetArrOfConstPtrs(flux[lev]), GetArrOfPtrs(flux[lev-1]),
                                  refRatio(lev-1), Geom(lev-1));
    }

    for (int lev = 0; lev < nlevs; ++lev) {
        amrex::computeDivergence(conv[lev], GetArrOfConstPtrs(flux[lev]), Geom(lev));
        conv[lev].mult(-1.0_rt);
    }
}

// visc = nu lap(u), applied with a linear operator so that walls (Dirichlet)
// and coarse/fine boundaries are handled by AMReX.
void
IncFlow::ComputeViscousTerm (Vector<MultiFab>& visc, Vector<MultiFab>& vel)
{
    int const nlevs = finest_level + 1;

    // Only the operator is needed, no multigrid hierarchy.
    LPInfo info;
    info.setMaxCoarseningLevel(0);
    MLABecLaplacian op(LevelGeom(), LevelGrids(), LevelDmap(), info, {}, AMREX_SPACEDIM);
    op.setDomainBC(LinOpBC(LinOpBCType::Dirichlet), LinOpBC(LinOpBCType::Dirichlet));
    for (int lev = 0; lev < nlevs; ++lev) {
        // The ghost cells outside walls hold the wall velocity.
        op.setLevelBC(lev, &vel[lev]);
    }
    // alpha a u - beta div(b grad(u)) with alpha = 0, beta = b = 1
    op.setScalars(0.0_rt, 1.0_rt);
    for (int lev = 0; lev < nlevs; ++lev) {
        op.setBCoeffs(lev, 1.0_rt);
    }

    MLMG mlmg(op);
    mlmg.apply(GetVecOfPtrs(visc), GetVecOfPtrs(vel));
    for (int lev = 0; lev < nlevs; ++lev) {
        visc[lev].mult(-nu);
    }
}
