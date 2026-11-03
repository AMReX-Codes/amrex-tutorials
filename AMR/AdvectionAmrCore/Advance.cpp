#include "AmrCoreAdv.H"
#include "Kernels.H"

#include <AMReX_ParReduce.H>

using namespace amrex;

// Advance one level by dt_lev with Heun's method (second-order SSP
// Runge-Kutta), and record the fluxes in the flux registers.
void
AmrCoreAdv::AdvanceLevel (int lev, Real time, Real dt_lev)
{
    std::swap(phi_old[lev], phi_new[lev]);
    t_old[lev] = time;
    t_new[lev] = time + dt_lev;

    // Two ghost cells for the limited slopes of the upwind cells.
    MultiFab Sborder(grids[lev], dmap[lev], 1, 2);
    Array<MultiFab,AMREX_SPACEDIM> flux;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        flux[idim].define(amrex::convert(grids[lev], IntVect::TheDimensionVector(idim)),
                          dmap[lev], 1, 0);
    }

    auto const& pold = phi_old[lev].const_arrays();
    auto const& pnew = phi_new[lev].arrays();
    AMREX_D_TERM(auto const& fx = flux[0].const_arrays();,
                 auto const& fy = flux[1].const_arrays();,
                 auto const& fz = flux[2].const_arrays();)
    GpuArray<Real,AMREX_SPACEDIM> dtdx;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        dtdx[idim] = dt_lev * geom[lev].InvCellSize(idim);
    }

    // Stage 0: phi_new = phi_old - dt div F(phi_old)
    // Stage 1: phi_new = phi_old - dt div (F(phi_old) + F(phi_new))/2, where
    //          FillPatch at t_new uses the stage 0 result stored in phi_new.
    for (int stage = 0; stage < 2; ++stage) {
        Real const stage_time = time + stage*dt_lev;
        FillPatch(lev, stage_time, Sborder);
        ComputeFluxes(lev, Sborder, stage_time, flux, stage == 1);

        ParallelFor(phi_new[lev], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            pnew[b](i,j,k) = pold[b](i,j,k)
                - (AMREX_D_TERM(  dtdx[0] * (fx[b](i+1,j,k) - fx[b](i,j,k)),
                                + dtdx[1] * (fy[b](i,j+1,k) - fy[b](i,j,k)),
                                + dtdx[2] * (fz[b](i,j,k+1) - fz[b](i,j,k))));
        });
    }

    if (do_reflux) {
        // The registers hold flux * area * dt. The coarse level of an
        // interface subtracts its flux, and the fine level adds its fluxes.
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Real area = 1.0_rt;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (d != idim) { area *= geom[lev].CellSize(d); }
            }
            if (lev < finest_level) {
                flux_reg[lev+1]->CrseInit(flux[idim], idim, 0, 0, 1, -dt_lev*area);
            }
            if (lev > 0) {
                flux_reg[lev]->FineAdd(flux[idim], idim, 0, 0, 1, dt_lev*area);
            }
        }
    }

    // Sborder and flux are freed on return while kernels may still be running.
    Gpu::streamSynchronize();
}

// Upwind fluxes u*phi on all faces of level lev. With average = true, the
// result is averaged with the fluxes already in flux.
void
AmrCoreAdv::ComputeFluxes (int lev, MultiFab const& state, Real time,
                           Array<MultiFab,AMREX_SPACEDIM>& flux, bool average) const
{
    auto const& q = state.const_arrays();
    auto const problo = Geom(lev).ProbLoArray();
    auto const dx = Geom(lev).CellSizeArray();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        auto const& f = flux[idim].arrays();
        Dim3 const e = IntVect::TheDimensionVector(idim).dim3();
        ParallelFor(flux[idim], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            Real const u = adv::face_velocity(idim, i, j, problo, dx, time);
            Real const fnew = adv::upwind_flux(q[b], i, j, k, e, u);
            f[b](i,j,k) = average ? 0.5_rt*(f[b](i,j,k) + fnew) : fnew;
        });
    }
}

// Stable time step on level lev: cfl / max over cells of sum_d |u_d|/dx_d.
Real
AmrCoreAdv::EstTimeStep (int lev, Real time) const
{
    auto const problo = Geom(lev).ProbLoArray();
    auto const dx = Geom(lev).CellSizeArray();
    auto const dxinv = Geom(lev).InvCellSizeArray();
    Real rate = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, phi_new[lev], IntVect(0),
        [=] AMREX_GPU_DEVICE (int, int i, int j, int) -> GpuTuple<Real>
        {
            Real r = 0.0_rt;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                IntVect const e = IntVect::TheDimensionVector(idim);
                Real const ulo = adv::face_velocity(idim, i,      j,      problo, dx, time);
                Real const uhi = adv::face_velocity(idim, i+e[0], j+e[1], problo, dx, time);
                r += amrex::max(std::abs(ulo), std::abs(uhi)) * dxinv[idim];
            }
            return {r};
        });
    ParallelDescriptor::ReduceRealMax(rate);
    return (rate > 0.0_rt) ? cfl / rate : stop_time;
}
