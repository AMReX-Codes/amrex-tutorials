#include "CNS.H"
#include "CNS_hydro_K.H"

using namespace amrex;

Real
CNS::advance (Real time, Real dt, int iteration, int ncycle)
{
    BL_PROFILE("CNS::advance()");

    // The current new data become the old data.
    for (int i = 0; i < NUM_STATE_DATA_TYPE; ++i) {
        state[i].allocOldData();
        state[i].swapTimeLevels(dt);
    }

    // As the coarse side, this level adds its fluxes to the register owned by
    // the next finer level.  As the fine side, it adds to its own register.
    FluxRegister* fr_as_crse = nullptr;
    if (do_reflux && level < parent->finestLevel()) {
        fr_as_crse = getLevel(level+1).flux_reg.get();
        fr_as_crse->setVal(Real(0.0));
    }

    FluxRegister* fr_as_fine = nullptr;
    if (do_reflux && level > 0) {
        fr_as_fine = flux_reg.get();
    }

    // AmrLevel::RK fills the ghost cells of each stage, including
    // coarse/fine boundaries interpolated in space and time.
    RK(rk_order, State_Type, time, dt, iteration, ncycle,
       // Given state S, compute dSdt.  dtsub is this stage's weight times dt.
       [&] (int /*stage*/, MultiFab& dSdt, MultiFab const& S, Real /*t*/, Real dtsub)
       {
           compute_dSdt(S, dSdt, dtsub, fr_as_crse, fr_as_fine);
       },
       // After each stage, update the derived components.
       [&] (int /*stage*/, MultiFab& S) { computeTemp(S, 0); });

    return dt;
}

void
CNS::compute_dSdt (const MultiFab& S, MultiFab& dSdt, Real dt,
                   FluxRegister* fr_as_crse, FluxRegister* fr_as_fine)
{
    BL_PROFILE("CNS::compute_dSdt()");

    const auto dx = geom.CellSizeArray();
    const auto dxinv = geom.InvCellSizeArray();
    auto const lparm = parm;

    Array<MultiFab,AMREX_SPACEDIM> fluxes;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        fluxes[dir].define(amrex::convert(S.boxArray(),IntVect::TheDimensionVector(dir)),
                           S.DistributionMap(), NCONS, 0);
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto const& s = S.const_array(mfi);
        auto const& dsdt = dSdt.array(mfi);

        // Temporary data live in the asynchronous arena; they are freed safely
        // after the kernels using them finish.  Each tile computes all fluxes
        // of its faces, because faces on tile boundaries are shared.
        Array<FArrayBox,AMREX_SPACEDIM> flxtmp;
        GpuArray<Array4<Real>,AMREX_SPACEDIM> flx;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            flxtmp[dir].resize(amrex::surroundingNodes(bx,dir), NCONS, The_Async_Arena());
            flx[dir] = flxtmp[dir].array();
        }
        const Box& bxg2 = amrex::grow(bx,2);
        FArrayBox qtmp(bxg2, NPRIM, The_Async_Arena());
        auto const& q = qtmp.array();
        auto const& qc = qtmp.const_array();

        ParallelFor(bxg2, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cns_ctoprim(i, j, k, s, q, lparm);
        });

        FArrayBox slopetmp(amrex::grow(bx,1), NEQNS, The_Async_Arena());
        auto const& slope = slopetmp.array();
        auto const& slopec = slopetmp.const_array();

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            ParallelFor(amrex::grow(bx,dir,1), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cns_slope(dir, i, j, k, slope, qc);
            });
            auto const& f = flx[dir];
            ParallelFor(amrex::surroundingNodes(bx,dir),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                cns_flux(dir, i, j, k, f, slopec, qc, lparm);
            });
        }

        // Flux divergence
        ParallelFor(bx, NCONS, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            dsdt(i,j,k,n) = AMREX_D_TERM(
                  dxinv[0] * (flx[0](i,j,k,n) - flx[0](i+1,j,k,n)),
                + dxinv[1] * (flx[1](i,j,k,n) - flx[1](i,j+1,k,n)),
                + dxinv[2] * (flx[2](i,j,k,n) - flx[2](i,j,k+1,n)));
        });

        // Save the fluxes for refluxing. Each face is owned by one tile.
        if (fr_as_crse || fr_as_fine) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                auto const& fout = fluxes[dir].array(mfi);
                auto const& fin = flx[dir];
                ParallelFor(mfi.nodaltilebox(dir), NCONS,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    fout(i,j,k,n) = fin(i,j,k,n);
                });
            }
        }

        // Gravity acts in the last direction.
        if (lparm.gravity != Real(0.0)) {
            const Real g = lparm.gravity;
            constexpr int umg = UMX + AMREX_SPACEDIM - 1;
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dsdt(i,j,k,umg) += g * s(i,j,k,URHO);
                dsdt(i,j,k,UEDEN) += g * s(i,j,k,umg);
            });
        }
    }

    // Register the fluxes, scaled by dt and face area, for refluxing.
    const Real dv = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const Real area = dv / dx[dir];
        if (fr_as_crse) {
            fr_as_crse->CrseInit(fluxes[dir], dir, 0, 0, NCONS, -dt*area, FluxRegister::ADD);
        }
        if (fr_as_fine) {
            fr_as_fine->FineAdd(fluxes[dir], dir, 0, 0, NCONS, dt*area);
        }
    }
}
