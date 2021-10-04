#ifndef PARTICLE_PUSHER_H_
#define PARTICLE_PUSHER_H_

#include "Electrostatic_PIC_K.H"
#include "field_solver/FieldSolver.H"

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MacBndry.H>

using namespace amrex;

namespace FieldSolver {

void fixRHSForSolve (Vector<std::unique_ptr<MultiFab> >& rhs,
                     const Vector<const iMultiFab*>& masks,
                     const Vector<Geometry>& geom, const IntVect& ratio) {
    int num_levels = rhs.size();
    for (int lev = 0; lev < num_levels; ++lev) {
        MultiFab& fine_rhs = *rhs[lev];
        const FabArray<BaseFab<int> >& mask = *masks[lev];
        const BoxArray& fine_ba = fine_rhs.boxArray();
        const DistributionMapping& fine_dm = fine_rhs.DistributionMap();
        MultiFab fine_bndry_data(fine_ba, fine_dm, 1, 1);
        zeroOutBoundary(fine_rhs, fine_bndry_data, mask);
    }
}

void computePhi (ScalarMeshData& rhs, ScalarMeshData& phi,
                 Vector<BoxArray>& grids,
                 Vector<DistributionMapping>& dm,
                 Vector<Geometry>& geom,
                 Vector<const iMultiFab*>& masks) {

    int num_levels = rhs.size();

    Vector<std::unique_ptr<MultiFab> > tmp_rhs(num_levels);
    for (int lev = 0; lev < num_levels; ++lev) {
        tmp_rhs[lev].reset(new MultiFab(rhs[lev]->boxArray(), dm[lev], 1, 0));
        MultiFab::Copy(*tmp_rhs[lev], *rhs[lev], 0, 0, 1, 0);
    }

    IntVect ratio(D_DECL(2, 2, 2));
    fixRHSForSolve(tmp_rhs, masks, geom, ratio);

    int verbose = 2;
    Real rel_tol = 1.0e-12;
    Real abs_tol = 1.0e-14;

    Vector<Geometry>            level_geom(1);
    Vector<BoxArray>            level_grids(1);
    Vector<DistributionMapping> level_dm(1);
    Vector<MultiFab*>           level_phi(1);
    Vector<const MultiFab*>     level_rhs(1);

    for (int lev = 0; lev < num_levels; ++lev) {
        level_phi[0]   = phi[lev].get();
        level_rhs[0]   = tmp_rhs[lev].get();
        level_geom[0]  = geom[lev];
        level_grids[0] = grids[lev];
        level_dm[0]    = dm[lev];

        MLNodeLaplacian linop(level_geom, level_grids, level_dm);

        linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet)},
            {AMREX_D_DECL(LinOpBCType::Dirichlet,
                          LinOpBCType::Dirichlet,
                          LinOpBCType::Dirichlet)});

        linop.setLevelBC(0, nullptr);

        MultiFab sigma(level_grids[0], level_dm[0], 1, 0);
        sigma.setVal(1.0);
        linop.setSigma(0, sigma);

        MLMG mlmg(linop);
        mlmg.setMaxIter(100);
        mlmg.setMaxFmgIter(0);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(0);

        mlmg.solve(level_phi, level_rhs, rel_tol, abs_tol);

        if (lev < num_levels-1) {

            PhysBCFunctNoOp cphysbc, fphysbc;

            int lo_bc[] = {INT_DIR, INT_DIR};
            int hi_bc[] = {INT_DIR, INT_DIR};

            Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
            NodeBilinear mapper;

            amrex::InterpFromCoarseLevel(*phi[lev+1], 0.0, *phi[lev],
                                         0, 0, 1, geom[lev], geom[lev+1],
                                         cphysbc, 0, fphysbc, 0,
                                         IntVect(D_DECL(2, 2, 2)), &mapper, bcs, 0);
        }
    }

    for (int lev = 0; lev < num_levels; ++lev) {
        const Geometry& gm = geom[lev];
        phi[lev]->FillBoundary(gm.periodicity());
    }
}

void sumFineToCrseNodal (const MultiFab& fine, MultiFab& crse,
                         const Geometry& cgeom, const IntVect& ratio) {

    const BoxArray& fine_BA = fine.boxArray();
    const DistributionMapping& fine_dm = fine.DistributionMap();
    BoxArray coarsened_fine_BA = fine_BA;
    coarsened_fine_BA.coarsen(ratio);

    MultiFab coarsened_fine_data(coarsened_fine_BA, fine_dm, 1, 0);
    coarsened_fine_data.setVal(0.0);

    for (MFIter mfi(coarsened_fine_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto crse_data = coarsened_fine_data[mfi].array();
        const auto fine_data = fine[mfi].array();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                                   sum_fine_to_crse_nodal(i, j, k, crse_data, fine_data, ratio);
                               });
    }

    crse.ParallelCopy(coarsened_fine_data, cgeom.periodicity(), FabArrayBase::ADD);
}

void zeroOutBoundary (MultiFab& input_data,
                      MultiFab& bndry_data,
                      const iMultiFab& mask) {
    bndry_data.setVal(0.0, 1);
    for (MFIter mfi(input_data); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto input_arr = input_data[mfi].array();
        auto bndry_arr = bndry_data[mfi].array();
        auto mask_arr = mask[mfi].array();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                                   zero_out_bndry(i, j, k, input_arr, bndry_arr, mask_arr);
                               });
    }
    bndry_data.FillBoundary();
}

void getLevelMasks (Vector<iMultiFab*>& masks,
                    const Vector<BoxArray>& grids,
                    const Vector<DistributionMapping>& dmap,
                    const Vector<Geometry>& geom,
                    const int ncells) {

    int num_levels = grids.size();
    BL_ASSERT(num_levels == dmap.size());

    int covered = 0;
    int notcovered = 1;
    int physbnd = 1;
    int interior = 0;

    for (int lev = 0; lev < num_levels; ++lev) {
        BoxArray nba = grids[lev];
        nba.surroundingNodes();

        FabArray<BaseFab<int> > tmp_mask(nba, dmap[lev], 1, ncells);
        tmp_mask.BuildMask(geom[lev].Domain(), geom[lev].periodicity(),
                           covered, notcovered, physbnd, interior);
        masks[lev].reset(new FabArray<BaseFab<int> >(nba, dmap[lev], 1, 0));
        for (MFIter mfi(tmp_mask); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            const auto tmp_arr = tmp_mask[mfi].array();
            auto mask_arr = (*masks[lev])[mfi].array();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                                       build_mask(i, j, k, tmp_arr, mask_arr, ncells);
                                   });
        }
    }
}

void computeE (      VectorMeshData& E,
               const ScalarMeshData& phi,
               const Vector<Geometry>& geom) {

    const int num_levels = E.size();

    for (int lev = 0; lev < num_levels; ++lev) {
        const auto& gm = geom[lev];
        const auto dx = gm.CellSizeArray();
        for (MFIter mfi(*phi[lev]); mfi.isValid(); ++mfi) {
            Box bx = mfi.validbox();
            bx.grow(1);
            auto Ex_arr = (*E[lev][0])[mfi].array();
            auto Ey_arr = (*E[lev][1])[mfi].array();

            const auto phi_arr = (*phi[lev])[mfi].array();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                                       compute_E_nodal(i, j, k,
                                                       Ex_arr, Ey_arr,
                                                       phi_arr, dx);
                                   });
        }

        E[lev][0]->FillBoundary(gm.periodicity());
        E[lev][1]->FillBoundary(gm.periodicity());
    }
}

};

#endif
