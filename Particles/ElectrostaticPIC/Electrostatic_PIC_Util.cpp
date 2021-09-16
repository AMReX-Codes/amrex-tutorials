#include "Electrostatic_PIC_Util.H"
#include "Electrostatic_PIC_2D.H"

using namespace amrex;

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
                      const FabArray<BaseFab<int> >& mask) {
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

void getLevelMasks (Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks,
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
#if AMREX_SPACEDIM == 3
            auto Ez_arr = (*E[lev][2])[mfi].array();
#endif
            const auto phi_arr = (*phi[lev])[mfi].array();
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                                       compute_E_nodal(i, j, k, phi_arr,
                                                       Ex_arr, Ey_arr,
#if AMREX_SPACEDIM == 3
                                                       Ez_arr,
#endif
                                                       dx);
                                   });
        }

        E[lev][0]->FillBoundary(gm.periodicity());
        E[lev][1]->FillBoundary(gm.periodicity());
#if AMREX_SPACEDIM == 3
        E[lev][2]->FillBoundary(gm.periodicity());
#endif
        //        VisMF::Write(*E[lev][0], amrex::MultiFabFileFullPrefix(lev, "tmp", "Level_", "Ex"));
    }
}
