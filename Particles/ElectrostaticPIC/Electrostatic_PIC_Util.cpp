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
