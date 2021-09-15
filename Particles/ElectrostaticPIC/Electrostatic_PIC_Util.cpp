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
