#include <AMReX.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using amrex::Real;
using namespace amrex::literals;

AMREX_ENUM(Strategy, knapsack, sfc);

// The work per cell varies: cells inside a ball (the "blob") iterate
// work_ratio times, the others once. Boxes in the blob cost up to
// work_ratio times more.
struct Work
{
    amrex::GpuArray<Real,AMREX_SPACEDIM> dx;
    Real radius;
    int work_ratio;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void operator() (int i, int j, int k, amrex::Array4<Real const> const& in,
                     amrex::Array4<Real> const& out) const
    {
        amrex::IntVect const iv(AMREX_D_DECL(i,j,k));
        Real r2 = 0.0_rt;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real const x = (iv[d] + 0.5_rt) * dx[d] - 0.25_rt;
            r2 += x*x;
        }
        int const niter = (r2 < radius*radius) ? work_ratio : 1;
        Real v = in(i,j,k);
        for (int it = 0; it < niter; ++it) {
            v = 4.0_rt * v * (1.0_rt - v); // logistic map, stays in [0,1]
        }
        out(i,j,k) = v;
    }
};

// Run the kernel on every box and add the time each box took to its cost.
void run_and_time (amrex::MultiFab const& in, amrex::MultiFab& out,
                   amrex::LayoutData<Real>& cost, Work const& work)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(in); mfi.isValid(); ++mfi) { // no tiling: one timer per box
        amrex::Box const& bx = mfi.validbox();
        auto const& a = in.const_array(mfi);
        auto const& b = out.array(mfi);
        Real const t0 = amrex::second();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            work(i, j, k, a, b);
        });
        amrex::Gpu::streamSynchronize();
        cost[mfi] += amrex::second() - t0;
    }
}

// The load imbalance: the largest cost of a rank over the average.
Real imbalance (amrex::LayoutData<Real> const& cost)
{
    Real local = 0.0_rt;
    for (amrex::MFIter mfi(cost); mfi.isValid(); ++mfi) { local += cost[mfi]; }
    Real cmax = local;
    Real csum = local;
    amrex::ParallelAllReduce::Max(cmax, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(csum, amrex::ParallelDescriptor::Communicator());
    return cmax / (csum / amrex::ParallelDescriptor::NProcs());
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 64;
        int max_grid_size = 8;
        int work_ratio = 50;
        int nrepeat = 5;
        auto strategy = Strategy::knapsack;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("work_ratio", work_ratio);
            pp.query("nrepeat", nrepeat);
            pp.query_enum_case_insensitive("strategy", strategy);
        }

        amrex::Box const domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        amrex::Geometry geom(domain,
                             amrex::RealBox({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                            {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}),
                             amrex::CoordSys::cartesian, {AMREX_D_DECL(0, 0, 0)});
        Work const work{geom.CellSizeArray(), 0.2_rt, work_ratio};

        // The default mapping orders the boxes along a space-filling curve
        // and gives every rank about the same number of cells.
        amrex::DistributionMapping dm(ba);

        amrex::MultiFab phi(ba, dm, 1, 0);
        auto const& pa = phi.arrays();
        amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            pa[b](i,j,k) = 0.1_rt + 0.01_rt * Real((i + j + k) % 7);
        });
        amrex::Gpu::streamSynchronize(); // finish before the timing starts
        amrex::MultiFab out(ba, dm, 1, 0);

        // 1. Measure the cost of every box with the default mapping.
        amrex::LayoutData<Real> cost(ba, dm);
        for (amrex::MFIter mfi(cost); mfi.isValid(); ++mfi) { cost[mfi] = 0.0_rt; }
        for (int irep = 0; irep < nrepeat; ++irep) {
            run_and_time(phi, out, cost, work);
        }
        Real const imb_old = imbalance(cost);
        amrex::Print() << ba.size() << " boxes on " << amrex::ParallelDescriptor::NProcs()
                       << " ranks\nDefault mapping: max/avg cost per rank = " << imb_old << "\n";

        // 2. A new mapping from the measured costs. The efficiencies (average
        //    over the maximum cost per rank) are computed on the I/O rank.
        Real eff_cur = 0.0_rt;
        Real eff_new = 0.0_rt;
        amrex::DistributionMapping dm_new;
        if (strategy == Strategy::knapsack) {
            dm_new = amrex::DistributionMapping::makeKnapSack(cost, eff_cur, eff_new);
        } else {
            dm_new = amrex::DistributionMapping::makeSFC(cost, eff_cur, eff_new);
        }
        amrex::Print() << "Efficiency of the default mapping " << eff_cur << ", of the "
                       << amrex::getEnumNameString(strategy) << " mapping " << eff_new << "\n";

        // 3. Move the data to the new mapping and measure again.
        amrex::MultiFab phi_new(ba, dm_new, 1, 0);
        phi_new.ParallelCopy(phi);
        amrex::MultiFab out_new(ba, dm_new, 1, 0);
        amrex::LayoutData<Real> cost_new(ba, dm_new);
        for (amrex::MFIter mfi(cost_new); mfi.isValid(); ++mfi) { cost_new[mfi] = 0.0_rt; }
        for (int irep = 0; irep < nrepeat; ++irep) {
            run_and_time(phi_new, out_new, cost_new, work);
        }
        Real const imb_new = imbalance(cost_new);
        amrex::Print() << "New mapping: max/avg cost per rank = " << imb_new << "\n";

        // The result does not depend on the mapping.
        amrex::MultiFab check(ba, dm, 1, 0);
        check.ParallelCopy(out_new);
        amrex::MultiFab::Subtract(check, out, 0, 0, 1, 0);
        Real const diff = check.norm0();
        amrex::Print() << "Max difference between the results: " << diff << "\n";

        // The measured times are not part of the verdict: they are noisy.
        // The efficiencies come from the same costs by different arithmetic.
        // The heuristics do not promise a better mapping; allow a small loss.
        if (diff == 0.0_rt && eff_new >= 0.99_rt * eff_cur) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
