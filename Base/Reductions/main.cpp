#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <cmath>

using amrex::Long;
using amrex::Real;
using namespace amrex::literals;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int n_vector = 1000;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("n_vector", n_vector);
        }

        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        amrex::DistributionMapping dm(ba);
        amrex::MultiFab mf(ba, dm, 1, 0);

        // phi = i + j + k, so that the sum and the max have simple exact values.
        auto const& ma = mf.arrays();
        amrex::ParallelFor(mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            ma[b](i,j,k) = Real(AMREX_D_TERM(i, + j, + k));
        });

        // With n = n_cell and d = AMREX_SPACEDIM: sum = d n^(d-1) sum_i i.
        double const n = n_cell;
        double const d = AMREX_SPACEDIM;
        double const exact_sum = d * std::pow(n, d-1) * n * (n-1) / 2;
        double const exact_max = d * (n-1);

        Real sum_parreduce = 0; // kept for the check at the end

        // 1. ParReduce: several reductions over a MultiFab in one pass.
        //    The result is local to this rank; MPI reduction is up to us.
        {
            auto const& a = mf.const_arrays();
            auto r = amrex::ParReduce(
                amrex::TypeList<amrex::ReduceOpSum, amrex::ReduceOpMin, amrex::ReduceOpMax>{},
                amrex::TypeList<Real, Real, Real>{},
                mf, amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                    -> amrex::GpuTuple<Real, Real, Real>
                {
                    Real v = a[b](i,j,k);
                    return {v, v, v};
                });
            Real sum = amrex::get<0>(r);
            Real vmin = amrex::get<1>(r);
            Real vmax = amrex::get<2>(r);
            amrex::ParallelAllReduce::Sum(sum, amrex::ParallelDescriptor::Communicator());
            amrex::ParallelAllReduce::Min(vmin, amrex::ParallelDescriptor::Communicator());
            amrex::ParallelAllReduce::Max(vmax, amrex::ParallelDescriptor::Communicator());
            amrex::Print() << "ParReduce:\n"
                           << "  sum = " << sum << " (exact " << exact_sum << ")\n"
                           << "  min = " << vmin << " (exact 0)\n"
                           << "  max = " << vmax << " (exact " << exact_max << ")\n";
            sum_parreduce = sum;
        }

        // 2. ReduceOps + ReduceData in an MFIter loop: one reduction kernel
        //    per box, with a single synchronization in value().
        {
            amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_op;
            amrex::ReduceData<Real, Long> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box const& bx = mfi.tilebox();
                auto const& a = mf.const_array(mfi);
                reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    Real v = a(i,j,k);
                    Long even = (Long(v) % 2 == 0) ? 1 : 0;
                    return {v*v, even};
                });
            }
            ReduceTuple r = reduce_data.value(reduce_op);
            Real sumsq = amrex::get<0>(r);
            Long neven = amrex::get<1>(r);
            amrex::ParallelDescriptor::ReduceRealSum(sumsq);
            amrex::ParallelDescriptor::ReduceLongSum(neven);
            amrex::Print() << "ReduceOps/ReduceData:\n"
                           << "  sum of squares = " << sumsq << "\n"
                           << "  even cells = " << neven << "\n";
        }

        // 3. The same with Reducer, which bundles ReduceOps and ReduceData.
        {
            amrex::Reducer<amrex::TypeList<amrex::ReduceOpSum, amrex::ReduceOpSum>,
                           amrex::TypeList<Real, Long>> reducer;
            using Result_t = decltype(reducer)::Result_t;
            auto const& a = mf.const_arrays();
            reducer.eval(mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) -> Result_t
            {
                Real v = a[b](i,j,k);
                Long even = (Long(v) % 2 == 0) ? 1 : 0;
                return {v*v, even};
            });
            Result_t r = reducer.getResult();
            Real sumsq = amrex::get<0>(r);
            Long neven = amrex::get<1>(r);
            amrex::ParallelDescriptor::ReduceRealSum(sumsq);
            amrex::ParallelDescriptor::ReduceLongSum(neven);
            amrex::Print() << "Reducer:\n"
                           << "  sum of squares = " << sumsq << "\n"
                           << "  even cells = " << neven << "\n";
        }

        // 4. Built-in MultiFab reductions, already reduced over MPI ranks.
        Real sum_mf = mf.sum(0);
        amrex::Print() << "MultiFab member functions:\n"
                       << "  sum() = " << sum_mf << " (exact " << exact_sum << ")\n"
                       << "  min() = " << mf.min(0) << "\n"
                       << "  max() = " << mf.max(0) << " (exact " << exact_max << ")\n"
                       << "  norm0() = " << mf.norm0(0) << "\n"
                       << "  norm1() = " << mf.norm1(0) << "\n"
                       << "  norm2() = " << mf.norm2(0) << "\n";

        // 5. Reduce::Sum over a device vector holding 1, 2, ..., n_vector on
        //    every rank.
        {
            amrex::Gpu::DeviceVector<Real> v(n_vector);
            Real* p = v.data();
            amrex::ParallelFor(n_vector, [=] AMREX_GPU_DEVICE (int i)
            {
                p[i] = Real(i+1);
            });
            Real sum = amrex::Reduce::Sum(n_vector, p);
            amrex::ParallelAllReduce::Sum(sum, amrex::ParallelDescriptor::Communicator());
            double const nv = n_vector;
            amrex::Print() << "Reduce::Sum on a Gpu::DeviceVector:\n"
                           << "  sum = " << sum << " (exact "
                           << amrex::ParallelDescriptor::NProcs() * nv * (nv+1) / 2 << ")\n";
        }

        // Sanity check: the ParReduce and MultiFab sums against the exact value.
        double const tol = 1.e-4 * exact_sum;
        if (std::abs(sum_parreduce - exact_sum) <= tol && std::abs(sum_mf - exact_sum) <= tol) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
