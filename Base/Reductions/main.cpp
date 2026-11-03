#include <AMReX.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParReduce.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <type_traits>

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

        // phi = i + j + k, so that the reductions have simple exact values.
        auto const& ma = mf.arrays();
        amrex::ParallelFor(mf, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            ma[b](i,j,k) = Real(AMREX_D_TERM(i, + j, + k));
        });

        // Exact values. With n = n_cell, d = AMREX_SPACEDIM, s1 = sum_i i and
        // s2 = sum_i i^2 over 0 <= i < n:
        //   sum   = d n^(d-1) s1
        //   sumsq = d n^(d-1) s2 + d (d-1) n^(d-2) s1^2
        //   #even = (n^d + (n%2)^d) / 2   (cells with i+j+k even)
        double const n = n_cell;
        double const d = AMREX_SPACEDIM;
        double const s1 = n * (n-1) / 2;
        double const s2 = (n-1) * n * (2*n-1) / 6;
        double const exact_sum = d * std::pow(n, d-1) * s1;
        double const exact_sumsq = d * std::pow(n, d-1) * s2 + d * (d-1) * std::pow(n, d-2) * s1 * s1;
        double const exact_max = d * (n-1);
        Long const exact_even = (Long(std::pow(n, d)) + ((n_cell % 2) ? 1 : 0)) / 2;

        constexpr double rtol = std::is_same_v<Real,float> ? 1.e-3 : 1.e-12;
        auto check = [&] (char const* name, double value, double exact) {
            bool ok = std::abs(value - exact) <= rtol * std::abs(exact);
            amrex::Print() << "  " << name << " = " << value << (ok ? "" : "  WRONG")
                           << " (exact " << exact << ")\n";
            if (!ok) { status = 1; }
        };

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
            amrex::Print() << "ParReduce:\n";
            check("sum", sum, exact_sum);
            check("min", vmin, 0.0);
            check("max", vmax, exact_max);
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
            amrex::Print() << "ReduceOps/ReduceData:\n";
            check("sum of squares", sumsq, exact_sumsq);
            check("even cells", double(neven), double(exact_even));
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
            amrex::Print() << "Reducer:\n";
            check("sum of squares", sumsq, exact_sumsq);
            check("even cells", double(neven), double(exact_even));
        }

        // 4. Built-in MultiFab reductions, already reduced over MPI ranks.
        amrex::Print() << "MultiFab member functions:\n";
        check("sum()", mf.sum(0), exact_sum);
        check("min()", mf.min(0), 0.0);
        check("max()", mf.max(0), exact_max);
        check("norm0()", mf.norm0(0), exact_max);
        check("norm1()", mf.norm1(0), exact_sum);
        check("norm2()", mf.norm2(0), std::sqrt(exact_sumsq));

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
            amrex::Print() << "Reduce::Sum on a Gpu::DeviceVector:\n";
            double const nv = n_vector;
            check("sum", sum, amrex::ParallelDescriptor::NProcs() * nv * (nv+1) / 2);
        }

        amrex::Print() << (status == 0 ? "PASSED\n" : "FAILED\n");
    }
    amrex::Finalize();
    return status;
}
