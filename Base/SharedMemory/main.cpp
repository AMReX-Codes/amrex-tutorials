#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_OpenMP.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdint>

using amrex::Long;
using amrex::Real;
using namespace amrex::literals;

// Bin the cells by the value of phi in [0,1] and sum the signed integer
// field q over each bin, with global Long atomics and with private int
// sums in shared memory. Atomic adds to the same global address are
// serialized, so the global version slows down as the bins get fewer.

constexpr int max_q = 1000; // |q| <= max_q

// The bin of a value in [0,1].
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int bin_of (Real v, int nbins)
{
    return amrex::min(amrex::max(int(v * Real(nbins)), 0), nbins-1);
}

// 1. Portable: one atomic add into global memory per cell. On CPU, the
//    atomics are OpenMP atomics.
void binned_sum_atomic (amrex::MultiFab const& phi, amrex::iMultiFab const& q,
                        int nbins, Long* sums)
{
    auto const& pa = phi.const_arrays();
    auto const& qa = q.const_arrays();
    amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        amrex::HostDevice::Atomic::Add(sums + bin_of(pa[b](i,j,k), nbins), Long(qa[b](i,j,k)));
    });
    amrex::Gpu::streamSynchronize();
}

// 2. Each thread block (or each OpenMP thread on CPU) accumulates private
//    sums and adds them to the global ones at the end.
void binned_sum_private (amrex::MultiFab const& phi, amrex::iMultiFab const& q,
                         int nbins, Long* sums)
{
#ifdef AMREX_USE_GPU
    {
        constexpr int nthreads = 256;   // threads per block
        constexpr int max_nbins = 1024; // size of the static shared-memory array
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nbins <= max_nbins, "nbins must not exceed max_nbins");
        for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi) {
            amrex::Box const& bx = mfi.validbox();
            auto const& pa = phi.const_array(mfi);
            auto const& qa = q.const_array(mfi);
            amrex::BoxIndexer const indexer(bx); // linear index -> (i,j,k)
            auto const ncells = std::uint64_t(bx.numPts());
            auto const nblocks = int(std::min<std::uint64_t>(
                (ncells + nthreads - 1) / nthreads,
                amrex::Gpu::Device::maxBlocksPerLaunch()));
            // The private sums are int: with |q| <= max_q = 1000, a block
            // would need over two million cells to overflow.
#ifdef AMREX_USE_SYCL
            // SYCL: dynamic shared memory, requested at launch.
            std::size_t const shared_mem_bytes = nbins * sizeof(int);
            amrex::launch<nthreads>(nblocks, shared_mem_bytes, amrex::Gpu::gpuStream(),
            [=] AMREX_GPU_DEVICE (amrex::Gpu::Handler const& h) noexcept
            {
                auto* sh = static_cast<int*>(h.sharedMemory());
                int const tid = int(h.threadIdx());
                int const bdim = int(h.blockDim());
                for (int b = tid; b < nbins; b += bdim) { sh[b] = 0; }
                h.syncThreads();
                auto const stride = std::uint64_t(bdim) * h.gridDim();
                for (auto icell = std::uint64_t(h.globalIdx()); icell < ncells; icell += stride) {
                    auto const [i, j, k] = indexer(icell);
                    amrex::Gpu::Atomic::AddNoRet<int, sycl::access::address_space::local_space>
                        (sh + bin_of(pa(i,j,k), nbins), qa(i,j,k));
                }
                h.syncThreads();
                for (int b = tid; b < nbins; b += bdim) {
                    if (sh[b] != 0) { amrex::Gpu::Atomic::AddNoRet(sums + b, Long(sh[b])); }
                }
            });
#else
            // CUDA and HIP: both static and dynamic shared memory are
            // supported; here a static array.
            amrex::launch<nthreads>(nblocks, amrex::Gpu::gpuStream(),
            [=] AMREX_GPU_DEVICE () noexcept
            {
                __shared__ int sh[max_nbins];
                int const tid = int(threadIdx.x);
                int const bdim = int(blockDim.x);
                for (int b = tid; b < nbins; b += bdim) { sh[b] = 0; }
                __syncthreads();
                auto const stride = std::uint64_t(bdim) * gridDim.x;
                for (auto icell = std::uint64_t(bdim) * blockIdx.x + tid; icell < ncells; icell += stride) {
                    auto const [i, j, k] = indexer(icell);
                    amrex::Gpu::Atomic::AddNoRet(sh + bin_of(pa(i,j,k), nbins), qa(i,j,k));
                }
                __syncthreads();
                for (int b = tid; b < nbins; b += bdim) {
                    if (sh[b] != 0) { amrex::Gpu::Atomic::AddNoRet(sums + b, Long(sh[b])); }
                }
            });
#endif
        }
    }
#else
    {
        // Private sums per OpenMP thread, merged after the parallel region.
        int const nthreads = amrex::OpenMP::get_max_threads();
        amrex::Vector<amrex::Vector<Long>> local(nthreads, amrex::Vector<Long>(nbins, 0));
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        {
            auto& mine = local[amrex::OpenMP::get_thread_num()];
            for (amrex::MFIter mfi(phi, true); mfi.isValid(); ++mfi) {
                auto const& pa = phi.const_array(mfi);
                auto const& qa = q.const_array(mfi);
                amrex::Loop(mfi.tilebox(), [&] (int i, int j, int k)
                {
                    mine[bin_of(pa(i,j,k), nbins)] += qa(i,j,k);
                });
            }
        }
        for (auto const& l : local) {
            for (int b = 0; b < nbins; ++b) { sums[b] += l[b]; }
        } // end of merge
    }
#endif
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 64;
        int max_grid_size = 32;
        int nbins = 64;
        int nrepeat = 10;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nbins", nbins);
            pp.query("nrepeat", nrepeat);
        }

        amrex::Box const domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        amrex::DistributionMapping dm(ba);
        amrex::Geometry geom(domain,
                             amrex::RealBox({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                            {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}),
                             amrex::CoordSys::cartesian, {AMREX_D_DECL(1, 1, 1)});
        auto const dx = geom.CellSizeArray();

        // phi in [0,1] selects the bin; q, a signed integer with |q| <= max_q
        // and a nonzero sum over the domain, is summed.
        amrex::MultiFab phi(ba, dm, 1, 0);
        amrex::iMultiFab q(ba, dm, 1, 0);
        auto const& pa = phi.arrays();
        auto const& qa = q.arrays();
        amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            amrex::IntVect const iv(AMREX_D_DECL(i,j,k));
            Real s = 1.0_rt;
            Real c = 1.0_rt;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Real const x = 2.0_rt * amrex::Math::pi<Real>() * (iv[d]+0.5_rt)*dx[d];
                s *= std::sin(x);
                c *= std::cos(2.0_rt*x);
            }
            pa[b](i,j,k) = 0.5_rt * (1.0_rt + s);
            qa[b](i,j,k) = int(std::round(Real(max_q) * (c + 0.5_rt) / 1.5_rt));
        });

        // The sums accumulate over the nrepeat passes.
        amrex::Gpu::DeviceVector<Long> sums_atomic(nbins, 0);
        amrex::Gpu::DeviceVector<Long> sums_private(nbins, 0);

        for (int irep = 0; irep < nrepeat; ++irep) {
            BL_PROFILE("binned_sum_atomic");
            binned_sum_atomic(phi, q, nbins, sums_atomic.data());
        }
        for (int irep = 0; irep < nrepeat; ++irep) {
            BL_PROFILE("binned_sum_private");
            binned_sum_private(phi, q, nbins, sums_private.data());
        }

        // To the host, then summed over MPI ranks.
        amrex::Vector<Long> sa(nbins), sp(nbins);
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, sums_atomic.begin(), sums_atomic.end(), sa.begin());
        amrex::Gpu::copy(amrex::Gpu::deviceToHost, sums_private.begin(), sums_private.end(), sp.begin());
        amrex::ParallelDescriptor::ReduceLongSum(sa.data(), nbins);
        amrex::ParallelDescriptor::ReduceLongSum(sp.data(), nbins);

        // The total over all bins must equal nrepeat times the sum of q.
        Long const expected = Long(nrepeat) * q.sum(0);

        Long total = 0;
        bool same = true;
        for (int b = 0; b < nbins; ++b) {
            total += sp[b];
            same = same && (sa[b] == sp[b]);
        }
        amrex::Print() << "bin  sum\n";
        for (int b = 0; b < nbins; b += std::max(nbins/8, 1)) {
            amrex::Print() << b << "  " << sp[b] << "\n";
        }
        amrex::Print() << "The two results are " << (same ? "identical" : "different")
                       << "; total " << total << " (expected " << expected << ")\n";
        if (same && total == expected) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
