#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_Scan.H>

#include <algorithm>
#include <iterator>
#include <numeric>
#include <vector>

using amrex::Real;
using namespace amrex::literals;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n = 100000;
        Real threshold = 0.75_rt;
        {
            amrex::ParmParse pp;
            pp.query("n", n);
            pp.query("threshold", threshold);
        }

        // Gpu::DeviceVector lives in device memory (from The_Arena()) on GPU
        // and in host memory on CPU. Its elements must not be accessed on the
        // host in GPU builds; use its data() pointer in kernels.
        amrex::Gpu::DeviceVector<Real> values(n);
        Real* pv = values.data();

        // Uniform random numbers in [0,1); each GPU thread gets its own engine.
        amrex::ParallelForRNG(n, [=] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine)
        {
            pv[i] = amrex::Random(engine);
        });

        // Stream compaction: keep the values above threshold, in order. The
        // exclusive prefix sum of the keep flags is each kept value's output
        // position, and the total is the number of kept values.
        amrex::Gpu::DeviceVector<Real> kept(n);
        Real* pk = kept.data();
        int const nkept = amrex::Scan::PrefixSum<int>(n,
            [=] AMREX_GPU_DEVICE (int i) -> int { return pv[i] > threshold; },
            [=] AMREX_GPU_DEVICE (int i, int const& pos) {
                if (pv[i] > threshold) { pk[pos] = pv[i]; }
            },
            amrex::Scan::Type::exclusive, amrex::Scan::retSum);
        kept.resize(nkept);
        amrex::Print() << "Kept " << nkept << " of " << n << " values above " << threshold << "\n";

        // Copy to the host. Pinned (page-locked) host memory, from
        // The_Pinned_Arena(), makes asynchronous copies fast. The copy is
        // asynchronous on GPU, so wait for the stream before using the data.
        amrex::Gpu::PinnedVector<Real> h_values(n);
        amrex::Gpu::PinnedVector<Real> h_kept(nkept);
        amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, values.begin(), values.end(), h_values.begin());
        amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, kept.begin(), kept.end(), h_kept.begin());
        amrex::Gpu::streamSynchronize();

        // Reference compaction on the host
        std::vector<Real> ref;
        std::copy_if(h_values.begin(), h_values.end(), std::back_inserter(ref),
                     [=] (Real v) { return v > threshold; });
        bool const compaction_ok = std::equal(ref.begin(), ref.end(), h_kept.begin(), h_kept.end());
        amrex::Print() << "Compaction matches the host result: "
                       << (compaction_ok ? "yes" : "no") << "\n";

        // Scan::ExclusiveSum on raw device memory. The_Async_Arena() is for
        // temporaries: on GPU, free() is ordered on the stream, so there is
        // no need to synchronize before freeing.
        amrex::Gpu::HostVector<int> h_counts(n);
        for (int i = 0; i < n; ++i) { h_counts[i] = i % 4; }
        amrex::Gpu::DeviceVector<int> counts(n);
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, h_counts.begin(), h_counts.end(), counts.begin());

        auto* offsets = static_cast<int*>(amrex::The_Async_Arena()->alloc(sizeof(int)*n));
        int const total = amrex::Scan::ExclusiveSum(n, counts.data(), offsets);
        amrex::Gpu::PinnedVector<int> h_offsets(n);
        amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, offsets, offsets+n, h_offsets.begin());
        amrex::Gpu::streamSynchronize();
        amrex::The_Async_Arena()->free(offsets);

        std::vector<int> ref_offsets(n);
        std::exclusive_scan(h_counts.begin(), h_counts.end(), ref_offsets.begin(), 0);
        bool const scan_ok = std::equal(ref_offsets.begin(), ref_offsets.end(), h_offsets.begin())
            && total == std::accumulate(h_counts.begin(), h_counts.end(), 0);
        amrex::Print() << "ExclusiveSum total = " << total << ", matches the host result: "
                       << (scan_ok ? "yes" : "no") << "\n";

        if (compaction_ok && scan_ok) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
