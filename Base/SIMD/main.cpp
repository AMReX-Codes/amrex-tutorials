#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuLaunchFunctsSIMD.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#include <AMReX_SIMD.H>

#include <cmath>
#include <limits>
#include <type_traits>

using amrex::Real;
using namespace amrex::literals;

// Part 1: the same function, the signed square root of p - c, written in
// two ways.
//
// For illustration only: a call to a function the compiler cannot inline,
// with a branch inside, stops the vectorization of the loop that calls it.
AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
Real f_branch (Real p, Real c)
{
    if (p > c) {
        return std::sqrt(p - c);
    } else {
        return -std::sqrt(c - p);
    }
}

// Inlined and without branches, the loop vectorizes. The arithmetic is done
// unconditionally, on an argument made safe for every lane (abs), and ?:
// only picks one of the computed values.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real f_select (Real p, Real c)
{
    Real const r = std::sqrt(std::abs(p - c));
    return (p > c) ? r : -r;
}

// Part 2: a kernel for ParallelForSIMD. T is amrex::Real, or a SIMD
// register of several Reals when AMReX is built with SIMD support.
struct Push : amrex::simd::Vectorized<>
{
    Real dt;

    explicit Push (Real a_dt) : dt(a_dt) {}

    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (T& AMREX_RESTRICT x, T const& AMREX_RESTRICT v) const
    {
        x += v * dt;
    }
};

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 64;
        amrex::Long n = 1000000;
        int nrepeat = 20;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("n", n);
            pp.query("nrepeat", nrepeat);
        }

        // Part 1: implicit vectorization of ParallelFor kernels. On CPU,
        // ParallelFor marks its innermost loop with AMREX_PRAGMA_SIMD; whether
        // the compiler vectorizes it depends on the kernel.
        amrex::Box const domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(32);
        amrex::DistributionMapping dm(ba);
        amrex::MultiFab phi(ba, dm, 1, 0);
        auto const& pa = phi.arrays();
        Real const dx = 1.0_rt / Real(n_cell);
        amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            amrex::IntVect const iv(AMREX_D_DECL(i,j,k));
            Real r = 1.0_rt;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                r *= std::sin(2.0_rt * amrex::Math::pi<Real>() * (iv[d]+0.5_rt)*dx);
            }
            pa[b](i,j,k) = 0.5_rt * (1.0_rt + r);
        });
        Real const c = 0.5_rt;

        amrex::MultiFab out_branch(ba, dm, 1, 0);
        amrex::MultiFab out_select(ba, dm, 1, 0);
        auto const& ca = phi.const_arrays();
        auto const& oba = out_branch.arrays();
        auto const& osa = out_select.arrays();
        for (int irep = 0; irep < nrepeat; ++irep) {
            BL_PROFILE("branch_not_inlined");
            amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                oba[b](i,j,k) = f_branch(ca[b](i,j,k), c);
            });
            amrex::Gpu::streamSynchronize();
        }
        for (int irep = 0; irep < nrepeat; ++irep) {
            BL_PROFILE("select_inlined");
            amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                osa[b](i,j,k) = f_select(ca[b](i,j,k), c);
            });
            amrex::Gpu::streamSynchronize();
        }
        amrex::MultiFab::Subtract(out_select, out_branch, 0, 0, 1, 0);
        Real const diff1 = out_select.norm0();
        amrex::Print() << "Part 1: max difference between the two kernels = " << diff1 << "\n";

        // Part 2: explicit SIMD. ParallelForSIMD<Push> hands the kernel a
        // SIMDindex that covers simd_width consecutive elements; load_1d
        // returns a SIMD register for it, store_1d writes the register back.
        // Without AMReX_SIMD, or on GPU, it is a plain ParallelFor and the
        // same code runs one element at a time.
        amrex::Print() << "Part 2: SIMD width " << amrex::simd::native_simd_size_real
#ifdef AMREX_USE_SIMD
                       << " (AMReX_SIMD=ON)\n";
#else
                       << " (AMReX_SIMD=OFF)\n";
#endif
        amrex::Gpu::DeviceVector<Real> x(n), x_ref(n), v(n);
        Real* xp = x.data();
        Real* xrp = x_ref.data();
        Real* vp = v.data();
        amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (amrex::Long i)
        {
            Real const xi = Real(i) / Real(n);
            xp[i] = xi;
            xrp[i] = xi;
            vp[i] = std::sin(2.0_rt * amrex::Math::pi<Real>() * xi);
        });
        amrex::Gpu::streamSynchronize();
        Push const push(0.1_rt / Real(nrepeat));

        for (int irep = 0; irep < nrepeat; ++irep) {
            BL_PROFILE("push_simd");
            amrex::ParallelForSIMD<Push>(n, [=] AMREX_GPU_DEVICE (auto si)
            {
                decltype(auto) xi = amrex::simd::load_1d(xp, si);
                decltype(auto) vi = amrex::simd::load_1d(vp, si);
                using T = std::remove_reference_t<decltype(xi)>;
                constexpr auto method = &Push::operator()<T>;
                push(xi, vi);
                // Writes back only for arguments the method takes by
                // non-const reference, here x; a no-op without SIMD.
                amrex::simd::store_1d<method, 0>(xi, xp, si);
            });
            amrex::Gpu::streamSynchronize();
        }
        for (int irep = 0; irep < nrepeat; ++irep) {
            BL_PROFILE("push_scalar");
            amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (amrex::Long i)
            {
                push(xrp[i], vp[i]);
            });
            amrex::Gpu::streamSynchronize();
        }

        Real const diff2 = amrex::Reduce::Max<Real>(n, [=] AMREX_GPU_DEVICE (amrex::Long i)
        {
            return std::abs(xp[i] - xrp[i]);
        });
        amrex::Print() << "Part 2: max difference between the SIMD and the scalar push = "
                       << diff2 << "\n";

        // Both pairs of kernels round the same way, up to fused multiply-adds
        // or a fast-math sqrt the compiler may use in one version only.
        Real const eps = std::numeric_limits<Real>::epsilon();
        if (diff1 <= 4.0_rt * eps && diff2 <= Real(nrepeat) * 4.0_rt * eps) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
