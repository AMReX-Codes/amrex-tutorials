#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <limits>

using amrex::Real;
using namespace amrex::literals;

AMREX_ENUM(Order, second, fourth);

// Laplacian of order 2 or 4. The branch is resolved at compile time.
template <int order>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real laplacian (amrex::IntVect const& iv, amrex::Array4<Real const> const& p,
                amrex::GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{
    Real r = 0.0_rt;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        auto const e = amrex::IntVect::TheDimensionVector(d);
        if constexpr (order == 2) {
            r += (p(iv+e) - 2.0_rt*p(iv) + p(iv-e)) * dxinv[d] * dxinv[d];
        } else {
            r += (-p(iv+2*e) + 16.0_rt*p(iv+e) - 30.0_rt*p(iv) + 16.0_rt*p(iv-e) - p(iv-2*e))
                * (dxinv[d] * dxinv[d] / 12.0_rt);
        }
    }
    return r;
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int nrepeat = 5;
        auto order = Order::second;
        bool absorption = false;
        Real sigma = 0.5_rt;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nrepeat", nrepeat);
            pp.query_enum_case_insensitive("order", order);
            pp.query("absorption", absorption);
            pp.query("sigma", sigma);
        }
        amrex::Print() << "order = " << amrex::getEnumNameString(order)
                       << ", absorption = " << absorption << "\n";

        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        amrex::DistributionMapping dm(ba);
        amrex::Geometry geom(domain,
                             amrex::RealBox({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                            {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}),
                             amrex::CoordSys::cartesian, {AMREX_D_DECL(1, 1, 1)});
        auto const prob_lo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();
        auto const dxinv = geom.InvCellSizeArray();

        amrex::MultiFab phi(ba, dm, 1, 2); // 2 ghost cells for the 4th-order stencil
        auto const& pa = phi.arrays();
        amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            amrex::IntVect const iv(AMREX_D_DECL(i,j,k));
            Real r = 1.0_rt;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                r *= std::sin(2.0_rt * amrex::Math::pi<Real>() * (prob_lo[d] + (iv[d]+0.5_rt)*dx[d]));
            }
            pa[b](i,j,k) = r;
        });
        phi.FillBoundary(geom.periodicity());

        amrex::MultiFab out_runtime(ba, dm, 1, 0);
        amrex::MultiFab out_cto(ba, dm, 1, 0);

        // 1. Runtime branches: one kernel tests the options for every cell.
        for (int n = 0; n < nrepeat; ++n) {
            BL_PROFILE("runtime_branches");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(out_runtime, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box const& bx = mfi.tilebox();
                auto const& p = phi.const_array(mfi);
                auto const& out = out_runtime.array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    amrex::IntVect const iv(AMREX_D_DECL(i,j,k));
                    Real r = (order == Order::second) ? laplacian<2>(iv, p, dxinv)
                                                      : laplacian<4>(iv, p, dxinv);
                    if (absorption) { r -= sigma * p(iv); }
                    out(i,j,k) = r;
                });
            }
            amrex::Gpu::streamSynchronize();
        }

        // 2. Compile-time options: a kernel is compiled for each of the 2x2
        //    combinations of the listed values, and the one matching the
        //    runtime values is launched. The extra kernel arguments are
        //    std::integral_constant<int,...>, so the options are constants.
        for (int n = 0; n < nrepeat; ++n) {
            BL_PROFILE("compile_time_options");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(out_cto, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                amrex::Box const& bx = mfi.tilebox();
                auto const& p = phi.const_array(mfi);
                auto const& out = out_cto.array(mfi);
                amrex::ParallelFor(
                    amrex::TypeList<amrex::CompileTimeOptions<int(Order::second), int(Order::fourth)>,
                                    amrex::CompileTimeOptions<0, 1>>{},
                    {int(order), int(absorption)},
                    bx, [=] AMREX_GPU_DEVICE (int i, int j, int k, auto order_ct, auto absorption_ct)
                {
                    amrex::IntVect const iv(AMREX_D_DECL(i,j,k));
                    constexpr int ord = (order_ct.value == int(Order::second)) ? 2 : 4;
                    Real r = laplacian<ord>(iv, p, dxinv);
                    // A plain `if` on a compile-time constant is optimized out
                    // too. (nvcc does not allow `if constexpr` blocks to be
                    // the first to use a captured variable, here sigma.)
                    if (absorption_ct.value) { r -= sigma * p(iv); }
                    out(i,j,k) = r;
                });
            }
            amrex::Gpu::streamSynchronize();
        }

        amrex::MultiFab::Subtract(out_cto, out_runtime, 0, 0, 1, 0);
        Real const diff = out_cto.norm0();
        Real const scale = out_runtime.norm0();
        amrex::Print() << "max |result| = " << scale << ", max diff between the two kernels = "
                       << diff << "\n";

        Real const tol = 10.0_rt * std::numeric_limits<Real>::epsilon() * scale;
        if (diff > tol) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
