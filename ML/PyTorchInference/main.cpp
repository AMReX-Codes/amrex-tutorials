#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <torch/script.h>
#ifdef AMREX_USE_CUDA
#include <c10/cuda/CUDAGuard.h>
#include <c10/cuda/CUDAStream.h>
#elif defined(AMREX_USE_GPU)
#error "This tutorial supports CPU and CUDA builds only"
#endif

#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        std::string model_file = "model.pt";
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("model_file", model_file);
        }

        Box const domain(IntVect(0), IntVect(n_cell-1));
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);
        Geometry geom(domain, RealBox(AMREX_D_DECL(0._rt,0._rt,0._rt),
                                      AMREX_D_DECL(1._rt,1._rt,1._rt)),
                      CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});

        // The model maps nin input fields to nout output fields cell by cell.
        constexpr int nin = 2;
        constexpr int nout = 2;
        MultiFab phi_in(ba, dm, nin, 0);
        MultiFab phi_out(ba, dm, nout, 0);

        // Input fields: sin(2 pi x) and cos(2 pi (x+y+z))
        {
            auto const dx = geom.CellSizeArray();
            auto const& a = phi_in.arrays();
            ParallelFor(phi_in, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                IntVect const iv(AMREX_D_DECL(i,j,k));
                Real s = 0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { s += (iv[d] + 0.5_rt)*dx[d]; }
                Real const x = (i + 0.5_rt)*dx[0];
                a[b](i,j,k,0) = std::sin(2._rt*Math::pi<Real>()*x);
                a[b](i,j,k,1) = std::cos(2._rt*Math::pi<Real>()*s);
            });
        }

        // Tensors use amrex::Real and live where the MultiFab data live.
        constexpr auto dtype = std::is_same_v<Real,double> ? torch::kFloat64 : torch::kFloat32;
#ifdef AMREX_USE_CUDA
        torch::Device const device(torch::kCUDA,
                                   static_cast<c10::DeviceIndex>(Gpu::Device::deviceId()));
#else
        torch::Device const device(torch::kCPU);
#endif
        auto const options = torch::TensorOptions().dtype(dtype).device(device);

        // Load the TorchScript model written by make_model.py.
        torch::jit::script::Module model;
        try {
            model = torch::jit::load(model_file, device);
        } catch (c10::Error const& e) {
            amrex::Abort("Failed to load " + model_file + ": " + e.what());
        }
        model.to(dtype);
        model.eval();
        amrex::Print() << "Loaded " << model_file << " on " << device << "\n";

        // No autograd bookkeeping for inference
        c10::InferenceMode inference_guard;

        // One batch per box, so no tiling
        for (MFIter mfi(phi_in); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.validbox();
            auto const ncells = static_cast<std::int64_t>(bx.numPts());
            auto const& in = phi_in.const_array(mfi);
            auto const& out = phi_out.array(mfi);

#ifdef AMREX_USE_CUDA
            // Run torch on AMReX's current stream so that their kernels are ordered.
            c10::cuda::CUDAStreamGuard stream_guard(
                c10::cuda::getStreamFromExternal(Gpu::gpuStream(), device.index()));
#endif

            // Gather the cells of the box into an (ncells, nin) batch in GPU
            // memory for GPU builds.
            Gpu::DeviceVector<Real> batch(ncells*nin);
            Real* pbatch = batch.data();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                auto const m = bx.index(IntVect(AMREX_D_DECL(i,j,k)));
                for (int n = 0; n < nin; ++n) {
                    pbatch[m*nin+n] = in(i,j,k,n);
                }
            });

            // Wrap the batch in a tensor without copying, and evaluate the model.
            torch::Tensor input = torch::from_blob(pbatch, {ncells, nin}, options);
            torch::Tensor output;
            {
                BL_PROFILE("model.forward");
                output = model.forward({input}).toTensor().contiguous();
            }
            AMREX_ALWAYS_ASSERT(output.dim() == 2 && output.size(0) == ncells &&
                                output.size(1) == nout && output.scalar_type() == dtype);

            // Scatter the result back into the MultiFab.
            Real const* pout = output.const_data_ptr<Real>();
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                auto const m = bx.index(IntVect(AMREX_D_DECL(i,j,k)));
                for (int n = 0; n < nout; ++n) {
                    out(i,j,k,n) = pout[m*nout+n];
                }
            });

            // batch and output are freed at the end of this iteration.
            Gpu::streamSynchronize();
        }

        // Check against the function the model represents, y = tanh(W x + b),
        // with the weights set in make_model.py.
        MultiFab error(ba, dm, nout, 0);
        {
            auto const& a = phi_in.const_arrays();
            auto const& y = phi_out.const_arrays();
            auto const& e = error.arrays();
            ParallelFor(error, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                Real const x0 = a[b](i,j,k,0);
                Real const x1 = a[b](i,j,k,1);
                e[b](i,j,k,0) = y[b](i,j,k,0) - std::tanh( 1.0_rt*x0 - 0.5_rt*x1 + 0.5_rt);
                e[b](i,j,k,1) = y[b](i,j,k,1) - std::tanh(0.25_rt*x0 + 2.0_rt*x1 - 0.25_rt);
            });
        }
        Real const max_error = error.norm0(0, nout, IntVect(0));
        amrex::Print() << "Output range: [" << phi_out.min(0) << ", " << phi_out.max(0) << "] and ["
                       << phi_out.min(1) << ", " << phi_out.max(1) << "]\n"
                       << "Max error vs. tanh(W x + b): " << max_error << "\n";

        if (!(max_error < 1000*std::numeric_limits<Real>::epsilon())) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
