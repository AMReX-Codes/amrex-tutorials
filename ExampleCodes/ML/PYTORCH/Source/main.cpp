#include<torch/script.h> // One-stop header

#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

#include "myfunc.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void main_main ()
{

    BL_PROFILE_VAR("Init",Init);

    // store the current time so we can later compute total run time.
    Real strt_time = ParallelDescriptor::second();
    Real eval_time = 0.0;

    // **********************************
    // SIMULATION PARAMETERS

    // number of cells on each side of the domain
    int n_cell;

    // size of each box (or grid)
    int max_grid_size;

    // pytorch model file
    std::string model_filename;

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        // pp.get means we require the inputs file to have it
        // pp.query means we optionally need the inputs file to have it - but we must supply a default here
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Read in model file name
        pp.query("model_file", model_filename);
    }

    // **********************************
    // SIMULATION SETUP

    // make BoxArray and Geometry
    // ba will contain a list of boxes that cover the domain
    // geom contains information such as the physical domain size,
    //               number of points in the domain, and periodicity
    BoxArray ba;
    Geometry geom;

    // AMREX_D_DECL means "do the first X of these, where X is the dimensionality of the simulation"
    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "domain"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [0,1] in each direction.
    RealBox real_box({AMREX_D_DECL( 0., 0., 0.)},
                     {AMREX_D_DECL( 1., 1., 1.)});

    // not periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

    // This defines a Geometry object
    geom.define(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // Nghost = number of ghost cells for each array
    int Nghost = 0;

    // Nc_in = number of components for input array
    // Nc_out = number of components for output array
    int Nc_in = 1;
    int Nc_out = 2;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two phi multifabs; one will store the input state, the other the output.
    MultiFab phi_in(ba, dm, Nc_in, Nghost);
    MultiFab phi_out(ba, dm, Nc_out, Nghost);

    // **********************************
    // INITIALIZE DATA

    BL_PROFILE_VAR_STOP(Init);

    BL_PROFILE_VAR("InitData",InitData);

    // loop over boxes
    for (MFIter mfi(phi_in); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& phi_input = phi_in.array(mfi);

        // set phi
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            phi_input(i,j,k) = (i+j+k)/(n_cell-1.0);
        });
    }

    BL_PROFILE_VAR_STOP(InitData);

    BL_PROFILE_VAR("WriteInitPlot",WriteInitPlot);

    // Write a plotfile of the initial data
    const std::string& pltfile = "plt_inputs";
    WriteSingleLevelPlotfile(pltfile, phi_in, {"dt"}, geom, 0.0, 0);

    BL_PROFILE_VAR_STOP(WriteInitPlot);

    // **********************************
    // LOAD PYTORCH MODEL

    BL_PROFILE_VAR("LoadPytorch",LoadPytorch);

    // Load pytorch module via torch script
    torch::jit::script::Module module;
    try {
        // Deserialize the ScriptModule from a file using torch::jit::load().
        module = torch::jit::load(model_filename);
    }
    catch (const c10::Error& e) {
        amrex::Abort("Error loading the model\n");
    }

    Print() << "Model loaded.\n";

    // set pytorch data type (default is float or torch::kFloat32)
    auto dtype0 = torch::kFloat64;

#ifdef AMREX_USE_CUDA
    torch::Device device0(torch::kCUDA);
    module.to(device0);
    amrex::Print() << "Copying model to GPU." << std::endl;

    // set tensor options
    auto tensoropt = torch::TensorOptions().dtype(dtype0).device(device0);
#else
    auto tensoropt = torch::TensorOptions().dtype(dtype0);
#endif

    BL_PROFILE_VAR_STOP(LoadPytorch);

    // **********************************
    // EVALUATE MODEL

    BL_PROFILE_VAR("Eval",Eval);

    // loop over boxes
    for ( MFIter mfi(phi_in); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& phi_input = phi_in.array(mfi);
        const Array4<Real>& phi_output = phi_out.array(mfi);

        // retrieve smallend and size of box
        const IntVect bx_lo = bx.smallEnd();
        const IntVect nbox = bx.size();

        // compute total cells in the box
        int ncell = AMREX_SPACEDIM == 2 ?
            nbox[0] * nbox[1] : nbox[0] * nbox[1] * nbox[2];

        // create a temporary array to store MultiFab data
        // this is needed to create a tensor from a contiguous block of memory
        amrex::Gpu::ManagedVector<Real> aux(ncell*Nc_in);
        Real* AMREX_RESTRICT auxPtr = aux.dataPtr();

        // copy input multifab to torch tensor
        amrex::ParallelFor(bx, Nc_in, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
        {
            int ii = i - bx_lo[0];
            int jj = j - bx_lo[1];
            int index = jj*nbox[0] + ii;
#if AMREX_SPACEDIM == 3
            int kk = k - bx_lo[2];
            index += kk*nbox[0]*nbox[1];
#endif
            // array order is row-based [index][comp]
            auxPtr[index*Nc_in + n] = phi_input(i, j, k, n);
        });

        // create torch tensor from array
        at::Tensor inputs_torch = torch::from_blob(auxPtr, {ncell, Nc_in}, tensoropt);

        // store the current time so we can later compute total eval time.
        Real eval_t_start = ParallelDescriptor::second();

        // evaluate torch model
        at::Tensor outputs_torch = module.forward({inputs_torch}).toTensor();
        outputs_torch = outputs_torch.to(dtype0);

        // add eval time
        eval_time += ParallelDescriptor::second() - eval_t_start;

        // get accessor to tensor (read-only)
#ifdef AMREX_USE_CUDA
        auto outputs_torch_acc = outputs_torch.packed_accessor64<Real,2>();
#else
        auto outputs_torch_acc = outputs_torch.accessor<Real,2>();
#endif

        // copy tensor to output multifab
        amrex::ParallelFor(bx, Nc_out, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept
        {
            int ii = i - bx_lo[0];
            int jj = j - bx_lo[1];
            int index = jj*nbox[0] + ii;
#if AMREX_SPACEDIM == 3
            int kk = k - bx_lo[2];
            index += kk*nbox[0]*nbox[1];
#endif
            phi_output(i, j, k, n) = outputs_torch_acc[index][n];
        });
    }

    BL_PROFILE_VAR_STOP(Eval);

    BL_PROFILE_VAR("Post",Post);

    // Tell the I/O Processor to write out that we're done
    amrex::Print() << "Finish evaluating model.\n";

    // Write a plotfile of the current data
    const std::string& pltfile2 = "plt_outputs";
    WriteSingleLevelPlotfile(pltfile2, phi_out, {"X_0", "X_1"}, geom, 1.0, 1);

    // Call the timer again and compute the maximum difference between the start time
    // and stop time over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(stop_time);
    amrex::Print() << "Run time = " << stop_time << std::endl;
    ParallelDescriptor::ReduceRealMax(eval_time);
    amrex::Print() << "Eval time = " << eval_time << std::endl;

    BL_PROFILE_VAR_STOP(Post);
}
