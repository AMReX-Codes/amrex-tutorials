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

    // **********************************
    // SIMULATION PARAMETERS

    // number of cells on each side of the domain
    int n_cell;

    // size of each box (or grid)
    int max_grid_size;

    // pytorch model file
    std::string model_filename;

    // largest time step
    Real dt;

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

        // time step
        pp.get("dt", dt);
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

    // Ncomp = number of components for each array
    int Ncomp = 1;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two phi multifabs; one will store the input state, the other the output.
    MultiFab phi_in(ba, dm, Ncomp, Nghost);
    MultiFab phi_out(ba, dm, Ncomp+1, Nghost);

    // random seed
    int seed = 42;

    // **********************************
    // INITIALIZE DATA

    ResetRandomSeed(seed);

    // loop over boxes
    for (MFIter mfi(phi_in); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        const Array4<Real>& phi_input = phi_in.array(mfi);

        // set phi = random(0, dt)
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            phi_input(i,j,k) = amrex::Random() * dt;
        });
    }

    // Write a plotfile of the initial data 
    const std::string& pltfile = amrex::Concatenate("plt",0,5);
    WriteSingleLevelPlotfile(pltfile, phi_in, {"dt"}, geom, 0.0, 0);

    // **********************************
    // LOAD PYTORCH MODEL

    // Load pytorch module via torch script
    torch::jit::script::Module module;
    try {
	// Deserialize the ScriptModule from a file using torch::jit::load().
	module = torch::jit::load(model_filename);
    }
    catch (const c10::Error& e) {
	std::cerr << "error loading the model\n";
	return;
    }
    
    Print() << "Model loaded.\n";

    // **********************************
    // EVALUATE MODEL
    
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
	
	// create torch tensor
        at::Tensor t1 = torch::zeros({ncell, Ncomp});

#ifdef USE_AMREX_CUDA
        t1 = t1.to(torch::kCUDA);
#endif

	// copy input multifab to torch tensor
	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
	    int ii = i - bx_lo[0];
	    int jj = j - bx_lo[1];
	    int index = jj*nbox[0] + ii;
#if AMREX_SPACEDIM == 3
	    int kk = k - bx_lo[2];
	    index += kk*nbox[0]*nbox[1];
#endif
	    t1[index][0] = phi_input(i, j, k, 0);
        });

        // create torch data array
        std::vector<torch::jit::IValue> inputs_torch{t1};
        at::Tensor outputs_torch = module.forward(inputs_torch).toTensor();
	
#ifdef USE_AMREX_CUDA
        outputs_torch = outputs_torch.to(torch::kCUDA);
#endif

	// copy tensor to output multifab
	ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // const int index = AMREX_SPACEDIM == 2 ?
            //         i*nbox[1] + j : (i*nbox[1] + j)*nbox[2] + k;
	    int ii = i - bx_lo[0];
	    int jj = j - bx_lo[1];
	    int index = jj*nbox[0] + ii;
#if AMREX_SPACEDIM == 3
	    int kk = k - bx_lo[2];
	    index += kk*nbox[0]*nbox[1];
#endif
	    phi_output(i, j, k, 0) = outputs_torch[index][0].item<double>();
	    phi_output(i, j, k, 1) = outputs_torch[index][1].item<double>();
        });
    }

    // Tell the I/O Processor to write out that we're done
    amrex::Print() << "Finish evaluating model.\n";
    
    // Write a plotfile of the current data
    const std::string& pltfile2 = amrex::Concatenate("plt",1,5);
    WriteSingleLevelPlotfile(pltfile2, phi_out, {"y0", "y1"}, geom, dt, 1);
}
