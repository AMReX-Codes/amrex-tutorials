#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TimeIntegrator.H>
#include <AMReX_HypreSolver.H>

#include "myfunc.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

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

    if (AMREX_SPACEDIM != 2) {
        amrex::Abort("Only 2D supported; recompile with DIM=2");
    }
    
    // **********************************
    // SIMULATION PARAMETERS

    // number of cells on each side of the domain
    int n_cell;

    // size of each box (or grid)
    int max_grid_size;

    // total steps in simulation
    int nsteps;

    // how often to write a plotfile
    int plot_int;

    // time step
    Real dt;

    // use adaptive time step (dt used to set output times)
    bool adapt_dt = false;

    // adaptive time step relative and absolute tolerances
    Real reltol = 1.0e-4;
    Real abstol = 1.0e-9;

    // Advection and Diffusion Coefficients
    Real advCoeffx = 1.0;
    Real advCoeffy = 1.0;
    Real diffCoeffx = 1.0;
    Real diffCoeffy = 1.0;

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

        // Default nsteps to 10, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be written
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // time step
        pp.get("dt",dt);

        // use adaptive step sizes
        pp.query("adapt_dt",adapt_dt);

        // adaptive step tolerances
        pp.query("reltol",reltol);
        pp.query("abstol",abstol);

        pp.query("advCoeffx",advCoeffx);
        pp.query("advCoeffy",advCoeffy);
        pp.query("diffCoeffx",diffCoeffx);
        pp.query("diffCoeffy",diffCoeffy);
    }

    // **********************************
    // SIMULATION SETUP

    // AMREX_D_DECL means "do the first X of these, where X is the dimensionality of the simulation"
    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));

    // Make a single box that is the entire domain
    Box domain(dom_lo, dom_hi);

    // ba will contain a list of boxes that cover the domain
    // Initialize the boxarray "ba" from the single box "domain"
    BoxArray ba(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [0,1] in each direction.
    RealBox real_box({AMREX_D_DECL(-1.,-1.,-1.)},
                     {AMREX_D_DECL( 1., 1., 1.)});

    // periodic in all direction
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    // geom contains information such as the physical domain size,
    //               number of points in the domain, and periodicity
    // This defines a Geometry object
    Geometry geom(domain, real_box, CoordSys::cartesian, is_periodic);

    // extract dx from the geometry object
    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();
    GpuArray<Real,AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();
    GpuArray<Real,AMREX_SPACEDIM> prob_hi = geom.ProbHiArray();

    // Nghost = number of ghost cells for each array
    int Nghost = 1;

    // Ncomp = number of components for each array
    int Ncomp = 1;

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // allocate phi MultiFab
    MultiFab phi(ba, dm, Ncomp, Nghost);

    // time = starting time in the simulation
    Real time = 0.0;

    // **********************************
    // INITIALIZE DATA

    InitializeData(phi,dx,prob_lo,prob_hi,time,advCoeffx,advCoeffy);

    // Write a plotfile of the initial data if plot_int > 0
    if (plot_int > 0)
    {
        int step = 0;
        const std::string& pltfile = amrex::Concatenate("plt",step,5);
        WriteSingleLevelPlotfile(pltfile, phi, {"phi"}, geom, time, 0);
    }

    auto rhs_function = [&](MultiFab& S_rhs, MultiFab& S_data, const Real /* time */) {

        // fill periodic ghost cells
        S_data.FillBoundary(geom.periodicity());

        S_rhs.setVal(0.);
        
        ComputeDiffusion(S_rhs, S_data, diffCoeffx, diffCoeffy, dx);
        ComputeAdvection(S_rhs, S_data, advCoeffx, advCoeffy, dx);
    };

    auto rhs_im_function = [&](MultiFab& S_rhs, MultiFab& S_data, const Real /* time */) {

        // fill periodic ghost cells
        S_data.FillBoundary(geom.periodicity());

        S_rhs.setVal(0.);
        
        ComputeDiffusion(S_rhs, S_data, diffCoeffx, diffCoeffy, dx);
    };

    auto rhs_ex_function = [&](MultiFab& S_rhs, MultiFab& S_data, const Real /* time */) {

        // fill periodic ghost cells
        S_data.FillBoundary(geom.periodicity());

        S_rhs.setVal(0.);
        
        ComputeAdvection(S_rhs, S_data, advCoeffx, advCoeffy, dx);
    };

    constexpr int hypre_stencil_size = 2 * AMREX_SPACEDIM + 1;
    using HyprePreconditioner = HypreSolver<hypre_stencil_size>;

    std::unique_ptr<HyprePreconditioner> hypre_preconditioner;
    Real hypre_gamma = std::numeric_limits<Real>::quiet_NaN();

    auto build_hypre_preconditioner = [&](Real gamma)
    {
        const IntVect domain_lo = geom.Domain().smallEnd();
        const IntVect domain_hi = geom.Domain().bigEnd();
        const HYPRE_Real hx = static_cast<HYPRE_Real>(gamma * diffCoeffx / (dx[0] * dx[0]));
        const HYPRE_Real hy = static_cast<HYPRE_Real>(gamma * diffCoeffy / (dx[1] * dx[1]));
        const HYPRE_Real diag = static_cast<HYPRE_Real>(1.0) + 2.0 * (hx + hy);

        auto marker = [=] AMREX_GPU_DEVICE (int /* boxno */, int /* i */, int /* j */,
                                            int /* k */, int /* n */) -> bool
        {
            return true;
        };

        auto filler = [=] AMREX_GPU_DEVICE (int /* boxno */, int i, int j, int k, int n,
                                            Array4<HYPRE_Int const> const* gid,
                                            HYPRE_Int& ncols, HYPRE_Int* cols,
                                            HYPRE_Real* mat)
        {
            const int ilo = domain_lo[0];
            const int ihi = domain_hi[0];
            const int jlo = domain_lo[1];
            const int jhi = domain_hi[1];

            const int im = (i == ilo) ? ihi : i - 1;
            const int ip = (i == ihi) ? ilo : i + 1;
            const int jm = (j == jlo) ? jhi : j - 1;
            const int jp = (j == jhi) ? jlo : j + 1;

            ncols = 0;

            cols[ncols] = gid[n](im,j,k);
            mat [ncols] = -hx;
            ++ncols;

            cols[ncols] = gid[n](ip,j,k);
            mat [ncols] = -hx;
            ++ncols;

            cols[ncols] = gid[n](i,jm,k);
            mat [ncols] = -hy;
            ++ncols;

            cols[ncols] = gid[n](i,jp,k);
            mat [ncols] = -hy;
            ++ncols;

            cols[ncols] = gid[n](i,j,k);
            mat [ncols] = diag;
            ++ncols;
        };

        hypre_preconditioner = std::make_unique<HyprePreconditioner>(
            Vector<IndexType>{IndexType::TheCellType()},
            IntVect(1), geom, ba, dm, marker, filler);
        hypre_gamma = gamma;
    };

    auto precond_setup = [&](MultiFab& /* S_data */, MultiFab& /* S_rhs */, const Real /* time */,
                             bool jok, bool& jcur, const Real gamma)
    {
        const bool same_gamma = hypre_preconditioner &&
            std::abs(gamma - hypre_gamma) <=
            (10.0 * std::numeric_limits<Real>::epsilon() *
             std::max(Real(1.0), std::max(std::abs(gamma), std::abs(hypre_gamma))));

        if (!jok) {
            build_hypre_preconditioner(gamma);
            jcur = true;
        }
    };

    auto precond_solve = [&](MultiFab& S_soln, MultiFab& S_rhs, MultiFab& /* S_data */,
                             MultiFab& /* S_state_rhs */, const Real /* time */,
                             const Real /* gamma */, const Real /* delta */,
                             int /* lr */)
    {
        AMREX_ALWAYS_ASSERT(hypre_preconditioner != nullptr);
        S_soln.setVal(0.0);
        hypre_preconditioner->solve(Vector<MultiFab*>{&S_soln},
                                    Vector<MultiFab const*>{&S_rhs},
                                    0.0_rt, 0.0_rt, 1);
    };

    TimeIntegrator<MultiFab> integrator(phi, time);
    integrator.set_rhs(rhs_function);
    integrator.set_imex_rhs(rhs_im_function, rhs_ex_function);
    integrator.set_preconditioner(precond_setup, precond_solve);

    if (adapt_dt) {
        integrator.set_adaptive_step();
        integrator.set_tolerances(reltol, abstol);
    } else {
        integrator.set_time_step(dt);
    }

    Real evolution_start_time = ParallelDescriptor::second();

    for (int step = 1; step <= nsteps; ++step)
    {
        // Set time to evolve to
        time += dt;

        Real step_start_time = ParallelDescriptor::second();

        // Advance to output time
        integrator.evolve(phi, time);

        Real step_stop_time = ParallelDescriptor::second() - step_start_time;
        ParallelDescriptor::ReduceRealMax(step_stop_time);

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << step << " in " << step_stop_time << " seconds; dt = " << dt << " time = " << time << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && step%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",step,5);
            WriteSingleLevelPlotfile(pltfile, phi, {"phi"}, geom, time, step);
        }
    }

    Real evolution_stop_time = ParallelDescriptor::second() - evolution_start_time;
    ParallelDescriptor::ReduceRealMax(evolution_stop_time);
    amrex::Print() << "Total evolution time = " << evolution_stop_time << " seconds\n";

    // exact solution
    BoxArray ba_exact(domain);
    DistributionMapping dm_exact(ba_exact);
    MultiFab phi_exact(ba_exact, dm_exact, Ncomp, n_cell);
    InitializeData(phi_exact,dx,prob_lo,prob_hi,time,advCoeffx,advCoeffy);
    phi_exact.SumBoundary(geom.periodicity());

    {
        const std::string& pltfile = amrex::Concatenate("exact",nsteps,5);
        WriteSingleLevelPlotfile(pltfile, phi_exact, {"phi"}, geom, time, nsteps);
    }

    MultiFab phi_exact_dist(ba,dm,Ncomp,0);
    phi_exact_dist.ParallelCopy(phi_exact,0,0,1);
    
    MultiFab::Subtract(phi_exact_dist,phi,0,0,1,0);

    {
        const std::string& pltfile = amrex::Concatenate("diff",nsteps,5);
        WriteSingleLevelPlotfile(pltfile, phi_exact_dist, {"phi"}, geom, time, nsteps);
    }

    Real error = phi_exact_dist.norm1(0,geom.periodicity());
    amrex::Print() << "L1 error = " << error << std::endl;
}

void InitializeData(MultiFab& phi,
                    const GpuArray<Real,AMREX_SPACEDIM> dx,
                    const GpuArray<Real,AMREX_SPACEDIM> prob_lo,
                    const GpuArray<Real,AMREX_SPACEDIM> prob_hi,
                    const Real& time,
                    const Real& Ax,
                    const Real& Ay) {

    int ng = phi.nGrow();
    
    GpuArray<Real,AMREX_SPACEDIM> L;
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
        L[d] = prob_hi[d] - prob_lo[d];
    }

    // loop over boxes
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);
        const Array4<Real>& phi_array = phi.array(mfi);

        Real sigma = 0.1;
        Real a = 1.0/(sigma*sqrt(2*M_PI));
        Real b = -0.5/(sigma*sigma);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Real x = prob_lo[0] + (((Real) i) + 0.5) * dx[0];
            Real y = prob_lo[1] + (((Real) j) + 0.5) * dx[1];
            Real r = (x-Ax*time) * (x-Ax*time) + (y-Ay*time) * (y-Ay*time);
            phi_array(i,j,k) = a * std::exp(b * r);
        });
    }

}

void ComputeDiffusion(MultiFab& S_rhs,
                      MultiFab& S_data,
                      const Real& Dx,
                      const Real& Dy,
                      const GpuArray<Real,AMREX_SPACEDIM> dx) {

    for ( MFIter mfi(S_data,TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.tilebox();

        const Array4<const Real>& phi_array = S_data.array(mfi);
        const Array4<      Real>& rhs_array = S_rhs.array(mfi);

        // fill the right-hand-side for phi
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rhs_array(i,j,k) += Dx * ( (phi_array(i+1,j,k) - 2.*phi_array(i,j,k) + phi_array(i-1,j,k)) / (dx[0]*dx[0]) )
                              + Dy * ( (phi_array(i,j+1,k) - 2.*phi_array(i,j,k) + phi_array(i,j-1,k)) / (dx[1]*dx[1]) );
        });
    }
}


void ComputeAdvection(MultiFab& S_rhs,
                      MultiFab& S_data,
                      const Real& Ax,
                      const Real& Ay,
                      const GpuArray<Real,AMREX_SPACEDIM> dx) {

    Real dxInv = 1.0 / dx[0];
    Real dyInv = 1.0 / dx[1];
    Real sideCoeffx = Ax * dxInv;
    Real sideCoeffy = Ay * dyInv;

    for (MFIter mfi(S_data,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        const Array4<const Real>& phi_array = S_data.array(mfi);
        const Array4<      Real>& rhs_array = S_rhs.array(mfi);

        // x-direction
        if (Ax > 0)
        {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                rhs_array(i,j,k) -= sideCoeffx * (phi_array(i,j,k) - phi_array(i-1,j,k));
            });
        }
        else
        {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                rhs_array(i,j,k) -= sideCoeffx * (phi_array(i+1,j,k) - phi_array(i,j,k));
            });
        }

        // y-direction
        if (Ay > 0)
        {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                rhs_array(i,j,k) -= sideCoeffy * (phi_array(i,j,k) - phi_array(i,j-1,k));
            });
        }
        else
        {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                rhs_array(i,j,k) -= sideCoeffy * (phi_array(i,j+1,k) - phi_array(i,j,k));
            });
        }
    }
}
