#include <AMReX_Particles.H>
#include <AMReX.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_TagBox.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_WriteEBSurface.H>

#include "Indexing.H"
#include "FluidParticleContainer.H"

using namespace amrex;

extern void make_eb_cylinder(const Geometry& geom);
extern void define_velocity(const Real time, const Geometry& geo, Array<MultiFab,AMREX_SPACEDIM>& vel_out, const MultiFab& phi);
extern void mac_project_velocity(Array<MultiFab,AMREX_SPACEDIM>& vel_out, const Geometry& geom, int use_hypre); 

Real est_time_step(const Real current_dt, const Geometry& geom, Array<MultiFab,AMREX_SPACEDIM>& vel, const Real cfl)
{
    const Real change_max = 1.1;

    Real dt_est = std::numeric_limits<Real>::max();

    const Real* dx =  geom.CellSize();

    const Vector<std::string> coord_dir {AMREX_D_DECL("x", "y", "z")};

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        Real est = vel[idim].norm0(0,0,false);
        // amrex::Print() << "Max vel in " << coord_dir[idim] << "-direction is " << est << std::endl;
        dt_est = amrex::min(dt_est, dx[idim]/est);
    }

    ParallelDescriptor::ReduceRealMin(dt_est);

    // Apply our CFL
    dt_est *= cfl;

    // Do not grow the timestep by more than a ratio of change_max
    dt_est = std::min(dt_est, current_dt * change_max);

    return dt_est;
}

void write_plotfile(int step, Real time, const Geometry& geom, MultiFab& plotmf, 
                    FluidParticleContainer& pc, int write_ascii)
{
    // Copy processor id into the component of plotfile_mf immediately following velocities
    int proc_comp = AMREX_SPACEDIM; 
    for (MFIter mfi(plotmf); mfi.isValid(); ++mfi)
       plotmf[mfi].setVal<RunOn::Host>(ParallelDescriptor::MyProc(),mfi.validbox(),proc_comp,1);

    std::stringstream sstream;
    sstream << "plt" << std::setw(5) << std::setfill('0') << step;
    std::string plotfile_name = sstream.str();
    
#if (AMREX_SPACEDIM == 2)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "xvel", "yvel", "proc", "phi" },
                                     geom, time, 0);
#elif (AMREX_SPACEDIM == 3)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "xvel", "yvel", "zvel", "proc", "phi" },
                                     geom, time, 0);
#endif

    pc.Checkpoint(plotfile_name, "particles", true); // Write particles to plotfile

    std::stringstream pstream;
    pstream << "part" << std::setw(5) << std::setfill('0') << step;
    const std::string ascii_filename = pstream.str();

    if (write_ascii) {
       pc.WriteAsciiFile(ascii_filename);
    }

}

int main (int argc, char* argv[])
{
    // Turn off amrex-related output
    amrex::SetVerbose(0);

    amrex::Initialize(argc, argv);

    Real strt_time = amrex::second();
    Real eb_strt_time;
    Real eb_stop_time;
    
    {
        int n_cell = 128;
        int max_grid_size = 32;
        int n_ppc = 4;
        int pic_interpolation = Interpolation::CIC;
        Real stop_time = 1000.0;
        int max_step = 100;
        int plot_int  = 1;
        int write_ascii  = 0;
        int write_initial_phi  = 0;
        int write_eb_geom      = 1;
        int use_hypre  = 0;
        Real phi_cutoff = 0.1;
        Real cfl = 0.7;

        Real dt = std::numeric_limits<Real>::max();

        // read parameters
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("n_ppc", n_ppc);
            pp.query("pic_interpolation", pic_interpolation);
            pp.query("stop_time", stop_time);
            pp.query("max_step", max_step);
            pp.query("plot_int", plot_int);
            pp.query("use_hypre", use_hypre);
            pp.query("write_ascii", write_ascii);
            pp.query("write_initial_phi", write_initial_phi);
            pp.query("write_eb_geom", write_eb_geom);
            pp.query("phi_cutoff", phi_cutoff);
            pp.query("cfl", cfl);
        }

#ifndef AMREX_USE_HYPRE
        if (use_hypre == 1) 
           amrex::Abort("Cant use hypre if we dont build with USE_HYPRE=TRUE");
#endif

        if (n_cell%8 != 0)
           amrex::Abort("n_cell must be a multiple of 8");

        int n_cell_x = n_cell;
        int n_cell_y = n_cell;
#if (AMREX_SPACEDIM == 3)
        int n_cell_z = n_cell/8;
#endif

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.0,1.0,0.125)});

            Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,1)};
            Geometry::Setup(&rb, 0, is_periodic.data());
            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)});
            geom.define(domain);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        const int ncomp_phi = 1;
        Vector<BCRec> phi_bc(ncomp_phi);
        for (int n = 0; n < ncomp_phi; ++n)
        {
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
            {
                // is_periodic overrides inputs in domain_(lo/hi)_bc_type
                if (geom.isPeriodic(i))
                {
                    // Set the BCs to interior Dirichlet if we are periodic
                    // in this dimension:
                    phi_bc[n].setLo(i, BCType::int_dir);
                    phi_bc[n].setHi(i, BCType::int_dir);
                }
                else
                {
                    // Use first order extrapolation to enforce Neumann BCs with 0 gradient
                    // at the boundaries if we are not periodic in this dimension:
                    phi_bc[n].setLo(i, BCType::foextrap);
                    phi_bc[n].setHi(i, BCType::foextrap);
                }
            }
        }

        Array<MultiFab,AMREX_SPACEDIM> vel;
        MultiFab plotfile_mf;
        MultiFab phi_mf;

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible

        eb_strt_time = amrex::second();
        make_eb_cylinder(geom);
        eb_stop_time = amrex::second() - eb_strt_time;

        std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > factory =
           makeEBFabFactory(geom, grids, dmap, {4, 4, 2}, EBSupport::full);
        const EBFArrayBoxFactory* ebfact = &(static_cast<amrex::EBFArrayBoxFactory const&>(*factory));

        // Velocities are face-centered
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].define (amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 1, MFInfo(), *factory);
        }

        // store plotfile variables; velocity, processor id, and phi (the EB writer appends volfrac)
        plotfile_mf.define(grids, dmap, AMREX_SPACEDIM+2, 0, MFInfo(), *factory);
        
        // make a separate phi MultiFab for the particle-mesh operations because we need a ghost cell
        phi_mf.define(grids, dmap, 1, 1, MFInfo(), *factory);

        // Get volume fraction for the embedded geometry
        // volume fraction = 0 for cells covered by the embedded geometry
        // volume fraction = 1 for cells not covered by the embedded geometry
        // 0 < volume fraction < 1 for cells cut by the embedded geometry
        const MultiFab& vol_mf = ebfact->getVolFrac();

        // Initialize a gaussian density profile in the domain for cells
        // not covered by the embedded geometry by multiplying phi
        // by the volume fraction.
        for (MFIter mfi(phi_mf); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            Array4<Real> phi = phi_mf[mfi].array();
            Array4<const Real> vol = vol_mf[mfi].array();
            const auto plo = geom.ProbLoArray();
            const auto dx = geom.CellSizeArray();

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real x = plo[0] + (0.5+i) * dx[0]; 
                Real y = plo[1] + (0.5+j) * dx[1];

                Real r2 = (pow(x-0.5, 2) + pow(y-0.75,2)) / 0.01;

                Real gauss = std::exp(-r2);
                gauss = gauss >= phi_cutoff ? gauss : 0.0;

                phi(i,j,k) = gauss * vol(i,j,k);
            });
        }

        // Fill periodic BCs and interior ghost cells
        phi_mf.FillBoundary(geom.periodicity());

        // Fill non-periodic domain BCs (e.g. Neumann in x & y dimensions)
        FillDomainBoundary(phi_mf, geom, phi_bc);

        // Set cells covered by the embedded boundary to phi = 0 so
        // when we interpolate to particles, the covered cells do not contribute.
        EB_set_covered(phi_mf,0.0);

        if (write_initial_phi) {
            const std::string pfname = "initial_phi";
            WriteSingleLevelPlotfile(pfname, phi_mf, {"phi"}, geom, 0.0, 0);
        }

        // Initialize Particles
        FluidParticleContainer FPC(geom, dmap, grids);

        // Initialize n_ppc randomly located particles per cell.
        // Particles are weighted by interpolated density field phi.
        // Only creates particles in regions not covered by the embedded geometry.
        FPC.InitParticles(phi_mf, vol_mf, phi_cutoff, n_ppc, pic_interpolation);

        FPC.DepositToMesh(phi_mf, pic_interpolation);
        EB_set_covered(phi_mf,-1.0);

        if (write_initial_phi) {
            const std::string pfname = "initial_phi_after_deposit";
            WriteSingleLevelPlotfile(pfname, phi_mf, {"phi"}, geom, 0.0, 0);
        }

        // set initial velocity to u=(1,0,0)
        AMREX_D_TERM(vel[0].setVal(1.0);,
                     vel[1].setVal(0.0);,
                     vel[2].setVal(0.0););

#if (AMREX_SPACEDIM == 3)
        if (write_eb_geom) 
        {
            amrex::Print() << "Writing EB surface" << std::endl;
            WriteEBSurface (grids, dmap, geom, ebfact);
        }
#endif

        Real time = 0.0;

        // Write out the initial data
        {
            amrex::Print() << "Creating the initial velocity field " << std::endl;
            define_velocity(time,geom,vel,phi_mf);
            mac_project_velocity(vel,geom,use_hypre);
            EB_average_face_to_cellcenter(plotfile_mf,0,amrex::GetArrOfConstPtrs(vel));

            // copy initial deposited phi into the plotfile
            MultiFab::Copy(plotfile_mf, phi_mf, 0, Idx::phi, 1, 0);

            amrex::Print() << "Writing the initial data into plt00000\n" << std::endl;
            write_plotfile(0, time, geom, plotfile_mf, FPC, write_ascii);
        }

        // This computes the first dt
        dt = est_time_step(dt, geom, vel, cfl);

        int nstep = 0;

        // Sum phi to check conservation
        // Real sum_phi = phi_mf.sum();
        Real sum_phi = FPC.SumPhi();

        amrex::Print() << "Initial sum of phi is " << sum_phi << std::endl;

        for (int i = 0; i < max_step; i++)
        {
            if (time < stop_time)
            {
                amrex::Print() << "STEP " << i+1 << " starts at TIME = " << time
                               << " DT = " << dt << std::endl;

                dt = amrex::min(dt, stop_time - time);

                Real t_nph = time + 0.5 * dt;

                define_velocity(t_nph,geom,vel,phi_mf);
                mac_project_velocity(vel,geom,use_hypre);

                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    vel[idim].FillBoundary(geom.periodicity());
                }

                // Step Particles
                FPC.AdvectWithUmac(vel.data(), 0, dt);

                // Redistribute Particles across MPI ranks with their new positions
                FPC.Redistribute();

                // Deposit Particles to the grid to update phi
                FPC.DepositToMesh(phi_mf, pic_interpolation);
                EB_set_covered(phi_mf,-1.0);

                // Increment time
                time += dt;
                nstep++;

                // Write to a plotfile
                if (plot_int > 0 && (i+1)%plot_int == 0)
                {
                   average_face_to_cellcenter(plotfile_mf,0,amrex::GetArrOfConstPtrs(vel));

                   // copy phi into the plotfile
                   MultiFab::Copy(plotfile_mf, phi_mf, 0, Idx::phi, 1, 0);

                   write_plotfile(i+1, time, geom, plotfile_mf, FPC, write_ascii);
                }

                // Sum phi to check conservation
                // sum_phi = phi_mf.sum();
                sum_phi = FPC.SumPhi();

                amrex::Print() << "STEP " << i+1 << " ends   at TIME = " << time
                               << " DT = " << dt << " Sum(Phi) = " << sum_phi << "\n" << std::endl;

                // Compute lagged dt for next time step based on this half-time velocity
                dt = est_time_step(dt, geom, vel, cfl);

            } else {

                // Copy velocity into plotfile
                average_face_to_cellcenter(plotfile_mf,0,amrex::GetArrOfConstPtrs(vel));

                // Write to a plotfile
                write_plotfile(i+1, time, geom, plotfile_mf, FPC, write_ascii);
                break;
            }
        }
    }

    Real stop_time = amrex::second() - strt_time;
    amrex::Print() << "\nTime to create EB geometry " << eb_stop_time << std::endl;
    amrex::Print() << "Total run time             " << stop_time << std::endl;

    amrex::Finalize();
}
