#include <AMReX_Particles.H>
#include <AMReX.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MacProjector.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_TagBox.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_WriteEBSurface.H>

#include <MyParticleContainer.H>

using namespace amrex;

void write_plotfile(int step_counter, const auto& geom, const auto& plotmf, auto& pc, int write_ascii)
{
    std::stringstream sstream;
    sstream << "plt" << std::setw(5) << std::setfill('0') << step_counter;
    std::string plotfile_name = sstream.str();

    // amrex::Print() << "Writing " << plotfile_name << std::endl;    
    
    if (step_counter == 0)
    {
#if (AMREX_SPACEDIM == 2)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "proc" ,"xvel", "yvel" },
                                     geom, 0.0, 0);
#elif (AMREX_SPACEDIM == 3)
       EB_WriteSingleLevelPlotfile(plotfile_name, plotmf,
                                   { "proc", "xvel", "yvel", "zvel" },
                                     geom, 0.0, 0);
#endif
    }

    pc.Checkpoint(plotfile_name, "particles", true); //Write Tracer particles to plotfile

    std::stringstream pstream;
    pstream << "part" << std::setw(5) << std::setfill('0') << step_counter;
    const std::string ascii_filename = pstream.str();

    if (write_ascii)
       pc.WriteAsciiFile(ascii_filename);

}

int main (int argc, char* argv[])
{
    // Turn off amrex-related output
    amrex::SetVerbose(0);

    amrex::Initialize(argc, argv);

    Real strt_time = amrex::second();
    
    {
        int mg_verbose = 0;
        int cg_verbose = 0;
        int n_cell = 128;
        int max_grid_size = 32;
        std::string particle_file = "";
        Real max_time = 1.0;
        int max_steps = 100;
        int plot_int  = 1;
        int write_ascii  = 0;
        int use_hypre  = 0;
        Real time_step = 0.01;

        Real obstacle_radius = 0.10;
        Real particle_radius = 0.02;

        amrex::Vector<int> ob_id;

        // read parameters
        {
            ParmParse pp;
            pp.query("mg_verbose", mg_verbose);
            pp.query("cg_verbose", cg_verbose);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("particle_file", particle_file);
            pp.query("max_time", max_time);
            pp.query("max_steps", max_steps);
            pp.query("plot_int", plot_int);
            pp.query("use_hypre", use_hypre);
            pp.query("write_ascii", write_ascii);
            pp.query("time_step", time_step);

            pp.queryarr("obstacles", ob_id);

            pp.query("obstacle_radius", obstacle_radius);
            pp.query("particle_radius", particle_radius);
        }

#ifndef AMREX_USE_HYPRE
        if (use_hypre == 1) 
           amrex::Abort("Cant use hypre if we dont build with USE_HYPRE=TRUE");
#endif

        if (n_cell%8 != 0)
           amrex::Abort("n_cell must be a multiple of 8");

        int n_cell_x = 2*n_cell;
        int n_cell_y =   n_cell;
        int n_cell_z =   n_cell/8;
        int num_obstacles;

        if (ob_id.empty())
        {
           amrex::Print() << " **************************************************** "     << std::endl;
           amrex::Print() << " You didn't specify any obstacles -- please try again " << std::endl;
           amrex::Print() << " ****************************************************\n "     << std::endl;
           exit(0);

        } else {

           num_obstacles = ob_id.size();

           if (num_obstacles > 9)
           {
              amrex::Print() << " **************************************************** "     << std::endl;
              amrex::Print() << " We only have 9 possible obstacles " << std::endl;
              amrex::Print() << " You specified too many -- please try again " << std::endl;
              amrex::Print() << " ****************************************************\n "     << std::endl;
              exit(0);
           } 

           for (int i = 0; i < num_obstacles; i++) 
              if (ob_id[i] < 0 || ob_id[i] > 8)
              {
                 amrex::Print() << " **************************************************** "     << std::endl;
                 amrex::Print() << " The obstacles must be identified using integers from 0 through 8 (inclusive) " << std::endl;
                 amrex::Print() << " You specified an invalid obstacle -- please try again " << std::endl;
                 amrex::Print() << " ****************************************************\n "     << std::endl;
                 exit(0);
              }

           amrex::Print() << " \n********************************************************************" << std::endl; 
           amrex::Print() << " You specified " << num_obstacles << " objects in the domain: ";
              for (int i = 0; i < num_obstacles; i++) 
                  amrex::Print() << ob_id[i] << " ";
             amrex::Print() << std::endl;
           amrex::Print() << " ********************************************************************" << std::endl; 
        } 

        Real zlen = 0.125;

        Geometry geom;
        BoxArray grids;
        DistributionMapping dmap;
        {
            RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(2.0,1.0,zlen)});

            Array<int,AMREX_SPACEDIM> isp{AMREX_D_DECL(0,1,1)};
            Geometry::Setup(&rb, 0, isp.data());
            Box domain(IntVect{AMREX_D_DECL(0,0,0)},
                       IntVect{AMREX_D_DECL(n_cell_x-1,n_cell_y-1,n_cell_z-1)});
            geom.define(domain);

            grids.define(domain);
            grids.maxSize(max_grid_size);

            dmap.define(grids);
        }

        Array<MultiFab,AMREX_SPACEDIM> vel;
        Array<MultiFab,AMREX_SPACEDIM> beta;
        MultiFab plotfile_mf;

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 100;    // typically a huge number so MG coarsens as much as possible

        amrex::Vector<amrex::RealArray> obstacle_center = {
            {AMREX_D_DECL(0.3,0.2,0.5)},
            {AMREX_D_DECL(0.3,0.5,0.5)},
            {AMREX_D_DECL(0.3,0.8,0.5)},
            {AMREX_D_DECL(0.7,0.25,0.5)},
            {AMREX_D_DECL(0.7,0.60,0.5)},
            {AMREX_D_DECL(0.7,0.85,0.5)},
            {AMREX_D_DECL(1.1,0.2,0.5)},
            {AMREX_D_DECL(1.1,0.5,0.5)},
            {AMREX_D_DECL(1.1,0.8,0.5)}};

        int direction =  2;
        Real height   = -1.0;  // Putting a negative number for height means it extends beyond the domain

        // The "false" below is the boolean that determines if the fluid is inside ("true") or 
        //     outside ("false") the object(s)

        Array<EB2::CylinderIF,9> obstacles{
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 0], false),
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 1], false),
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 2], false),
            EB2::CylinderIF(0.9*obstacle_radius, height, direction, obstacle_center[ 3], false),
            EB2::CylinderIF(0.9*obstacle_radius, height, direction, obstacle_center[ 4], false),
            EB2::CylinderIF(0.9*obstacle_radius, height, direction, obstacle_center[ 5], false),
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 6], false),
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 7], false),
            EB2::CylinderIF(obstacle_radius, height, direction, obstacle_center[ 8], false)};

        switch(num_obstacles) {

           case 1:
              {
              auto gshop1 = EB2::makeShop(obstacles[ob_id[0]]);
              EB2::Build(gshop1, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 2:
              {
              amrex::Print() << "Objects " << ob_id[0] << " " << ob_id[1] << std::endl;
              auto all2 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]]);
              auto gshop2  = EB2::makeShop(all2);
              EB2::Build(gshop2, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 3:
              {
              auto all3 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto gshop3  = EB2::makeShop(all3);
              EB2::Build(gshop3, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 4:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto all     = EB2::makeUnion(group_1,obstacles[ob_id[3]]);
              auto gshop4  = EB2::makeShop(all);
              EB2::Build(gshop4, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 5:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto group_2 = EB2::makeUnion(obstacles[ob_id[3]],obstacles[ob_id[4]]);
              auto all     = EB2::makeUnion(group_1,group_2);
              auto gshop5  = EB2::makeShop(all);
              EB2::Build(gshop5, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 6:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto group_2 = EB2::makeUnion(obstacles[ob_id[3]],obstacles[ob_id[4]],obstacles[ob_id[5]]);
              auto all     = EB2::makeUnion(group_1,group_2);
              auto gshop6  = EB2::makeShop(all);
              EB2::Build(gshop6, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 7:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto group_2 = EB2::makeUnion(obstacles[ob_id[3]],obstacles[ob_id[4]],obstacles[ob_id[5]]);
              auto group_3 = obstacles[ob_id[6]];
              auto all     = EB2::makeUnion(group_1,group_2,group_3);
              auto gshop7  = EB2::makeShop(all);
              EB2::Build(gshop7, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 8:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto group_2 = EB2::makeUnion(obstacles[ob_id[3]],obstacles[ob_id[4]],obstacles[ob_id[5]]);
              auto group_3 = EB2::makeUnion(obstacles[ob_id[6]],obstacles[ob_id[7]]);
              auto all     = EB2::makeUnion(group_1,group_2,group_3);
              auto gshop8  = EB2::makeShop(all);
              EB2::Build(gshop8, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           case 9:
              {
              auto group_1 = EB2::makeUnion(obstacles[ob_id[0]],obstacles[ob_id[1]],obstacles[ob_id[2]]);
              auto group_2 = EB2::makeUnion(obstacles[ob_id[3]],obstacles[ob_id[4]],obstacles[ob_id[5]]);
              auto group_3 = EB2::makeUnion(obstacles[ob_id[6]],obstacles[ob_id[7]],obstacles[ob_id[8]]);
              auto all     = EB2::makeUnion(group_1,group_2,group_3);
              auto gshop9  = EB2::makeShop(all);
              EB2::Build(gshop9, geom, required_coarsening_level, max_coarsening_level);
              break;
              }

           default:;
        }
   
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom);
   
        // options are basic, volume, or full
        EBSupport ebs = EBSupport::full;
  
        // number of ghost cells for each of the 3 EBSupport types
        Vector<int> ng_ebs = {2,2,2};
 
        // This object provides access to the EB database in the format of basic AMReX objects
        // such as BaseFab, FArrayBox, FabArray, and MultiFab
        EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, ebs);
        
//      const FabFactory<FArrayBox>& test_factory = (num_obstacles > 0) ? 
//          EBFArrayBoxFactory(eb_level, geom, grids, dmap, ng_ebs, ebs): FArrayBoxFactory();

	// Velocities and Beta are face-centered
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].define (amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 1, MFInfo(), factory);
            beta[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 0, MFInfo(), factory);
            beta[idim].setVal(1.0);
        }

        // store plotfile variables; velocity and processor id
        plotfile_mf.define(grids, dmap, AMREX_SPACEDIM+1, 0, MFInfo(), factory);

        // Initialize Particles
        MyParticleContainer MyPC(geom, dmap, grids);
        MyPC.InitParticles(particle_file,zlen);

        // set initial velocity to u=(1,0,0)
        AMREX_D_TERM(vel[0].setVal(1.0);,
                     vel[1].setVal(0.0);,
                     vel[2].setVal(0.0););

        LPInfo lp_info;

        // If we want to use hypre to solve the full problem we need to not coarsen inside AMReX
        if (use_hypre) 
            lp_info.setMaxCoarseningLevel(0);

        MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // mac velocity
			     MLMG::Location::FaceCenter,       // velocity located on face centers
                             {amrex::GetArrOfConstPtrs(beta)}, // beta
			     MLMG::Location::FaceCenter,       // beta located on face centers
			     MLMG::Location::CellCenter,       // location of velocity divergence
                             {geom},
                             lp_info);                          // structure for passing info to the operator

        // Set bottom-solver to use hypre instead of native BiCGStab 
        if (use_hypre) 
           macproj.getMLMG().setBottomSolver(MLMG::BottomSolver::hypre);

        macproj.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                          LinOpBCType::Periodic,
                                          LinOpBCType::Periodic)},
            {AMREX_D_DECL(LinOpBCType::Dirichlet,
                          LinOpBCType::Periodic,
                          LinOpBCType::Periodic)});

        macproj.setVerbose(mg_verbose);
        macproj.getMLMG().setBottomVerbose(cg_verbose);

        Real reltol = 1.e-8;
        Real abstol = 1.e-12;

        amrex::Print() << " \n********************************************************************" << std::endl; 
        amrex::Print() << " First let's project the initial velocity to find " << std::endl;
        amrex::Print() << "   the flow field around the obstacles ... " << std::endl;
        amrex::Print() << "******************************************************************** \n" << std::endl; 

        macproj.project(reltol, abstol);

        amrex::Print() << " \n********************************************************************" << std::endl; 
        amrex::Print() << " Done!  Now let's advect the particles ... " << std::endl;
        amrex::Print() << "******************************************************************** \n" << std::endl; 

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            vel[idim].FillBoundary(geom.periodicity());
        }

        // copy processor id into plotfile_mf
        for (MFIter mfi(vel[0]); mfi.isValid(); ++mfi)
            plotfile_mf[mfi].setVal(ParallelDescriptor::MyProc());

        // copy velocity into plotfile
        average_face_to_cellcenter(plotfile_mf,1,amrex::GetArrOfConstPtrs(vel));

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

        // get max velocity on grid vx_max, vy_max
        Real vx_max = vel[0].max(0, 0);
        Real vy_max = vel[1].max(0, 0);
        // calculate dt_x = dx/vx_max and dt_y
        Real dt_x = geom.CellSize(0)/vx_max;
        Real dt_y = geom.CellSize(1)/vy_max;
        // dt_limit = min(dt_x, dt_y)
        Real dt_limit = std::min(dt_x, dt_y);
        time_step = 0.1 * dt_limit;

        if (AMREX_SPACEDIM > 2)
           vel[2].setVal(0.0);

        Real x = MyPC.FindWinner(0);

        amrex::Print() << " \n********************************************************************" << std::endl; 
        amrex::Print() << "At time = 0, the lead particles is at x = " << x << std::endl;
        amrex::Print() << "********************************************************************\n " << std::endl; 

        Real time = 0.0;
        for (int i = 0; i < max_steps; i++)
        {
            if (time < max_time) {
                time_step = std::min(time_step, max_time - time);

                // Step Particles
                MyPC.AdvectWithUmac(vel.data(), 0, time_step);

                MyPC.Redistribute();

                // Write to a plotfile
                if (i%plot_int == 0)
                   write_plotfile(i, geom, plotfile_mf, MyPC, write_ascii);

                // Increment time
                time += time_step;

                // Find the maximum particle position "x" to determine the winning particle
                using ParticleType = MyParticleContainer::ParticleType;

                // This finds the particle with the maximum "x"
                x = MyPC.FindWinner(0);

                if (i%100 == 0)
                {
                   amrex::Print() << "Timestep " << i << ", Time = " << time << " and leading particle now at " << x << std::endl;
                }

                if (x > 1.99) 
                {
                   amrex::Print() << " \n********************************************************************" << std::endl; 
                   amrex::Print() << "We have a winner...and the winning time is " << time << std::endl;
                   amrex::Print() << "********************************************************************\n " << std::endl; 
                   write_plotfile(i, geom, plotfile_mf, MyPC, write_ascii);
                   break;
                }

            } else {
                // Write to a plotfile
                write_plotfile(i, geom, plotfile_mf, MyPC, write_ascii);
                break;
            }
        }

        Print() << "Writing EB surface" << std::endl;
        WriteEBSurface (grids, dmap, geom, &factory);
    }
  
    Real stop_time = amrex::second() - strt_time;
    amrex::Print() << "Total run time " << stop_time << std::endl;

    amrex::Finalize();
}
