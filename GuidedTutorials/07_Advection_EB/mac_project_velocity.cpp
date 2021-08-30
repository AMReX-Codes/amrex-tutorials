#include <AMReX_MacProjector.H>

void mac_project_velocity(amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>& vel, const amrex::Geometry& geom, int use_hypre)
{
    using namespace amrex; 

    LPInfo lp_info;

    // If we want to use hypre to solve the full problem we need to not coarsen inside AMReX
    if (use_hypre) 
        lp_info.setMaxCoarseningLevel(0);

    Array<MultiFab,AMREX_SPACEDIM> beta;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        beta[idim].define(vel[idim].boxArray(), vel[idim].DistributionMap(), 1, 0, MFInfo(), vel[idim].Factory());
        beta[idim].setVal(1.0);
    }

    MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // mac velocity
                         MLMG::Location::FaceCenter,       // velocity located on face centers
                         {amrex::GetArrOfConstPtrs(beta)}, // beta
                         MLMG::Location::FaceCenter,       // beta located on face centers
		         MLMG::Location::CellCenter,       // location of mac_phi
                         {geom},
                         lp_info);                          // structure for passing info to the operator

        // Set bottom-solver to use hypre instead of native BiCGStab 
        if (use_hypre) 
           macproj.getMLMG().setBottomSolver(MLMG::BottomSolver::hypre);

        macproj.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                          LinOpBCType::Neumann,
                                          LinOpBCType::Periodic)},
                            {AMREX_D_DECL(LinOpBCType::Neumann,
                                          LinOpBCType::Neumann,
                                          LinOpBCType::Periodic)});

        Real reltol = 1.e-8;
        Real abstol = 1.e-12;

        macproj.project(reltol, abstol);
}
