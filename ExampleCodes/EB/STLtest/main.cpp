#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB_utils.H>
#include <AMReX_EB_STL_utils.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        int nghost = 1;
        int max_grid_size=64;

        Vector<Real> plo;
        Vector<Real> phi;
        Vector<int>  ncells;
        Real dx[3];

        amrex::ParmParse pp;
        pp.getarr("prob_lo",plo);
        pp.getarr("prob_hi",phi);
        pp.getarr("ncells",ncells);
        pp.query("max_grid_size",max_grid_size);

        RealBox real_box({AMREX_D_DECL(plo[0], plo[1], plo[2])},
                {AMREX_D_DECL(phi[0], phi[1], phi[2])});

        Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

        IntVect domain_lo(AMREX_D_DECL(0,0,0));
        IntVect domain_hi(AMREX_D_DECL(ncells[0]-1,ncells[1]-1,ncells[2]-1));

        dx[0]=(phi[0]-plo[0])/ncells[0];
        dx[1]=(phi[1]-plo[1])/ncells[1];
        dx[2]=(phi[2]-plo[2])/ncells[2];

        Box domain(domain_lo, domain_hi);
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);

        Geometry geom(domain,real_box,CoordSys::cartesian,is_periodic);
        DistributionMapping dm(ba);

        int required_coarsening_level = 0; // typically the same as the max AMR level index
        int max_coarsening_level = 0;    // typically a huge number so MG coarsens as much as possible
        // build a simple geometry using the "eb2." parameters in the inputs file
        EB2::Build(geom, required_coarsening_level, max_coarsening_level);

        std::string pltfile;
        auto const& factory = makeEBFabFactory(geom, ba, dm, {1,1,1}, EBSupport::full);
        MultiFab const& vfrc = factory->getVolFrac();
        MultiFab object(ba,dm,1,0); 
        object.setVal(1.0);  
        
        MultiFab::Subtract(object, vfrc, 0,0,1,0);
        amrex::WriteMLMF("plt", {&object}, {geom});

    }

    amrex::Finalize();
}
