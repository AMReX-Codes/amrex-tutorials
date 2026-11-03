#include "AmrCoreAdv.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <type_traits>
#include <vector>

using namespace amrex;

namespace {
// Defaults for the parameters AmrCore requires, so that the code runs
// without an inputs file. Values from inputs or the command line win.
void add_default_parameters ()
{
    {
        ParmParse pp("amr");
        std::vector<int> n_cell{AMREX_D_DECL(32, 32, 8)};
        pp.queryAdd("n_cell", n_cell);
        int max_level = 2;
        pp.queryAdd("max_level", max_level);
    }
    {
        ParmParse pp("geometry");
        std::vector<Real> prob_hi{AMREX_D_DECL(1.0_rt, 1.0_rt, 0.25_rt)};
        pp.queryAdd("prob_hi", prob_hi);
        std::vector<int> is_periodic(AMREX_SPACEDIM, 1);
        pp.queryAdd("is_periodic", is_periodic);
    }
}
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        add_default_parameters();

        // Reads amr.* and geometry.* and builds the Geometry of all levels.
        AmrCoreAdv amr_core_adv;

        amr_core_adv.InitData();
        Real const mass0 = amr_core_adv.TotalMass();

        amr_core_adv.Evolve();
        Real const mass = amr_core_adv.TotalMass();

        // Reflux and conservative interpolation keep the total mass constant.
        Real const rel_err = std::abs(mass - mass0) / mass0;
        amrex::Print() << "Relative change in mass: " << rel_err << "\n";
        Real const tol = std::is_same_v<Real,float> ? 1.e-5_rt : 1.e-12_rt;
        if (rel_err > tol) {
            amrex::Print() << "FAILED: mass is not conserved\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
