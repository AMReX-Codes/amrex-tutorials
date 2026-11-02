#include "IncFlow.H"
#include "prob.H"

#include <AMReX.H>
#include <AMReX_Print.H>

#include <functional>

using namespace amrex;

int main (int argc, char* argv[])
{
    // prob_add_defaults adds problem defaults to ParmParse after it has read
    // the inputs file and command line.
    amrex::Initialize(argc, argv, std::function<void()>(prob_add_defaults));
    int status = 0;
    {
        // Reads amr.* and geometry.* and builds the Geometry of all levels.
        IncFlow incflow;

        incflow.InitData();
        incflow.Evolve();

        // Compare with the exact solution, or check that the result is sane.
        status = prob_check(incflow);
        amrex::Print() << (status == 0 ? "PASSED\n" : "FAILED\n");
    }
    amrex::Finalize();
    return status;
}
