#include "IncFlow.H"
#include "prob.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <type_traits>
#include <vector>

using namespace amrex;

void prob_add_defaults ()
{
    // queryAdd keeps a value given in the inputs file or on the command line
    // and adds the default otherwise.
    ParmParse ppg("geometry");
    std::vector<Real> prob_hi{AMREX_D_DECL(1.0_rt, 1.0_rt, 0.25_rt)};
    ppg.queryAdd("prob_hi", prob_hi);
    std::vector<int> is_periodic(AMREX_SPACEDIM, 1);
    ppg.queryAdd("is_periodic", is_periodic);

    ParmParse ppa("amr");
    // Coarser in 3D to keep the test short
    std::vector<int> n_cell{AMREX_D_DECL((AMREX_SPACEDIM == 2) ? 64 : 32,
                                         (AMREX_SPACEDIM == 2) ? 64 : 32, 8)};
    ppa.queryAdd("n_cell", n_cell);
    int max_level = 1;
    ppa.queryAdd("max_level", max_level);
}

int prob_check (IncFlow const& incflow)
{
    // There is no exact solution. On a periodic domain, the conservative
    // scheme keeps the mean velocity (zero here) to round-off, also across
    // coarse/fine boundaries.
    MultiFab const& vel = incflow.Velocity(0);
    Real umax = 0.0_rt;
    Real mean = 0.0_rt;
    auto const npts = Real(vel.boxArray().numPts());
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        umax = std::max(umax, vel.norminf(n, 0));
        mean = std::max(mean, std::abs(vel.sum(n)) / npts);
    }
    amrex::Print() << "Max |velocity| = " << umax << ", max |mean velocity| = " << mean << "\n";

    Real const tol = std::is_same_v<Real,float> ? 1.e-5_rt : 1.e-12_rt;
    return (std::isfinite(umax) && mean < tol) ? 0 : 1;
}
