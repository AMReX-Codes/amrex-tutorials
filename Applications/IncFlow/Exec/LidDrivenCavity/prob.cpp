#include "IncFlow.H"
#include "prob.H"

#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace amrex;

void prob_add_defaults ()
{
    // queryAdd keeps a value given in the inputs file or on the command line
    // and adds the default otherwise.
    ParmParse ppg("geometry");
    std::vector<Real> prob_hi(AMREX_SPACEDIM, 1.0_rt);
    ppg.queryAdd("prob_hi", prob_hi);
    std::vector<int> is_periodic(AMREX_SPACEDIM, 0); // walls everywhere
    ppg.queryAdd("is_periodic", is_periodic);

    ParmParse ppa("amr");
    std::vector<int> n_cell(AMREX_SPACEDIM, (AMREX_SPACEDIM == 2) ? 32 : 16);
    ppa.queryAdd("n_cell", n_cell);
    int max_level = 1;
    ppa.queryAdd("max_level", max_level);

    // The lid
    ParmParse ppl("yhi");
    std::vector<Real> lid_velocity{AMREX_D_DECL(1.0_rt, 0.0_rt, 0.0_rt)};
    ppl.queryAdd("velocity", lid_velocity);
}

int prob_check (IncFlow const& incflow)
{
    // There is no exact solution. The velocity must stay finite and, at this
    // Reynolds number, not exceed the lid velocity by much.
    MultiFab const& vel = incflow.Velocity(0);
    Real umax = 0.0_rt;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        umax = std::max(umax, vel.norminf(n, 0));
    }
    Real const ke = incflow.KineticEnergy();
    amrex::Print() << "Max |velocity| = " << umax << ", kinetic energy = " << ke << "\n";

    return (std::isfinite(ke) && ke > 0.0_rt && umax < 1.05_rt) ? 0 : 1;
}
