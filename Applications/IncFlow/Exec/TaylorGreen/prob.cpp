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
    std::vector<Real> prob_hi{AMREX_D_DECL(1.0_rt, 1.0_rt, 0.25_rt)};
    ppg.queryAdd("prob_hi", prob_hi);
    std::vector<int> is_periodic(AMREX_SPACEDIM, 1);
    ppg.queryAdd("is_periodic", is_periodic);

    ParmParse ppa("amr");
    std::vector<int> n_cell{AMREX_D_DECL(32, 32, 8)};
    ppa.queryAdd("n_cell", n_cell);
    int max_level = 0;
    ppa.queryAdd("max_level", max_level);
}

int prob_check (IncFlow const& incflow)
{
    // Level 0 covers the domain and holds the average of finer levels.
    MultiFab const& vel = incflow.Velocity(0);
    auto const geomdata = incflow.Geom(0).data();
    Real const decay = std::exp(-8.0_rt * Math::pi<Real>() * Math::pi<Real>()
                                * incflow.Viscosity() * incflow.Time());

    MultiFab err(vel.boxArray(), vel.DistributionMap(), AMREX_SPACEDIM, 0);
    auto const& v = vel.const_arrays();
    auto const& e = err.arrays();
    ParallelFor(err, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        // The initial velocity times the decay factor
        prob_initdata(i, j, k, e[b], geomdata);
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            e[b](i,j,k,n) = v[b](i,j,k,n) - decay * e[b](i,j,k,n);
        }
    });

    Real rel_err = 0.0_rt;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        rel_err = std::max(rel_err, err.norminf(n, 0) / decay);
    }
    amrex::Print() << "Relative max error in velocity: " << rel_err << "\n";

    return (rel_err < 0.01_rt) ? 0 : 1;
}
