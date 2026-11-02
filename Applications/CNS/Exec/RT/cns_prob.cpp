#include "CNS.H"

#include <AMReX_ParmParse.H>

#include <string>
#include <vector>

using namespace amrex;

void cns_prob_defaults ()
{
    // queryAdd keeps a value given in the inputs file or on the command line
    // and adds the default otherwise.
    ParmParse pp;
    Real stop_time = Real(1.9);
    pp.queryAdd("stop_time", stop_time);

    // A tall box, periodic in the horizontal directions
    ParmParse ppg("geometry");
    std::vector<Real> prob_hi{AMREX_D_DECL(Real(0.5),Real(1.0),Real(1.0))};
    std::vector<int> is_periodic{AMREX_D_DECL(1,0,0)};
#if (AMREX_SPACEDIM == 3)
    prob_hi[1] = Real(0.5);
    is_periodic[1] = 1;
#endif
    if (!ppg.contains("prob_extent")) {
        ppg.queryAdd("prob_hi", prob_hi);
    }
    ppg.queryAdd("is_periodic", is_periodic);

    ParmParse ppa("amr");
    std::vector<int> n_cell{AMREX_D_DECL(32,64,64)};
#if (AMREX_SPACEDIM == 3)
    n_cell = {16, 16, 32};
#endif
    ppa.queryAdd("n_cell", n_cell);
    int max_level = 1;
    ppa.queryAdd("max_level", max_level);

    // Slip walls at the bottom and top
    ParmParse ppc("cns");
    std::vector<int> lo_bc(AMREX_SPACEDIM, 0);
    std::vector<int> hi_bc(AMREX_SPACEDIM, 0);
    lo_bc.back() = 4;
    hi_bc.back() = 4;
    ppc.queryAdd("lo_bc", lo_bc);
    ppc.queryAdd("hi_bc", hi_bc);
    Real gravity = Real(-1.0);
    ppc.queryAdd("gravity", gravity);
    Real eos_gamma = Real(5.0)/Real(3.0);
    ppc.queryAdd("eos_gamma", eos_gamma);
    int refine_max_dengrad_lev = 10;
    ppc.queryAdd("refine_max_dengrad_lev", refine_max_dengrad_lev);
    Real refine_dengrad = Real(0.01);
    ppc.queryAdd("refine_dengrad", refine_dengrad);
}

void cns_prob_read_params (ProbParm& prob_parm)
{
    ParmParse pp("prob");
    pp.query("rho_1", prob_parm.rho_1);
    pp.query("rho_2", prob_parm.rho_2);
    pp.query("p0_base", prob_parm.p0_base);
}
