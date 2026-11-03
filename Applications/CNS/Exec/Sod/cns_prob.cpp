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
    Real stop_time = Real(0.2);
    pp.queryAdd("stop_time", stop_time);

    // Unit domain, periodic except in x
    ParmParse ppg("geometry");
    if (!ppg.contains("prob_extent")) {
        std::vector<Real> prob_hi(AMREX_SPACEDIM, Real(1.0));
        ppg.queryAdd("prob_hi", prob_hi);
    }
    std::vector<int> is_periodic{AMREX_D_DECL(0,1,1)};
    ppg.queryAdd("is_periodic", is_periodic);

    ParmParse ppa("amr");
    std::vector<int> n_cell{AMREX_D_DECL(64,8,8)};
    ppa.queryAdd("n_cell", n_cell);
    int max_level = 1;
    ppa.queryAdd("max_level", max_level);
    std::vector<std::string> derive_plot_vars{"pressure", "velocity"};
    ppa.queryAdd("derive_plot_vars", derive_plot_vars);

    // Slip walls at the ends of the tube, so that no mass leaves the domain
    ParmParse ppc("cns");
    std::vector<int> lo_bc{AMREX_D_DECL(4,0,0)};
    std::vector<int> hi_bc{AMREX_D_DECL(4,0,0)};
    ppc.queryAdd("lo_bc", lo_bc);
    ppc.queryAdd("hi_bc", hi_bc);
    int refine_max_dengrad_lev = 10;
    ppc.queryAdd("refine_max_dengrad_lev", refine_max_dengrad_lev);
    Real refine_dengrad = Real(0.05);
    ppc.queryAdd("refine_dengrad", refine_dengrad);
}

void cns_prob_read_params (ProbParm& prob_parm)
{
    ParmParse pp("prob");
    pp.query("p_l", prob_parm.p_l);
    pp.query("p_r", prob_parm.p_r);
    pp.query("rho_l", prob_parm.rho_l);
    pp.query("rho_r", prob_parm.rho_r);
    pp.query("u_l", prob_parm.u_l);
    pp.query("u_r", prob_parm.u_r);
    pp.query("x0", prob_parm.x0);
}
