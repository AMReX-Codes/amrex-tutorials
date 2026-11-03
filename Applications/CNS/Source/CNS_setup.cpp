#include "CNS.H"
#include "CNS_derive.H"

#include <AMReX_LevelBld.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

namespace {

// Map the physical boundary types (cns.lo_bc, cns.hi_bc) to the BCType of
// each state component.  Index: Interior, Inflow, Outflow, Symmetry,
// SlipWall, NoSlipWall.
constexpr int scalar_bc[] = {
    BCType::int_dir, BCType::ext_dir, BCType::foextrap,
    BCType::reflect_even, BCType::reflect_even, BCType::reflect_even
};
constexpr int norm_vel_bc[] = {
    BCType::int_dir, BCType::ext_dir, BCType::foextrap,
    BCType::reflect_odd, BCType::reflect_odd, BCType::reflect_odd
};
constexpr int tang_vel_bc[] = {
    BCType::int_dir, BCType::ext_dir, BCType::foextrap,
    BCType::reflect_even, BCType::reflect_even, BCType::reflect_odd
};

BCRec component_bc (int comp, BCRec const& phys_bc)
{
    BCRec bc;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        const int* table = scalar_bc;
        if (comp >= UMX && comp <= UMZ) {
            table = (comp == UMX + dir) ? norm_vel_bc : tang_vel_bc;
        }
        bc.setLo(dir, table[phys_bc.lo(dir)]);
        bc.setHi(dir, table[phys_bc.hi(dir)]);
    }
    return bc;
}

// Fills ghost cells outside the domain with BCType::ext_dir (inflow).  The
// other BC types (periodic, extrapolation, reflection) are handled by
// GpuBndryFuncFab.  None of the problems here use inflow.
struct CnsFillExtDir
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& /*iv*/, Array4<Real> const& /*dest*/,
                     const int /*dcomp*/, const int /*numcomp*/,
                     GeometryData const& /*geom*/, const Real /*time*/,
                     const BCRec* /*bcr*/, const int /*bcomp*/,
                     const int /*orig_comp*/) const
    {}
};

}

void
CNS::variableSetUp ()
{
    read_params();
    cns_prob_read_params(prob_parm);

    // One cell-centered StateData with NCONS components and NUM_GROW ghost
    // cells, interpolated conservatively from coarse to fine.
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, NUM_GROW, NCONS,
                           &cell_cons_interp, /*extrap=*/ false,
                           /*store_in_checkpoint=*/ true);

    Vector<std::string> names{"density", "xmom", "ymom", "zmom", "rho_E", "rho_e", "Temp"};
    Vector<BCRec> bcs(NCONS);
    for (int n = 0; n < NCONS; ++n) {
        bcs[n] = component_bc(n, phys_bc);
    }

    StateDescriptor::BndryFunc bndryfunc(GpuBndryFuncFab<CnsFillExtDir>(CnsFillExtDir{}));
    bndryfunc.setRunOnGPU(true);  // the BC function launches GPU kernels itself

    desc_lst.setComponent(State_Type, URHO, names, bcs, bndryfunc);

    // Derived quantities, computed on demand for plotfiles
    // (amr.derive_plot_vars) and tagging.
    derive_lst.add("pressure", IndexType::TheCellType(), 1, cns_derpres,
                   &DeriveRec::TheSameBox);
    derive_lst.addComponent("pressure", desc_lst, State_Type, UEINT, 1);

    derive_lst.add("velocity", IndexType::TheCellType(), AMREX_SPACEDIM,
                   {AMREX_D_DECL("x_velocity","y_velocity","z_velocity")}, cns_dervel,
                   &DeriveRec::TheSameBox);
    derive_lst.addComponent("velocity", desc_lst, State_Type, URHO, 1+AMREX_SPACEDIM);
}

void
CNS::variableCleanUp ()
{
    desc_lst.clear();
    derive_lst.clear();
    error_tags.clear();
}

// LevelBld tells Amr how to build a CNS level.
class CNSBld
    :
    public LevelBld
{
    void variableSetUp () override { CNS::variableSetUp(); }
    void variableCleanUp () override { CNS::variableCleanUp(); }
    AmrLevel* operator() () override { return new CNS; }
    AmrLevel* operator() (Amr& papa, int lev, const Geometry& level_geom,
                          const BoxArray& ba, const DistributionMapping& dm,
                          Real time) override
    {
        return new CNS(papa, lev, level_geom, ba, dm, time);
    }
};

LevelBld*
getLevelBld ()
{
    static CNSBld cns_bld;
    return &cns_bld;
}
