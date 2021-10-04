#include "ElectrostaticParticleContainer.H"
#include "Electrostatic_PIC_Util.H"
#include "Electrostatic_PIC_K.H"

#include "diagnostics/FieldIO.H"
#include "field_solver/FieldSolver.H"

#include <AMReX.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>

#include <iostream>
#include <iomanip>
#include <random>
#include <cassert>

using namespace amrex;

void run_espic ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    run_espic();

    amrex::Finalize();
}

void run_espic ()
{
    // inputs parameters
    int max_level, n_cell, max_grid_size, particle_output_int, n_buffer, max_step, is_periodic[AMREX_SPACEDIM];
    Real dt;
    {
        ParmParse pp;
        pp.get("max_level", max_level);
        pp.get("n_cell", n_cell);
        pp.get("n_buffer", n_buffer);
        pp.get("max_grid_size", max_grid_size);
        pp.get("max_step", max_step);
        pp.get("particle_output_int", particle_output_int);
        pp.get("dt", dt);
    }

    // Example assumes at most one level of refinement
    assert(max_level < 2);
    int num_levels = max_level + 1;

    // Hard code ref ratio to 2
    Vector<int> rr(num_levels-1);
    for (int lev = 1; lev < num_levels; lev++)
        rr[lev-1] = 2;

    // Hard code physical problem domain
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n,-20.0e-6);
        real_box.setHi(n, 20.0e-6);
    }

    // This sets the boundary conditions to be doubly or triply periodic
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        is_periodic[i] = 0;
    }

    IntVect dom_lo(IntVect(D_DECL(0,0,0)));
    IntVect dom_hi(IntVect(D_DECL(n_cell-1, n_cell-1, n_cell-1)));
    Box domain(dom_lo, dom_hi);

    // make Geometry for each level
    Vector<Geometry> geom(num_levels);
    geom[0].define(domain,&real_box,CoordSys::cartesian,is_periodic);
    for (int lev = 1; lev < num_levels; ++lev) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_periodic);
    }

    // make grids for each level
    Vector<BoxArray> grids(num_levels);
    grids[0].define(domain);
    if (num_levels > 1) {
        int n_fine = n_cell*rr[0];
        IntVect refined_lo(D_DECL(3*n_fine/8,3*n_fine/8,3*n_fine/8));
        IntVect refined_hi(D_DECL(5*n_fine/8-1,5*n_fine/8-1,5*n_fine/8-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        grids[1].define(refined_patch);
    }

    for (int lev = 0; lev < num_levels; lev++) {
        grids[lev].maxSize(max_grid_size);
    }

    int Ncomp  = 1;

    Vector<DistributionMapping> dm(num_levels);
    Vector<std::unique_ptr<MultiFab> > phi(num_levels);
    Vector<std::unique_ptr<MultiFab> > rhs(num_levels);
    Vector<std::array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> > eField(num_levels);

    for (int lev = 0; lev < num_levels; ++lev) {
        BoxArray nba = grids[lev];
        nba.surroundingNodes();
        dm[lev].define(grids[lev]);

        rhs[lev].reset(new MultiFab(nba, dm[lev], Ncomp, 1));
        phi[lev].reset(new MultiFab(nba, dm[lev], Ncomp, 2));

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            eField[lev][idim].reset(new MultiFab(nba, dm[lev], Ncomp, 1));
            eField[lev][idim]->setVal(0.0);
        }

        rhs[lev]->setVal(0.0);
        phi[lev]->setVal(0.0);
    }

    Vector<std::unique_ptr<FabArray<BaseFab<int> > > > masks(num_levels);
    getLevelMasks(masks, grids, dm, geom);

    Vector<std::unique_ptr<FabArray<BaseFab<int> > > > gather_masks(num_levels);
    getLevelMasks(gather_masks, grids, dm, geom, n_buffer + 1); // convert from num nodes to num cells

    ElectrostaticParticleContainer myPC(geom, dm, grids, rr);

    myPC.InitParticles();

    for (int step = 0; step <= max_step; ++step) {

        myPC.DepositCharge(rhs);

        computePhi(rhs, phi, grids, dm, geom, masks);

        computeE(eField, phi, geom);

        myPC.FieldGather(eField, gather_masks);

        if (step % particle_output_int == 0) myPC.writeParticles(step);

        //        WritePlotFile(rhs, phi, eField, myPC, geom, step);

        myPC.Evolve(eField, rhs, dt);

        myPC.Redistribute();
    }
}
