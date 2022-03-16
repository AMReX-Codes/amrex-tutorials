#include "diagnostics/FieldIO.H"
#include "particles/ElectrostaticParticleContainer.H"

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_InterpBndryData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Vector.H>

#include <memory>
#include <string>

using namespace amrex;

void WritePlotFile (const Vector<const MultiFab*>& rhs,
                    const Vector<const MultiFab*>& phi,
                    const Vector<std::array<const MultiFab*, AMREX_SPACEDIM> >& E,
                    const ElectrostaticParticleContainer& pc,
                    const Vector<Geometry>& geom,
                    int nstep)
{
    int num_output_comp = 2 + AMREX_SPACEDIM;
    int num_levels = rhs.size();
    IntVect cc_flag = IntVect::TheZeroVector();
    Vector<std::unique_ptr<MultiFab> > output_cc(num_levels);
    for (int lev = 0; lev < num_levels; ++lev) {
        const BoxArray& nodal_ba = rhs[lev]->boxArray();
        output_cc[lev].reset(new MultiFab(amrex::convert(nodal_ba, cc_flag),
                                          rhs[lev]->DistributionMap(), num_output_comp, 0));
        amrex::average_node_to_cellcenter(*output_cc[lev], 0, *rhs[lev],  0, 1);
        amrex::average_node_to_cellcenter(*output_cc[lev], 1, *phi[lev],  0, 1);
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            amrex::average_node_to_cellcenter(*output_cc[lev], 2+i, *E[lev][i], 0, 1);
        }
    }

    Vector<std::string> varnames;
    varnames.push_back("rhs");
    varnames.push_back("phi");
    varnames.push_back("Ex");
    varnames.push_back("Ey");

    Vector<std::string> particle_varnames;
    particle_varnames.push_back("weight");
    particle_varnames.push_back("vx");
    particle_varnames.push_back("vy");
    particle_varnames.push_back("Ex");
    particle_varnames.push_back("Ey");

    Vector<int> level_steps;
    level_steps.push_back(0);
    level_steps.push_back(0);

    int output_levs = num_levels;

    Vector<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputRR[lev] = IntVect(AMREX_D_DECL(2, 2, 2));
    }

    const std::string& pltfile = amrex::Concatenate("plt", nstep, 5);
    WriteMultiLevelPlotfile(pltfile, output_levs, GetVecOfConstPtrs(output_cc),
                            varnames, geom, 0.0, level_steps, outputRR);

    pc.Checkpoint(pltfile, "particle0", true, particle_varnames);
}
