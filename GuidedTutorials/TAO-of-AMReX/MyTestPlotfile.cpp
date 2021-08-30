
#include "MyTest.H"
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void
MyTest::writePlotfile ()
{
    const int ncomp = 6;
    Vector<std::string> varname = {"solution", "rhs", "exact_solution", "error", "adjoint", "adjoint_rhs"};

    const int nlevels = max_level+1;

    MultiFab plotmf;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        plotmf.define(grids, dmap, ncomp, 0);
        MultiFab::Copy(plotmf, solution      , 0, 0, 1, 0);
        MultiFab::Copy(plotmf, rhs           , 0, 1, 1, 0);
        MultiFab::Copy(plotmf, exact_solution, 0, 2, 1, 0);
        MultiFab::Copy(plotmf, solution      , 0, 3, 1, 0);
        MultiFab::Copy(plotmf, adjoint       , 0, 4, 1, 0);
        MultiFab::Copy(plotmf, adjoint_rhs   , 0, 5, 1, 0);
        MultiFab::Subtract(plotmf, plotmf, 2, 3, 1, 0); // error = soln - exact
    }

    std::string filename = get_iteration_filename("plt");
    WriteSingleLevelPlotfile(filename, plotmf, varname, geom, 0.0, 1);
}

void
MyTest::writeMinimalPlotfile ()
{
    const int ncomp = 2;
    Vector<std::string> varname = {"solution", "rhs"};

    const int nlevels = max_level+1;

    MultiFab plotmf;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        plotmf.define(grids, dmap, ncomp, 0);
        MultiFab::Copy(plotmf, solution      , 0, 0, 1, 0);
        MultiFab::Copy(plotmf, rhs           , 0, 1, 1, 0);
    }

    std::string filename = get_iteration_filename("plt");
    WriteSingleLevelPlotfile(filename, plotmf, varname, geom, 0.0, 1);
}

