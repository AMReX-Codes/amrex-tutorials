#include <AMReX.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <cmath>
#include <type_traits>

using namespace amrex;

namespace {

// Initial condition: a Gaussian that straddles the coarse/fine interface.
void init_phi (MultiFab& phi, Geometry const& geom)
{
    auto const& a = phi.arrays();
    auto const problo = geom.ProbLoArray();
    auto const dx = geom.CellSizeArray();
    ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        IntVect iv(AMREX_D_DECL(i,j,k));
        Real r2 = 0.0_rt;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real x = problo[d] + (iv[d] + 0.5_rt) * dx[d];
            r2 += (x - 0.35_rt) * (x - 0.35_rt);
        }
        a[b](i,j,k) = std::exp(-r2 / 0.01_rt);
    });
}

// Diffusive flux F = -D grad(phi) on every face of the level.
void compute_flux (MultiFab const& phi, Array<MultiFab,AMREX_SPACEDIM>& flux,
                   Geometry const& geom, Real diff_coef)
{
    auto const& p = phi.const_arrays();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        auto const& f = flux[d].arrays();
        Real const fac = -diff_coef * geom.InvCellSize(d);
        Dim3 const e = IntVect::TheDimensionVector(d).dim3();
        ParallelFor(flux[d], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            f[b](i,j,k) = fac * (p[b](i,j,k) - p[b](i-e.x,j-e.y,k-e.z));
        });
    }
}

// Conservative update phi -= dt * div(F).
void update (MultiFab& phi, Array<MultiFab,AMREX_SPACEDIM> const& flux,
             Geometry const& geom, Real dt)
{
    auto const& p = phi.arrays();
    AMREX_D_TERM(auto const& fx = flux[0].const_arrays();,
                 auto const& fy = flux[1].const_arrays();,
                 auto const& fz = flux[2].const_arrays();)
    auto const dxinv = geom.InvCellSizeArray();
    ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        p[b](i,j,k) -= dt * (AMREX_D_TERM(  dxinv[0] * (fx[b](i+1,j,k) - fx[b](i,j,k)),
                                          + dxinv[1] * (fy[b](i,j+1,k) - fy[b](i,j,k)),
                                          + dxinv[2] * (fz[b](i,j,k+1) - fz[b](i,j,k))));
    });
}

Real total_mass (MultiFab const& phi_crse, Geometry const& geom)
{
    // Level 0 covers the domain, and average_down keeps it consistent with level 1.
    return phi_crse.sum(0) * AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int nsteps = 100;
        int do_reflux = 1;
        int plot = 0;
        Real diff_coef = 1.0_rt;
        Real cfl = 0.9_rt;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nsteps", nsteps);
            pp.query("do_reflux", do_reflux);
            pp.query("plot", plot);
            pp.query("diff_coef", diff_coef);
            pp.query("cfl", cfl);
        }
        AMREX_ALWAYS_ASSERT(n_cell % 4 == 0);

        // Two levels on the unit square/cube with non-periodic walls.
        constexpr int nlevs = 2;
        IntVect const ratio(2);
        RealBox const rb({AMREX_D_DECL(0.0_rt,0.0_rt,0.0_rt)},
                         {AMREX_D_DECL(1.0_rt,1.0_rt,1.0_rt)});
        Array<int,AMREX_SPACEDIM> const is_periodic{AMREX_D_DECL(0,0,0)};
        Box const domain0(IntVect(0), IntVect(n_cell-1));

        Vector<Geometry> geom(nlevs);
        geom[0].define(domain0, rb, CoordSys::cartesian, is_periodic);
        geom[1].define(amrex::refine(domain0, ratio), rb, CoordSys::cartesian, is_periodic);

        // Level 0 covers the domain. Level 1 refines the middle half of it.
        Vector<BoxArray> grids(nlevs);
        grids[0].define(domain0);
        grids[1].define(amrex::refine(Box(IntVect(n_cell/4), IntVect(3*n_cell/4-1)), ratio));
        Vector<DistributionMapping> dmap(nlevs);
        for (int lev = 0; lev < nlevs; ++lev) {
            grids[lev].maxSize(max_grid_size);
            dmap[lev].define(grids[lev]);
        }

        Vector<MultiFab> phi(nlevs);
        for (int lev = 0; lev < nlevs; ++lev) {
            phi[lev].define(grids[lev], dmap[lev], 1, 0);
            init_phi(phi[lev], geom[lev]);
        }
        amrex::average_down(phi[1], phi[0], 0, 1, ratio);

        // Homogeneous Neumann walls: fill ghost cells by first-order extrapolation.
        Vector<BCRec> bcs(1);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            bcs[0].setLo(d, BCType::foextrap);
            bcs[0].setHi(d, BCType::foextrap);
        }
        // GpuBndryFuncFab fills the BCType-based ghost cells; FabFillNoOp is the
        // user hook for external Dirichlet (BCType::ext_dir), unused here.
        using BndryFunc = GpuBndryFuncFab<FabFillNoOp>;
        Vector<PhysBCFunct<BndryFunc>> physbc;
        for (int lev = 0; lev < nlevs; ++lev) {
            physbc.emplace_back(geom[lev], bcs, BndryFunc(FabFillNoOp{}));
        }

        // The flux register lives on the fine level's coarse/fine boundary.
        FluxRegister flux_reg(grids[1], dmap[1], ratio, 1, 1);

        // No subcycling: both levels use the stable time step of the fine level.
        Real const dxf = geom[1].CellSize(0);
        Real const dt = cfl * dxf * dxf / (2.0_rt * AMREX_SPACEDIM * diff_coef);

        Real const mass0 = total_mass(phi[0], geom[0]);
        amrex::Print() << "Initial mass: " << mass0 << "\n";

        Vector<MultiFab> phi_ghost(nlevs);
        Vector<Array<MultiFab,AMREX_SPACEDIM>> flux(nlevs);
        for (int lev = 0; lev < nlevs; ++lev) {
            phi_ghost[lev].define(grids[lev], dmap[lev], 1, 1);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                flux[lev][d].define(amrex::convert(grids[lev], IntVect::TheDimensionVector(d)),
                                    dmap[lev], 1, 0);
            }
        }

        Real time = 0.0_rt;
        for (int step = 1; step <= nsteps; ++step)
        {
            // Level 0 ghost cells: periodic/interior copies plus physical boundaries.
            FillPatchSingleLevel(phi_ghost[0], time, {&phi[0]}, {time}, 0, 0, 1,
                                 geom[0], physbc[0], 0);

            // Level 1 ghost cells: copy from level 1 where possible, otherwise
            // interpolate conservatively from level 0.
            FillPatchTwoLevels(phi_ghost[1], time, {&phi[0]}, {time}, {&phi[1]}, {time},
                               0, 0, 1, geom[0], geom[1], physbc[0], 0, physbc[1], 0,
                               ratio, &cell_cons_interp, bcs, 0);

            for (int lev = 0; lev < nlevs; ++lev) {
                compute_flux(phi_ghost[lev], flux[lev], geom[lev], diff_coef);
                update(phi[lev], flux[lev], geom[lev], dt);
            }

            if (do_reflux) {
                // Register holds dt * area * (sum of fine fluxes - coarse flux) on
                // the coarse faces of the coarse/fine boundary.
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    Real const area_c = AMREX_D_TERM(1.0_rt, *geom[0].CellSize((d+1)%AMREX_SPACEDIM),
                                                             *geom[0].CellSize((d+2)%AMREX_SPACEDIM));
                    Real const area_f = AMREX_D_TERM(1.0_rt, *geom[1].CellSize((d+1)%AMREX_SPACEDIM),
                                                             *geom[1].CellSize((d+2)%AMREX_SPACEDIM));
                    flux_reg.CrseInit(flux[0][d], d, 0, 0, 1, -dt*area_c);
                    flux_reg.FineAdd (flux[1][d], d, 0, 0, 1,  dt*area_f);
                }
                // Correct the coarse cells next to the fine level.
                flux_reg.Reflux(phi[0], 1.0_rt, 0, 0, 1, geom[0]);
            }

            // Covered coarse cells get the average of the fine cells.
            amrex::average_down(phi[1], phi[0], 0, 1, ratio);

            time += dt;
        }

        Real const mass = total_mass(phi[0], geom[0]);
        Real const rel_err = std::abs(mass - mass0) / mass0;
        amrex::Print() << "Final mass:   " << mass << " after " << nsteps
                       << " steps, time = " << time << "\n"
                       << "Relative change in mass: " << rel_err << "\n";

        if (plot) {
            WriteMultiLevelPlotfile("plt_static", nlevs, GetVecOfConstPtrs(phi), {"phi"},
                                    geom, time, Vector<int>(nlevs, nsteps),
                                    Vector<IntVect>(nlevs-1, ratio));
        }

        Real const tol = std::is_same_v<Real,float> ? 1.e-5_rt : 1.e-12_rt;
        if (rel_err > tol) {
            amrex::Print() << "FAILED: mass is not conserved\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
