#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <type_traits>

using namespace amrex;

namespace {

// Exact solution: zero on the domain boundary of [0,1]^D
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_solution (AMREX_D_DECL(Real x, Real y, Real z))
{
    constexpr Real pi = Math::pi<Real>();
    return AMREX_D_TERM(std::sin(pi*x), *std::sin(pi*y), *std::sin(pi*z));
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int verbose = 1;
        Real radius = 0.25_rt;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("verbose", verbose);
            pp.query("radius", radius);
        }

        const Box domain(IntVect(0), IntVect(n_cell-1));
        const RealBox rb({AMREX_D_DECL(0.0_rt,0.0_rt,0.0_rt)},
                         {AMREX_D_DECL(1.0_rt,1.0_rt,1.0_rt)});
        const Geometry geom(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
        BoxArray grids(domain);
        grids.maxSize(max_grid_size);
        const DistributionMapping dmap(grids);

        // A spherical obstacle with the fluid outside.  Multigrid needs
        // the EB on coarsened levels, so allow as many as possible.
        EB2::SphereIF sphere(radius, {AMREX_D_DECL(0.5_rt,0.5_rt,0.5_rt)}, false);
        const int required_coarsening_level = 0;
        const int max_coarsening_level = 30;
        EB2::Build(EB2::makeShop(sphere), geom, required_coarsening_level,
                   max_coarsening_level);

        EBFArrayBoxFactory factory(EB2::IndexSpace::top().getLevel(geom), geom,
                                   grids, dmap, {2,2,2}, EBSupport::full);
        auto const& ebdata = factory.getEBDataArrays();

        // MultiFabs built with the EB factory
        MultiFab phi   (grids, dmap, 1, 1, MFInfo(), factory);
        MultiFab rhs   (grids, dmap, 1, 0, MFInfo(), factory);
        MultiFab phi_eb(grids, dmap, 1, 0, MFInfo(), factory);

        // The right-hand side -lap(phi) = D pi^2 phi_exact is evaluated at
        // cell centroids, and the EB Dirichlet value at the EB centroids.
        // Both centroids are in [-0.5,0.5] relative to the cell center.
        const auto problo = geom.ProbLoArray();
        const auto dx = geom.CellSizeArray();
        const Real fac = AMREX_SPACEDIM * Math::pi<Real>() * Math::pi<Real>();
        auto const& rhs_a = rhs.arrays();
        auto const& phieb_a = phi_eb.arrays();
        ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
        {
            auto flag = ebdata.get<EBData_t::cellflag>(b,i,j,k);
            IntVect iv(AMREX_D_DECL(i,j,k));
            if (flag.isCovered()) {
                rhs_a[b](i,j,k) = 0;
                phieb_a[b](i,j,k) = 0;
            } else if (flag.isRegular()) {
                GpuArray<Real,AMREX_SPACEDIM> x;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    x[d] = problo[d] + (Real(iv[d])+Real(0.5))*dx[d];
                }
                rhs_a[b](i,j,k) = fac * exact_solution(AMREX_D_DECL(x[0],x[1],x[2]));
                phieb_a[b](i,j,k) = 0;
            } else {
                GpuArray<Real,AMREX_SPACEDIM> xc, xb;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    Real c = ebdata.get<EBData_t::centroid>(b,i,j,k,d);
                    Real bc = ebdata.get<EBData_t::bndrycent>(b,i,j,k,d);
                    xc[d] = problo[d] + (Real(iv[d])+Real(0.5)+c)*dx[d];
                    xb[d] = problo[d] + (Real(iv[d])+Real(0.5)+bc)*dx[d];
                }
                rhs_a[b](i,j,k) = fac * exact_solution(AMREX_D_DECL(xc[0],xc[1],xc[2]));
                phieb_a[b](i,j,k) = exact_solution(AMREX_D_DECL(xb[0],xb[1],xb[2]));
            }
        });

        // Operator: (a alpha - b div beta grad) phi = rhs with a = 0, b = 1,
        // beta = 1, i.e., -lap(phi) = rhs.
        MLEBABecLap mlebabec({geom}, {grids}, {dmap}, LPInfo(), {&factory});

        // Homogeneous Dirichlet on the domain boundary.  setDomainBC must be
        // called first; a nullptr in setLevelBC means zero boundary values.
        mlebabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                           LinOpBCType::Dirichlet,
                                           LinOpBCType::Dirichlet)},
                             {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                           LinOpBCType::Dirichlet,
                                           LinOpBCType::Dirichlet)});
        mlebabec.setLevelBC(0, nullptr);

        mlebabec.setScalars(0.0_rt, 1.0_rt);
        mlebabec.setACoeffs(0, 0.0_rt);
        mlebabec.setBCoeffs(0, 1.0_rt);

        // Inhomogeneous Dirichlet on the EB: values in phi_eb, beta = 1 there.
        // (The default EB boundary condition is homogeneous Neumann.)
        mlebabec.setEBDirichlet(0, phi_eb, 1.0_rt);

        MLMG mlmg(mlebabec);
        mlmg.setVerbose(verbose);
        const Real tol_rel = std::is_same_v<Real,float> ? 1.e-5_rt : 1.e-11_rt;
        const Real tol_abs = 0.0_rt;
        phi.setVal(0.0_rt);
        mlmg.solve({&phi}, {&rhs}, tol_rel, tol_abs);

        // Error at cell centers of the uncovered cells
        MultiFab plt(grids, dmap, 3, 0, MFInfo(), factory);
        auto const& plt_a = plt.arrays();
        auto const& phi_a = phi.const_arrays();
        ParallelFor(plt, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
        {
            AMREX_D_TERM(Real x = problo[0] + (Real(i)+Real(0.5))*dx[0];,
                         Real y = problo[1] + (Real(j)+Real(0.5))*dx[1];,
                         Real z = problo[2] + (Real(k)+Real(0.5))*dx[2];)
            bool covered = ebdata.get<EBData_t::cellflag>(b,i,j,k).isCovered();
            Real exact = covered ? Real(0) : exact_solution(AMREX_D_DECL(x,y,z));
            plt_a[b](i,j,k,0) = phi_a[b](i,j,k);
            plt_a[b](i,j,k,1) = exact;
            plt_a[b](i,j,k,2) = covered ? Real(0) : phi_a[b](i,j,k) - exact;
        });
        const Real max_error = plt.norminf(2, 0);

        // The scheme is second order, so the tolerance scales with dx^2.
        const Real tol = 10.0_rt * dx[0] * dx[0];
        amrex::Print() << "Max error: " << max_error << " (tolerance " << tol << ")\n";
        if (max_error > tol) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }

        EB_WriteSingleLevelPlotfile("plt", plt, {"phi", "exact", "error"}, geom, 0.0, 0);
    }
    amrex::Finalize();
    return status;
}
