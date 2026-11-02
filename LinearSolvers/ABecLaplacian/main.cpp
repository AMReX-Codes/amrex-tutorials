#include <AMReX.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <type_traits>

using namespace amrex;

AMREX_ENUM(BCKind, dirichlet, neumann, periodic);

namespace {

using RealVect_t = GpuArray<Real,AMREX_SPACEDIM>;

// Manufactured problem for  a*alpha*phi - b*div(beta*grad(phi)) = rhs
// on the unit square/cube.  beta is about beta_in inside a sphere of
// radius 0.15 and 1 outside, with a smooth but steep transition.  The fine
// level covers the sphere.
struct Problem
{
    Real a = 1;
    Real b = 1;
    Real beta_in = 10;
    Real steepness = 20;

    static constexpr Real r0 = Real(0.15);

    // phi = prod_d sin(2 pi x_d + pi/4)
    [[nodiscard]] AMREX_GPU_HOST_DEVICE Real phi (RealVect_t const& x) const
    {
        Real r = 1;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { r *= std::sin(k()*x[d] + q()); }
        return r;
    }

    // d(phi)/d(x_dir)
    [[nodiscard]] AMREX_GPU_HOST_DEVICE Real dphi (RealVect_t const& x, int dir) const
    {
        Real r = 1;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            r *= (d == dir) ? k()*std::cos(k()*x[d] + q()) : std::sin(k()*x[d] + q());
        }
        return r;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE Real alpha (RealVect_t const& x) const
    {
        return Real(1) + Real(0.5)*std::sin(k()*x[0]);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE Real radius (RealVect_t const& x) const
    {
        Real r2 = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { r2 += (x[d]-Real(0.5))*(x[d]-Real(0.5)); }
        return std::sqrt(r2);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE Real beta (RealVect_t const& x) const
    {
        return Real(1) + Real(0.5)*(beta_in-Real(1))
            * (Real(1) - std::tanh(steepness*(radius(x)-r0)));
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE Real dbeta (RealVect_t const& x, int dir) const
    {
        Real const r = radius(x);
        if (r == Real(0)) { return Real(0); }
        Real const c = std::cosh(steepness*(r-r0));
        return -Real(0.5)*(beta_in-Real(1)) * steepness / (c*c) * (x[dir]-Real(0.5)) / r;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE Real rhs (RealVect_t const& x) const
    {
        // div(beta grad phi) = beta*lap(phi) + grad(beta).grad(phi)
        Real div = -Real(AMREX_SPACEDIM)*k()*k()*beta(x)*phi(x);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { div += dbeta(x,d)*dphi(x,d); }
        return a*alpha(x)*phi(x) - b*div;
    }

private:
    [[nodiscard]] AMREX_GPU_HOST_DEVICE static Real k () { return Real(2)*Math::pi<Real>(); }
    [[nodiscard]] AMREX_GPU_HOST_DEVICE static Real q () { return Real(0.25)*Math::pi<Real>(); }
};

// Physical location of index iv for a MultiFab of the given index type.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
RealVect_t location (IntVect const& iv, IndexType const& ixt,
                     RealVect_t const& problo, RealVect_t const& dx)
{
    RealVect_t x;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real const offset = ixt.cellCentered(d) ? Real(0.5) : Real(0);
        x[d] = problo[d] + (Real(iv[d]) + offset) * dx[d];
    }
    return x;
}

// Set phi to zero in valid cells (the initial guess) and store the boundary
// data in the ghost cells outside the domain: the value of phi on the face
// for Dirichlet, and d(phi)/dx_d on the face for Neumann.
void set_boundary_data (MultiFab& phi, Geometry const& geom, Problem const& prob,
                        Array<BCKind,AMREX_SPACEDIM> const& bc)
{
    Box const domain = geom.Domain();
    auto const problo = geom.ProbLoArray();
    auto const probhi = geom.ProbHiArray();
    auto const dx = geom.CellSizeArray();
    GpuArray<int,AMREX_SPACEDIM> is_neumann{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { is_neumann[d] = (bc[d] == BCKind::neumann); }

    auto const& phi_a = phi.arrays();
    ParallelFor(phi, phi.nGrowVect(), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        IntVect const iv(AMREX_D_DECL(i,j,k));
        RealVect_t x = location(iv, IndexType::TheCellType(), problo, dx);
        int nout = 0, dir = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (iv[d] < domain.smallEnd(d)) {
                x[d] = problo[d];
                ++nout; dir = d;
            } else if (iv[d] > domain.bigEnd(d)) {
                x[d] = probhi[d];
                ++nout; dir = d;
            }
        }
        if (nout == 1) {
            phi_a[b](i,j,k) = is_neumann[dir] ? prob.dphi(x,dir) : prob.phi(x);
        } else {
            phi_a[b](i,j,k) = Real(0);
        }
    });
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int max_level = 1;
        Array<BCKind,AMREX_SPACEDIM> bc{AMREX_D_DECL(BCKind::dirichlet,
                                                     BCKind::dirichlet,
                                                     BCKind::dirichlet)};
        Problem prob;
        int verbose = 1;
        int max_iter = 100;
        int linop_maxorder = 2;
        Real tol_rel = std::is_same_v<Real,float> ? Real(1.e-5) : Real(1.e-10);
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("max_level", max_level);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                pp.query_enum_case_insensitive("bc", bc[d], d);
            }
            pp.query("a", prob.a);
            pp.query("b", prob.b);
            pp.query("beta_in", prob.beta_in);

            ParmParse ppm("mlmg");
            ppm.query("verbose", verbose);
            ppm.query("max_iter", max_iter);
            ppm.query("linop_maxorder", linop_maxorder);
            ppm.query("tol_rel", tol_rel);
        }

        // Fixed AMR hierarchy: each finer level refines the middle half of
        // the level below by a factor of two.
        int const nlevels = max_level + 1;
        int const ref_ratio = 2;
        Vector<Geometry> geom(nlevels);
        Vector<BoxArray> grids(nlevels);
        Vector<DistributionMapping> dmap(nlevels);

        RealBox real_box({AMREX_D_DECL(0._rt,0._rt,0._rt)}, {AMREX_D_DECL(1._rt,1._rt,1._rt)});
        Array<int,AMREX_SPACEDIM> is_periodic{};
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { is_periodic[d] = (bc[d] == BCKind::periodic); }

        Box domain(IntVect(0), IntVect(n_cell-1));
        Box fine_region = domain;
        for (int lev = 0; lev < nlevels; ++lev) {
            geom[lev].define(domain, real_box, CoordSys::cartesian, is_periodic);
            grids[lev].define(fine_region);
            grids[lev].maxSize(max_grid_size);
            dmap[lev].define(grids[lev]);
            fine_region.grow(-fine_region.length(0)/4).refine(ref_ratio);
            domain.refine(ref_ratio);
        }

        // Cell-centered phi (with one ghost cell for the boundary data), rhs,
        // exact solution and alpha; face-centered beta.
        Vector<MultiFab> phi(nlevels), rhs(nlevels), exact(nlevels), acoef(nlevels);
        Vector<Array<MultiFab,AMREX_SPACEDIM>> bcoef(nlevels);
        for (int lev = 0; lev < nlevels; ++lev) {
            phi  [lev].define(grids[lev], dmap[lev], 1, 1);
            rhs  [lev].define(grids[lev], dmap[lev], 1, 0);
            exact[lev].define(grids[lev], dmap[lev], 1, 0);
            acoef[lev].define(grids[lev], dmap[lev], 1, 0);

            auto const problo = geom[lev].ProbLoArray();
            auto const dx = geom[lev].CellSizeArray();

            auto const& rhs_a = rhs[lev].arrays();
            auto const& exact_a = exact[lev].arrays();
            auto const& acoef_a = acoef[lev].arrays();
            ParallelFor(rhs[lev], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto x = location(IntVect(AMREX_D_DECL(i,j,k)), IndexType::TheCellType(),
                                  problo, dx);
                rhs_a  [b](i,j,k) = prob.rhs(x);
                exact_a[b](i,j,k) = prob.phi(x);
                acoef_a[b](i,j,k) = prob.alpha(x);
            });

            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                bcoef[lev][d].define(amrex::convert(grids[lev], IntVect::TheDimensionVector(d)),
                                     dmap[lev], 1, 0);
                IndexType const ixt = bcoef[lev][d].ixType();
                auto const& bcoef_a = bcoef[lev][d].arrays();
                ParallelFor(bcoef[lev][d], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
                {
                    bcoef_a[b](i,j,k) = prob.beta(location(IntVect(AMREX_D_DECL(i,j,k)), ixt,
                                                           problo, dx));
                });
            }

            set_boundary_data(phi[lev], geom[lev], prob, bc);
        }

        // The operator for all levels of the composite solve.
        MLABecLaplacian mlabec(geom, grids, dmap);

        // Order of the stencil at Dirichlet and coarse/fine boundaries.
        mlabec.setMaxOrder(linop_maxorder);

        Array<LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (bc[d] == BCKind::dirichlet) {
                lobc[d] = hibc[d] = LinOpBCType::Dirichlet;
            } else if (bc[d] == BCKind::neumann) {
                lobc[d] = hibc[d] = LinOpBCType::inhomogNeumann;
            } else {
                lobc[d] = hibc[d] = LinOpBCType::Periodic;
            }
        }
        mlabec.setDomainBC(lobc, hibc);
        for (int lev = 0; lev < nlevels; ++lev) {
            mlabec.setLevelBC(lev, &phi[lev]);
        }

        mlabec.setScalars(prob.a, prob.b);
        for (int lev = 0; lev < nlevels; ++lev) {
            mlabec.setACoeffs(lev, acoef[lev]);
            mlabec.setBCoeffs(lev, GetArrOfConstPtrs(bcoef[lev]));
        }

        MLMG mlmg(mlabec);
        mlmg.setVerbose(verbose);
        mlmg.setMaxIter(max_iter);
        mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), tol_rel, Real(0));

        // MLMG averages the fine solution down onto the covered coarse cells,
        // so we can compare with the exact solution everywhere.  The error is
        // O(dx^2).  Without Dirichlet boundaries, the mean of phi is only
        // controlled by the a*alpha term, and the error is larger.
        Real const h = geom[0].CellSize(0);
        Real const tol = Real(100) * h * h;
        for (int lev = 0; lev < nlevels; ++lev) {
            MultiFab::Subtract(exact[lev], phi[lev], 0, 0, 1, 0);
            Real const err = exact[lev].norm0();
            amrex::Print() << "Level " << lev << ": " << grids[lev].size() << " boxes, "
                           << grids[lev].numPts() << " cells, max error = " << err << "\n";
            if (err > tol) { status = 1; }
        }

        if (status) {
            amrex::Print() << "FAILED: error larger than " << tol << "\n";
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
