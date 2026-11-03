#include <AMReX.H>
#include <AMReX_AlgMG.H>
#include <AMReX_Algebra.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <exception>
#include <type_traits>

using namespace amrex;
using namespace amrex::literals;

namespace {

// -div(beta grad phi) = f on the unit domain with phi = 0 on the boundary.
// beta is a checkerboard of blocks of `block` cells alternating between 1
// and `jump`. Geometric multigrid converges slowly or not at all for large
// jumps, because coarse grids cannot represent the coefficient.
struct Coef
{
    int block = 4;
    Real jump = 1.e4;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE Real operator() (IntVect const& cell) const
    {
        int parity = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { parity += cell[d] / block; }
        return (parity % 2) ? jump : Real(1);
    }
};

// Face coefficient: harmonic mean of the two cells. A face on the domain
// boundary takes the coefficient of the cell inside.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real face_coef (Coef const& coef, IntVect const& lo_cell, IntVect const& hi_cell,
                Box const& domain)
{
    bool const lo_in = domain.contains(lo_cell);
    bool const hi_in = domain.contains(hi_cell);
    if (!lo_in) { return coef(hi_cell); }
    if (!hi_in) { return coef(lo_cell); }
    Real const b0 = coef(lo_cell);
    Real const b1 = coef(hi_cell);
    return Real(2)*b0*b1/(b0+b1);
}

// The exact solution, prod_d sin(pi x_d), at the center of a cell.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_phi (IntVect const& cell, Real dx)
{
    Real r = 1;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        r *= std::sin(Math::pi<Real>() * (Real(cell[d]) + Real(0.5)) * dx);
    }
    return r;
}

constexpr bool is_float = std::is_same_v<Real,float>;

// The error is bounded by the condition number times the tolerance; this
// loose bound is enough to catch a wrong solution.
Real error_tolerance (Real reltol) { return Real(1000) * reltol; }

// Part 1: assemble the matrix as an SpMatrix and solve with AlgMG.
bool solve_standalone (Box const& domain, Real dx, Coef const& coef, Real reltol)
{
    // One row per cell, numbered like Box::index (x fastest). AlgPartition
    // splits the rows evenly among the MPI ranks.
    AlgPartition const partition(domain.numPts());
    AlgVector<Real> x(partition), b(partition), exact(partition);
    Long const nlocal = partition.numLocalRows();
    Long const row_begin = partition.globalRowBegin();
    BoxIndexer const indexer(domain);

    auto* pexact = exact.data();
    ParallelFor(nlocal, [=] AMREX_GPU_DEVICE (Long lrow)
    {
        pexact[lrow] = exact_phi(indexer.intVect(lrow+row_begin), dx);
    });

    // CSR arrays of the local rows with 2*dim+1 slots per row. On the
    // domain boundary a neighbor is missing; its slot gets column -1, and
    // CsrValid{false} tells define to drop it. The Dirichlet boundary face
    // adds 2 beta / dx^2 to the diagonal (phi = 0 on the face).
    int const nslots = 2*AMREX_SPACEDIM+1;
    Gpu::DeviceVector<Real> vals(nlocal*nslots);
    Gpu::DeviceVector<Long> cols(nlocal*nslots);
    Gpu::DeviceVector<Long> offsets(nlocal+1);
    auto* pv = vals.data();
    auto* pc = cols.data();
    auto* po = offsets.data();
    Real const dxinv2 = Real(1)/(dx*dx);
    ParallelFor(nlocal+1, [=] AMREX_GPU_DEVICE (Long lrow)
    {
        po[lrow] = lrow*nslots;
        if (lrow == nlocal) { return; }
        IntVect const cell = indexer.intVect(lrow+row_begin);
        Long* col = pc + lrow*nslots;
        Real* val = pv + lrow*nslots;
        Real diag = 0;
        int s = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            for (int side = -1; side <= 1; side += 2) {
                IntVect nbr = cell;
                nbr[d] += side;
                Real const bf = (side < 0) ? face_coef(coef, nbr, cell, domain)
                                           : face_coef(coef, cell, nbr, domain);
                if (domain.contains(nbr)) {
                    col[s] = domain.index(nbr);
                    val[s] = -bf*dxinv2;
                    diag += bf*dxinv2;
                } else {
                    col[s] = -1;
                    val[s] = 0;
                    diag += Real(2)*bf*dxinv2;
                }
                ++s;
            }
        }
        col[s] = lrow+row_begin;
        val[s] = diag;
    });
    Gpu::streamSynchronize();

    SpMatrix<Real> A;
    A.define(partition, pv, pc, nlocal*nslots, po, CsrSorted{false}, CsrValid{false});

    // The right-hand side is A times the exact solution, so the solver
    // should recover it to within the tolerance.
    SpMV(b, A, exact);

    AlgMG<Real> amg(A); // A must outlive amg
    {
        ParmParse pp("algmg");
        int verbose = 1;
        pp.query("verbose", verbose);
        amg.setVerbose(verbose);
        // V-cycles accelerated by conjugate gradients; the matrix is
        // symmetric positive definite. Also: none, bicgstab, gmres.
        AlgMGKrylovSolver krylov = AlgMGKrylovSolver::pcg;
        pp.query_enum_case_insensitive("krylov", krylov);
        amg.setKrylovSolver(krylov);
        AlgMGSmoother smoother = AlgMGSmoother::chebyshev;
        pp.query_enum_case_insensitive("smoother", smoother);
        amg.setSmoother(smoother);
    }
    amg.setRelTol(reltol);
    amg.setThrowException(true);

    x.setVal(0); // initial guess
    bool converged = true;
    try {
        amg.solve(x, b);
    } catch (std::exception const& e) {
        amrex::Print() << "AlgMG failed: " << e.what() << "\n";
        converged = false;
    }

    // Check the residual and the error independently of the solver.
    AlgVector<Real> r(partition);
    SpMV(r, A, x);
    amrex::Axpy(r, Real(-1), b);
    Real const rel_res = r.norm2() / b.norm2();
    amrex::Axpy(x, Real(-1), exact);
    Real const err = x.norminf();

    amrex::Print() << "AlgMG standalone: " << amg.numLevels() << " levels, "
                   << amg.getNumIters() << " iterations, relative residual "
                   << rel_res << ", max error " << err << "\n";
    return converged && rel_res <= Real(2)*reltol && err <= error_tolerance(reltol);
}

// Part 2: the same problem with MLABecLaplacian and MLMG, using AlgMG as
// MLMG's bottom solver and as the multigrid type.
bool solve_mlmg (Box const& domain, Real dx, Coef const& coef, Real reltol,
                 int max_grid_size, int max_iter)
{
    Geometry const geom(domain, RealBox({AMREX_D_DECL(0._rt,0._rt,0._rt)},
                                        {AMREX_D_DECL(1._rt,1._rt,1._rt)}),
                        CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping const dm(ba);

    MultiFab phi(ba, dm, 1, 1), rhs(ba, dm, 1, 0), exact(ba, dm, 1, 1);
    auto const& ea = exact.arrays();
    ParallelFor(exact, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        ea[b](i,j,k) = exact_phi(IntVect(AMREX_D_DECL(i,j,k)), dx);
    });

    // The same face coefficients as in the matrix of part 1.
    Array<MultiFab,AMREX_SPACEDIM> bcoef;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bcoef[d].define(amrex::convert(ba, IntVect::TheDimensionVector(d)), dm, 1, 0);
        auto const& ba_d = bcoef[d].arrays();
        ParallelFor(bcoef[d], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            IntVect const hi_cell(AMREX_D_DECL(i,j,k));
            IntVect const lo_cell = hi_cell - IntVect::TheDimensionVector(d);
            ba_d[b](i,j,k) = face_coef(coef, lo_cell, hi_cell, domain);
        });
    }

    // alpha*a - beta*div(b grad) with alpha = 0 and beta = 1. Max order 2
    // gives the same boundary stencil as the matrix of part 1.
    MLABecLaplacian mlabec({geom}, {ba}, {dm});
    mlabec.setMaxOrder(2);
    mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet, LinOpBCType::Dirichlet,
                                     LinOpBCType::Dirichlet)},
                       {AMREX_D_DECL(LinOpBCType::Dirichlet, LinOpBCType::Dirichlet,
                                     LinOpBCType::Dirichlet)});
    mlabec.setLevelBC(0, nullptr); // phi = 0 on the boundary
    mlabec.setScalars(0.0_rt, 1.0_rt);
    mlabec.setBCoeffs(0, GetArrOfConstPtrs(bcoef));

    MLMG mlmg(mlabec);
    mlmg.setMaxIter(max_iter);
    mlmg.setThrowException(true); // throw instead of abort if it fails
    {
        ParmParse pp("mlmg");
        int verbose = 1;
        pp.query("verbose", verbose);
        mlmg.setVerbose(verbose);
    }
    // Called when MLMG builds its AlgMG solver, to change AlgMG options.
    // Tolerances and verbosity come from MLMG's bottom-solver settings.
    mlmg.setAlgMGOptions([] (AlgMG<Real>& amg)
    {
        amg.setKrylovSolver(AlgMGKrylovSolver::pcg);
    });

    mlmg.apply({&rhs}, {&exact}); // rhs = A * exact

    struct Variant {
        char const* name;
        MultigridType type;
        BottomSolver bottom;
        bool must_converge;
    };
    Variant const variants[] = {
        // Geometric multigrid down to the default bottom solver
        {"geometric", MultigridType::geometric, BottomSolver::Default, false},
        // Geometric multigrid with AlgMG solving the bottom level. This does
        // not help here, because the geometric coarse levels above the
        // bottom are the problem.
        {"geometric + AlgMG bottom", MultigridType::geometric, BottomSolver::algmg, false},
        // AlgMG on the whole level in every MLMG iteration
        {"algebraic", MultigridType::algebraic, BottomSolver::Default, true},
        // Geometric first, AlgMG if the geometric iterations stall
        {"hybrid", MultigridType::hybrid, BottomSolver::Default, true}
    };

    bool passed = true;
    for (auto const& v : variants) {
        mlmg.setMultigridType(v.type);
        mlmg.setBottomSolver(v.bottom);
        phi.setVal(0); // initial guess
        try {
            mlmg.solve({&phi}, {&rhs}, reltol, 0.0_rt);
            MultiFab::Subtract(phi, exact, 0, 0, 1, 0);
            Real const err = phi.norminf(0, 0);
            amrex::Print() << "MLMG " << v.name << ": " << mlmg.getNumIters()
                           << " iterations, max error " << err;
            if (v.type == MultigridType::hybrid) {
                amrex::Print() << (mlmg.usedAlgMG() ? ", switched to AlgMG" : ", no switch");
            }
            amrex::Print() << "\n";
            if (err > error_tolerance(reltol)) { passed = false; }
        } catch (std::exception const&) {
            amrex::Print() << "MLMG " << v.name << ": did not converge after "
                           << mlmg.getNumIters() << " iterations\n";
            if (v.must_converge) { passed = false; }
        }
    }
    return passed;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int max_iter = 100;
        Coef coef;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("max_iter", max_iter);
            pp.query("block", coef.block);
            pp.query("jump", coef.jump);
        }
        Box const domain(IntVect(0), IntVect(n_cell-1));
        Real const dx = 1.0_rt / Real(n_cell);
        // Relative tolerance of both AlgMG and MLMG
        Real reltol = is_float ? 1.e-4_rt : 1.e-10_rt;
        ParmParse("mlmg").query("tol_rel", reltol);

        amrex::Print() << "Checkerboard of " << coef.block << "-cell blocks, jump "
                       << coef.jump << ", " << n_cell << "^" << AMREX_SPACEDIM << " cells\n";

        bool const ok1 = solve_standalone(domain, dx, coef, reltol);
        bool const ok2 = solve_mlmg(domain, dx, coef, reltol, max_grid_size, max_iter);

        if (ok1 && ok2) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
