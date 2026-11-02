#include "MultiCompNodeLap.H"

#include <cmath>
#include <utility>

using namespace amrex;

namespace {

// Invert the n x n row-major matrix a with Gauss-Jordan elimination.
template <int N>
GpuArray<Real,N*N> invert (GpuArray<Real,N*N> a, int n)
{
    GpuArray<Real,N*N> inv{};
    for (int i = 0; i < n; ++i) { inv[i*n+i] = 1; }
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int r = col+1; r < n; ++r) {
            if (std::abs(a[r*n+col]) > std::abs(a[piv*n+col])) { piv = r; }
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a[piv*n+col] != 0, "MultiCompNodeLap: singular C");
        for (int c = 0; c < n; ++c) {
            std::swap(a[col*n+c], a[piv*n+c]);
            std::swap(inv[col*n+c], inv[piv*n+c]);
        }
        Real const p = a[col*n+col];
        for (int c = 0; c < n; ++c) { a[col*n+c] /= p; inv[col*n+c] /= p; }
        for (int r = 0; r < n; ++r) {
            if (r == col) { continue; }
            Real const f = a[r*n+col];
            for (int c = 0; c < n; ++c) {
                a[r*n+c] -= f*a[col*n+c];
                inv[r*n+c] -= f*inv[col*n+c];
            }
        }
    }
    return inv;
}

}

MultiCompNodeLap::MultiCompNodeLap (Vector<Geometry> const& a_geom,
                                    Vector<BoxArray> const& a_grids,
                                    Vector<DistributionMapping> const& a_dmap,
                                    Vector<Real> const& a_coeff,
                                    LPInfo const& a_info)
{
    m_ncomp = static_cast<int>(std::lround(std::sqrt(Real(a_coeff.size()))));
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_ncomp >= 1 && m_ncomp <= max_ncomp &&
                                     m_ncomp*m_ncomp == static_cast<int>(a_coeff.size()),
                                     "MultiCompNodeLap: coeff must have N*N entries, N <= 4");
    for (int i = 0; i < m_ncomp*m_ncomp; ++i) { m_coeff[i] = a_coeff[i]; }
    m_coeff_inv = invert<max_ncomp>(m_coeff, m_ncomp);

    // getNComp() must be valid before define, which allocates internal data.
    MLNodeLinOp::define(a_geom, a_grids, a_dmap, a_info);
}

void
MultiCompNodeLap::prepareForSolve ()
{
    MLNodeLinOp::prepareForSolve();
    buildMasks(); // e.g., m_dirichlet_mask, which is 1 on Dirichlet nodes
}

void
MultiCompNodeLap::Fapply (int amrlev, int mglev, MultiFab& out, MultiFab const& in) const
{
    int const ncomp = m_ncomp;
    auto const coeff = m_coeff;
    auto const dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    auto const& out_a = out.arrays();
    auto const& in_a = in.const_arrays();
    auto const& msk_a = m_dirichlet_mask[amrlev][mglev]->const_arrays();

    ParallelFor(out, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        if (msk_a[b](i,j,k)) { // Dirichlet node
            for (int n = 0; n < ncomp; ++n) { out_a[b](i,j,k,n) = Real(0); }
            return;
        }
        auto const& x = in_a[b];
        GpuArray<Real,max_ncomp> lap{};
        for (int m = 0; m < ncomp; ++m) {
            lap[m] = AMREX_D_TERM(
                (x(i-1,j,k,m) - Real(2)*x(i,j,k,m) + x(i+1,j,k,m)) * dxinv[0]*dxinv[0],
              + (x(i,j-1,k,m) - Real(2)*x(i,j,k,m) + x(i,j+1,k,m)) * dxinv[1]*dxinv[1],
              + (x(i,j,k-1,m) - Real(2)*x(i,j,k,m) + x(i,j,k+1,m)) * dxinv[2]*dxinv[2]);
        }
        for (int n = 0; n < ncomp; ++n) {
            Real r = 0;
            for (int m = 0; m < ncomp; ++m) { r += coeff[n*ncomp+m] * lap[m]; }
            out_a[b](i,j,k,n) = r;
        }
    });
}

void
MultiCompNodeLap::Fsmooth (int amrlev, int mglev, MultiFab& sol, MultiFab const& rhs) const
{
    // Jacobi with the diagonal blocks D = C * diag(lap):
    //     sol += omega * D^{-1} (rhs - L(sol))
    MultiFab Lsol(sol.boxArray(), sol.DistributionMap(), m_ncomp, 0);
    Fapply(amrlev, mglev, Lsol, sol);

    int const ncomp = m_ncomp;
    auto const cinv = m_coeff_inv;
    auto const dxinv = m_geom[amrlev][mglev].InvCellSizeArray();
    Real diag_lap = 0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { diag_lap -= Real(2)*dxinv[d]*dxinv[d]; }
    // Damping factor with good smoothing for the Laplacian: 2/3, 4/5, 6/7 in 1D-3D
    Real const omega = Real(2*AMREX_SPACEDIM)/Real(2*AMREX_SPACEDIM+1);

    auto const& sol_a = sol.arrays();
    auto const& rhs_a = rhs.const_arrays();
    auto const& Lsol_a = Lsol.const_arrays();
    auto const& msk_a = m_dirichlet_mask[amrlev][mglev]->const_arrays();

    ParallelFor(sol, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        if (msk_a[b](i,j,k)) { return; }
        GpuArray<Real,max_ncomp> res{};
        for (int m = 0; m < ncomp; ++m) {
            res[m] = rhs_a[b](i,j,k,m) - Lsol_a[b](i,j,k,m);
        }
        for (int n = 0; n < ncomp; ++n) {
            Real z = 0;
            for (int m = 0; m < ncomp; ++m) { z += cinv[n*ncomp+m] * res[m]; }
            sol_a[b](i,j,k,n) += omega * z / diag_lap;
        }
    });

    // Make nodes shared by several boxes consistent.
    nodalSync(amrlev, mglev, sol);
}

void
MultiCompNodeLap::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    // Fill the ghost nodes of the fine residual.
    applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

    // crse and fine may have different layouts, e.g., after agglomeration.
    bool const need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        cfine.define(amrex::coarsen(fine.boxArray(), 2), fine.DistributionMap(), m_ncomp, 0);
    }
    MultiFab& cmf = need_parallel_copy ? cfine : crse;

    // Fine nodes with valid data: the nodal domain, plus one layer of ghost
    // nodes in periodic directions.
    Box const fdomain = amrex::surroundingNodes(m_geom[amrlev][cmglev-1].growPeriodicDomain(1));

    int const ncomp = m_ncomp;
    auto const& c_a = cmf.arrays();
    auto const& f_a = fine.const_arrays();

    // Weights 1/4, 1/2, 1/4 in each direction
    constexpr int dj = (AMREX_SPACEDIM > 1) ? 1 : 0;
    constexpr int dk = (AMREX_SPACEDIM > 2) ? 1 : 0;
    ParallelFor(cmf, IntVect(0), ncomp, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
    {
        Real r = 0;
        for (int kk = -dk; kk <= dk; ++kk) {
        for (int jj = -dj; jj <= dj; ++jj) {
        for (int ii = -1; ii <= 1; ++ii) {
            if (fdomain.contains(2*i+ii, 2*j+jj, 2*k+kk)) {
                Real w = Real(2-std::abs(ii)) * Real(0.25);
                if (dj) { w *= Real(2-std::abs(jj)) * Real(0.25); }
                if (dk) { w *= Real(2-std::abs(kk)) * Real(0.25); }
                r += w * f_a[b](2*i+ii, 2*j+jj, 2*k+kk, n);
            }
        }}}
        c_a[b](i,j,k,n) = r;
    });

    if (need_parallel_copy) {
        crse.ParallelCopy(cfine);
    }
    setDirichletNodesToZero(amrlev, cmglev, crse);
}

void
MultiCompNodeLap::interpolation (int amrlev, int fmglev, MultiFab& fine, MultiFab const& crse) const
{
    bool const need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        cfine.define(amrex::coarsen(fine.boxArray(), 2), fine.DistributionMap(), m_ncomp, 0);
        cfine.ParallelCopy(crse);
    }
    MultiFab const& cmf = need_parallel_copy ? cfine : crse;

    auto const& f_a = fine.arrays();
    auto const& c_a = cmf.const_arrays();
    auto const& msk_a = m_dirichlet_mask[amrlev][fmglev]->const_arrays();

    // A fine node coincides with a coarse node, or lies between 2, 4 or 8 of
    // them.  i>>1 and (i+1)>>1 are equal for even i and are the two
    // neighbors for odd i, so averaging over them handles all cases.
    ParallelFor(fine, IntVect(0), m_ncomp,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
    {
        if (msk_a[b](i,j,k)) { return; }
        int const ic[2] = {i>>1, (i+1)>>1};
        int const jc[2] = {j>>1, (j+1)>>1};
        int const kc[2] = {k>>1, (k+1)>>1};
        Real r = 0;
        for (int kk = 0; kk < 2; ++kk) {
        for (int jj = 0; jj < 2; ++jj) {
        for (int ii = 0; ii < 2; ++ii) {
            r += c_a[b](ic[ii],jc[jj],kc[kk],n);
        }}}
        f_a[b](i,j,k,n) += Real(0.125) * r;
    });
}
