#ifndef MULTI_COMP_NODE_LAP_H_
#define MULTI_COMP_NODE_LAP_H_

#include <AMReX_MLNodeLinOp.H>

#include <string>

/**
 * A user-defined nodal linear operator with N coupled components,
 *
 *     L(phi)_n = sum_m C_nm lap(phi_m),
 *
 * where C is a constant, invertible N x N matrix.  lap is the standard
 * (2D+1)-point Laplacian on the nodes.
 *
 * MLNodeLinOp provides the multigrid machinery, ghost node exchange,
 * Dirichlet node masks and dot products.  A derived class implements
 * the operator (Fapply), a smoother (Fsmooth), and the transfers between
 * multigrid levels (restriction and interpolation).  This class only
 * supports single-level solves.  Non-periodic boundaries must be Dirichlet.
 */
class MultiCompNodeLap
    : public amrex::MLNodeLinOp
{
public:
    static constexpr int max_ncomp = 4;

    /**
     * \param a_grids cell-centered BoxArrays; the unknowns live on their nodes.
     * \param a_coeff the matrix C in row-major order, with N*N entries.
     */
    MultiCompNodeLap (amrex::Vector<amrex::Geometry> const& a_geom,
                      amrex::Vector<amrex::BoxArray> const& a_grids,
                      amrex::Vector<amrex::DistributionMapping> const& a_dmap,
                      amrex::Vector<amrex::Real> const& a_coeff,
                      amrex::LPInfo const& a_info = amrex::LPInfo());

    [[nodiscard]] std::string name () const override { return "MultiCompNodeLap"; }

    [[nodiscard]] int getNComp () const override { return m_ncomp; }

    // out = L(in).  The ghost nodes of in have been filled.
    void Fapply (int amrlev, int mglev, amrex::MultiFab& out,
                 amrex::MultiFab const& in) const override;

    // One sweep of block Jacobi on L(sol) = rhs.
    void Fsmooth (int amrlev, int mglev, amrex::MultiFab& sol,
                  amrex::MultiFab const& rhs) const override;

    // Full-weighting restriction of the residual to MG level cmglev.
    void restriction (int amrlev, int cmglev, amrex::MultiFab& crse,
                      amrex::MultiFab& fine) const override;

    // fine += linear interpolation of the coarse correction.
    void interpolation (int amrlev, int fmglev, amrex::MultiFab& fine,
                        amrex::MultiFab const& crse) const override;

    void prepareForSolve () override;

private:
    int m_ncomp = 0;
    amrex::GpuArray<amrex::Real,max_ncomp*max_ncomp> m_coeff{};     // C
    amrex::GpuArray<amrex::Real,max_ncomp*max_ncomp> m_coeff_inv{}; // C^{-1}
};

#endif
