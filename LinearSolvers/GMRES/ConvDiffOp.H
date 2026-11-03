#ifndef CONV_DIFF_OP_H_
#define CONV_DIFF_OP_H_

#include <AMReX_GMRES.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>

/**
 * Matrix-free convection-diffusion operator A(u) = lap(u) - c.grad(u)
 * with homogeneous Dirichlet boundaries, discretized with central
 * differences on cell centers.  A is not symmetric, so we solve with GMRES.
 *
 * The class provides the vector operations and the matrix-vector product
 * that amrex::GMRES<MultiFab,ConvDiffOp> needs.  The preconditioner is one
 * MLMG V-cycle of the Laplacian alone.
 */
class ConvDiffOp
{
public:
    using RT = amrex::Real; // required by amrex::GMRES

    ConvDiffOp (amrex::Geometry const& geom, amrex::BoxArray const& ba,
                amrex::DistributionMapping const& dm,
                amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& velocity,
                amrex::MLMG& precond_mlmg)
        : m_geom(geom), m_ba(ba), m_dm(dm), m_vel(velocity), m_mlmg(&precond_mlmg)
    {}

    void usePrecond (bool flag) { m_use_precond = flag; }

    // Vectors for the right-hand side b and the unknowns x of A x = b.  x
    // has a ghost cell, so that apply can fill it.
    [[nodiscard]] amrex::MultiFab makeVecRHS () const
    {
        return amrex::MultiFab(m_ba, m_dm, 1, 0);
    }

    [[nodiscard]] amrex::MultiFab makeVecLHS () const
    {
        amrex::MultiFab mf(m_ba, m_dm, 1, 1);
        mf.setBndry(0);
        return mf;
    }

    // Linear algebra
    [[nodiscard]] RT norm2 (amrex::MultiFab const& mf) const { return mf.norm2(); }

    [[nodiscard]] RT dotProduct (amrex::MultiFab const& a, amrex::MultiFab const& b) const
    {
        return amrex::MultiFab::Dot(a, 0, b, 0, 1, 0);
    }

    static void scale (amrex::MultiFab& mf, RT fac) { mf.mult(fac); }

    static void setToZero (amrex::MultiFab& mf) { mf.setVal(0); }

    static void assign (amrex::MultiFab& lhs, amrex::MultiFab const& rhs)
    {
        amrex::MultiFab::Copy(lhs, rhs, 0, 0, 1, 0);
    }

    // lhs += a*rhs
    static void increment (amrex::MultiFab& lhs, amrex::MultiFab const& rhs, RT a)
    {
        amrex::MultiFab::Saxpy(lhs, a, rhs, 0, 0, 1, 0);
    }

    // lhs = a*rhs_a + b*rhs_b
    static void linComb (amrex::MultiFab& lhs, RT a, amrex::MultiFab const& rhs_a,
                         RT b, amrex::MultiFab const& rhs_b)
    {
        amrex::MultiFab::LinComb(lhs, a, rhs_a, 0, b, rhs_b, 0, 0, 1, 0);
    }

    // Ax = A(x).  x was made by makeVecLHS.
    void apply (amrex::MultiFab& Ax, amrex::MultiFab const& x) const;

    // Right preconditioning: lhs = P^{-1}(rhs), where P approximates A.
    void precond (amrex::MultiFab& lhs, amrex::MultiFab const& rhs) const;

private:
    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> m_vel;
    amrex::MLMG* m_mlmg;
    bool m_use_precond = true;
};

#endif
