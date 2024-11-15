#ifndef AMREX_GMRES_MLMG_H_
#define AMREX_GMRES_MLMG_H_
#include <AMReX_Config.H>

#include <AMReX_GMRES.H>
#include <utility>

namespace amrex {

/**
 * \brief Solve using GMRES with multigrid as preconditioner
 *
 * The linear system to solve is provided by MLMG, which is also being used
 * as the preconditioner.
 *
 */
template <typename MF>
class GMRESPOISSONT
{
public:
    using RT = amrex::Real; // typename MF::RT; // double or float
    using GM = GMRES<MF,GMRESPOISSONT<MF>>;

    explicit GMRESPOISSONT (const BoxArray& ba, const DistributionMapping& dm, const Geometry& geom);

    /**
     * \brief Solve the linear system
     *
     * \param a_sol     unknowns, i.e., x in A x = b.
     * \param a_rhs     RHS, i.e., b in A x = b.
     * \param a_tol_rel relative tolerance.
     * \param a_tol_abs absolute tolerance.
     */
    void solve (MF& a_sol, MF const& a_rhs, RT a_tol_rel, RT a_tol_abs);

    //! Sets verbosity.
    void setVerbose (int v) { m_gmres.setVerbose(v); }

    //! Get the GMRES object.
    GM& getGMRES () { return m_gmres; }

    //! Make MultiFab without ghost cells
    MF makeVecRHS () const;

    //! Make MultiFab with ghost cells and set ghost cells to zero
    MF makeVecLHS () const;

    RT norm2 (MF const& mf) const;

    static void scale (MF& mf, RT scale_factor);

    RT dotProduct (MF const& mf1, MF const& mf2) const;

    //! lhs = 0
    static void setToZero (MF& lhs);

    //! lhs = rhs
    static void assign (MF& lhs, MF const& rhs);

    //! lhs += a*rhs
    static void increment (MF& lhs, MF const& rhs, RT a);

    //! lhs = a*rhs_a + b*rhs_b
    static void linComb (MF& lhs, RT a, MF const& rhs_a, RT b, MF const& rhs_b);

    //! lhs = L(rhs)
    void apply (MF& lhs, MF& rhs) const;

    void precond (MF& lhs, MF const& rhs) const;

    //! Control whether or not to use MLMG as preconditioner.
    bool usePrecond (bool new_flag) { return std::exchange(m_use_precond, new_flag); }

private:
    GM m_gmres;
    BoxArray m_ba;
    DistributionMapping m_dm;
    Geometry m_geom;
    bool m_use_precond = false;
};

template <typename MF>
GMRESPOISSONT<MF>::GMRESPOISSONT (const BoxArray& ba, const DistributionMapping& dm, const Geometry& geom)
    : m_ba(ba), m_dm(dm), m_geom(geom)
{
    m_gmres.define(*this);
}

template <typename MF>
auto GMRESPOISSONT<MF>::makeVecRHS () const -> MF
{
    return MultiFab(m_ba, m_dm, 1, 0);
}

template <typename MF>
auto GMRESPOISSONT<MF>::makeVecLHS () const -> MF
{
    MF mf(m_ba, m_dm, 1, 1);
    mf.setVal(0.);
    return mf;
}

template <typename MF>
auto GMRESPOISSONT<MF>::norm2 (MF const& mf) const -> RT
{
    auto r = MultiFab::Dot(mf,0,1,0);
    return std::sqrt(r);
}

template <typename MF>
void GMRESPOISSONT<MF>::scale (MF& mf, RT scale_factor)
{
    Scale(mf, scale_factor, 0, nComp(mf), 0);
}

template <typename MF>
auto GMRESPOISSONT<MF>::dotProduct (MF const& mf1, MF const& mf2) const -> RT
{
    return MultiFab::Dot(mf1,0,mf2,0,1,0);
}

template <typename MF>
void GMRESPOISSONT<MF>::setToZero (MF& lhs)
{
    setVal(lhs, RT(0.0));
}

template <typename MF>
void GMRESPOISSONT<MF>::assign (MF& lhs, MF const& rhs)
{
    LocalCopy(lhs, rhs, 0, 0, nComp(lhs), IntVect(0));
}

template <typename MF>
void GMRESPOISSONT<MF>::increment (MF& lhs, MF const& rhs, RT a)
{
    Saxpy(lhs, a, rhs, 0, 0, nComp(lhs), IntVect(0));
}

template <typename MF>
void GMRESPOISSONT<MF>::linComb (MF& lhs, RT a, MF const& rhs_a, RT b, MF const& rhs_b)
{
    LinComb(lhs, a, rhs_a, 0, b, rhs_b, 0, 0, nComp(lhs), IntVect(0));
}

template <typename MF>
void GMRESPOISSONT<MF>::apply (MF& lhs, MF& rhs) const
{
    // apply matrix to rhs for output lhs
    rhs.FillBoundary(m_geom.periodicity());
    
    const GpuArray<Real, AMREX_SPACEDIM> dx = m_geom.CellSizeArray();

    for ( MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        
        const Box& bx = mfi.tilebox();

        const Array4<const Real> & rhs_p = rhs.array(mfi);
        const Array4<      Real> & lhs_p = lhs.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            lhs_p(i,j,k) = ( rhs_p(i+1,j,k) - 2.*rhs_p(i,j,k) + rhs_p(i-1,j,k) ) / (dx[0]*dx[0])
                         + ( rhs_p(i,j+1,k) - 2.*rhs_p(i,j,k) + rhs_p(i,j-1,k) ) / (dx[1]*dx[1])
#if (AMREX_SPACEDIM == 3)
                         + ( rhs_p(i,j,k+1) - 2.*rhs_p(i,j,k) + rhs_p(i,j,k-1) ) / (dx[2]*dx[2])
#endif
                ;
        });

    }

    
}

template <typename MF>
void GMRESPOISSONT<MF>::precond (MF& lhs, MF const& rhs) const
{
    if (m_use_precond) {

    } else {
        LocalCopy(lhs, rhs, 0, 0, nComp(lhs), IntVect(0));
    }
}

template <typename MF>
void GMRESPOISSONT<MF>::solve (MF& a_sol, MF const& a_rhs, RT a_tol_rel, RT a_tol_abs)
{
    /*
    auto res = makeVecRHS();
    m_mlmg->apply({&res}, {&a_sol}); // res = L(sol)
    increment(res, a_rhs, RT(-1)); // res = L(sol) - rhs
    auto cor = makeVecLHS();
    m_linop->setDirichletNodesToZero(0,0,res);
    m_gmres.solve(cor, res, a_tol_rel, a_tol_abs); // L(cor) = res
    increment(a_sol, cor, RT(-1)); // sol = sol - cor
    */
    auto res = makeVecRHS();
    apply(res,a_sol);
    increment(res, a_rhs, RT(-1)); // res = L(sol) - rhs
    auto cor = makeVecLHS();
    m_gmres.solve(cor, res, a_tol_rel, a_tol_abs); // L(cor) = res
    increment(a_sol, cor, RT(-1)); // sol = sol - cor
}

using GMRESPOISSON = GMRESPOISSONT<MultiFab>;

}

#endif
