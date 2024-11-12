#ifndef AMREX_GMRES_POISSON_H_
#define AMREX_GMRES_POISSON_H_
#include <AMReX_Config.H>

#include <AMReX_GMRES.H>

#include <utility>

namespace amrex {

/**
 * \brief Solve Poisson equation using null preconditioner
 */
template <typename MF>
class GMRESPoissonT
{
public:
    using GM = GMRES<MF,GMRESPoissonT<MF>>;

    explicit GMRESPoissonT ();

    /**
     * \brief Solve the linear system
     *
     * \param a_sol     unknowns, i.e., x in A x = b.
     * \param a_rhs     RHS, i.e., b in A x = b.
     * \param a_tol_rel relative tolerance.
     * \param a_tol_abs absolute tolerance.
     */
    void solve (MF& a_sol, MF const& a_rhs, amrex::Real a_tol_rel, amrex::Real a_tol_abs);

    //! Sets verbosity.
    void setVerbose (int v) { m_gmres.setVerbose(v); }

    //! Sets the max number of iterations
    void setMaxIters (int niters) { m_gmres.setMaxIters(niters); }

    //! Gets the number of iterations.
    [[nodiscard]] int getNumIters () const { return m_gmres.getNumIters(); }

    //! Gets the 2-norm of the residual.
    [[nodiscard]] amrex::Real getResidualNorm () const { return m_gmres.getResidualNorm(); }

    //! Get the GMRES object.
    GM& getGMRES () { return m_gmres; }

    //! Make MultiFab without ghost cells
    MF makeVecRHS () const;

    //! Make MultiFab with ghost cells and set ghost cells to zero
    MF makeVecLHS () const;

    amrex::Real norm2 (MF const& mf) const;

    static void scale (MF& mf, amrex::Real scale_factor);

    amrex::Real dotProduct (MF const& mf1, MF const& mf2) const;

    //! lhs = 0
    static void setToZero (MF& lhs);

    //! lhs = rhs
    static void assign (MF& lhs, MF const& rhs);

    //! lhs += a*rhs
    static void increment (MF& lhs, MF const& rhs, amrex::Real a);

    //! lhs = a*rhs_a + b*rhs_b
    static void linComb (MF& lhs, amrex::Real a, MF const& rhs_a, amrex::Real b, MF const& rhs_b);

    //! lhs = L(rhs)
    void apply (MF& lhs, MF const& rhs) const;

    void precond (MF& lhs, MF const& rhs) const;

private:
    GM m_gmres;
};

template <typename MF>
GMRESPoissonT<MF>::GMRESPoissonT ()
{

}

template <typename MF>
auto GMRESPoissonT<MF>::makeVecRHS () const -> MF
{
//    return m_linop->make(0, 0, IntVect(0));
    MultiFab x;
    return x;
}

template <typename MF>
auto GMRESPoissonT<MF>::makeVecLHS () const -> MF
{
//    auto mf = m_linop->make(0, 0, IntVect(1));
//    setBndry(mf, amrex::Real(0), 0, nComp(mf));
//    return mf;
    MultiFab x;
    return x;
}

template <typename MF>
auto GMRESPoissonT<MF>::norm2 (MF const& mf) const -> amrex::Real
{
//    auto r = m_linop->xdoty(0, 0, mf, mf, false);
//    return std::sqrt(r);
    return 0.;
}

template <typename MF>
void GMRESPoissonT<MF>::scale (MF& mf, amrex::Real scale_factor)
{
    Scale(mf, scale_factor, 0, nComp(mf), 0);
}

template <typename MF>
auto GMRESPoissonT<MF>::dotProduct (MF const& mf1, MF const& mf2) const -> amrex::Real
{
//    return m_linop->xdoty(0, 0, mf1, mf2, false);
    return 0.;
}

template <typename MF>
void GMRESPoissonT<MF>::setToZero (MF& lhs)
{
    setVal(lhs, amrex::Real(0.0));
}

template <typename MF>
void GMRESPoissonT<MF>::assign (MF& lhs, MF const& rhs)
{
    LocalCopy(lhs, rhs, 0, 0, nComp(lhs), IntVect(0));
}

template <typename MF>
void GMRESPoissonT<MF>::increment (MF& lhs, MF const& rhs, amrex::Real a)
{
    Saxpy(lhs, a, rhs, 0, 0, nComp(lhs), IntVect(0));
}

template <typename MF>
void GMRESPoissonT<MF>::linComb (MF& lhs, amrex::Real a, MF const& rhs_a, amrex::Real b, MF const& rhs_b)
{
    LinComb(lhs, a, rhs_a, 0, b, rhs_b, 0, 0, nComp(lhs), IntVect(0));
}

template <typename MF>
void GMRESPoissonT<MF>::apply (MF& lhs, MF const& rhs) const
{
    // fixme
}

template <typename MF>
void GMRESPoissonT<MF>::precond (MF& lhs, MF const& rhs) const
{
    LocalCopy(lhs, rhs, 0, 0, nComp(lhs), IntVect(0));
}

template <typename MF>
void GMRESPoissonT<MF>::solve (MF& a_sol, MF const& a_rhs, amrex::Real a_tol_rel, amrex::Real a_tol_abs)
{
    /*
    auto res = makeVecRHS();
    auto cor = makeVecLHS();
    m_gmres.solve(cor, res, a_tol_rel, a_tol_abs); // L(cor) = res
    */
}

using GMRESPoisson = GMRESPoissonT<MultiFab>;

}

#endif
