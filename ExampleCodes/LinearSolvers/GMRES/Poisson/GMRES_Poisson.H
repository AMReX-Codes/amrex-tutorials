#ifndef GMRES_POISSON_H_
#define GMRES_POISSON_H_

#include <AMReX_GMRES.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

/**
 * \brief Solve Poisson's equation using amrex GMRES class
 *
 * Refer to comments in amrex/Src/LinearSolvers/AMReX_GMRES.H
 * for details on function implementation requirements
 */
class GMRESPOISSON
{
public:
    using RT = amrex::Real; // double or float
    using GM = amrex::GMRES<amrex::MultiFab,GMRESPOISSON>;

    GMRESPOISSON (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::Geometry& geom);

    /**
     * \brief Solve the linear system
     *
     * \param a_sol     unknowns, i.e., x in A x = b.
     * \param a_rhs     RHS, i.e., b in A x = b.
     * \param a_tol_rel relative tolerance.
     * \param a_tol_abs absolute tolerance.
     */
    void solve (amrex::MultiFab& a_sol, amrex::MultiFab const& a_rhs, RT a_tol_rel, RT a_tol_abs);

    //! Sets verbosity.
    void setVerbose (int v) { m_gmres.setVerbose(v); }

    //! Get the GMRES object.
    GM& getGMRES () { return m_gmres; }

    //! Make MultiFab without ghost cells
    amrex::MultiFab makeVecRHS () const;

    //! Make MultiFab with ghost cells and set ghost cells to zero
    amrex::MultiFab makeVecLHS () const;

    RT norm2 (amrex::MultiFab const& mf) const;

    static void scale (amrex::MultiFab& mf, RT scale_factor);

    RT dotProduct (amrex::MultiFab const& mf1, amrex::MultiFab const& mf2) const;

    //! lhs = 0
    static void setToZero (amrex::MultiFab& lhs);

    //! lhs = rhs
    static void assign (amrex::MultiFab& lhs, amrex::MultiFab const& rhs);

    //! lhs += a*rhs
    static void increment (amrex::MultiFab& lhs, amrex::MultiFab const& rhs, RT a);

    //! lhs = a*rhs_a + b*rhs_b
    static void linComb (amrex::MultiFab& lhs, RT a, amrex::MultiFab const& rhs_a, RT b, amrex::MultiFab const& rhs_b);

    //! lhs = L(rhs)
    void apply (amrex::MultiFab& lhs, amrex::MultiFab& rhs) const;

    void precond (amrex::MultiFab& lhs, amrex::MultiFab const& rhs) const;

    //! Control whether or not to use MLMG as preconditioner.
    bool usePrecond (bool new_flag) { return std::exchange(m_use_precond, new_flag); }

private:
    GM m_gmres;
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    amrex::Geometry m_geom;
    bool m_use_precond;
};

#endif
