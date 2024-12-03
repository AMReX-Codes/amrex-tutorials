#ifndef AMREX_GMRES_POISSON_H_
#define AMREX_GMRES_POISSON_H_
#include <AMReX_Config.H>

#include <AMReX_GMRES.H>
#include <utility>

namespace amrex {

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
    using GM = GMRES<MultiFab,GMRESPOISSON>;

    explicit GMRESPOISSON (const BoxArray& ba, const DistributionMapping& dm, const Geometry& geom);

    /**
     * \brief Solve the linear system
     *
     * \param a_sol     unknowns, i.e., x in A x = b.
     * \param a_rhs     RHS, i.e., b in A x = b.
     * \param a_tol_rel relative tolerance.
     * \param a_tol_abs absolute tolerance.
     */
    void solve (MultiFab& a_sol, MultiFab const& a_rhs, RT a_tol_rel, RT a_tol_abs);

    //! Sets verbosity.
    void setVerbose (int v) { m_gmres.setVerbose(v); }

    //! Get the GMRES object.
    GM& getGMRES () { return m_gmres; }

    //! Make MultiFab without ghost cells
    MultiFab makeVecRHS () const;

    //! Make MultiFab with ghost cells and set ghost cells to zero
    MultiFab makeVecLHS () const;

    RT norm2 (MultiFab const& mf) const;

    static void scale (MultiFab& mf, RT scale_factor);

    RT dotProduct (MultiFab const& mf1, MultiFab const& mf2) const;

    //! lhs = 0
    static void setToZero (MultiFab& lhs);

    //! lhs = rhs
    static void assign (MultiFab& lhs, MultiFab const& rhs);

    //! lhs += a*rhs
    static void increment (MultiFab& lhs, MultiFab const& rhs, RT a);

    //! lhs = a*rhs_a + b*rhs_b
    static void linComb (MultiFab& lhs, RT a, MultiFab const& rhs_a, RT b, MultiFab const& rhs_b);

    //! lhs = L(rhs)
    void apply (MultiFab& lhs, MultiFab& rhs) const;

    void precond (MultiFab& lhs, MultiFab const& rhs) const;

    //! Control whether or not to use MLMG as preconditioner.
    bool usePrecond (bool new_flag) { return std::exchange(m_use_precond, new_flag); }

private:
    GM m_gmres;
    BoxArray m_ba;
    DistributionMapping m_dm;
    Geometry m_geom;
    bool m_use_precond;
};

GMRESPOISSON::GMRESPOISSON (const BoxArray& ba, const DistributionMapping& dm, const Geometry& geom)
    : m_ba(ba), m_dm(dm), m_geom(geom)
{
    m_gmres.define(*this);
}

auto GMRESPOISSON::makeVecRHS () const -> MultiFab
{
    return MultiFab(m_ba, m_dm, 1, 0);
}

auto GMRESPOISSON::makeVecLHS () const -> MultiFab
{
    return MultiFab(m_ba, m_dm, 1, 1);
}

auto GMRESPOISSON::norm2 (MultiFab const& mf) const -> RT
{
    return mf.norm2();
}

void GMRESPOISSON::scale (MultiFab& mf, RT scale_factor)
{
    mf.mult(scale_factor);
}

auto GMRESPOISSON::dotProduct (MultiFab const& mf1, MultiFab const& mf2) const -> RT
{
    return MultiFab::Dot(mf1,0,mf2,0,1,0);
}

void GMRESPOISSON::setToZero (MultiFab& lhs)
{
    lhs.setVal(0.);
}

void GMRESPOISSON::assign (MultiFab& lhs, MultiFab const& rhs)
{
    MultiFab::Copy(lhs,rhs,0,0,1,0);
}

void GMRESPOISSON::increment (MultiFab& lhs, MultiFab const& rhs, RT a)
{
    MultiFab::Saxpy(lhs,a,rhs,0,0,1,0);
}

void GMRESPOISSON::linComb (MultiFab& lhs, RT a, MultiFab const& rhs_a, RT b, MultiFab const& rhs_b)
{
    MultiFab::LinComb(lhs,a,rhs_a,0,b,rhs_b,0,0,1,0);
}

void GMRESPOISSON::apply (MultiFab& lhs, MultiFab& rhs) const
{
    // apply matrix to rhs for output lhs
    rhs.FillBoundary(m_geom.periodicity());

    const GpuArray<Real, AMREX_SPACEDIM> dx = m_geom.CellSizeArray();

    for ( MFIter mfi(lhs,TilingIfNotGPU()); mfi.isValid(); ++mfi ) {

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

void GMRESPOISSON::precond (MultiFab& lhs, MultiFab const& rhs) const
{
    if (m_use_precond) {

        // in the preconditioner we use right-preconditioning to solve
        // lhs = P^inv(rhs), where P^inv approximates A^inv
        // Here we use Jacobi iterations to represent P^inv with an initial guess of lhs=0

        const GpuArray<Real, AMREX_SPACEDIM> dx = m_geom.CellSizeArray();

        amrex::Real fac = 0.;
        for (int d=0; d<AMREX_SPACEDIM; ++d) { fac -= 2./(dx[d]*dx[d]); }

        MultiFab tmp(m_ba, m_dm, 1, 1);
        auto const& tmp_ma = tmp.const_arrays();
        auto const& rhs_ma = rhs.const_arrays();
        auto const& lhs_ma = lhs.arrays();

        lhs.setVal(0.);

        const int niters = 8;
        for (int iter = 0; iter < niters; ++iter) {

            MultiFab::Copy(tmp, lhs, 0, 0, 1, 0);
            tmp.FillBoundary(m_geom.periodicity());

            ParallelFor(lhs, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
            {
                auto const& tmp_ptr = tmp_ma[b];
                auto ax = ( tmp_ptr(i+1,j,k) - 2.*tmp_ptr(i,j,k) + tmp_ptr(i-1,j,k) ) / (dx[0]*dx[0])
                        + ( tmp_ptr(i,j+1,k) - 2.*tmp_ptr(i,j,k) + tmp_ptr(i,j-1,k) ) / (dx[1]*dx[1])
#if (AMREX_SPACEDIM == 3)
                        + ( tmp_ptr(i,j,k+1) - 2.*tmp_ptr(i,j,k) + tmp_ptr(i,j,k-1) ) / (dx[2]*dx[2])
#endif
                    ;
                auto res = rhs_ma[b](i,j,k) - ax;

                lhs_ma[b](i,j,k) += res / fac * Real(2./3.); // 2/3: weighted jacobi

            });
            Gpu::streamSynchronize();

        }
    } else {
        MultiFab::Copy(lhs,rhs,0,0,1,0);
    }
}

void GMRESPOISSON::solve (MultiFab& a_sol, MultiFab const& a_rhs, RT a_tol_rel, RT a_tol_abs)
{
    m_gmres.solve(a_sol, a_rhs, a_tol_rel, a_tol_abs);
}

}

#endif
