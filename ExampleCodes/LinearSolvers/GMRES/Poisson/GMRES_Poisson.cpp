#include "GMRES_Poisson.H"

using namespace amrex;

GMRESPOISSON::GMRESPOISSON (const BoxArray& ba, const DistributionMapping& dm, const Geometry& geom)
    : m_ba(ba), m_dm(dm), m_geom(geom)
{
    m_gmres.define(*this);
}

MultiFab GMRESPOISSON::makeVecRHS () const
{
    return MultiFab(m_ba, m_dm, 1, 0);
}

MultiFab GMRESPOISSON::makeVecLHS () const
{
    return MultiFab(m_ba, m_dm, 1, 1);
}

Real GMRESPOISSON::norm2 (MultiFab const& mf) const
{
    return mf.norm2();
}

void GMRESPOISSON::scale (MultiFab& mf, RT scale_factor)
{
    mf.mult(scale_factor);
}

Real GMRESPOISSON::dotProduct (MultiFab const& mf1, MultiFab const& mf2) const
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
