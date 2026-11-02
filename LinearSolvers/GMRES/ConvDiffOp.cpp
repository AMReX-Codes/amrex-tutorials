#include "ConvDiffOp.H"

using namespace amrex;

void ConvDiffOp::apply (MultiFab& Ax, MultiFab const& x) const
{
    // GMRES allows us to fill the ghost cells of x.
    const_cast<MultiFab&>(x).FillBoundary(m_geom.periodicity());

    auto const dxinv = m_geom.InvCellSizeArray();
    auto const vel = m_vel;
    Box const domain = m_geom.Domain();
    auto const& x_a = x.const_arrays();
    auto const& ax_a = Ax.arrays();
    ParallelFor(Ax, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        auto const& xa = x_a[b];
        IntVect const iv(AMREX_D_DECL(i,j,k));
        Real const xc = xa(iv);
        Real lap = 0, adv = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            IntVect const ivm = iv - IntVect::TheDimensionVector(d);
            IntVect const ivp = iv + IntVect::TheDimensionVector(d);
            // Dirichlet u = 0 on the domain face: the ghost value is -xc.
            Real const xm = (ivm[d] < domain.smallEnd(d)) ? -xc : xa(ivm);
            Real const xp = (ivp[d] > domain.bigEnd(d))   ? -xc : xa(ivp);
            lap += (xm - Real(2)*xc + xp) * dxinv[d]*dxinv[d];
            adv += vel[d] * (xp - xm) * Real(0.5)*dxinv[d];
        }
        ax_a[b](i,j,k) = lap - adv;
    });
}

void ConvDiffOp::precond (MultiFab& lhs, MultiFab const& rhs) const
{
    if (m_use_precond) {
        // Approximately solve lap(lhs) = rhs with one V-cycle, starting
        // from zero.  A fixed number of cycles keeps the preconditioner a
        // linear operator, as GMRES requires.
        lhs.setVal(0);
        m_mlmg->setPrecondIter(1);
        m_mlmg->precond({&lhs}, {&rhs}, Real(0), Real(0));
    } else {
        MultiFab::Copy(lhs, rhs, 0, 0, 1, 0);
    }
}
