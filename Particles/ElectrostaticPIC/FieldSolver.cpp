#include "FieldSolver.H"

#include <type_traits>

using namespace amrex::literals;

FieldSolver::FieldSolver (const amrex::Geometry& geom, const amrex::BoxArray& ba,
                          const amrex::DistributionMapping& dm, amrex::Real epsilon0)
    : m_geom(geom),
      m_epsilon0(epsilon0),
      m_linop({geom}, {ba}, {dm}),
      m_mlmg(m_linop)
{
    // All-periodic domain. The operator is singular; MLMG removes the mean
    // of the right-hand side to make the problem solvable.
    m_linop.setDomainBC({AMREX_D_DECL(amrex::LinOpBCType::Periodic,
                                      amrex::LinOpBCType::Periodic,
                                      amrex::LinOpBCType::Periodic)},
                        {AMREX_D_DECL(amrex::LinOpBCType::Periodic,
                                      amrex::LinOpBCType::Periodic,
                                      amrex::LinOpBCType::Periodic)});
    m_linop.setLevelBC(0, nullptr);
    m_mlmg.setVerbose(0);
}

void
FieldSolver::solvePhi (amrex::MultiFab& phi, const amrex::MultiFab& rho)
{
    // MLPoisson solves lap(phi) = rhs.
    amrex::MultiFab rhs(rho.boxArray(), rho.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(rhs, rho, 0, 0, 1, 0);
    rhs.mult(-1.0_rt/m_epsilon0);

    constexpr amrex::Real rel_tol = std::is_same_v<amrex::Real, float> ? 1.e-4_rt : 1.e-8_rt;
    const amrex::Real abs_tol = 0.0_rt;
    m_mlmg.solve({&phi}, {&rhs}, rel_tol, abs_tol);

    phi.FillBoundary(m_geom.periodicity());
}

void
FieldSolver::computeE (amrex::MultiFab& efield, const amrex::MultiFab& phi) const
{
    const auto dxi = m_geom.InvCellSizeArray();
    auto const& e = efield.arrays();
    auto const& p = phi.const_arrays();
    amrex::ParallelFor(efield,
        [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
        {
            amrex::ignore_unused(j, k);
            // Centered differences of phi at cell centers.
            AMREX_D_TERM(
                e[b](i,j,k,0) = -0.5_rt*dxi[0] * (p[b](i+1,j,k) - p[b](i-1,j,k));,
                e[b](i,j,k,1) = -0.5_rt*dxi[1] * (p[b](i,j+1,k) - p[b](i,j-1,k));,
                e[b](i,j,k,2) = -0.5_rt*dxi[2] * (p[b](i,j,k+1) - p[b](i,j,k-1));)
        });
    // Ghost cells are needed for interpolation to particles near grid edges.
    efield.FillBoundary(m_geom.periodicity());
}

amrex::Real
FieldSolver::fieldEnergy (const amrex::MultiFab& efield) const
{
    amrex::Real e2 = 0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        e2 += amrex::MultiFab::Dot(efield, d, efield, d, 1, 0);
    }
    const auto dx = m_geom.CellSizeArray();
    return 0.5_rt * m_epsilon0 * e2 * AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
}
