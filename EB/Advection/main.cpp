#include <AMReX.H>
#include <AMReX_BCRec.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB_Redistribution.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <limits>
#include <string>
#include <type_traits>

using namespace amrex;

namespace {

// Ghost cells of the state and of dS/dt needed by the redistribution
constexpr int nghost = 4;

// Area fraction of the low face of cell (i,j,k) in direction d
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real face_area (EBDataArrays const& eb, int b, int d, int i, int j, int k)
{
    IntVect iv(AMREX_D_DECL(i,j,k));
    iv[d] -= 1;
    const Dim3 lo = iv.dim3();
    auto flag_lo = eb.get<EBData_t::cellflag>(b, lo.x, lo.y, lo.z);
    auto flag_hi = eb.get<EBData_t::cellflag>(b,i,j,k);
    if (flag_lo.isRegular() && flag_hi.isRegular()) { return Real(1.); }
    if (flag_lo.isCovered() || flag_hi.isCovered()) { return Real(0.); }
    // A cut cell is next to this face, so the area fractions exist.
    if (d == 0) { return eb.get<EBData_t::apx>(b,i,j,k); }
#if (AMREX_SPACEDIM == 3)
    if (d == 2) { return eb.get<EBData_t::apz>(b,i,j,k); }
#endif
    return eb.get<EBData_t::apy>(b,i,j,k);
}

// Monotonized central slope
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real mc_slope (Real sm, Real s0, Real sp)
{
    const Real dl = s0 - sm;
    const Real dr = sp - s0;
    if (dl*dr <= Real(0.)) { return Real(0.); }
    const Real dc = Real(0.5)*(sp - sm);
    const Real lim = Real(2.)*amrex::min(std::abs(dl), std::abs(dr));
    return std::copysign(amrex::min(std::abs(dc), lim), dc);
}

// Make the face velocity discretely divergence free in the cut-cell sense,
// sum over faces of (area fraction * u) = 0, with no flow through the EB.
// This is a MAC projection: solve lap(phi) = div(u), then u -= grad(phi).
void project_velocity (Array<MultiFab,AMREX_SPACEDIM>& u, Geometry const& geom,
                       EBFArrayBoxFactory const& factory, int verbose)
{
    BoxArray const& ba = factory.boxArray();
    DistributionMapping const& dm = factory.DistributionMap();

    MultiFab divu(ba, dm, 1, 0, MFInfo(), factory);
    MultiFab phi (ba, dm, 1, 1, MFInfo(), factory);
    phi.setVal(0.0_rt);

    // u lives on face centroids, so no interpolation is needed.
    EB_computeDivergence(divu, GetArrOfConstPtrs(u), geom, true);
    amrex::Print() << "max |div u| before projection: " << divu.norminf() << "\n";
    divu.mult(-1.0_rt); // rhs of -lap(phi) = -div(u)

    MLEBABecLap op({geom}, {ba}, {dm}, LPInfo(), {&factory});
    op.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic, LinOpBCType::Periodic,
                                 LinOpBCType::Periodic)},
                   {AMREX_D_DECL(LinOpBCType::Periodic, LinOpBCType::Periodic,
                                 LinOpBCType::Periodic)});
    op.setLevelBC(0, nullptr);
    op.setScalars(0.0_rt, 1.0_rt);
    op.setBCoeffs(0, 1.0_rt); // homogeneous Neumann on the EB by default

    MLMG mlmg(op);
    mlmg.setVerbose(verbose);
    const Real tol_rel = std::is_same_v<Real,float> ? 1.e-5_rt : 1.e-12_rt;
    mlmg.solve({&phi}, {&divu}, tol_rel, 0.0_rt);

    // The fluxes are -grad(phi) on face centroids.
    Array<MultiFab,AMREX_SPACEDIM> gradphi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        gradphi[d].define(u[d].boxArray(), dm, 1, 0, MFInfo(), factory);
    }
    mlmg.getFluxes({GetArrOfPtrs(gradphi)}, MLMG::Location::FaceCentroid);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        MultiFab::Add(u[d], gradphi[d], 0, 0, 1, 0);
    }

    EB_computeDivergence(divu, GetArrOfConstPtrs(u), geom, true);
    amrex::Print() << "max |div u| after projection:  " << divu.norminf() << "\n";
}

// Volume flux through each face: area fraction times velocity.  Faces with
// zero area carry no flux, so the EB is a wall.
void volume_flux (Array<MultiFab,AMREX_SPACEDIM>& ua,
                  Array<MultiFab,AMREX_SPACEDIM> const& u,
                  EBFArrayBoxFactory const& factory)
{
    auto const& ebdata = factory.getEBDataArrays();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        auto const& uaa = ua[d].arrays();
        auto const& ud = u[d].const_arrays();
        ParallelFor(ua[d], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
        {
            uaa[b](i,j,k) = face_area(ebdata, b, d, i, j, k) * ud[b](i,j,k);
        });
    }
}

// dS/dt = -div(u S) on cut cells, followed by redistribution
class Advection
{
public:
    Advection (Geometry const& geom, EBFArrayBoxFactory const& factory,
               Array<MultiFab,AMREX_SPACEDIM> const& u, int ncomp,
               std::string redistribution_type)
        : m_geom(geom), m_factory(factory), m_ncomp(ncomp),
          m_redistribution_type(std::move(redistribution_type))
    {
        BoxArray const& ba = factory.boxArray();
        DistributionMapping const& dm = factory.DistributionMap();

        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            m_ua[d].define(u[d].boxArray(), dm, 1, 0, MFInfo(), factory);
            m_flux[d].define(u[d].boxArray(), dm, ncomp, 0, MFInfo(), factory);
        }
        volume_flux(m_ua, u, factory);
        m_dSdt.define(ba, dm, ncomp, nghost, MFInfo(), factory);

        // Boundary conditions used by state redistribution: all periodic
        Vector<BCRec> bcrec(ncomp);
        for (auto& bc : bcrec) {
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                bc.setLo(d, BCType::int_dir);
                bc.setHi(d, BCType::int_dir);
            }
        }
        m_bcrec.resize(ncomp);
        Gpu::copy(Gpu::hostToDevice, bcrec.begin(), bcrec.end(), m_bcrec.begin());
    }

    // Largest stable time step for a given CFL number
    [[nodiscard]] Real estTimeStep (Real cfl) const
    {
        const auto dxinv = m_geom.InvCellSizeArray();
        Real rate = 0.0_rt;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            rate += m_ua[d].norminf() * dxinv[d];
        }
        return cfl / rate;
    }

    // Compute the redistributed dS/dt.  S must have nghost filled ghost cells.
    void compute_dSdt (MultiFab& dSdt, MultiFab const& S, Real dt)
    {
        const int ncomp = m_ncomp;
        auto const& ebdata = m_factory.getEBDataArrays();

        // Upwind fluxes with limited slopes; first order next to the EB.
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            auto const& flux = m_flux[d].arrays();
            auto const& ua = m_ua[d].const_arrays();
            auto const& s = S.const_arrays();
            ParallelFor(m_flux[d], IntVect(0), ncomp,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
            {
                const Real vel = ua[b](i,j,k);
                IntVect iv(AMREX_D_DECL(i,j,k));
                IntVect const e = IntVect::TheDimensionVector(d);
                // Upwind cell and its slope in direction d
                IntVect const c = (vel >= Real(0.)) ? iv - e : iv;
                bool regular = true;
                for (int m = -1; m <= 1; ++m) {
                    const Dim3 cm = (c + m*e).dim3();
                    regular = regular && ebdata.get<EBData_t::cellflag>
                        (b, cm.x, cm.y, cm.z).isRegular();
                }
                auto const& sb = s[b];
                const Real s0 = sb(c,n);
                const Real slope = regular ? mc_slope(sb(c-e,n), s0, sb(c+e,n)) : Real(0.);
                const Real sface = (vel >= Real(0.)) ? s0 + Real(0.5)*slope
                                                     : s0 - Real(0.5)*slope;
                flux[b](i,j,k,n) = vel * sface;
            });
        }

        // Conservative divergence, (sum of a u S) / (vfrac dx), with vfrac = 1
        // in regular cells
        {
            const auto dxinv = m_geom.InvCellSizeArray();
            auto const& rhs = m_dSdt.arrays();
            AMREX_D_TERM(auto const& fx = m_flux[0].const_arrays();,
                         auto const& fy = m_flux[1].const_arrays();,
                         auto const& fz = m_flux[2].const_arrays();)
            ParallelFor(m_dSdt, IntVect(0), ncomp,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
            {
                const Real vf = ebdata.get<EBData_t::volfrac>(b,i,j,k);
                if (vf > Real(0.)) {
                    rhs[b](i,j,k,n) = -(AMREX_D_TERM(
                          dxinv[0]*(fx[b](i+1,j,k,n) - fx[b](i,j,k,n)),
                        + dxinv[1]*(fy[b](i,j+1,k,n) - fy[b](i,j,k,n)),
                        + dxinv[2]*(fz[b](i,j,k+1,n) - fz[b](i,j,k,n)))) / vf;
                } else {
                    rhs[b](i,j,k,n) = Real(0.);
                }
            });
        }
        m_dSdt.FillBoundary(m_geom.periodicity());

        // Redistribution moves part of the update of small cut cells to their
        // neighbors so that the time step is not limited by the small cells.
        auto const& flags = m_factory.getMultiEBCellFlagFab();
        auto const& vfrac = m_factory.getVolFrac();
        auto const& area = m_factory.getAreaFrac();
        auto const& fcent = m_factory.getFaceCent();
        auto const& ccent = m_factory.getCentroid();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dSdt, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            auto const& out = dSdt.array(mfi);
            auto const& in = m_dSdt.array(mfi);
            const FabType type = flags[mfi].getType(amrex::grow(bx,nghost));
            if (type == FabType::singlevalued) {
                // Scratch space: weights for FluxRedist, a temporary state for
                // StateRedist.  It is local to the tile.
                FArrayBox scratch(amrex::grow(bx,nghost), ncomp, The_Async_Arena());
                scratch.setVal<RunOn::Device>(1.0_rt);
                ApplyRedistribution(bx, ncomp, out, in, S.const_array(mfi),
                                    scratch.array(), flags[mfi].const_array(),
                                    AMREX_D_DECL(area[0]->const_array(mfi),
                                                 area[1]->const_array(mfi),
                                                 area[2]->const_array(mfi)),
                                    vfrac.const_array(mfi),
                                    AMREX_D_DECL(fcent[0]->const_array(mfi),
                                                 fcent[1]->const_array(mfi),
                                                 fcent[2]->const_array(mfi)),
                                    ccent.const_array(mfi), m_bcrec.data(), m_geom,
                                    dt, m_redistribution_type);
            } else {
                // All regular or all covered: nothing to redistribute
                ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    out(i,j,k,n) = in(i,j,k,n);
                });
            }
        }
    }

private:
    Geometry m_geom;
    EBFArrayBoxFactory const& m_factory;
    int m_ncomp;
    std::string m_redistribution_type;
    Array<MultiFab,AMREX_SPACEDIM> m_ua;
    Array<MultiFab,AMREX_SPACEDIM> m_flux;
    MultiFab m_dSdt; // before redistribution
    Gpu::DeviceVector<BCRec> m_bcrec;
};

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int max_step = 1000;
        Real stop_time = 0.1_rt;
        Real cfl = 0.5_rt;
        int verbose = 0;
        int plot_int = -1;
        std::string redistribution_type = "StateRedist";
        Real radius = 0.15_rt;
        Vector<Real> velocity{AMREX_D_DECL(1.0_rt, 0.25_rt, 0.0_rt)};
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("max_step", max_step);
            pp.query("stop_time", stop_time);
            pp.query("cfl", cfl);
            pp.query("verbose", verbose);
            pp.query("plot_int", plot_int);
            pp.query("redistribution_type", redistribution_type);
            pp.query("radius", radius);
            pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
        }
        AMREX_ALWAYS_ASSERT(redistribution_type == "StateRedist" ||
                            redistribution_type == "FluxRedist" ||
                            redistribution_type == "NoRedist");

        // Periodic unit square/cube
        const Box domain(IntVect(0), IntVect(n_cell-1));
        const RealBox rb({AMREX_D_DECL(0.0_rt,0.0_rt,0.0_rt)},
                         {AMREX_D_DECL(1.0_rt,1.0_rt,1.0_rt)});
        const Geometry geom(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        const DistributionMapping dm(ba);

        // A spherical obstacle.  The MAC projection uses multigrid, which
        // needs coarsened EB levels.
        EB2::SphereIF sphere(radius, {AMREX_D_DECL(0.5_rt,0.5_rt,0.5_rt)}, false);
        EB2::Build(EB2::makeShop(sphere), geom, 0, 30);
        EBFArrayBoxFactory factory(EB2::IndexSpace::top().getLevel(geom), geom, ba, dm,
                                   {nghost+1, nghost+1, nghost+1}, EBSupport::full);

        // Uniform flow, projected to flow around the obstacle
        Array<MultiFab,AMREX_SPACEDIM> u;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            u[d].define(amrex::convert(ba, IntVect::TheDimensionVector(d)), dm, 1, 0,
                        MFInfo(), factory);
            u[d].setVal(velocity[d]);
        }
        project_velocity(u, geom, factory, verbose);

        // Component 0: a Gaussian blob; component 1: a constant, which a
        // divergence-free velocity must leave unchanged.
        const int ncomp = 2;
        MultiFab S   (ba, dm, ncomp, nghost, MFInfo(), factory);
        MultiFab Snew(ba, dm, ncomp, nghost, MFInfo(), factory);
        MultiFab dSdt(ba, dm, ncomp, 0, MFInfo(), factory);
        {
            const auto problo = geom.ProbLoArray();
            const auto dx = geom.CellSizeArray();
            auto const& sa = S.arrays();
            ParallelFor(S, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
            {
                AMREX_D_TERM(Real x = problo[0] + (Real(i)+Real(0.5))*dx[0] - Real(0.25);,
                             Real y = problo[1] + (Real(j)+Real(0.5))*dx[1] - Real(0.45);,
                             Real z = problo[2] + (Real(k)+Real(0.5))*dx[2] - Real(0.5);)
                Real r2 = AMREX_D_TERM(x*x, + y*y, + z*z);
                sa[b](i,j,k,0) = std::exp(-r2/Real(2.*0.08*0.08));
                sa[b](i,j,k,1) = Real(1.);
            });
        }
        S.FillBoundary(geom.periodicity());

        MultiFab const& vfrac = factory.getVolFrac();
        const Real dv = AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
        auto total_mass = [&] () { return MultiFab::Dot(vfrac, 0, S, 0, 1, 0) * dv; };

        auto write_plotfile = [&] (int step, Real time)
        {
            MultiFab plt(ba, dm, ncomp+AMREX_SPACEDIM, 0, MFInfo(), factory);
            MultiFab::Copy(plt, S, 0, 0, ncomp, 0);
            EB_average_face_to_cellcenter(plt, ncomp, GetArrOfConstPtrs(u));
            EB_WriteSingleLevelPlotfile(amrex::Concatenate("plt", step, 5), plt,
                                        {"blob", "constant", AMREX_D_DECL("u","v","w")},
                                        geom, time, step);
        };

        Advection advection(geom, factory, u, ncomp, redistribution_type);
        const Real dt = advection.estTimeStep(cfl);
        const Real mass0 = total_mass();
        amrex::Print() << "Redistribution: " << redistribution_type
                       << ", dt = " << dt << ", initial mass = " << mass0 << "\n";

        if (plot_int > 0) { write_plotfile(0, 0.0_rt); }

        // Second-order SSP Runge-Kutta
        Real time = 0.0_rt;
        int step = 0;
        while (step < max_step && time < stop_time*(1.0_rt-1.e-6_rt))
        {
            const Real dt_step = amrex::min(dt, stop_time - time);
            auto const& sa = S.const_arrays();
            auto const& sn = Snew.arrays();
            auto const& rhs = dSdt.const_arrays();

            // Stage 1: Snew = S + dt L(S)
            advection.compute_dSdt(dSdt, S, dt_step);
            ParallelFor(S, IntVect(0), ncomp,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
            {
                sn[b](i,j,k,n) = sa[b](i,j,k,n) + dt_step*rhs[b](i,j,k,n);
            });
            EB_set_covered(Snew, 0, ncomp, 0, 0.0_rt);
            Snew.FillBoundary(geom.periodicity());

            // Stage 2: S = (S + Snew + dt L(Snew)) / 2
            advection.compute_dSdt(dSdt, Snew, dt_step);
            auto const& so = S.arrays();
            ParallelFor(S, IntVect(0), ncomp,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
            {
                so[b](i,j,k,n) = Real(0.5)*(so[b](i,j,k,n) + sn[b](i,j,k,n)
                                            + dt_step*rhs[b](i,j,k,n));
            });
            EB_set_covered(S, 0, ncomp, 0, 0.0_rt);
            S.FillBoundary(geom.periodicity());

            time += dt_step;
            ++step;
            if (verbose > 0) {
                amrex::Print() << "Step " << step << " time " << time
                               << " mass " << total_mass() << "\n";
            }
            if (plot_int > 0 && step % plot_int == 0) { write_plotfile(step, time); }
        }

        // Self checks: no NaNs, mass conservation, a constant stays constant
        const Real mass = total_mass();
        const Real mass_err = std::abs(mass - mass0) / mass0;
        auto const& sa = S.const_arrays();
        auto const& vf = vfrac.const_arrays();
        Real const_err = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, S, IntVect(0),
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Real>
            {
                return (vf[b](i,j,k) > Real(0.)) ? std::abs(sa[b](i,j,k,1) - Real(1.))
                                                 : Real(0.);
            });
        ParallelDescriptor::ReduceRealMax(const_err);
        const bool has_nan = S.contains_nan(0, ncomp, 0);

        const Real eps = std::numeric_limits<Real>::epsilon();
        const Real mass_tol = 1.e3_rt * eps;
        const Real const_tol = std::is_same_v<Real,float> ? 1.e-3_rt : 1.e-8_rt;
        amrex::Print() << "Steps: " << step << ", time = " << time << "\n"
                       << "Blob min and max: " << S.min(0) << " " << S.max(0) << "\n"
                       << "Relative change in mass: " << mass_err
                       << " (tolerance " << mass_tol << ")\n"
                       << "Max deviation of the constant: " << const_err
                       << " (tolerance " << const_tol << ")\n";
        if (has_nan || mass_err > mass_tol || const_err > const_tol) {
            amrex::Print() << "FAILED" << (has_nan ? " (NaN)" : "") << "\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }

        write_plotfile(step, time);
    }
    amrex::Finalize();
    return status;
}
