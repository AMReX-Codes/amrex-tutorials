#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB_utils.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleInterpolators.H>
#include <AMReX_ParticleMesh.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>
#include <AMReX_TracerParticle_mod_K.H>

#include <cmath>
#include <type_traits>

using namespace amrex;
using namespace amrex::literals;

namespace {

// Real components: the position (the first AMREX_SPACEDIM components of a
// pure SoA container), the initial position, and the weight the particle
// carries. The integer component counts how often a particle was reflected
// off the EB.
struct RealIdx
{
    enum {
        AMREX_D_DECL(x, y, z),
        AMREX_D_DECL(x0, y0, z0),
        w,
        ncomps
    };
};

struct IntIdx
{
    enum { reflected, ncomps };
};

using PC = ParticleContainerPureSoA<RealIdx::ncomps, IntIdx::ncomps>;
using ConstPTD = PC::ConstPTDType;
using ConstParticle = PC::ParticleType::ConstType;

// The spherical obstacle and the initial scalar, a Gaussian blob upstream
struct Problem
{
    Real radius = 0.15_rt;
    GpuArray<Real,AMREX_SPACEDIM> center{AMREX_D_DECL(0.5_rt, 0.5_rt, 0.5_rt)};
    GpuArray<Real,AMREX_SPACEDIM> blob{AMREX_D_DECL(0.25_rt, 0.45_rt, 0.5_rt)};
    Real blob_width = 0.08_rt;

    // Exact distance inside the sphere (> 0 inside), used only for checking
    [[nodiscard]] AMREX_GPU_HOST_DEVICE Real depth (RealVect const& x) const {
        Real r2 = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { r2 += (x[d]-center[d])*(x[d]-center[d]); }
        return radius - std::sqrt(r2);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE Real phi0 (RealVect const& x) const {
        Real r2 = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { r2 += (x[d]-blob[d])*(x[d]-blob[d]); }
        return std::exp(-r2/(2.0_rt*blob_width*blob_width));
    }
};

// Make a uniform velocity flow around the EB with a MAC projection, as in
// EB/Advection, and set it to zero on covered faces. The velocity is kept at
// face centers, where mac_interpolate expects it; EB_computeDivergence
// interpolates it to the face centroids for the cut-cell divergence.
void make_velocity (Array<MultiFab,AMREX_SPACEDIM>& u, Vector<Real> const& u0,
                    Geometry const& geom, EBFArrayBoxFactory const& factory)
{
    BoxArray const& ba = factory.boxArray();
    DistributionMapping const& dm = factory.DistributionMap();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { u[d].setVal(u0[d]); }

    MultiFab divu(ba, dm, 1, 0, MFInfo(), factory);
    MultiFab phi (ba, dm, 1, 1, MFInfo(), factory);
    phi.setVal(0.0_rt);
    EB_computeDivergence(divu, GetArrOfConstPtrs(u), geom, false);
    divu.mult(-1.0_rt); // rhs of -lap(phi) = -div(u)

    MLEBABecLap op({geom}, {ba}, {dm}, LPInfo(), {&factory});
    op.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic, LinOpBCType::Periodic,
                                 LinOpBCType::Periodic)},
                   {AMREX_D_DECL(LinOpBCType::Periodic, LinOpBCType::Periodic,
                                 LinOpBCType::Periodic)});
    op.setLevelBC(0, nullptr);
    op.setScalars(0.0_rt, 1.0_rt);
    op.setBCoeffs(0, 1.0_rt); // homogeneous Neumann on the EB
    MLMG mlmg(op);
    mlmg.setVerbose(0);
    const Real tol_rel = std::is_same_v<Real,float> ? 1.e-5_rt : 1.e-12_rt;
    mlmg.solve({&phi}, {&divu}, tol_rel, 0.0_rt);

    // The fluxes are -grad(phi): u -= grad(phi)
    Array<MultiFab,AMREX_SPACEDIM> gradphi;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        gradphi[d].define(u[d].boxArray(), dm, 1, 0, MFInfo(), factory);
    }
    mlmg.getFluxes({GetArrOfPtrs(gradphi)}, MLMG::Location::FaceCenter);

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        MultiFab::Add(u[d], gradphi[d], 0, 0, 1, 0);
    }
    EB_set_covered_faces(GetArrOfPtrs(u), 0.0_rt);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        // Ghost faces are needed to interpolate at particles near box edges.
        u[d].FillBoundary(geom.periodicity());
    }
}

// nppc particles at random positions in each cell outside the sphere. The
// weight is phi0 times cell volume / nppc, so the deposited density is phi0.
void init_particles (PC& pc, int nppc, Problem const& prob)
{
    const auto plo = pc.Geom(0).ProbLoArray();
    const auto dx = pc.Geom(0).CellSizeArray();
    const Real wcell = AMREX_D_TERM(dx[0], *dx[1], *dx[2]) / Real(nppc);
    const int myproc = ParallelDescriptor::MyProc();

    for (MFIter mfi = pc.MakeMFIter(0); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const Long np = bx.numPts() * nppc;
        const Long id0 = PC::ParticleType::NextID();
        PC::ParticleType::NextID(id0 + np);

        auto& tile = pc.DefineAndReturnParticleTile(0, mfi);
        tile.resize(np);
        auto ptd = tile.getParticleTileData();
        ParallelForRNG(np, [=] AMREX_GPU_DEVICE (Long ip, RandomEngine const& engine) noexcept
        {
            const IntVect iv = bx.atOffset(ip / nppc);
            ptd.idcpu(ip) = SetParticleIDandCPU(id0 + ip, myproc);
            RealVect x;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                x[d] = plo[d] + (Real(iv[d]) + Random(engine)) * dx[d];
                ptd.pos(d, ip) = ParticleReal(x[d]);
                ptd.rdata(RealIdx::x0+d)[ip] = ParticleReal(x[d]);
            }
            ptd.rdata(RealIdx::w)[ip] = ParticleReal(prob.phi0(x) * wcell);
            ptd.idata(IntIdx::reflected)[ip] = 0;
            // Redistribute removes invalid particles.
            if (prob.depth(x) > 0.0_rt) { ptd.id(ip).make_invalid(); }
        });
    }
    pc.Redistribute();
}

// A position, for interpolating at a point that is not a particle
struct Position
{
    ParticleReal x[AMREX_SPACEDIM];
    [[nodiscard]] AMREX_GPU_HOST_DEVICE ParticleReal pos (int d) const { return x[d]; }
};

// Signed distance s (positive in the fluid) and its gradient at x, by
// multilinear interpolation of the nodal values of the cell containing x
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real signed_distance (Array4<Real const> const& sd, ParticleReal const* x,
                      GpuArray<Real,AMREX_SPACEDIM> const& plo,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxi, Real* grad)
{
    int iv[3] = {0, 0, 0};
    Real f[AMREX_SPACEDIM];
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const Real l = (Real(x[d]) - plo[d]) * dxi[d];
        iv[d] = int(std::floor(l));
        f[d] = l - Real(iv[d]);
        grad[d] = 0.0_rt;
    }
    Real s = 0.0_rt;
    for (int c = 0; c < (1 << AMREX_SPACEDIM); ++c) {
        int n[3] = {0, 0, 0};
        Real w = 1.0_rt;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            n[d] = (c >> d) & 1;
            w *= n[d] ? f[d] : 1.0_rt - f[d];
        }
        const Real v = sd(iv[0]+n[0], iv[1]+n[1], iv[2]+n[2]);
        s += w * v;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            Real dw = (n[d] ? 1.0_rt : -1.0_rt) * dxi[d];
            for (int e = 0; e < AMREX_SPACEDIM; ++e) {
                if (e != d) { dw *= n[e] ? f[e] : 1.0_rt - f[e]; }
            }
            grad[d] += dw * v;
        }
    }
    return s;
}

// One midpoint (RK2) step with the face velocity u interpolated to the
// particles. The interpolation does not know about the EB, so a few
// particles near it end up inside the body; reflect them back across the
// EB using the signed distance sd.
void advect (PC& pc, Array<MultiFab,AMREX_SPACEDIM> const& u, MultiFab const& sd, Real dt)
{
    const auto plo = pc.Geom(0).ProbLoArray();
    const auto dxi = pc.Geom(0).InvCellSizeArray();
    for (PC::ParIterType pti(pc, 0); pti.isValid(); ++pti) {
        GpuArray<Array4<Real const>,AMREX_SPACEDIM> ua{
            AMREX_D_DECL(u[0].const_array(pti), u[1].const_array(pti), u[2].const_array(pti))};
        auto const& sda = sd.const_array(pti);
        auto ptd = pti.GetParticleTile().getParticleTileData();
        ParallelFor(pti.numParticles(), [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleReal v[AMREX_SPACEDIM];
            mac_interpolate(ptd[i], plo, dxi, ua, v);
            Position p{};
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                p.x[d] = ptd.pos(d, i) + ParticleReal(0.5_rt*dt) * v[d];
            }
            mac_interpolate(p, plo, dxi, ua, v);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                p.x[d] = ptd.pos(d, i) + ParticleReal(dt) * v[d];
            }

            // Inside the body (s < 0): move 2|s| along the normal grad(s)/|grad(s)|.
            Real grad[AMREX_SPACEDIM];
            const Real s = signed_distance(sda, p.x, plo, dxi, grad);
            if (s < 0.0_rt) {
                Real gnorm = 0.0_rt;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { gnorm += grad[d]*grad[d]; }
                gnorm = std::sqrt(gnorm);
                if (gnorm > 0.0_rt) {
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        p.x[d] += ParticleReal(-2.0_rt * s * grad[d] / gnorm);
                    }
                    ptd.idata(IntIdx::reflected)[i] += 1;
                }
            }
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { ptd.pos(d, i) = p.x[d]; }
        });
    }
    // Particles move less than a cell per step, so only neighboring ranks
    // communicate (local = 1). Redistribute also wraps periodic positions.
    pc.Redistribute(0, -1, 0, 1);
}

// Deposit the weights to the mesh with cloud-in-cell (linear) weights.
// ParticleToMesh adds the ghost cell contributions to the valid cells.
void deposit (PC const& pc, MultiFab& phi)
{
    const auto plo = pc.Geom(0).ProbLoArray();
    const auto dxi = pc.Geom(0).InvCellSizeArray();
    const Real inv_vol = AMREX_D_TERM(dxi[0], *dxi[1], *dxi[2]);
    ParticleToMesh(pc, phi, 0,
        [=] AMREX_GPU_DEVICE (ConstPTD const& ptd, int i, Array4<Real> const& phia) noexcept
        {
            ConstParticle const p = ptd[i];
            ParticleInterpolator::Linear interp(p, plo, dxi);
            interp.ParticleToMesh(p, phia, 0, 0, 1,
                [=] (ConstParticle const&, int) noexcept
                {
                    return ptd.rdata(RealIdx::w)[i] * inv_vol;
                });
        });
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int nppc = 4;
        Real stop_time = 0.3_rt;
        Real cfl = 0.5_rt;
        int plot_int = -1;
        Problem prob;
        Vector<Real> velocity{AMREX_D_DECL(1.0_rt, 0.25_rt, 0.0_rt)};
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nppc", nppc);
            pp.query("stop_time", stop_time);
            pp.query("cfl", cfl);
            pp.query("plot_int", plot_int);
            pp.query("radius", prob.radius);
            pp.queryarr("velocity", velocity, 0, AMREX_SPACEDIM);
        }

        // Periodic unit square/cube with a spherical obstacle in the middle
        const Box domain(IntVect(0), IntVect(n_cell-1));
        const RealBox rb({AMREX_D_DECL(0.0_rt,0.0_rt,0.0_rt)},
                         {AMREX_D_DECL(1.0_rt,1.0_rt,1.0_rt)});
        const Geometry geom(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(1,1,1)});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        const DistributionMapping dm(ba);

        EB2::SphereIF sphere(prob.radius, {AMREX_D_DECL(prob.center[0], prob.center[1],
                                                        prob.center[2])}, false);
        EB2::Build(EB2::makeShop(sphere), geom, 0, 30); // coarse levels for multigrid
        EBFArrayBoxFactory factory(EB2::IndexSpace::top().getLevel(geom), geom, ba, dm,
                                   {3, 3, 3}, EBSupport::full);

        // Face velocity with two ghost faces for the interpolation
        Array<MultiFab,AMREX_SPACEDIM> u;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            u[d].define(amrex::convert(ba, IntVect::TheDimensionVector(d)), dm, 1, 2,
                        MFInfo(), factory);
        }
        make_velocity(u, velocity, geom, factory);

        // Nodal signed distance to the EB, positive in the fluid. It is
        // accurate within a few cells of the EB, which is where it is used.
        MultiFab sd(amrex::convert(ba, IntVect(1)), dm, 1, 2, MFInfo(), factory);
        FillSignedDistance(sd, true);
        sd.FillBoundary(geom.periodicity());

        PC pc(geom, dm, ba);
        pc.SetSoACompileTimeNames({AMREX_D_DECL("x", "y", "z"),
                                   AMREX_D_DECL("x0", "y0", "z0"), "w"}, {"reflected"});
        init_particles(pc, nppc, prob);
        const Long np = pc.TotalNumberOfParticles();

        Real total_weight = ReduceSum(pc, [=] AMREX_GPU_HOST_DEVICE (ConstPTD const& ptd, int i)
                                      { return Real(ptd.rdata(RealIdx::w)[i]); });
        ParallelDescriptor::ReduceRealSum(total_weight);
        amrex::Print() << np << " particles, total weight " << total_weight << "\n";

        // phi: the particle weights deposited on the mesh
        const Real dv = AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
        MultiFab phi(ba, dm, 1, 1, MFInfo(), factory);
        deposit(pc, phi);

        auto write_plotfile = [&] (int step, Real time)
        {
            MultiFab plt(ba, dm, 1+AMREX_SPACEDIM, 0, MFInfo(), factory);
            MultiFab::Copy(plt, phi, 0, 0, 1, 0);
            EB_average_face_to_cellcenter(plt, 1, GetArrOfConstPtrs(u));
            const std::string name = amrex::Concatenate("plt", step, 5);
            EB_WriteSingleLevelPlotfile(name, plt, {"phi", AMREX_D_DECL("u","v","w")},
                                        geom, time, step);
            pc.WritePlotFile(name, "particles");
        };

        // Exact maximum depth of particles inside the sphere, in cells
        const Real dx = geom.CellSize(0);
        auto max_depth = [&] () {
            Real r = ReduceMax(pc, [=] AMREX_GPU_HOST_DEVICE (ConstPTD const& ptd, int i) -> Real
            {
                RealVect x;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { x[d] = Real(ptd.pos(d, i)); }
                return std::max(0.0_rt, prob.depth(x)) / dx;
            });
            ParallelDescriptor::ReduceRealMax(r);
            return r;
        };

        // Move forward to stop_time, then back with the velocity reversed.
        // The flow is steady, so the particles should return to where they
        // started, up to the time stepping error.
        Real umax = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { umax = std::max(umax, u[d].norm0()); }
        const int nsteps = std::max(1, int(std::ceil(stop_time * umax / (cfl * dx))));
        const Real dt = stop_time / Real(nsteps);
        amrex::Print() << nsteps << " steps forward and back with dt = " << dt << "\n";

        Real depth = 0;
        int step = 0;
        if (plot_int > 0) { write_plotfile(step, 0.0_rt); }
        for (int dir : {1, -1}) {
            if (dir < 0) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { u[d].mult(-1.0_rt, 0, 1, 2); }
            }
            for (int n = 0; n < nsteps; ++n) {
                advect(pc, u, sd, dt);
                ++step;
                depth = std::max(depth, max_depth());
                if (plot_int > 0 && step % plot_int == 0) {
                    deposit(pc, phi);
                    write_plotfile(step, (dir > 0) ? Real(n+1)*dt : stop_time - Real(n+1)*dt);
                }
            }
        }
        deposit(pc, phi);

        // Distance (periodic) of each particle that was never reflected from
        // its initial position: the mean and the maximum, in cells
        const GpuArray<Real,AMREX_SPACEDIM> prob_len{AMREX_D_DECL(geom.ProbLength(0),
                                                                 geom.ProbLength(1),
                                                                 geom.ProbLength(2))};
        auto distance = [=] AMREX_GPU_HOST_DEVICE (ConstPTD const& ptd, int i) -> Real
        {
            if (ptd.idata(IntIdx::reflected)[i] > 0) { return 0.0_rt; }
            Real r2 = 0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                Real dd = Real(ptd.pos(d, i) - ptd.rdata(RealIdx::x0+d)[i]);
                dd -= prob_len[d] * std::round(dd / prob_len[d]);
                r2 += dd*dd;
            }
            return std::sqrt(r2) / dx;
        };
        Real max_disp = ReduceMax(pc, distance);
        Real mean_disp = ReduceSum(pc, distance);
        Long nreflected = ReduceSum(pc, [=] AMREX_GPU_HOST_DEVICE (ConstPTD const& ptd, int i) -> Long
        {
            return (ptd.idata(IntIdx::reflected)[i] > 0) ? 1 : 0;
        });
        ParallelDescriptor::ReduceRealMax(max_disp);
        ParallelDescriptor::ReduceRealSum(mean_disp);
        ParallelDescriptor::ReduceLongSum(nreflected);
        if (nreflected < np) { mean_disp /= Real(np - nreflected); }

        const Real weight_err = std::abs(phi.sum(0)*dv - total_weight) / total_weight;
        amrex::Print() << "Back at t = 0:\n"
                       << "  particles reflected off the EB: " << nreflected << "\n"
                       << "  max depth inside the sphere after any step: " << depth
                       << " cells\n"
                       << "  distance of the other particles from their start: mean "
                       << mean_disp << ", max " << max_disp << " cells\n"
                       << "  relative error of the deposited weight: " << weight_err << "\n";

        // The reflection uses the interpolated signed distance, so a particle
        // can stay a small fraction of a cell inside the sphere.
        constexpr bool is_single = std::is_same_v<Real,float>
                                || std::is_same_v<ParticleReal,float>;
        const Real weight_tol = is_single ? 1.e-4_rt : 1.e-12_rt;
        if (pc.TotalNumberOfParticles() != np || weight_err > weight_tol ||
            depth > 0.2_rt || nreflected >= np || mean_disp > 1.e-3_rt) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
