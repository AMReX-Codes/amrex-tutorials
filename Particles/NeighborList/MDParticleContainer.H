#ifndef MD_PARTICLE_CONTAINER_H_
#define MD_PARTICLE_CONTAINER_H_

#include <AMReX_NeighborParticles.H>
#include <AMReX_Particles.H>

#include <array>

// Positions first, then velocity, acceleration and potential energy (unit mass).
struct RealIdx
{
    enum {
        AMREX_D_DECL(x, y, z),
        AMREX_D_DECL(vx, vy, vz),
        AMREX_D_DECL(ax, ay, az),
        pe,
        ncomps
    };
};

// Parameters of the soft-sphere potential U(r) = k/2 (sigma - r)^2 for r < sigma.
struct MDParams
{
    amrex::Real sigma = amrex::Real(0.4);       // particle diameter (interaction range)
    amrex::Real stiffness = amrex::Real(100.0);
    amrex::Real r_list = amrex::Real(1.0);      // neighbor list radius, sigma + skin
};

// Pairs closer than r_list go into the neighbor list.
struct CheckPair
{
    amrex::ParticleReal r2_list;

    template <class P1, class P2>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator() (const P1& p1, const P2& p2) const noexcept
    {
        amrex::ParticleReal r2 = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            const amrex::ParticleReal dx = p1.pos(d) - p2.pos(d);
            r2 += dx*dx;
        }
        return r2 <= r2_list;
    }
};

struct MDDiagnostics
{
    amrex::Long np = 0;
    amrex::Real kinetic = 0;
    amrex::Real potential = 0;
    std::array<amrex::Real,AMREX_SPACEDIM> momentum{};
    amrex::Real abs_momentum = 0;   // sum of |p| over particles and directions
    amrex::Real max_speed = 0;
};

// A neighbor particle container keeps copies of nearby particles from other
// tiles (and periodic images) next to the real particles of each tile.
class MDParticleContainer
    : public amrex::NeighborParticleContainerPureSoA<RealIdx::ncomps, 0>
{
public:
    MDParticleContainer (const amrex::Geometry& geom,
                         const amrex::DistributionMapping& dm,
                         const amrex::BoxArray& ba,
                         int num_neighbor_cells)
        : amrex::NeighborParticleContainerPureSoA<RealIdx::ncomps, 0>(geom, dm, ba,
                                                                      num_neighbor_cells)
    {}

    // nppc^AMREX_SPACEDIM particles per cell on a lattice, with Gaussian velocities.
    void initParticles (int nppc, amrex::Real v_thermal);

    // Acceleration and potential energy from the neighbor list.
    void computeForces (MDParams const& params);

    // v += a dt_kick, then x += v dt_drift.
    void kickDrift (amrex::Real dt_kick, amrex::Real dt_drift);

    [[nodiscard]] MDDiagnostics diagnostics () const;
};

#endif
