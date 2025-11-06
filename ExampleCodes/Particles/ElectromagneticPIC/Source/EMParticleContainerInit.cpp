#include "EMParticleContainer.H"
#include "Constants.H"

using namespace amrex;

namespace
{
    AMREX_GPU_HOST_DEVICE void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
    {
        int nx = nppc[0];
        int ny = nppc[1];
        int nz = nppc[2];

        int ix_part = i_part/(ny * nz);
        int iy_part = (i_part % (ny * nz)) % ny;
        int iz_part = (i_part % (ny * nz)) / ny;

        r[0] = (0.5+ix_part)/nx;
        r[1] = (0.5+iy_part)/ny;
        r[2] = (0.5+iz_part)/nz;
    }

    AMREX_GPU_HOST_DEVICE void get_gaussian_random_momentum(Real* u, Real u_mean, Real u_std,
                                                            amrex::RandomEngine const& engine) {
        Real ux_th = amrex::RandomNormal(0.0, u_std, engine);
        Real uy_th = amrex::RandomNormal(0.0, u_std, engine);
        Real uz_th = amrex::RandomNormal(0.0, u_std, engine);

        u[0] = u_mean + ux_th;
        u[1] = u_mean + uy_th;
        u[2] = u_mean + uz_th;
    }
}

EMParticleContainer::
EMParticleContainer(const Geometry            & a_geom,
                    const DistributionMapping & a_dmap,
                    const BoxArray            & a_ba,
                    const int                   a_species_id,
                    const Real                  a_charge,
                    const Real                  a_mass)
    : ParticleContainer<0, 0, PIdx::nattribs, 0>(a_geom, a_dmap, a_ba),
    m_species_id(a_species_id), m_charge(a_charge), m_mass(a_mass)
{}

void
EMParticleContainer::
InitParticles(const IntVect& a_num_particles_per_cell,
              const Real     a_thermal_momentum_std,
              const Real     a_thermal_momentum_mean,
              const Real     a_density,
              const RealBox& a_bounds,
              const int      a_problem)
{
    BL_PROFILE("EMParticleContainer::InitParticles");

    const int lev = 0;
    const auto dx = Geom(lev).CellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();

    const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0],
                                      *a_num_particles_per_cell[1],
                                      *a_num_particles_per_cell[2]);
    const Real scale_fac = dx[0]*dx[1]*dx[2]/num_ppc;

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();

        Gpu::DeviceVector<unsigned int> counts(tile_box.numPts(), 0);
        Array4<unsigned int> counts_arr {counts.dataPtr(), amrex::begin(tile_box),
                                         amrex::end(tile_box), 1};

        Gpu::DeviceVector<unsigned int> offsets(tile_box.numPts());
        Array4<unsigned int> offsets_arr {offsets.dataPtr(), amrex::begin(tile_box),
                                          amrex::end(tile_box), 1};

        amrex::ParallelFor(tile_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int i_part=0; i_part<num_ppc;i_part++)
            {
                Real r[3];

                get_position_unit_cell(r, a_num_particles_per_cell, i_part);

                Real x = plo[0] + (i + r[0])*dx[0];
                Real y = plo[1] + (j + r[1])*dx[1];
                Real z = plo[2] + (k + r[2])*dx[2];

                if (x >= a_bounds.hi(0) || x < a_bounds.lo(0) ||
                    y >= a_bounds.hi(1) || y < a_bounds.lo(1) ||
                    z >= a_bounds.hi(2) || z < a_bounds.lo(2) ) continue;

                counts_arr(i, j, k) += 1;
            }
        });

        int num_to_add = Scan::ExclusiveSum(counts.size(), counts.data(), offsets.data());

        auto& particle_tile = DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());

        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + num_to_add;
        particle_tile.resize(new_size);

        if (num_to_add == 0) continue;

        auto ptd = particle_tile.getParticleTileData()

        int procID = ParallelDescriptor::MyProc();

        amrex::ParallelForRNG(tile_box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, amrex::RandomEngine const& engine) noexcept
        {
            int pidx = offsets_arr(i, j, k) + old_size;

            for (int i_part=0; i_part<num_ppc;i_part++)
            {
                Real r[3];
                Real u[3];

                get_position_unit_cell(r, a_num_particles_per_cell, i_part);

                Real x = plo[0] + (i + r[0])*dx[0];
                Real y = plo[1] + (j + r[1])*dx[1];
                Real z = plo[2] + (k + r[2])*dx[2];

                if (a_problem == 0) {
                    get_gaussian_random_momentum(u, a_thermal_momentum_mean,
                                                 a_thermal_momentum_std,
                                                 engine);
                }
                else if (a_problem == 1 ) {
                    u[0] = 0.01;
                    u[1] = 0.0;
                    u[2] = 0.0;
                } else {
                    amrex::Abort("problem type not valid");
                }

                if (x >= a_bounds.hi(0) || x < a_bounds.lo(0) ||
                    y >= a_bounds.hi(1) || y < a_bounds.lo(1) ||
                    z >= a_bounds.hi(2) || z < a_bounds.lo(2) ) continue;

                ptd.id(pidx) = pidx + 1;
                ptd.cpu(pidx) = procID;
                ptd.pos(0, pidx) = x;
                ptd.pos(1, pidx) = y;
                ptd.pos(2, pidx) = z;

                ptd.rdata(PIdx::ux  )[pidx] = u[0] * PhysConst::c;
                ptd.rdata(PIdx::uy  )[pidx] = u[1] * PhysConst::c;
                ptd.rdata(PIdx::uz  )[pidx] = u[2] * PhysConst::c;
                ptd.rdata(PIdx::w   )[pidx] = a_density * scale_fac;
                ptd.rdata(PIdx::Ex  )[pidx] = 0.0;
                ptd.rdata(PIdx::Ey  )[pidx] = 0.0;
                ptd.rdata(PIdx::Ez  )[pidx] = 0.0;
                ptd.rdata(PIdx::Bx  )[pidx] = 0.0;
                ptd.rdata(PIdx::By  )[pidx] = 0.0;
                ptd.rdata(PIdx::Bz  )[pidx] = 0.0;
                ptd.rdata(PIdx::ginv)[pidx] = 0.0;

                ++pidx;
            }
        });
    }

    AMREX_ASSERT(OK());
}
