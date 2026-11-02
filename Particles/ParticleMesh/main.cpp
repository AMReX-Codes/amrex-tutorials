#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleInterpolators.H>
#include <AMReX_ParticleMesh.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <string>
#include <type_traits>

using namespace amrex::literals;

// Positions first, then the particle mass and the field gathered from the mesh.
struct RealIdx
{
    enum {
        AMREX_D_DECL(x, y, z),
        mass,
        AMREX_D_DECL(Ex, Ey, Ez),
        ncomps
    };
};

using PC = amrex::ParticleContainerPureSoA<RealIdx::ncomps, 0>;
using PTD = PC::PTDType;
using ConstPTD = PC::ConstPTDType;
using ConstParticle = PC::ParticleType::ConstType;

// A linear vector field, E_d(x) = (d+1) * (1 + x + 2y + 3z). Linear
// interpolation reproduces it exactly.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real linearField (int d, AMREX_D_DECL(amrex::Real x, amrex::Real y, amrex::Real z))
{
    return amrex::Real(d+1) * (1.0_rt + AMREX_D_TERM(x, + 2.0_rt*y, + 3.0_rt*z));
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int particle_max_grid_size = 8;
        amrex::Long num_particles = 100000;
        bool write_plotfile = true;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("particle_max_grid_size", particle_max_grid_size);
            pp.query("num_particles", num_particles);
            pp.query("write_plotfile", write_plotfile);
        }

        const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        const amrex::RealBox real_box({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                      {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
        const amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
        const amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        // Mesh data and particles may live on different grids ("dual grid").
        // ParticleToMesh and MeshToParticle handle the communication.
        amrex::BoxArray mesh_ba(domain);
        mesh_ba.maxSize(max_grid_size);
        const amrex::DistributionMapping mesh_dm(mesh_ba);

        amrex::BoxArray particle_ba(domain);
        particle_ba.maxSize(particle_max_grid_size);
        const amrex::DistributionMapping particle_dm(particle_ba);

        PC pc(geom, particle_dm, particle_ba);

        // Randomly placed particles, all with the same mass. For pure SoA
        // particles, the real_array_data includes (ignored) position slots.
        const amrex::Real particle_mass = 1.0_rt / amrex::Real(num_particles);
        PC::ParticleInitData pdata = {{}, {},
            {AMREX_D_DECL(0.0, 0.0, 0.0), particle_mass, AMREX_D_DECL(0.0, 0.0, 0.0)}, {}};
        pc.InitRandom(num_particles, 451, pdata);

        const amrex::Long np_total = pc.TotalNumberOfParticles();
        amrex::Real total_mass = amrex::ReduceSum(pc,
            [=] AMREX_GPU_HOST_DEVICE (const ConstPTD& ptd, int i) noexcept -> amrex::Real
            {
                return ptd.rdata(RealIdx::mass)[i];
            });
        amrex::ParallelDescriptor::ReduceRealSum(total_mass);
        amrex::Print() << "Number of particles: " << np_total << "\n";

        const auto plo = geom.ProbLoArray();
        const auto dx = geom.CellSizeArray();
        const auto dxi = geom.InvCellSizeArray();
        const amrex::Real inv_cell_volume = AMREX_D_TERM(dxi[0], *dxi[1], *dxi[2]);

        // Deposit the mass density with cloud-in-cell (linear) weights. A
        // particle touches 2^AMREX_SPACEDIM cells, so rho needs one ghost cell;
        // ParticleToMesh sums ghost cell contributions into the valid cells.
        amrex::MultiFab rho(mesh_ba, mesh_dm, 1, 1);
        amrex::ParticleToMesh(pc, rho, 0,
            [=] AMREX_GPU_DEVICE (const ConstPTD& ptd, int i,
                                  amrex::Array4<amrex::Real> const& rho_arr) noexcept
            {
                const ConstParticle p = ptd[i];
                amrex::ParticleInterpolator::Linear interp(p, plo, dxi);
                // The callable returns the value to deposit (uses atomic adds).
                interp.ParticleToMesh(p, rho_arr, 0, 0, 1,
                    [=] (ConstParticle const&, int) noexcept
                    {
                        return ptd.rdata(RealIdx::mass)[i] * inv_cell_volume;
                    });
            });

        const amrex::Real deposited_mass = rho.sum(0) / inv_cell_volume;

        // Fill a vector field on the mesh. This linear field is not periodic,
        // so we set the ghost cells directly instead of calling FillBoundary.
        amrex::MultiFab efield(mesh_ba, mesh_dm, AMREX_SPACEDIM, 1);
        auto const& ef = efield.arrays();
        amrex::ParallelFor(efield, efield.nGrowVect(), AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
            {
                amrex::ignore_unused(j, k);
                ef[b](i,j,k,n) = linearField(n, AMREX_D_DECL(plo[0] + (i+0.5_rt)*dx[0],
                                                              plo[1] + (j+0.5_rt)*dx[1],
                                                              plo[2] + (k+0.5_rt)*dx[2]));
            });

        // Gather the field to the particles with the same linear weights.
        amrex::MeshToParticle(pc, efield, 0,
            [=] AMREX_GPU_DEVICE (const PTD& ptd, int i,
                                  amrex::Array4<const amrex::Real> const& e_arr) noexcept
            {
                auto p = ptd[i];
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    ptd.rdata(RealIdx::Ex+d)[i] = 0;
                }
                amrex::ParticleInterpolator::Linear interp(p, plo, dxi);
                // The first callable reads the mesh, the second updates the
                // particle with each weighted contribution.
                interp.MeshToParticle(p, e_arr, 0, 0, AMREX_SPACEDIM,
                    [=] (amrex::Array4<const amrex::Real> const& arr,
                         int ii, int jj, int kk, int comp) noexcept
                    {
                        return arr(ii, jj, kk, comp);
                    },
                    [=] (PC::ParticleType&, int comp, amrex::Real val) noexcept
                    {
                        ptd.rdata(RealIdx::Ex+comp)[i] += static_cast<amrex::ParticleReal>(val);
                    });
            });

        // Largest error of the gathered field.
        amrex::Real gather_error = amrex::ReduceMax(pc,
            [=] AMREX_GPU_HOST_DEVICE (const ConstPTD& ptd, int i) noexcept -> amrex::Real
            {
                amrex::Real err = 0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    const amrex::Real exact = linearField(d, AMREX_D_DECL(ptd.pos(0,i),
                                                                          ptd.pos(1,i),
                                                                          ptd.pos(2,i)));
                    err = amrex::max(err, std::abs(ptd.rdata(RealIdx::Ex+d)[i] - exact)
                                     / std::abs(exact));
                }
                return err;
            });
        amrex::ParallelDescriptor::ReduceRealMax(gather_error);

        const amrex::Real mass_error = std::abs(deposited_mass - total_mass) / total_mass;
        amrex::Print() << "Total particle mass: " << total_mass << "\n"
                       << "Deposited mass:      " << deposited_mass << "\n"
                       << "Relative mass error:   " << mass_error << "\n"
                       << "Relative gather error: " << gather_error << "\n";

        if (write_plotfile) {
            amrex::MultiFab plotmf(mesh_ba, mesh_dm, 1+AMREX_SPACEDIM, 0);
            amrex::MultiFab::Copy(plotmf, rho, 0, 0, 1, 0);
            amrex::MultiFab::Copy(plotmf, efield, 0, 1, AMREX_SPACEDIM, 0);
            amrex::WriteSingleLevelPlotfile("plt00000", plotmf,
                                            {"rho", AMREX_D_DECL("Ex", "Ey", "Ez")}, geom, 0.0, 0);
            const amrex::Vector<std::string> real_names{"mass", AMREX_D_DECL("Ex", "Ey", "Ez")};
            pc.WritePlotFile("plt00000", "particles", real_names, {});
        }

        constexpr bool is_single = std::is_same_v<amrex::Real, float>
                                || std::is_same_v<amrex::ParticleReal, float>;
        // Summing many cells in a different order than the particles gives
        // a slightly larger round-off error for the mass.
        constexpr amrex::Real mass_tol = is_single ? 1.e-3_rt : 1.e-10_rt;
        constexpr amrex::Real gather_tol = is_single ? 1.e-4_rt : 1.e-12_rt;
        if (mass_error > mass_tol || gather_error > gather_tol) {
            amrex::Print() << "FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "PASSED\n";
        }
    }
    amrex::Finalize();
    return status;
}
