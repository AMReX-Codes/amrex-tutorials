#include <AMReX.H>
#include <AMReX_DenseBins.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_ParticleTransformation.H>
#include <AMReX_ParticleUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <cmath>

using namespace amrex::literals;

struct RealIdx
{
    enum {
        AMREX_D_DECL(x, y, z),
        AMREX_D_DECL(vx, vy, vz),
        ncomps
    };
};

using PC = amrex::ParticleContainerPureSoA<RealIdx::ncomps, 0>;
using ParticleType = PC::ParticleType;
using PTD = PC::PTDType;
using ConstPTD = PC::ConstPTDType;

// Squared speed of particle i; works for both PTD and ConstPTD.
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::ParticleReal speed2 (P const& ptd, int i)
{
    amrex::ParticleReal v2 = 0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        v2 += ptd.rdata(RealIdx::vx+d)[i] * ptd.rdata(RealIdx::vx+d)[i];
    }
    return v2;
}

// Random positions and velocities, so particles start in no particular order.
void initParticles (PC& pc, int nppc)
{
    const int lev = 0;
    const auto plo = pc.Geom(lev).ProbLoArray();
    const auto dx = pc.Geom(lev).CellSizeArray();
    const int myproc = amrex::ParallelDescriptor::MyProc();

    for (amrex::MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const amrex::Box bx = mfi.tilebox();
        const amrex::Long np = bx.numPts() * nppc;
        const amrex::Long id0 = ParticleType::NextID();
        ParticleType::NextID(id0 + np);

        auto& tile = pc.DefineAndReturnParticleTile(lev, mfi);
        tile.resize(np);
        auto ptd = tile.getParticleTileData();
        const amrex::IntVect lo = bx.smallEnd();
        const amrex::IntVect len = bx.length();

        amrex::ParallelForRNG(np,
        [=] AMREX_GPU_DEVICE (amrex::Long ip, amrex::RandomEngine const& engine) noexcept
        {
            ptd.idcpu(ip) = amrex::SetParticleIDandCPU(id0 + ip, myproc);
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                // Stay strictly inside the tile so that no Redistribute is needed.
                const amrex::Real r = amrex::Random(engine) * 0.999_rt + 0.0005_rt;
                ptd.pos(d, ip) = static_cast<amrex::ParticleReal>(
                    plo[d] + (lo[d] + r*len[d]) * dx[d]);
                ptd.rdata(RealIdx::vx+d)[ip] = static_cast<amrex::ParticleReal>(
                    amrex::RandomNormal(0.0_rt, 1.0_rt, engine));
            }
        });
    }
}

// Count, over all particles, how often the bin index of particle i+1 is
// smaller than that of particle i. Zero means the tiles are sorted by bin.
amrex::Long countUnsorted (PC& pc, amrex::IntVect const& bin_size)
{
    const int lev = 0;
    const auto plo = pc.Geom(lev).ProbLoArray();
    const auto dxi = pc.Geom(lev).InvCellSizeArray();
    const auto domain = pc.Geom(lev).Domain();

    amrex::ReduceOps<amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Long> reduce_data(reduce_op);
    for (PC::ParIterType pti(pc, lev); pti.isValid(); ++pti)
    {
        const amrex::Box box = pti.validbox();
        const auto ptd = pti.GetParticleTile().getConstParticleTileData();
        auto bin = [=] AMREX_GPU_DEVICE (int i) -> int
        {
            amrex::Box tbx;
            const amrex::IntVect iv = amrex::getParticleCell(ptd, i, plo, dxi, domain);
            return amrex::getTileIndex(iv, box, true, bin_size, tbx);
        };
        reduce_op.eval(pti.numParticles()-1, reduce_data,
            [=] AMREX_GPU_DEVICE (int i) -> amrex::GpuTuple<amrex::Long>
            {
                return { bin(i) > bin(i+1) ? 1 : 0 };
            });
    }
    amrex::Long n = amrex::get<0>(reduce_data.value(reduce_op));
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

// Sum of all velocity components over all particles.
amrex::Real velocitySum (PC const& pc)
{
    amrex::Real vsum = amrex::ReduceSum(pc,
        [=] AMREX_GPU_HOST_DEVICE (const ConstPTD& ptd, int i) noexcept -> amrex::Real
        {
            return AMREX_D_TERM(ptd.rdata(RealIdx::vx)[i],
                                + ptd.rdata(RealIdx::vy)[i],
                                + ptd.rdata(RealIdx::vz)[i]);
        });
    amrex::ParallelDescriptor::ReduceRealSum(vsum);
    return vsum;
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int nppc = 4;
        int bin_size = 4;
        amrex::Real vcut = 1.5_rt;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nppc", nppc);
            pp.query("bin_size", bin_size);
            pp.query("vcut", vcut);
        }

        const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        const amrex::RealBox real_box({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                      {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)});
        const amrex::Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
        const amrex::Geometry geom(domain, real_box, amrex::CoordSys::cartesian, is_periodic);

        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        const amrex::DistributionMapping dm(ba);

        PC pc(geom, dm, ba);
        initParticles(pc, nppc);
        const amrex::Long np_total = pc.TotalNumberOfParticles();
        amrex::Print() << "Number of particles: " << np_total << "\n";

        const int lev = 0;
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();

        // 1. DenseBins: bin the particles of each tile by cell without moving
        //    them. The permutation lists particle indices bin by bin, and
        //    offsets[b] is where bin b starts. Here we count particles per cell
        //    and compare with ParticleContainer::Increment.
        amrex::MultiFab count_bins(ba, dm, 1, 0);
        amrex::MultiFab count_ref(ba, dm, 1, 0);
        count_bins.setVal(0.0);
        count_ref.setVal(0.0);
        pc.Increment(count_ref, lev);

        amrex::DenseBins<ConstPTD> bins;
        for (PC::ParIterType pti(pc, lev); pti.isValid(); ++pti)
        {
            const amrex::Box bx = pti.validbox();
            const auto ptd = pti.GetParticleTile().getConstParticleTileData();
            const amrex::IntVect lo = bx.smallEnd();
            const auto domain_box = geom.Domain();
            // Map each particle to its cell, relative to the lower corner of bx.
            bins.build(pti.numParticles(), ptd, bx,
                [=] AMREX_GPU_DEVICE (const ConstPTD& p, int i) noexcept -> amrex::IntVect
                {
                    return amrex::getParticleCell(p, i, plo, dxi, domain_box) - lo;
                });

            const auto* offsets = bins.offsetsPtr();
            const auto* perm = bins.permutationPtr();
            auto const& cnt = count_bins.array(pti);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int b = static_cast<int>(bx.index(amrex::IntVect(AMREX_D_DECL(i, j, k))));
                // The particles in cell (i,j,k) are perm[offsets[b]:offsets[b+1]].
                int n = 0;
                for (int m = offsets[b]; m < offsets[b+1]; ++m) {
                    const amrex::IntVect iv = amrex::getParticleCell(ptd, perm[m], plo, dxi, domain_box);
                    if (iv == amrex::IntVect(AMREX_D_DECL(i, j, k))) { ++n; }
                }
                cnt(i,j,k) = amrex::Real(n);
            });
            amrex::Gpu::streamSynchronize(); // bins is reused for the next tile
        }
        amrex::MultiFab::Subtract(count_bins, count_ref, 0, 0, 1, 0);
        const amrex::Real bin_error = count_bins.norminf(0);
        amrex::Print() << "DenseBins: max difference in particles per cell = " << bin_error << "\n";

        // 2. Sort the particles in memory, by cell and by groups of cells.
        const amrex::Long unsorted0 = countUnsorted(pc, amrex::IntVect(1));
        pc.SortParticlesByCell();
        const amrex::Long unsorted1 = countUnsorted(pc, amrex::IntVect(1));
        pc.SortParticlesByBin(amrex::IntVect(bin_size));
        const amrex::Long unsorted2 = countUnsorted(pc, amrex::IntVect(bin_size));
        amrex::Print() << "Out-of-order neighbors: initial " << unsorted0
                       << ", after SortParticlesByCell " << unsorted1
                       << ", after SortParticlesByBin " << unsorted2 << "\n";

        // 3. Partition each tile so that fast particles come first.
        const amrex::ParticleReal vcut2 = vcut*vcut;
        amrex::Long nfast = 0;
        amrex::Long partition_errors = 0;
        for (PC::ParIterType pti(pc, lev); pti.isValid(); ++pti)
        {
            auto& tile = pti.GetParticleTile();
            const int nleft = amrex::partitionParticles(tile,
                [=] AMREX_GPU_DEVICE (const PTD& ptd, int i) noexcept
                {
                    return speed2(ptd, i) > vcut2;
                });
            nfast += nleft;

            const auto ptd = tile.getConstParticleTileData();
            partition_errors += amrex::Reduce::Sum<amrex::Long>(tile.numParticles(),
                [=] AMREX_GPU_DEVICE (int i) noexcept -> amrex::Long
                {
                    const bool fast = speed2(ptd, i) > vcut2;
                    return (fast != (i < nleft)) ? 1 : 0;
                });
        }
        amrex::ParallelDescriptor::ReduceLongSum(nfast);
        amrex::ParallelDescriptor::ReduceLongSum(partition_errors);
        amrex::Print() << "partitionParticles: " << nfast << " fast particles, "
                       << partition_errors << " misplaced\n";

        // 4. Copy the slow particles to a second container (filterParticles), and
        //    copy all particles with reversed velocities to a third one
        //    (transformParticles). The destination tiles must be large enough.
        PC slow_pc(geom, dm, ba);
        PC reversed_pc(geom, dm, ba);
        for (PC::ParIterType pti(pc, lev); pti.isValid(); ++pti)
        {
            const auto& src = pti.GetParticleTile();
            const int np = src.numParticles();

            auto& slow = slow_pc.DefineAndReturnParticleTile(lev, pti);
            slow.resize(np);
            const int nslow = amrex::filterParticles(slow, src,
                [=] AMREX_GPU_DEVICE (const ConstPTD& ptd, int i) noexcept
                {
                    return speed2(ptd, i) <= vcut2;
                });
            slow.resize(nslow);

            auto& reversed = reversed_pc.DefineAndReturnParticleTile(lev, pti);
            reversed.resize(np);
            amrex::transformParticles(reversed, src,
                [=] AMREX_GPU_HOST_DEVICE (const PTD& dst, const ConstPTD& s, int si, int di) noexcept
                {
                    dst.idcpu(di) = s.idcpu(si);
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        dst.pos(d, di) = s.pos(d, si);
                        dst.rdata(RealIdx::vx+d)[di] = -s.rdata(RealIdx::vx+d)[si];
                    }
                });
        }

        const amrex::Long nslow = slow_pc.TotalNumberOfParticles();
        amrex::Print() << "filterParticles: " << nslow << " slow particles\n";

        const amrex::Real v0 = velocitySum(pc);
        const amrex::Real v1 = velocitySum(reversed_pc);
        amrex::Print() << "transformParticles: velocity sum " << v0 << " -> " << v1 << "\n";

        const bool ok = bin_error == 0
            && unsorted0 > 0 && unsorted1 == 0 && unsorted2 == 0
            && partition_errors == 0
            && nfast + nslow == np_total
            && reversed_pc.TotalNumberOfParticles() == np_total
            && std::abs(v0 + v1) <= 1.e-3_rt * std::abs(v0) + 1.e-6_rt
            && slow_pc.OK() && reversed_pc.OK();
        if (ok) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
