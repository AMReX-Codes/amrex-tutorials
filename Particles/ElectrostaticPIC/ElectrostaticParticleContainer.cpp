#include <iomanip>

#include "ElectrostaticParticleContainer.H"
#include "AMReX_PlotFileUtil.H"

#include "Particle_Pusher.H"
#include "Electrostatic_PIC_Util.H"
#include "Electrostatic_PIC_K.H"

using namespace amrex;

void ElectrostaticParticleContainer::InitParticles() {

    if ( ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber() ) {

        ParticleType p;

        p.id()   = ParticleType::NextID();
        p.cpu()  = ParallelDescriptor::MyProc();

        p.pos(0) = -2.5e-6;
        p.pos(1) =  0.0;

        std::array<ParticleReal,PIdx::nattribs> attribs;
        attribs[PIdx::w]  = 1.0;
        attribs[PIdx::vx] = 0.0;
        attribs[PIdx::vy] = 0.0;
        attribs[PIdx::Ex] = 0.0;
        attribs[PIdx::Ey] = 0.0;

        // Add to level 0, grid 0, and tile 0
        // Redistribute() will move it to the proper place.
        std::pair<int,int> key {0,0};
        auto& particle_tile = GetParticles(0)[key];

        particle_tile.push_back(p);
        particle_tile.push_back_real(attribs);

    }

    Redistribute();
}

void
ElectrostaticParticleContainer::DepositCharge(ScalarMeshData& rho) {

    int num_levels = rho.size();
    int finest_level = num_levels - 1;

    // each level deposits it's own particles
    const int ng = rho[0]->nGrow();
    for (int lev = 0; lev < num_levels; ++lev) {
        rho[lev]->setVal(0.0, ng);
        const auto& gm = m_gdb->Geom(lev);
        auto plo = gm.ProbLoArray();
        auto dxi = gm.InvCellSizeArray();
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            const Long np  = pti.numParticles();
            const auto& wp = pti.GetAttribs(PIdx::w);
            const auto& particles = pti.GetArrayOfStructs();

            amrex::Real q = this->charge;
            auto rhoarr = (*rho[lev])[pti].array();
            const auto wp_ptr = wp.data();
            const auto p_ptr = particles().data();
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept {
                          deposit_cic(p_ptr[i], wp_ptr[i], q, rhoarr, plo, dxi);
                                   });
        }

        rho[lev]->SumBoundary(gm.periodicity());
    }

    // now we average down fine to crse
    IntVect ratio(D_DECL(2, 2, 2));  // FIXME
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        sumFineToCrseNodal(*rho[lev+1], *rho[lev], m_gdb->Geom(lev), ratio);
    }

    for (int lev = 0; lev < num_levels; ++lev) {
        rho[lev]->mult(-1.0/PhysConst::ep0, ng);
    }
}

void
ElectrostaticParticleContainer::
FieldGather(const VectorMeshData& E,
            const Vector<std::unique_ptr<FabArray<BaseFab<int> > > >& masks) {

    const int num_levels = E.size();
    const int ng = E[0][0]->nGrow();

    if (num_levels == 1) {
        const int lev = 0;
        const auto& gm = m_gdb->Geom(lev);
        auto plo = gm.ProbLoArray();
        auto dxi = gm.InvCellSizeArray();
        AMREX_ASSERT(OnSameGrids(lev, *E[lev][0]));

        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            auto& particles = pti.GetArrayOfStructs();
            auto p_ptr = particles().data();
            int nstride = particles.dataShape().first;
            const Long np  = pti.numParticles();

            auto& attribs = pti.GetAttribs();
            auto Ex_p = attribs[PIdx::Ex].data();
            auto Ey_p = attribs[PIdx::Ey].data();

            const auto& exarr = (*E[lev][0])[pti].array();
            const auto& eyarr = (*E[lev][1])[pti].array();

            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept {
                                       interpolate_cic(p_ptr[i], Ex_p[i], Ey_p[i],
                                                       exarr, eyarr, plo, dxi);
                                   });
        }
        return;
    }

    const BoxArray& fine_BA = E[1][0]->boxArray();
    const DistributionMapping& fine_dm = E[1][0]->DistributionMap();
    BoxArray coarsened_fine_BA = fine_BA;
    coarsened_fine_BA.coarsen(IntVect(D_DECL(2,2,2)));

    MultiFab coarse_Ex(coarsened_fine_BA, fine_dm, 1, 1);
    MultiFab coarse_Ey(coarsened_fine_BA, fine_dm, 1, 1);

    coarse_Ex.ParallelCopy(*E[0][0], 0, 0, 1, 1, 1);
    coarse_Ey.ParallelCopy(*E[0][1], 0, 0, 1, 1, 1);

    for (int lev = 0; lev < num_levels; ++lev) {
        const auto& gm = Geom(lev);

        BL_ASSERT(OnSameGrids(lev, *E[lev][0]));

        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const Long np  = pti.numParticles();

            auto& attribs = pti.GetAttribs();
            auto Ex_p = attribs[PIdx::Ex].data();
            auto Ey_p = attribs[PIdx::Ey].data();

            const auto exarr = (*E[lev][0])[pti].array();
            const auto eyarr = (*E[lev][1])[pti].array();

            auto p_ptr = particles().data();
            auto ploarr = gm.ProbLoArray();
            auto dxi = gm.InvCellSizeArray();

            if (lev == 0) {
                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept {
                                           interpolate_cic(p_ptr[i], Ex_p[i], Ey_p[i],
                                                           exarr, eyarr, ploarr, dxi);
                                       });
            } else {
                const auto& cgm = Geom(lev-1);
                auto cdxi = cgm.InvCellSizeArray();

                const auto cexarr = coarse_Ex[pti].array();
                const auto ceyarr = coarse_Ey[pti].array();

                const auto maskarr = (*masks[lev])[pti].array();

                amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept {
                    interpolate_cic_two_levels(p_ptr[i], Ex_p[i], Ey_p[i],
                                               exarr, eyarr, cexarr, ceyarr, maskarr,
                                               ploarr, dxi, cdxi, lev);
                                       });
            }
        }
    }
}

void ElectrostaticParticleContainer:: Evolve (const VectorMeshData& E, ScalarMeshData& rho,
                                              const Real& dt) {

    const int num_levels = E.size();

    for (int lev = 0; lev < num_levels; ++lev) {

        const auto& gm = m_gdb->Geom(lev);
        const RealBox& prob_domain = gm.ProbDomain();

        BL_ASSERT(OnSameGrids(lev, *rho[lev]));
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {

            // Particle structs
            auto& particles = pti.GetArrayOfStructs();
            const Long np  = pti.numParticles();

            // Particle attributes
            auto& attribs = pti.GetAttribs();
            auto vxp = attribs[PIdx::vx].data();
            auto vyp = attribs[PIdx::vy].data();

            auto Exp = attribs[PIdx::Ex].data();
            auto Eyp = attribs[PIdx::Ey].data();

            auto p_ptr = particles().data();
            const auto plo = gm.ProbLoArray();
            const auto phi = gm.ProbHiArray();
            amrex::Real q = this->charge;
            amrex::Real m = this->mass;
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept {
                                       push_leapfrog(p_ptr[i].pos(0), p_ptr[i].pos(1),
                                                     vxp[i], vyp[i], Exp[i], Eyp[i],
                                                     q, m, dt, plo, phi);
                                   });
        }
    }
}

void ElectrostaticParticleContainer::pushX (const Real& dt) {
    for (int lev = 0; lev <= finestLevel(); ++lev) {
        const auto& gm = m_gdb->Geom(lev);
        const RealBox& prob_domain = gm.ProbDomain();
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const Long np  = pti.numParticles();

            auto& attribs = pti.GetAttribs();
            auto vxp = attribs[PIdx::vx].data();
            auto vyp = attribs[PIdx::vy].data();

            auto p_ptr = particles().data();
            auto plo = gm.ProbLoArray();
            auto phi = gm.ProbHiArray();
            amrex::Real q = this->charge;
            amrex::Real m = this->mass;
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept {
                                       push_leapfrog_positions(p_ptr[i].pos(0), p_ptr[i].pos(1),
                                                               vxp[i], vyp[i], dt, plo, phi);
                                   });

        }
    }
}

void ElectrostaticParticleContainer::writeParticles(int n) {
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}
