#include <iostream>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include "EMParticleContainer.H"
#include "Evolve.H"
#include "NodalFlags.H"
#include "Constants.H"
#include "IO.H"

using namespace amrex;

struct TestParams
{
    IntVect ncell;      // num cells in domain
    IntVect nppc;       // number of particles per cell in each dim
    int max_grid_size;
    int nsteps;
    int problem_type;
    Real cfl;
    bool write_plot;
    bool write_particles;
};

enum ProblemType {UniformPlasma = 0, Langmuir};

void check_solution(const MultiFab& jx, const Geometry& geom, const Real time);

void test_em_pic(const TestParams& parms)
{
    BL_PROFILE("test_em_pic");
    BL_PROFILE_VAR_NS("evolve_time", blp_evolve);
    BL_PROFILE_VAR_NS("init_time", blp_init);

    BL_PROFILE_VAR_START(blp_init);

    RealBox real_box({AMREX_D_DECL(-20e-6, -20e-6, -20e-6)},
                     {AMREX_D_DECL( 20e-6,  20e-6,  20e-6)});

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.ncell[0]-1,parms.ncell[1]-1,parms.ncell[2]-1));
    const Box domain(domain_lo, domain_hi);

    int coord = 0;
    IntArray is_periodic {AMREX_D_DECL(1,1,1)};
    Geometry geom(domain, real_box, coord, is_periodic);

    BoxArray ba(domain);
    ba.maxSize(parms.max_grid_size);
    DistributionMapping dm(ba);

    const int ng = 1;

    MultiFab Bx(amrex::convert(ba, YeeGrid::Bx_nodal_flag), dm, 1, ng);
    MultiFab By(amrex::convert(ba, YeeGrid::By_nodal_flag), dm, 1, ng);
    MultiFab Bz(amrex::convert(ba, YeeGrid::Bz_nodal_flag), dm, 1, ng);

    MultiFab Ex(amrex::convert(ba, YeeGrid::Ex_nodal_flag), dm, 1, ng);
    MultiFab Ey(amrex::convert(ba, YeeGrid::Ey_nodal_flag), dm, 1, ng);
    MultiFab Ez(amrex::convert(ba, YeeGrid::Ez_nodal_flag), dm, 1, ng);

    MultiFab jx(amrex::convert(ba, YeeGrid::jx_nodal_flag), dm, 1, ng);
    MultiFab jy(amrex::convert(ba, YeeGrid::jy_nodal_flag), dm, 1, ng);
    MultiFab jz(amrex::convert(ba, YeeGrid::jz_nodal_flag), dm, 1, ng);

    Ex.setVal(0.0); Ey.setVal(0.0); Ez.setVal(0.0);
    Bx.setVal(0.0); By.setVal(0.0); Bz.setVal(0.0);
    jx.setVal(0.0); jy.setVal(0.0); jz.setVal(0.0);

    amrex::Print() << "Initializing particles... ";

    int num_species;
    Vector<std::unique_ptr<EMParticleContainer> > particles(2);

    if (parms.problem_type == UniformPlasma)
    {
        num_species = 2;

        EMParticleContainer* electrons;
        EMParticleContainer* protons;

        electrons = new EMParticleContainer(geom, dm, ba, 0, -PhysConst::q_e, PhysConst::m_e);
        protons   = new EMParticleContainer(geom, dm, ba, 1,  PhysConst::q_e, PhysConst::m_p);

        electrons->InitParticles(parms.nppc, 0.01, 10.0, 1e25, real_box, parms.problem_type);
        protons  ->InitParticles(parms.nppc, 0.01, 10.0, 1e25, real_box, parms.problem_type);

        particles[0].reset(electrons);
        particles[1].reset(protons);
    }
    else if (parms.problem_type == Langmuir)
    {
        num_species = 1;

        EMParticleContainer* electrons;

        electrons = new EMParticleContainer(geom, dm, ba, 0, -PhysConst::q_e, PhysConst::m_e);

        electrons->InitParticles(parms.nppc, 0.01, 10.0, 1e25, real_box, parms.problem_type);

        particles[0].reset(electrons);
    }

    amrex::Print() << "Done. " << std::endl;

    for (int i = 0; i < num_species; ++i) particles[i]->OK();

    amrex::Print() << "Starting main PIC loop... " << std::endl;

    int nsteps = parms.nsteps;
    const Real dt = compute_dt(geom, parms.cfl);
    bool synchronized = true;

    BL_PROFILE_VAR_STOP(blp_init);

    BL_PROFILE_VAR_START(blp_evolve);

    Real time = 0.0;
    for (int step = 0; step < nsteps; ++step)
    {
        amrex::Print() << "    Time step: " <<  step << std::endl;

        if (synchronized)
        {
            for (int i = 0; i < num_species; ++i)
            {
                // u(time/dt = -0.5). For readability only since fields null.
                particles[i]->PushParticleMomenta(Ex, Ey, Ez, Bx, By, Bz, -0.5*dt);
            }
            synchronized = false;
        }
        else
        {
            fill_boundary_electric_field(Ex, Ey, Ez, geom);
            // B(time/dt = step).
            evolve_magnetic_field(Ex, Ey, Ez, Bx, By, Bz, geom, 0.5*dt);
            fill_boundary_magnetic_field(Bx, By, Bz, geom);
        }

        jx.setVal(0.0); jy.setVal(0.0), jz.setVal(0.0);
        for (int i = 0; i < num_species; ++i) {
            // u(time/dt = step + 0.5). x(time/dt = step + 1). j(time/dt = step + 0.5).
            particles[i]->PushAndDeposeParticles(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, dt);
        }
        jx.SumBoundary(geom.periodicity());
        jy.SumBoundary(geom.periodicity());
        jz.SumBoundary(geom.periodicity());

        // B(time/dt = step + 0.5).
        evolve_magnetic_field(Ex, Ey, Ez, Bx, By, Bz, geom, 0.5*dt);
        fill_boundary_magnetic_field(Bx, By, Bz, geom);

        // E(time/dt = step + 1).
        evolve_electric_field(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, geom, dt);

        if (step == nsteps - 1)
        {
            // B(time/dt = step + 1).
            evolve_magnetic_field(Ex, Ey, Ez, Bx, By, Bz, geom, 0.5*dt);
            for (int i = 0; i < num_species; ++i)
            {
                // u(time/dt = step + 1).
                particles[i]->PushParticleMomenta(Ex, Ey, Ez, Bx, By, Bz, 0.5*dt);
            }
            synchronized = true;
        }

        for (int i = 0; i < num_species; ++i)
        {
            particles[i]->RedistributeLocal();
        }

        if (parms.problem_type == Langmuir)
        {
            check_solution(jx, geom, time + 0.5*dt);
        }

        time += dt;
    }

    amrex::Print() << "Done. " << std::endl;

    BL_PROFILE_VAR_STOP(blp_evolve);

    if (parms.problem_type == UniformPlasma)
    {
        amrex::Print() << "Not computing error - no exact solution" << std::endl;
    }

    if (parms.write_plot)
    {
        if (parms.write_particles)
        {
            // Since, with ParaView, the particle file will be indexed over
            // the naturals, starting from zero, we do the same with the mesh
            // data. This is useful when the user executes the code in
            // succession with nsteps = 1, 2, 3, 4, 5... and wishes to
            // visualise the fields in tandem with the particles.
            // Note that the current is half a step behind.
            WritePlotFile(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, geom, nsteps - 1, nsteps - 1);
            if (parms.problem_type == UniformPlasma)
            {
                WriteParticleFile(*particles[0], "electrons", nsteps - 1);
                WriteParticleFile(*particles[1], "protons",   nsteps - 1);
            }
            else if (parms.problem_type == Langmuir)
            {
                WriteParticleFile(*particles[0], "electrons", nsteps - 1);
            }
        }
        else
        {
            WritePlotFile(Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, geom, time, nsteps);
        }
    }
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {

    amrex::InitRandom(451);

    ParmParse pp;
    TestParams parms;

    pp.get("ncell", parms.ncell);
    pp.get("nppc",  parms.nppc);
    pp.get("max_grid_size", parms.max_grid_size);
    pp.get("nsteps", parms.nsteps);
    pp.get("write_plot", parms.write_plot);
    pp.get("write_particles", parms.write_particles);
    pp.get("cfl", parms.cfl);

    std::string problem_name;
    pp.get("problem_type", problem_name);
    if (problem_name == "UniformPlasma")
    {
        parms.problem_type = UniformPlasma;
    }
    else if (problem_name == "Langmuir")
    {
        parms.problem_type = Langmuir;
    }
    else
    {
        amrex::Abort("Problem not recognized");
    }

    test_em_pic(parms);

    }

    amrex::Finalize();
}
