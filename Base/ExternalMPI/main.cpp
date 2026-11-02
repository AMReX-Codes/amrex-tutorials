#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <mpi.h>

#include <string>

int main (int argc, char* argv[])
{
    // The application initializes MPI itself.
    MPI_Init(&argc, &argv);
    int world_rank = 0;
    int world_size = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Split the ranks into two groups (even and odd world ranks). Each group
    // runs its own AMReX on its own communicator.
    int const group = world_rank % 2;
    MPI_Comm group_comm;
    MPI_Comm_split(MPI_COMM_WORLD, group, world_rank, &group_comm);
    int group_size = 0;
    MPI_Comm_size(group_comm, &group_size);

    // The application handles its own command line; here the first argument
    // is an inputs file, unless it looks like `name=value`.
    std::string inputs_file;
    if (argc > 1 && std::string(argv[1]).find('=') == std::string::npos) {
        inputs_file = argv[1];
    }

    // Called by amrex::Initialize after the (empty) ParmParse database is
    // created, but before AMReX reads its own parameters.
    auto add_parameters = [&] ()
    {
        amrex::ParmParse pp_amrex("amrex");
        pp_amrex.add("fpe_trap_invalid", 1); // trap NaNs
        pp_amrex.add("verbose", 0);
        // On GPU, AMReX reserves most of the device memory up front, shared
        // among the ranks of its own communicator. It cannot see the other
        // group's ranks on the same GPU, so do not reserve memory up front.
        pp_amrex.add("the_arena_init_size", 0);

        if (!inputs_file.empty()) {
            amrex::ParmParse::addfile(inputs_file);
        }

        amrex::ParmParse pp;
        pp.add("group", group);
    };

    // build_parm_parse = false: do not parse argc/argv into ParmParse.
    // AMReX uses a duplicate of group_comm.
    amrex::Initialize(argc, argv, false, group_comm, add_parameters);
    int status = 0;
    {
        int n_cell = 16;
        int my_group = -1;
        amrex::ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.get("group", my_group);

        int const nprocs = amrex::ParallelDescriptor::NProcs();
        amrex::Print() << "Group " << my_group << ": AMReX runs on " << nprocs
                       << " of " << world_size << " MPI ranks\n";

        // A MultiFab distributed over this group's ranks only
        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(n_cell/2);
        amrex::MultiFab mf(ba, amrex::DistributionMapping(ba), 1, 0);
        mf.setVal(amrex::Real(my_group + 1));
        amrex::Real const sum = mf.sum(0); // reduced over this group's ranks
        amrex::Real const exact = amrex::Real(my_group + 1) * amrex::Real(domain.numPts());
        amrex::Print() << "Group " << my_group << ": sum = " << sum
                       << " (exact " << exact << ")\n";

        if (nprocs != group_size || my_group != group || sum != exact) {
            amrex::Print() << "Group " << my_group << ": FAILED\n";
            status = 1;
        } else {
            amrex::Print() << "Group " << my_group << ": PASSED\n";
        }
    }
    amrex::Finalize();

    // The application created these, so it frees them.
    MPI_Comm_free(&group_comm);
    MPI_Finalize();
    return status;
}
