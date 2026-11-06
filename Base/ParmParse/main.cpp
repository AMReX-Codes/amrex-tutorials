#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <sstream>
#include <string>
#include <vector>

// An enum class whose enumerator names can be read from the inputs.
AMREX_ENUM(TimeStepper, ForwardEuler, RK2, RK3);

int main (int argc, char* argv[])
{
    // Builds the ParmParse database from the inputs file (argv[1]) and the
    // remaining `name=value` command-line arguments, which take precedence.
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        using amrex::Real;
        using namespace amrex::literals;

        // Unknown symbols in expressions are also looked up with the
        // "my_constants." prefix.
        amrex::ParmParse pp("", "my_constants");

        // query() leaves the variable unchanged if the name is not found and
        // returns whether it was found. get() aborts if it is not found.
        int nsteps = 100;
        bool found = pp.query("nsteps", nsteps);
        amrex::Print() << "nsteps        = " << nsteps << (found ? " (inputs)" : " (default)") << "\n";

        std::string title = "ParmParse tutorial";
        found = pp.query("title", title);
        amrex::Print() << "title         = " << title << (found ? " (inputs)" : " (default)") << "\n";

        // queryarr only overwrites as many entries as are given and keeps the
        // rest, so read into an empty vector to replace a default array.
        std::vector<int> n_cell{64, 32, 16};
        std::vector<int> n_cell_in;
        found = pp.queryarr("n_cell", n_cell_in);
        if (found) { n_cell = n_cell_in; }
        amrex::Print() << "n_cell        =";
        for (int n : n_cell) { amrex::Print() << " " << n; }
        amrex::Print() << (found ? " (inputs)" : " (default)") << "\n";

        // Math expressions
        int max_grid_size = 32;
        pp.query("max_grid_size", max_grid_size);
        Real dt = 0.015625_rt;
        pp.query("dt", dt);
        Real dx = 0.0625_rt;
        pp.query("dx", dx);
        // An expression containing spaces needs queryWithParser.
        Real omega = 3.141592653589793_rt;
        pp.queryWithParser("omega", omega);
        amrex::Print() << "max_grid_size = " << max_grid_size << "\n"
                       << "dt            = " << dt << "\n"
                       << "dx            = " << dx << "\n"
                       << "omega         = " << omega << "\n";

        // Enum
        auto stepper = TimeStepper::RK2;
        pp.query_enum_case_insensitive("time_stepper", stepper);
        amrex::Print() << "time_stepper  = " << amrex::getEnumNameString(stepper) << "\n";

        // 2D table
        std::vector<std::vector<double>> butcher_a{{0.0, 0.0}, {1.0, 0.0}};
        pp.querytable("butcher_a", butcher_a);
        amrex::Print() << "butcher_a     =";
        for (auto const& row : butcher_a) {
            amrex::Print() << " {";
            for (double a : row) { amrex::Print() << " " << a; }
            amrex::Print() << " }";
        }
        amrex::Print() << "\n";

        // A prefix selects the entries of a [physics] table, i.e., physics.*
        amrex::ParmParse pp_phys("physics");
        Real diffusivity = 0.1_rt;
        found = pp_phys.query("diffusivity", diffusivity);
        amrex::Print() << "physics.diffusivity = " << diffusivity << (found ? " (inputs)" : " (default)") << "\n";
        std::vector<Real> gravity{0.0_rt, -9.81_rt};
        std::vector<Real> gravity_in;
        found = pp_phys.queryarr("gravity", gravity_in);
        if (found) { gravity = gravity_in; }
        amrex::Print() << "physics.gravity     =";
        for (Real g : gravity) { amrex::Print() << " " << g; }
        amrex::Print() << "\n";

        // queryAdd() also adds the default to the database if the name is
        // not found, so later queries and the printed table see it.
        amrex::ParmParse pp_out("output");
        int plot_int = 10;
        pp_out.queryAdd("plot_int", plot_int);
        std::string plot_file = "plt";
        pp_out.queryAdd("plot_file", plot_file);
        std::string plot_file_2;
        pp_out.get("plot_file", plot_file_2); // always found now
        amrex::Print() << "output.plot_int  = " << plot_int << "\n"
                       << "output.plot_file = " << plot_file_2 << "\n";

        // Inputs that were never queried are often typos. Also see
        // prettyPrintTable and prettyPrintUsedInputs, whose output is
        // valid ParmParse input.
        std::ostringstream unused;
        amrex::ParmParse::prettyPrintUnusedInputs(unused);
        amrex::Print() << "\nUnused inputs:\n" << unused.str() << "\n";

        // Sanity check: the array and the expression from the shipped inputs.
        bool ok = n_cell == std::vector<int>{64, 32, 16} && dx == 0.0625_rt;
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
