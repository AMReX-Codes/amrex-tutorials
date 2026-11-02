#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

// An enum class whose enumerator names can be read from the inputs.
AMREX_ENUM(TimeStepper, ForwardEuler, RK2, RK3);

namespace {

template <typename T>
std::string to_string (std::vector<T> const& v)
{
    std::ostringstream os;
    for (auto const& x : v) { os << x << " "; }
    return os.str();
}

char const* source (bool found) { return found ? "(inputs)" : "(default)"; }

}

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
        amrex::Print() << "nsteps        = " << nsteps << " " << source(found) << "\n";

        std::string title = "ParmParse tutorial";
        found = pp.query("title", title);
        amrex::Print() << "title         = " << title << " " << source(found) << "\n";

        std::vector<int> n_cell{64, 32, 16};
        found = pp.queryarr("n_cell", n_cell);
        amrex::Print() << "n_cell        = " << to_string(n_cell) << source(found) << "\n";

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
        amrex::Print() << "butcher_a     = ";
        for (auto const& row : butcher_a) { amrex::Print() << "{ " << to_string(row) << "} "; }
        amrex::Print() << "\n";

        // A prefix selects the entries of a [physics] table, i.e., physics.*
        amrex::ParmParse pp_phys("physics");
        Real diffusivity = 0.1_rt;
        found = pp_phys.query("diffusivity", diffusivity);
        std::vector<Real> gravity{0.0_rt, -9.81_rt};
        pp_phys.queryarr("gravity", gravity);
        amrex::Print() << "physics.diffusivity = " << diffusivity << " " << source(found) << "\n"
                       << "physics.gravity     = " << to_string(gravity) << "\n";

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

        // Check the values that the shipped inputs (or the defaults) give.
        auto const eps = 1.e-6_rt;
        bool ok = n_cell == std::vector<int>{64, 32, 16}
            && max_grid_size == 32
            && dt == 0.015625_rt
            && dx == 0.0625_rt
            && std::abs(omega - 3.141592653589793_rt) < eps
            && stepper == TimeStepper::RK2
            && butcher_a == std::vector<std::vector<double>>{{0.0, 0.0}, {1.0, 0.0}}
            && std::abs(diffusivity - 0.1_rt) < eps
            && gravity.size() == 2 && std::abs(gravity[1] + 9.81_rt) < eps
            && plot_int == 10
            && plot_file_2 == "plt";
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
