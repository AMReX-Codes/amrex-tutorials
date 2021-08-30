#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

void test_parameters ();

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    test_parameters();

    amrex::Finalize();
}

void test_parameters ()
{
    {
        amrex::ParmParse pp;
        int i;
        bool b;
        std::vector<amrex::Real> ra;
        pp.get("an_int_scalar", i);
        pp.get("a_bool_scalar",b);
        pp.getarr("a_real_array", ra);
        amrex::Print() << "an_int_scalar = " << i << "\n"
                       << "a_bool_scalar = " << b << "\n";
        amrex::Print() << "a_real_array = ";
        for (auto x : ra) {
            amrex::Print() << x << " ";
        }
        amrex::Print() << "\n";
    }

    {
        amrex::ParmParse pp("a_prefix");
        std::vector<int> ia;
        amrex::Real r;
        std::string s;
        pp.getarr("an_int_array", ia);
        pp.get("a_real_scalar", r);
        pp.get("a_string", s);
        amrex::Print() << "an_int_array = ";
        for (auto x : ia) {
            amrex::Print() << x << " ";
        }
        amrex::Print() << "\n";
        amrex::Print() << "a_prefix.a_real_scalar = " << r << "\n"
                       << "a_prefix.a_string = " << s << "\n";
    }
}
