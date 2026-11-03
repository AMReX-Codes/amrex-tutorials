#include <AMReX.H>
#include <AMReX_IParser.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_Print.H>

#include <cmath>
#include <string>
#include <type_traits>

using amrex::Real;
using namespace amrex::literals;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        Real time = 0.5_rt;
        std::string const default_function = "A * sin(2*pi*x) * cos(2*pi*y) * exp(-t)";
        std::string const default_tag = "max(i, max(j, k)) < n/2";
        std::string function = default_function;
        std::string tag = default_tag;
        double A = 2.0;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("time", time);
            pp.query("ic.function", function);
            pp.query("ic.tag", tag);

            // Constants used by the expressions. queryAdd puts the defaults
            // into the database, where makeParser can find them.
            amrex::ParmParse ppc("my_constants");
            ppc.queryAdd("A", A);
            double pi = amrex::Math::pi<double>();
            ppc.queryAdd("pi", pi);
        }

        amrex::Box domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        amrex::DistributionMapping dm(ba);
        amrex::Geometry geom(domain,
                             amrex::RealBox({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                            {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}),
                             amrex::CoordSys::cartesian, {AMREX_D_DECL(1, 1, 1)});
        auto const prob_lo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();

        // 1. Parser built by hand: set the constants, then the variables in
        //    the order they are passed to the compiled executor.
        amrex::Parser parser(function);
        parser.setConstant("pi", amrex::Math::pi<double>());
        parser.setConstant("A", A);
        parser.registerVariables({"x", "y", "z", "t"});
        auto const f1 = parser.compile<4>(); // usable on host and device

        // 2. ParmParse::makeParser looks up unknown symbols in the database,
        //    here with the "my_constants" prefix.
        amrex::ParmParse pp("", "my_constants");
        amrex::Parser parser2 = pp.makeParser(function, {"x", "y", "z", "t"});
        auto const f2 = parser2.compile<4>();

        // 3. IParser for integer expressions, here of cell indices.
        amrex::IParser iparser(tag);
        iparser.setConstant("n", n_cell);
        iparser.registerVariables({"i", "j", "k"});
        auto const ftag = iparser.compile<3>();

        amrex::Print() << "ic.function = " << function << "\n"
                       << "ic.tag      = " << tag << "\n";

        // Evaluate on the mesh; the executors are captured by value. Also
        // compare with a hard-coded C++ version of the default expressions.
        amrex::MultiFab phi(ba, dm, 1, 0);
        amrex::iMultiFab tagged(ba, dm, 1, 0);
        amrex::MultiFab err(ba, dm, 3, 0);
        auto const& pa = phi.arrays();
        auto const& ta = tagged.arrays();
        auto const& ea = err.arrays();
        Real const A_r = Real(A);
        amrex::ParallelFor(phi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            amrex::IntVect const iv(AMREX_D_DECL(i,j,k));
            amrex::GpuArray<Real,3> x{0.0_rt, 0.0_rt, 0.0_rt}; // y = z = 0 in lower dimensions
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                x[d] = prob_lo[d] + (iv[d] + 0.5_rt) * dx[d];
            }
            pa[b](i,j,k) = Real(f1(x[0], x[1], x[2], time));
            ta[b](i,j,k) = int(ftag(i, j, k));

            Real const twopi = 2.0_rt * amrex::Math::pi<Real>();
            Real const ref = A_r * std::sin(twopi*x[0]) * std::cos(twopi*x[1]) * std::exp(-time);
            int const n = n_cell;
            int const ref_tag = amrex::max(i, amrex::max(j, k)) < n/2;
            ea[b](i,j,k,0) = std::abs(pa[b](i,j,k) - ref);
            ea[b](i,j,k,1) = std::abs(pa[b](i,j,k) - Real(f2(x[0], x[1], x[2], time)));
            ea[b](i,j,k,2) = Real(std::abs(ta[b](i,j,k) - ref_tag));
        });

        amrex::Print() << "max |phi| = " << phi.norm0() << ", tagged cells = " << tagged.sum(0) << "\n";

        Real const tol = std::is_same_v<Real,float> ? 1.e-5_rt : 1.e-12_rt;
        Real const err_ref = err.norm0(0);
        Real const err_pp = err.norm0(1);
        Real const err_tag = err.norm0(2);
        amrex::Print() << "Parser vs. makeParser: max diff = " << err_pp << "\n";
        bool ok = err_pp <= tol;
        if (function == default_function) {
            amrex::Print() << "Parser vs. C++: max diff = " << err_ref << "\n";
            ok = ok && (err_ref <= tol * Real(A));
        }
        if (tag == default_tag) {
            amrex::Print() << "IParser vs. C++: max diff = " << err_tag << "\n";
            ok = ok && (err_tag == 0.0_rt);
        }
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
