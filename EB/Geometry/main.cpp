#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#if (AMREX_SPACEDIM == 3)
#include <AMReX_WriteEBSurface.H>
#endif

#include <cmath>
#include <string>
#include <vector>

using namespace amrex;

namespace {

constexpr Real pi = Math::pi<Real>();

// Volume (area in 2D) and surface area (perimeter in 2D) of a ball
Real ball_volume (Real r)
{
    return (AMREX_SPACEDIM == 2) ? pi*r*r : Real(4./3.)*pi*r*r*r;
}

Real ball_surface (Real r)
{
    return (AMREX_SPACEDIM == 2) ? Real(2.)*pi*r : Real(4.)*pi*r*r;
}

// Build the EB index space for geometry `name`.  The implicit functions
// follow the convention f < 0 in the fluid, f > 0 in the body, and the
// CSG operations (union, intersection, difference) act on the bodies.
// Returns the expected fluid volume, or a negative value if unknown.
Real build_geometry (std::string const& name, Geometry const& geom)
{
    const Real domain_volume = AMREX_D_TERM(geom.ProbLength(0), *geom.ProbLength(1),
                                            *geom.ProbLength(2));
    // A single-level code without multigrid needs no coarsened EB levels.
    const int required_coarsening_level = 0;
    const int max_coarsening_level = 0;

    if (name == "sphere")
    {
        // Fluid inside a sphere (a circle in 2D)
        ParmParse pp("sphere");
        Real radius = 0.6_rt;
        pp.query("radius", radius);

        EB2::SphereIF sphere(radius, {AMREX_D_DECL(0.0_rt,0.0_rt,0.0_rt)}, true);
        EB2::Build(EB2::makeShop(sphere), geom, required_coarsening_level,
                   max_coarsening_level);
        return ball_volume(radius);
    }
    else if (name == "csg")
    {
        // Body 1: a box with a cylindrical hole, rotated and translated
        const Real a = 0.35_rt;  // half width of the box
        const Real rc = 0.15_rt; // radius of the hole
        EB2::BoxIF box({AMREX_D_DECL(-a,-a,-a)}, {AMREX_D_DECL(a,a,a)}, false);
#if (AMREX_SPACEDIM == 3)
        EB2::CylinderIF hole(rc, 2, {0.0_rt,0.0_rt,0.0_rt}, false); // along z
#else
        EB2::SphereIF hole(rc, {0.0_rt,0.0_rt}, false); // a disk in 2D
#endif
        auto block = EB2::translate(EB2::rotate(EB2::makeDifference(box, hole),
                                                pi/6, 2),
                                    {AMREX_D_DECL(-0.35_rt,-0.35_rt,0.0_rt)});

        // Body 2: a hemisphere, the intersection of a ball and a half space
        const Real rs = 0.25_rt;
        EB2::SphereIF ball(rs, {AMREX_D_DECL(0.5_rt,0.5_rt,0.0_rt)}, false);
        EB2::PlaneIF half_space({AMREX_D_DECL(0.5_rt,0.5_rt,0.0_rt)},
                                {AMREX_D_DECL(0.0_rt,1.0_rt,0.0_rt)});
        auto dome = EB2::makeIntersection(ball, half_space);

        auto bodies = EB2::makeUnion(block, dome);
        EB2::Build(EB2::makeShop(bodies), geom, required_coarsening_level,
                   max_coarsening_level);

        const Real block_volume = std::pow(2*a, AMREX_SPACEDIM)
            - pi*rc*rc*std::pow(2*a, AMREX_SPACEDIM-2);
        return domain_volume - block_volume - Real(0.5)*ball_volume(rs);
    }
    else if (name == "parser")
    {
        // Implicit function given as a runtime math expression
        ParmParse pp("parser");
        std::string function = "(x/0.8)^2 + (y/0.5)^2 + (z/0.4)^2 - 1";
        const bool user_function = pp.query("function", function);

        Parser parser(function);
        parser.registerVariables({"x","y","z"});
        EB2::ParserIF pif(parser.compile<3>());
        // The GeometryShop keeps the Parser alive while the EB is built.
        EB2::Build(EB2::makeShop(pif, parser), geom, required_coarsening_level,
                   max_coarsening_level);
        // Ellipse/ellipsoid with the fluid inside
        const Real exact = (AMREX_SPACEDIM == 2) ? pi*0.8_rt*0.5_rt
                                                 : Real(4./3.)*pi*0.8_rt*0.5_rt*0.4_rt;
        return user_function ? Real(-1.) : exact;
    }
    else if (name == "stl")
    {
        // Without a GeometryShop, EB2::Build reads the eb2.* parameters,
        // e.g., eb2.geom_type = stl and eb2.stl_file.
        EB2::Build(geom, required_coarsening_level, max_coarsening_level);

        // octahedron.stl is |x|+|y|+|z| <= 0.6 with the body inside.
        ParmParse pp("eb2");
        std::string geom_type;
        std::string stl_file;
        Real stl_scale = 1.0_rt;
        pp.query("geom_type", geom_type);
        pp.query("stl_file", stl_file);
        pp.query("stl_scale", stl_scale);
        if (geom_type == "stl" && stl_file == "octahedron.stl") {
            const Real r = 0.6_rt*stl_scale;
            // The 2D cut through z = 0 is a diamond.
            const Real body = (AMREX_SPACEDIM == 2) ? 2*r*r : Real(4./3.)*r*r*r;
            return domain_volume - body;
        }
        return Real(-1.);
    }
    else
    {
        amrex::Abort("Unknown geometry " + name);
        return Real(-1.);
    }
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        std::vector<std::string> geometries{"sphere", "csg", "parser"};
        std::string surface_geometry;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.queryarr("geometries", geometries);
            surface_geometry = geometries[0];
            pp.query("surface_geometry", surface_geometry);
        }

        // [-1,1]^D with n_cell cells in each direction
        const Box domain(IntVect(0), IntVect(n_cell-1));
        const RealBox rb({AMREX_D_DECL(-1.0_rt,-1.0_rt,-1.0_rt)},
                         {AMREX_D_DECL( 1.0_rt, 1.0_rt, 1.0_rt)});
        const Geometry geom(domain, rb, CoordSys::cartesian,
                            {AMREX_D_DECL(0,0,0)});
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        const DistributionMapping dm(ba);

        const auto dx = geom.CellSizeArray();
        const Real dv = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
#if (AMREX_SPACEDIM == 2)
        const Real face_area_x = dx[1];
#else
        const Real face_area_x = dx[1]*dx[2];
#endif
        // The cells whose low x-face lies on the plane x = 0
        const int i_mid = n_cell/2;

        // Relative tolerance for the second-order EB volumes and areas
        const Real tol = 0.02_rt * (32.0_rt/Real(n_cell)) * (32.0_rt/Real(n_cell));

        for (auto const& name : geometries)
        {
            amrex::Print() << "\nGeometry: " << name << "\n";

            const Real expected = build_geometry(name, geom);

            // The factory gives access to the EB data on the grids: cell
            // flags, volume fractions, area fractions, centroids, ...
            auto factory = makeEBFabFactory(geom, ba, dm, {1,1,1}, EBSupport::full);
            MultiFab const& vfrac = factory->getVolFrac();
            auto const& ebdata = factory->getEBDataArrays();

            // Count cell types, sum fluid volume and EB area, and measure
            // the fluid cross section at x = 0 with the x-area fractions.
            auto r = ParReduce(TypeList<ReduceOpSum,ReduceOpSum,ReduceOpSum,
                                        ReduceOpSum,ReduceOpSum,ReduceOpSum>{},
                               TypeList<Long,Long,Long,Real,Real,Real>{},
                               vfrac, IntVect(0),
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
                -> GpuTuple<Long,Long,Long,Real,Real,Real>
            {
                auto flag = ebdata.get<EBData_t::cellflag>(b,i,j,k);
                Long nreg = 0, ncut = 0, ncov = 0;
                Real eb_area = 0, xsec = 0;
                if (flag.isRegular()) {
                    nreg = 1;
                    xsec = (i == i_mid) ? Real(1.) : Real(0.);
                } else if (flag.isCovered()) {
                    ncov = 1;
                } else {
                    ncut = 1;
                    eb_area = ebdata.get<EBData_t::bndryarea>(b,i,j,k);
                    xsec = (i == i_mid) ? ebdata.get<EBData_t::apx>(b,i,j,k) : Real(0.);
                }
                return {nreg, ncut, ncov,
                        ebdata.get<EBData_t::volfrac>(b,i,j,k), eb_area, xsec};
            });
            Long ncells[3] = {amrex::get<0>(r), amrex::get<1>(r), amrex::get<2>(r)};
            Real sums[3] = {amrex::get<3>(r), amrex::get<4>(r), amrex::get<5>(r)};
            ParallelDescriptor::ReduceLongSum(ncells, 3);
            ParallelDescriptor::ReduceRealSum(sums, 3);

            const Real volume = sums[0] * dv;
            // bndryarea is in units of the cell face area (dx is uniform here).
            const Real eb_area = sums[1] * std::pow(dx[0], AMREX_SPACEDIM-1);
            const Real xsection = sums[2] * face_area_x;

            amrex::Print() << "  regular, cut, covered cells: " << ncells[0] << ", "
                           << ncells[1] << ", " << ncells[2] << "\n"
                           << "  fluid volume:        " << volume << "\n"
                           << "  EB surface area:     " << eb_area << "\n"
                           << "  fluid area at x = 0: " << xsection << "\n";

            if (expected >= 0) {
                const Real err = std::abs(volume - expected) / expected;
                amrex::Print() << "  expected volume:     " << expected
                               << ", relative error " << err << "\n";
                if (err > tol) { status = 1; }
            }
            if (name == "sphere") {
                Real radius = 0.6_rt;
                ParmParse("sphere").query("radius", radius);
                const Real area_exact = ball_surface(radius);
                // Cross section of the ball through its center: a disk in 3D
                // or a chord in 2D.
                const Real xsec_exact = (AMREX_SPACEDIM == 2) ? 2*radius : pi*radius*radius;
                const Real err_area = std::abs(eb_area - area_exact) / area_exact;
                const Real err_xsec = std::abs(xsection - xsec_exact) / xsec_exact;
                amrex::Print() << "  expected EB area:    " << area_exact
                               << ", relative error " << err_area << "\n"
                               << "  expected area at x = 0: " << xsec_exact
                               << ", relative error " << err_xsec << "\n";
                if (err_area > tol || err_xsec > tol) { status = 1; }
            }

            // A MultiFab built with the EB factory carries the EB data.
            // EB_WriteSingleLevelPlotfile adds the volume fraction as "vfrac".
            MultiFab cell_type(ba, dm, 1, 0, MFInfo(), *factory);
            auto const& ct = cell_type.arrays();
            ParallelFor(cell_type, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
            {
                auto flag = ebdata.get<EBData_t::cellflag>(b,i,j,k);
                ct[b](i,j,k) = flag.isCovered() ? Real(0.) : (flag.isRegular() ? Real(2.) : Real(1.));
            });
            EB_WriteSingleLevelPlotfile("plt_" + name, cell_type, {"cell_type"}, geom, 0.0, 0);

#if (AMREX_SPACEDIM == 3)
            // Triangulated EB surface for ParaView: eb.pvtp and eb_*.vtp
            if (name == surface_geometry) {
                WriteEBSurface(ba, dm, geom, factory.get());
            }
#endif
        }

        amrex::Print() << "\n" << (status == 0 ? "PASSED" : "FAILED") << "\n";
    }
    amrex::Finalize();
    return status;
}
