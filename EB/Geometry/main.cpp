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

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        std::vector<std::string> geometries{"sphere", "csg", "parser"};
        std::string surface_geometry;
        Real radius = 0.6_rt; // of the sphere geometry
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            // queryarr only overwrites as many entries as are given, so read
            // into an empty vector to replace the default list.
            std::vector<std::string> g;
            if (pp.queryarr("geometries", g)) { geometries = g; }
            surface_geometry = geometries[0];
            pp.query("surface_geometry", surface_geometry);
            ParmParse("sphere").query("radius", radius);
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

        const Real pi = Math::pi<Real>();
        const auto dx = geom.CellSizeArray();
        const Real dv = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

        // A single-level code without multigrid needs no coarsened EB levels.
        const int required_coarsening_level = 0;
        const int max_coarsening_level = 0;

        // Relative tolerance for the second-order EB volume and area
        const Real tol = 0.02_rt;

        for (auto const& name : geometries)
        {
            amrex::Print() << "\nGeometry: " << name << "\n";

            // Build the EB index space. The implicit functions follow the
            // convention f < 0 in the fluid, f > 0 in the body, and the CSG
            // operations (union, intersection, difference) act on the bodies.
            if (name == "sphere")
            {
                // Fluid inside a sphere (a circle in 2D)
                EB2::SphereIF sphere(radius, {AMREX_D_DECL(0.0_rt,0.0_rt,0.0_rt)}, true);
                EB2::Build(EB2::makeShop(sphere), geom, required_coarsening_level,
                           max_coarsening_level);
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
            }
            else if (name == "parser")
            {
                // Implicit function given as a runtime math expression; the
                // default is an ellipsoid with the fluid inside.
                std::string function = "(x/0.8)^2 + (y/0.5)^2 + (z/0.4)^2 - 1";
                ParmParse("parser").query("function", function);

                Parser parser(function);
                parser.registerVariables({"x","y","z"});
                EB2::ParserIF pif(parser.compile<3>());
                // The GeometryShop keeps the Parser alive while the EB is built.
                EB2::Build(EB2::makeShop(pif, parser), geom, required_coarsening_level,
                           max_coarsening_level);
            }
            else if (name == "stl")
            {
                // Without a GeometryShop, EB2::Build reads the eb2.* parameters,
                // e.g., eb2.geom_type = stl and eb2.stl_file. octahedron.stl is
                // |x|+|y|+|z| <= 0.6 with the body inside.
                EB2::Build(geom, required_coarsening_level, max_coarsening_level);
            }
            else
            {
                amrex::Abort("Unknown geometry " + name);
            }

            // The factory gives access to the EB data on the grids: cell
            // flags, volume fractions, area fractions, centroids, ...
            auto factory = makeEBFabFactory(geom, ba, dm, {1,1,1}, EBSupport::full);
            MultiFab const& vfrac = factory->getVolFrac();
            auto const& ebdata = factory->getEBDataArrays();

            // Count the regular, cut and covered cells
            auto r = ParReduce(TypeList<ReduceOpSum,ReduceOpSum,ReduceOpSum>{},
                               TypeList<Long,Long,Long>{}, vfrac, IntVect(0),
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
                -> GpuTuple<Long,Long,Long>
            {
                auto flag = ebdata.get<EBData_t::cellflag>(b,i,j,k);
                if (flag.isRegular()) {
                    return {1, 0, 0};
                } else if (flag.isSingleValued()) {
                    return {0, 1, 0};
                } else {
                    return {0, 0, 1};
                }
            });
            Long ncells[3] = {amrex::get<0>(r), amrex::get<1>(r), amrex::get<2>(r)};
            ParallelDescriptor::ReduceLongSum(ncells, 3);

            // Sum the volume fractions and the EB areas of the cut cells.
            // bndryarea is in units of the cell face area (dx is uniform here).
            const Real volume = vfrac.sum(0) * dv;
            Real eb_area = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{},
                                     vfrac, IntVect(0),
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Real>
            {
                auto flag = ebdata.get<EBData_t::cellflag>(b,i,j,k);
                return flag.isSingleValued() ? ebdata.get<EBData_t::bndryarea>(b,i,j,k)
                                             : Real(0.);
            });
            ParallelDescriptor::ReduceRealSum(eb_area);
            eb_area *= std::pow(dx[0], AMREX_SPACEDIM-1);

            amrex::Print() << "  regular, cut, covered cells: " << ncells[0] << ", "
                           << ncells[1] << ", " << ncells[2] << "\n"
                           << "  fluid volume:    " << volume << "\n"
                           << "  EB surface area: " << eb_area << "\n";

            if (name == "sphere") {
                // Volume and surface of a ball (area and perimeter of a disk in 2D)
                const Real volume_exact = (AMREX_SPACEDIM == 2) ? pi*radius*radius
                                                                : Real(4./3.)*pi*radius*radius*radius;
                const Real area_exact = (AMREX_SPACEDIM == 2) ? Real(2.)*pi*radius
                                                              : Real(4.)*pi*radius*radius;
                const Real err_volume = std::abs(volume - volume_exact) / volume_exact;
                const Real err_area = std::abs(eb_area - area_exact) / area_exact;
                amrex::Print() << "  expected volume: " << volume_exact
                               << ", relative error " << err_volume << "\n"
                               << "  expected area:   " << area_exact
                               << ", relative error " << err_area << "\n";
                if (err_volume > tol || err_area > tol) { status = 1; }
            }

            // A MultiFab built with the EB factory carries the EB data.
            // EB_WriteSingleLevelPlotfile adds the volume fraction as "vfrac".
            MultiFab cell_type(ba, dm, 1, 0, MFInfo(), *factory);
            auto const& ct = cell_type.arrays();
            ParallelFor(cell_type, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
            {
                auto flag = ebdata.get<EBData_t::cellflag>(b,i,j,k);
                if (flag.isRegular()) {
                    ct[b](i,j,k) = Real(2.);
                } else if (flag.isSingleValued()) {
                    ct[b](i,j,k) = Real(1.);
                } else {
                    ct[b](i,j,k) = Real(0.);
                }
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
