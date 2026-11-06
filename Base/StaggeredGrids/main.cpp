#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <cmath>
#include <limits>

using amrex::Real;
using namespace amrex::literals;

// A linear function, reproduced exactly by averaging.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real linear (amrex::GpuArray<Real,AMREX_SPACEDIM> const& x)
{
    return 1.0_rt + AMREX_D_TERM(x[0], + 2.0_rt*x[1], + 3.0_rt*x[2]);
}

// Velocity component d at x: a shear flow, u_d depends on the next
// coordinate only, so its discrete divergence is exactly zero.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real shear (int d, amrex::GpuArray<Real,AMREX_SPACEDIM> const& x)
{
    int const d1 = (d + 1) % AMREX_SPACEDIM;
    return Real(d + 1) * std::sin(2.0_rt * amrex::Math::pi<Real>() * x[d1]);
}

// Position of index (i,j,k) for data of the given index type: a node at i,
// a cell center at i + 1/2.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<Real,AMREX_SPACEDIM>
position (int i, int j, int k, amrex::IndexType itype,
          amrex::GpuArray<Real,AMREX_SPACEDIM> const& plo,
          amrex::GpuArray<Real,AMREX_SPACEDIM> const& dx)
{
    int const idx[3] = {i, j, k};
    amrex::GpuArray<Real,AMREX_SPACEDIM> x;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        Real const offset = itype.nodeCentered(d) ? 0.0_rt : 0.5_rt;
        x[d] = plo[d] + (Real(idx[d]) + offset) * dx[d];
    }
    return x;
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 32;
        int max_grid_size = 16;
        bool write_plotfile = true;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("write_plotfile", write_plotfile);
        }

        amrex::Box const domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain); // cell-centered
        ba.maxSize(max_grid_size);
        amrex::DistributionMapping dm(ba);
        amrex::Geometry geom(domain,
                             amrex::RealBox({AMREX_D_DECL(0.0_rt, 0.0_rt, 0.0_rt)},
                                            {AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}),
                             amrex::CoordSys::cartesian, {AMREX_D_DECL(1, 1, 1)});
        auto const plo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();

        // Data on cells, on the faces normal to each direction, and on nodes.
        // The BoxArray of a MultiFab carries the index type: convert() shifts
        // it from cell-centered to face- or node-centered.
        amrex::MultiFab phi(ba, dm, 1, 1);
        amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> umac;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            umac[d].define(amrex::convert(ba, amrex::IntVect::TheDimensionVector(d)), dm, 1, 0);
        }
        amrex::MultiFab psi(amrex::convert(ba, amrex::IntVect::TheNodeVector()), dm, 1, 0);
        amrex::Print() << "Index types: phi " << phi.ixType() << ", umac[0] " << umac[0].ixType()
                       << ", psi " << psi.ixType() << "\n";

        // 1. Face velocity. The faces of a cell-centered tile are
        //    mfi.nodaltilebox(d): every face of a box belongs to exactly one
        //    tile. ParallelFor takes a box and a function per direction and
        //    launches them together.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(phi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            AMREX_D_TERM(amrex::Array4<Real> const& u0 = umac[0].array(mfi);,
                         amrex::Array4<Real> const& u1 = umac[1].array(mfi);,
                         amrex::Array4<Real> const& u2 = umac[2].array(mfi););
            amrex::ParallelFor(AMREX_D_DECL(mfi.nodaltilebox(0), mfi.nodaltilebox(1), mfi.nodaltilebox(2)),
                AMREX_D_DECL(
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    amrex::IndexType const xface(amrex::IntVect::TheDimensionVector(0));
                    u0(i,j,k) = shear(0, position(i, j, k, xface, plo, dx));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    amrex::IndexType const yface(amrex::IntVect::TheDimensionVector(1));
                    u1(i,j,k) = shear(1, position(i, j, k, yface, plo, dx));
                },
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    amrex::IndexType const zface(amrex::IntVect::TheDimensionVector(2));
                    u2(i,j,k) = shear(2, position(i, j, k, zface, plo, dx));
                }));
        }
        // A face on a box boundary is stored in both boxes. OverrideSync makes
        // the copies agree; not needed here, the value depends only on the index.
        // FillBoundaryAndSync does this and fills ghost faces in one call.
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            umac[d].OverrideSync(geom.periodicity());
        }

        // The divergence of the face velocity is cell-centered, and the
        // velocity averaged to cell centers has AMREX_SPACEDIM components.
        amrex::MultiFab divu(ba, dm, 1, 0);
        amrex::computeDivergence(divu, amrex::GetArrOfConstPtrs(umac), geom);
        amrex::MultiFab ucc(ba, dm, AMREX_SPACEDIM, 0);
        amrex::average_face_to_cellcenter(ucc, amrex::GetArrOfConstPtrs(umac), geom);

        // The exact velocity at the cell centers, minus the averaged one.
        amrex::MultiFab exact_ucc(ba, dm, AMREX_SPACEDIM, 0);
        {
            auto const& ea = exact_ucc.arrays();
            amrex::IndexType const itype = exact_ucc.ixType();
            amrex::ParallelFor(exact_ucc, amrex::IntVect(0), AMREX_SPACEDIM,
                [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
            {
                ea[b](i,j,k,n) = shear(n, position(i, j, k, itype, plo, dx));
            });
        }
        amrex::MultiFab::Subtract(exact_ucc, ucc, 0, 0, AMREX_SPACEDIM, 0);
        Real const err_div = divu.norm0();
        Real const err_ucc = exact_ucc.norm0(0, AMREX_SPACEDIM, amrex::IntVect(0));
        amrex::Print() << "Max |div u| = " << err_div
                       << ", max error of the cell-averaged velocity = " << err_ucc << "\n";

        // 2. A cell-centered scalar averaged to the faces. The linear function
        //    is not periodic, so the ghost cells are set directly.
        {
            auto const& pa = phi.arrays();
            amrex::IndexType const itype = phi.ixType();
            amrex::ParallelFor(phi, phi.nGrowVect(), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
            {
                pa[b](i,j,k) = linear(position(i, j, k, itype, plo, dx));
            });
        }
        amrex::Array<amrex::MultiFab,AMREX_SPACEDIM> phi_face;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            phi_face[d].define(umac[d].boxArray(), dm, 1, 0);
        }
        amrex::average_cellcenter_to_face(amrex::GetArrOfPtrs(phi_face), phi, geom);
        Real err_face = 0.0_rt;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            // The exact values on the faces, minus the averaged ones.
            amrex::MultiFab exact_face(phi_face[d].boxArray(), dm, 1, 0);
            auto const& ea = exact_face.arrays();
            amrex::IndexType const itype = exact_face.ixType();
            amrex::ParallelFor(exact_face, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
            {
                ea[b](i,j,k) = linear(position(i, j, k, itype, plo, dx));
            });
            amrex::MultiFab::Subtract(exact_face, phi_face[d], 0, 0, 1, 0);
            err_face = amrex::max(err_face, exact_face.norm0());
        }
        amrex::Print() << "Max error of the face-averaged scalar = " << err_face << "\n";

        // 3. Nodal data averaged to cell centers.
        {
            auto const& na = psi.arrays();
            amrex::IndexType const itype = psi.ixType();
            amrex::ParallelFor(psi, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
            {
                na[b](i,j,k) = linear(position(i, j, k, itype, plo, dx));
            });
        }
        amrex::MultiFab psi_cc(ba, dm, 1, 0);
        amrex::average_node_to_cellcenter(psi_cc, 0, psi, 0, 1);
        // phi holds the same linear function at the cell centers.
        amrex::MultiFab::Subtract(phi, psi_cc, 0, 0, 1, 0);
        Real const err_node = phi.norm0();
        amrex::Print() << "Max error of the cell-averaged nodal data = " << err_node << "\n";

        if (write_plotfile) {
            // Plotfile viewers usually expect cell-centered data: face and
            // nodal data are written after averaging.
            amrex::MultiFab out(ba, dm, 2 + AMREX_SPACEDIM, 0);
            amrex::MultiFab::Copy(out, divu, 0, 0, 1, 0);
            amrex::MultiFab::Copy(out, ucc, 0, 1, AMREX_SPACEDIM, 0);
            amrex::MultiFab::Copy(out, psi_cc, 0, 1 + AMREX_SPACEDIM, 1, 0);
            amrex::WriteSingleLevelPlotfile("plt_staggered", out,
                                            {"divu", AMREX_D_DECL("u", "v", "w"), "psi"},
                                            geom, 0.0_rt, 0);
        }

        Real const tol = 100.0_rt * std::numeric_limits<Real>::epsilon();
        if (err_div == 0.0_rt && err_ucc < tol && err_face < tol && err_node < tol) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
