#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_TagParallelFor.H>

using amrex::Real;
using namespace amrex::literals;

// The valid data: a function of the cell index and the component.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real exact (int i, int j, int k, int n)
{
    return Real(i + 2*j + 3*k) + 100.0_rt * Real(n);
}

// Fill ghost cell (i,j,k) outside a domain face with the value of its mirror
// image inside the domain (a reflect_even boundary condition).
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void reflect (amrex::Array4<Real> const& a, int i, int j, int k, int n,
              amrex::Orientation face, amrex::IntVect const& dlo,
              amrex::IntVect const& dhi)
{
    int const d = face.coordDir();
    int idx[3] = {i, j, k};
    idx[d] = face.isLow() ? 2*dlo[d] - 1 - idx[d] : 2*dhi[d] + 1 - idx[d];
    a(i,j,k,n) = a(idx[0],idx[1],idx[2],n);
}

// Does the valid box touch the domain face?
bool touches (amrex::Box const& vbx, amrex::Box const& domain, amrex::Orientation face)
{
    int const d = face.coordDir();
    return face.isLow() ? vbx.smallEnd(d) == domain.smallEnd(d)
                        : vbx.bigEnd(d)   == domain.bigEnd(d);
}

// One tag per (local box, domain face it touches): the Array4 of the box,
// the thin box of ghost cells outside the face, and the face.
amrex::Vector<amrex::Array4BoxOrientationTag<Real>>
make_boundary_tags (amrex::MultiFab& mf, amrex::Box const& domain, int nghost)
{
    amrex::Vector<amrex::Array4BoxOrientationTag<Real>> tags;
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        amrex::Box const& vbx = mfi.validbox();
        for (amrex::OrientationIter oit; oit; ++oit) {
            amrex::Orientation const face = oit();
            if (touches(vbx, domain, face)) {
                tags.push_back({mf.array(mfi), amrex::adjCell(vbx, face, nghost), face});
            }
        }
    }
    return tags;
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int status = 0;
    {
        int n_cell = 64;
        int max_grid_size = 8;
        int nghost = 2;
        int ncomp = 2;
        int nrepeat = 10;
        {
            amrex::ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("nghost", nghost);
            pp.query("ncomp", ncomp);
            pp.query("nrepeat", nrepeat);
        }

        amrex::Box const domain(amrex::IntVect(0), amrex::IntVect(n_cell-1));
        amrex::BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        amrex::DistributionMapping dm(ba);
        amrex::IntVect const dlo = domain.smallEnd();
        amrex::IntVect const dhi = domain.bigEnd();

        // Valid cells hold the exact function, ghost cells zero. One copy
        // per method.
        amrex::MultiFab phi_loop(ba, dm, ncomp, nghost);
        phi_loop.setVal(0.0_rt);
        auto const& pa = phi_loop.arrays();
        amrex::ParallelFor(phi_loop, amrex::IntVect(0), ncomp,
            [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n)
        {
            pa[b](i,j,k,n) = exact(i,j,k,n);
        });
        amrex::MultiFab phi_tags(ba, dm, ncomp, nghost);
        amrex::MultiFab phi_cached(ba, dm, ncomp, nghost);
        amrex::MultiFab::Copy(phi_tags, phi_loop, 0, 0, ncomp, nghost);
        amrex::MultiFab::Copy(phi_cached, phi_loop, 0, 0, ncomp, nghost);

        // 1. One launch per box and face: the baseline. Boxes not touching
        //    a domain face launch nothing; the others launch up to
        //    2*AMREX_SPACEDIM small kernels each. Not recommended on GPU
        //    for many small boxes, since no launch fills the device.
        for (int irep = 0; irep < nrepeat; ++irep) {
            BL_PROFILE("per_face_launches");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(phi_loop); mfi.isValid(); ++mfi) {
                amrex::Box const& vbx = mfi.validbox();
                auto const& a = phi_loop.array(mfi);
                for (amrex::OrientationIter oit; oit; ++oit) {
                    amrex::Orientation const face = oit();
                    if (touches(vbx, domain, face)) {
                        amrex::ParallelFor(amrex::adjCell(vbx, face, nghost), ncomp,
                            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
                        {
                            reflect(a, i, j, k, n, face, dlo, dhi);
                        });
                    }
                }
            }
            amrex::Gpu::streamSynchronize();
        }

        // 2. One launch for all tags. The kernel receives the tag of its
        //    cell. ParallelFor builds a TagVector (copies the tags to the
        //    device) and synchronizes on every call, so no explicit
        //    synchronization is needed here.
        auto const tags = make_boundary_tags(phi_tags, domain, nghost);
        for (int irep = 0; irep < nrepeat; ++irep) {
            BL_PROFILE("tags_rebuilt");
            amrex::ParallelFor(tags, ncomp,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                      amrex::Array4BoxOrientationTag<Real> const& tag)
            {
                reflect(tag.fab, i, j, k, n, tag.face, dlo, dhi);
            });
        }

        // 3. The TagVector is built once and reused; it only needs to be
        //    rebuilt when the boxes or the MultiFab change.
        amrex::TagVector<amrex::Array4BoxOrientationTag<Real>> tv(
            make_boundary_tags(phi_cached, domain, nghost));
        for (int irep = 0; irep < nrepeat; ++irep) {
            BL_PROFILE("tags_cached");
            amrex::ParallelFor(tv, ncomp,
                [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                      amrex::Array4BoxOrientationTag<Real> const& tag)
            {
                reflect(tag.fab, i, j, k, n, tag.face, dlo, dhi);
            });
            amrex::Gpu::streamSynchronize();
        }

        amrex::Long ntags = tags.size();
        amrex::ParallelDescriptor::ReduceLongSum(ntags);
        amrex::Print() << ba.size() << " boxes, " << ntags << " boundary faces\n";

        // Methods 2 and 3 must agree with the baseline, ghost cells included.
        amrex::MultiFab::Subtract(phi_tags, phi_loop, 0, 0, ncomp, nghost);
        amrex::MultiFab::Subtract(phi_cached, phi_loop, 0, 0, ncomp, nghost);
        Real const d2 = phi_tags.norm0(0, ncomp, amrex::IntVect(nghost));
        Real const d3 = phi_cached.norm0(0, ncomp, amrex::IntVect(nghost));
        amrex::Print() << "Max difference from the per-face launches: tags rebuilt "
                       << d2 << ", tags cached " << d3 << "\n";
        if (d2 == 0.0_rt && d3 == 0.0_rt) {
            amrex::Print() << "PASSED\n";
        } else {
            amrex::Print() << "FAILED\n";
            status = 1;
        }
    }
    amrex::Finalize();
    return status;
}
