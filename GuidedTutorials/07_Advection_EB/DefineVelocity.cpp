#include <AMReX_MultiFabUtil.H>

#include "face_velocity.H"

using namespace amrex;

void
define_velocity (const Real time, const Geometry& geom, Array<MultiFab,AMREX_SPACEDIM>& vel_out, const MultiFab& phi)
{
    const auto      dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();
    const Box& domain  = geom.Domain();

    const auto domlo = amrex::lbound(domain);
    const auto domhi = amrex::ubound(domain);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        for (MFIter mfi(phi,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            GpuArray<Box, AMREX_SPACEDIM> nbx;
            AMREX_D_TERM(nbx[0] = mfi.tilebox(IntVect{AMREX_D_DECL(1,0,0)});,   // x-face-based tilebox
                         nbx[1] = mfi.tilebox(IntVect{AMREX_D_DECL(0,1,0)});,   // y-face-based tilebox
                         nbx[2] = mfi.tilebox(IntVect{AMREX_D_DECL(0,0,1)}););  // z-face-based tilebox

            AMREX_D_TERM(const Box& ngbxx = amrex::grow(mfi.nodaltilebox(0),1);,
                         const Box& ngbxy = amrex::grow(mfi.nodaltilebox(1),1);,
                         const Box& ngbxz = amrex::grow(mfi.nodaltilebox(2),1););

            GpuArray<Array4<Real>, AMREX_SPACEDIM> vel{ AMREX_D_DECL( vel_out[0].array(mfi),
                                                                      vel_out[1].array(mfi),
                                                                      vel_out[2].array(mfi)) };

            const Box& psibox = Box(IntVect(AMREX_D_DECL(std::min(ngbxx.smallEnd(0)-1, ngbxy.smallEnd(0)-1),
                                                         std::min(ngbxx.smallEnd(1)-1, ngbxy.smallEnd(0)-1),
                                                         0)),
                                    IntVect(AMREX_D_DECL(std::max(ngbxx.bigEnd(0),   ngbxy.bigEnd(0)+1),
                                                         std::max(ngbxx.bigEnd(1)+1, ngbxy.bigEnd(1)),
                                                         0)));

            FArrayBox psifab(psibox, 1);
            Elixir psieli = psifab.elixir();
            Array4<Real> psi = psifab.array();
            GeometryData geomdata = geom.data();

            amrex::launch(psibox,
            [=] AMREX_GPU_DEVICE (const Box& tbx)
            {
                get_face_velocity_psi(tbx, time, psi, geomdata); 
            });

            AMREX_D_TERM(
                         amrex::ParallelFor(ngbxx,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_x(i, j, k, vel[0], psi, prob_lo, dx); 
                             if (i == domlo.x or i == domhi.x+1) vel[0](i,j,k) = 0.;
                         });,

                         amrex::ParallelFor(ngbxy,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_y(i, j, k, vel[1], psi, prob_lo, dx);
                             if (j == domlo.y or j == domhi.y+1) vel[1](i,j,k) = 0.;
                         });,

                         amrex::ParallelFor(ngbxz,
                         [=] AMREX_GPU_DEVICE (int i, int j, int k)
                         {
                             get_face_velocity_z(i, j, k, vel[2], psi, prob_lo, dx);
                             if (k == domlo.z or k == domhi.z+1) vel[2](i,j,k) = 0.;
                         });
                        );
        }
    }
}
