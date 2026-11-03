! Heat equation kernels. Each routine loops over the tiles of a MultiFab
! with an amrex_mfiter and passes the data of one tile to a plain loop.

module heat_module

  use amrex_base_module

  implicit none

  private

  public :: init_phi, advance, max_error, pi

  real(amrex_real), parameter :: pi = 3.14159265358979323846264338327950288_amrex_real

contains

  ! phi = 1 + prod_d sin(2 pi x_d)
  subroutine init_phi (phi, geom)
    type(amrex_multifab), intent(inout) :: phi
    type(amrex_geometry), intent(in) :: geom

    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: p(:,:,:,:)
    real(amrex_real) :: problo(3)

    ! amrex_problo is threadprivate, so read it outside the parallel region.
    problo = amrex_problo

    !$omp parallel private(mfi, bx, p)
    call amrex_mfiter_build(mfi, phi, tiling=.true.)
    do while (mfi%next())
       bx = mfi%tilebox()
       ! p(i,j,k,n) covers the whole box including ghost cells.
       p => phi%dataptr(mfi)
       call init_phi_tile(bx%lo, bx%hi, p, lbound(p), ubound(p), problo, geom%dx)
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel
  end subroutine init_phi

  ! phi_new = phi_old + dt * D * Laplacian(phi_old)
  subroutine advance (phi_old, phi_new, geom, dt, diffusion_coef)
    type(amrex_multifab), intent(inout) :: phi_old, phi_new
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: dt, diffusion_coef

    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: po(:,:,:,:), pn(:,:,:,:)
    real(amrex_real) :: dxinv(3)

    ! Fill ghost cells, here all from periodic images and neighbors.
    call phi_old%fill_boundary(geom)

    dxinv = 1.0_amrex_real / geom%dx

    !$omp parallel private(mfi, bx, po, pn)
    call amrex_mfiter_build(mfi, phi_new, tiling=.true.)
    do while (mfi%next())
       bx = mfi%tilebox()
       po => phi_old%dataptr(mfi)
       pn => phi_new%dataptr(mfi)
       ! phi_old and phi_new have the same boxes and ghost cells.
       call advance_tile(bx%lo, bx%hi, po, pn, lbound(po), ubound(po), &
                         dxinv, dt*diffusion_coef)
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel
  end subroutine advance

  ! Max over all cells of |phi - (1 + decay * prod_d sin(2 pi x_d))|
  function max_error (phi, geom, decay) result(err)
    type(amrex_multifab), intent(in) :: phi
    type(amrex_geometry), intent(in) :: geom
    real(amrex_real), intent(in) :: decay
    real(amrex_real) :: err

    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: p(:,:,:,:)
    real(amrex_real) :: problo(3)

    problo = amrex_problo
    err = 0.0_amrex_real

    !$omp parallel private(mfi, bx, p) reduction(max:err)
    call amrex_mfiter_build(mfi, phi, tiling=.true.)
    do while (mfi%next())
       bx = mfi%tilebox()
       p => phi%dataptr(mfi)
       err = max(err, error_tile(bx%lo, bx%hi, p, lbound(p), ubound(p), &
                                 problo, geom%dx, decay))
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

    ! Max over MPI ranks
    call amrex_parallel_reduce_max(err)
  end function max_error

  ! The tile routines take the tile box (lo, hi) and the bounds of the data
  ! (dlo, dhi). Indices beyond amrex_spacedim are 1 in both.

  pure function sine_mode (i, j, k, problo, dx) result(r)
    integer, intent(in) :: i, j, k
    real(amrex_real), intent(in) :: problo(3), dx(3)
    real(amrex_real) :: r
    integer :: iv(3), d
    iv = [i, j, k]
    r = 1.0_amrex_real
    do d = 1, amrex_spacedim
       r = r * sin(2.0_amrex_real * pi * (problo(d) + (iv(d) + 0.5_amrex_real) * dx(d)))
    end do
  end function sine_mode

  subroutine init_phi_tile (lo, hi, phi, dlo, dhi, problo, dx)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    real(amrex_real), intent(inout) :: phi(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3))
    real(amrex_real), intent(in) :: problo(3), dx(3)
    integer :: i, j, k
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             phi(i,j,k) = 1.0_amrex_real + sine_mode(i, j, k, problo, dx)
          end do
       end do
    end do
  end subroutine init_phi_tile

  subroutine advance_tile (lo, hi, pold, pnew, dlo, dhi, dxinv, fac)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    real(amrex_real), intent(in   ) :: pold(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3))
    real(amrex_real), intent(inout) :: pnew(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3))
    real(amrex_real), intent(in) :: dxinv(3), fac
    integer :: i, j, k
    real(amrex_real) :: lap
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             lap = (pold(i+1,j,k) - 2.0_amrex_real*pold(i,j,k) + pold(i-1,j,k)) * dxinv(1)**2
             if (amrex_spacedim >= 2) then
                lap = lap + (pold(i,j+1,k) - 2.0_amrex_real*pold(i,j,k) + pold(i,j-1,k)) * dxinv(2)**2
             end if
             if (amrex_spacedim == 3) then
                lap = lap + (pold(i,j,k+1) - 2.0_amrex_real*pold(i,j,k) + pold(i,j,k-1)) * dxinv(3)**2
             end if
             pnew(i,j,k) = pold(i,j,k) + fac * lap
          end do
       end do
    end do
  end subroutine advance_tile

  pure function error_tile (lo, hi, phi, dlo, dhi, problo, dx, decay) result(err)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    real(amrex_real), intent(in) :: phi(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3))
    real(amrex_real), intent(in) :: problo(3), dx(3), decay
    real(amrex_real) :: err
    integer :: i, j, k
    err = 0.0_amrex_real
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             err = max(err, abs(phi(i,j,k) - (1.0_amrex_real + decay * sine_mode(i, j, k, problo, dx))))
          end do
       end do
    end do
  end function error_tile

end module heat_module
