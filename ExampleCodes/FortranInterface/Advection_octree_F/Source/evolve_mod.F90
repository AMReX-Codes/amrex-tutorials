#include <AMReX_SPACE.H>

module evolve_module

  use amrex_amr_module
  use amrex_octree_module

  implicit none
  private

  public :: evolve

contains

  subroutine evolve ()
    use my_amr_module, only : stepno, max_step, stop_time, dtstep, plot_int
    use amr_data_module, only : phi_old, phi_new, t_new
    use plotfile_module, only : writeplotfile
    real(amrex_real) :: cur_time
    integer :: last_plot_file_step, step, lev, substep, finest_level

    cur_time = t_new(0)
    last_plot_file_step = 0;

    do step = stepno, max_step-1
       if (cur_time .ge. stop_time) exit

       if (amrex_parallel_ioprocessor()) then
          print *, ""
          print *, "STEP", step+1, "starts ..."
       end if

       call timestep(cur_time)

       cur_time = cur_time + dtstep

       if (amrex_parallel_ioprocessor()) then
          print *, "STEP", step+1, "end. TIME =", cur_time, "DT =", dtstep
       end if

       ! sync up time to avoid roundoff errors
       finest_level = amrex_get_finest_level()
       do lev = 0, finest_level
          t_new(lev) = cur_time
       end do

       if (plot_int .gt. 0 .and. mod(step+1,plot_int) .eq. 0) then
          last_plot_file_step = step+1
          call writeplotfile()
       end if

       if (cur_time .ge. stop_time - 1.e-6_amrex_real*dtstep) exit
    end do

    if (plot_int .gt. 0 .and. stepno .gt. last_plot_file_step) then
       call writeplotfile()
    end if

  end subroutine evolve

  subroutine timestep (time)
    use my_amr_module, only : regrid_int, stepno, dtstep
    use amr_data_module, only : t_old, t_new, phi_old, phi_new
    use averagedown_module, only : averagedown
    use compute_dt_module, only : compute_dt
    real(amrex_real), intent(in) :: time

    integer, save :: last_regrid_step = 0
    integer :: lev, finest_level

    if (regrid_int .gt. 0) then
       if (stepno .gt. last_regrid_step .and. mod(stepno, regrid_int) .eq. 0) then
          call amrex_regrid(0, time)
          last_regrid_step = stepno
       end if
    end if

    finest_level = amrex_get_finest_level()
    stepno = stepno+1

    call compute_dt()

    ! We need to update t_old(lev) and t_new(lev) before advance is called because of fillpath.
    do lev = 0, finest_level
       t_old(lev) = time
       t_new(lev) = time + dtstep
       ! swap phi_new(lev) and phi_old(lev) so they are consistent with t_new(lev) and t_old(lev)
       call amrex_multifab_swap(phi_old(lev), phi_new(lev))
    end do

    call advance(time, dtstep, stepno)

    call averagedown()
  end subroutine timestep

  ! update phi_new(lev)
  subroutine advance (time, dt, step)
    use my_amr_module, only : verbose
    use amr_data_module, only : phi_new, flux_reg
    use face_velocity_module, only : get_face_velocity
    use advect_module, only : advect
    use fillpatch_module, only : fillpatch
    integer, intent(in) :: step
    real(amrex_real), intent(in) :: time, dt

    integer, parameter :: ngrow = 3
    integer :: ncomp, idim, ilev, finest_level, igrd
    logical :: nodal(3)
    type(amrex_multifab), allocatable :: phiborder(:)
    type(amrex_octree_iter) :: oti
    type(amrex_box) :: bx, tbx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin,pout,pux,puy,puz,pfx,pfy,pfz
    type(amrex_fab) :: uface(amrex_spacedim)
    type(amrex_multifab), allocatable :: fluxes(:,:)

    if (verbose .gt. 0 .and. amrex_parallel_ioprocessor()) then
       write(*,'(A, 1X, I0, A, 1X, G0)') &
            "[Step", step, "] ADVANCE with dt =", dt
    end if

    ncomp = phi_new(0)%ncomp()
    finest_level = amrex_get_finest_level()

    allocate(fluxes(amrex_spacedim,0:finest_level))
    do ilev = 0, finest_level
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(fluxes(idim,ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, &
                                    ncomp, 0, nodal)
          call fluxes(idim,ilev)%setVal(0.0_amrex_real)
       end do
    end do

    allocate(phiborder(0:finest_level))
    do ilev = 0, finest_level
       call amrex_multifab_build(phiborder(ilev), phi_new(ilev)%ba, phi_new(ilev)%dm, ncomp, ngrow)
       call fillpatch(ilev, time, phiborder(ilev))
    end do

    !$omp parallel private(idim,oti,ilev,igrd,bx,tbx,pin,pout,pux,puy,puz,pfx,pfy,pfz,uface)
    do idim = 1, amrex_spacedim
       call uface(idim)%reset_omp_private()
    end do

    call amrex_octree_iter_build(oti)

    do while(oti%next())
       ilev = oti%level()
       igrd = oti%grid_index()
       bx   = oti%box()

       pin  => phiborder(ilev)%dataptr(igrd)
       pout => phi_new(ilev)%dataptr(igrd)

       pfx  => fluxes(1,ilev)%dataptr(igrd)
       pfy  => fluxes(2,ilev)%dataptr(igrd)
#if BL_SPACEDIM == 3
       pfz  => fluxes(3,ilev)%dataptr(igrd)
#endif

       do idim = 1, amrex_spacedim
          tbx = bx
          call tbx%nodalize(idim)
          call tbx%grow(1)
          call uface(idim)%resize(tbx,1)
       end do

       pux => uface(1)%dataptr()
       puy => uface(2)%dataptr()
#if BL_SPACEDIM == 3
       puz => uface(3)%dataptr()
#endif

       call get_face_velocity(time+0.5_amrex_real*dt, &
            pux, lbound(pux), ubound(pux), &
            puy, lbound(puy), ubound(puy), &
#if BL_SPACEDIM == 3
            puz, lbound(puz), ubound(puz), &
#endif
            amrex_geom(ilev)%dx, amrex_problo)

       call advect(time, bx%lo, bx%hi, &
            pin, lbound(pin), ubound(pin), &
            pout,lbound(pout),ubound(pout), &
            pux, lbound(pux), ubound(pux), &
            puy, lbound(puy), ubound(puy), &
#if BL_SPACEDIM == 3
            puz, lbound(puz), ubound(puz), &
#endif
            pfx, lbound(pfx), ubound(pfx), &
            pfy, lbound(pfy), ubound(pfy), &
#if BL_SPACEDIM == 3
            pfz, lbound(pfz), ubound(pfz), &
#endif
            amrex_geom(ilev)%dx, dt)
    end do

    call amrex_octree_iter_destroy(oti)

    do idim = 1, amrex_spacedim
       call amrex_fab_destroy(uface(idim))
    end do
    !$omp end parallel

    ! Note that the fluxes have already been scaled by dt and area.
    do ilev = 0, finest_level-1
       call flux_reg(ilev+1)%setval(0.0_amrex_real)
       call flux_reg(ilev+1)%crseadd(fluxes(:,ilev), -1.0_amrex_real)
       call flux_reg(ilev+1)%fineadd(fluxes(:,ilev+1), 1.0_amrex_real)
       call flux_reg(ilev+1)%reflux(phi_new(ilev), 1.0_amrex_real)
    end do

    do ilev = 0, finest_level
       do idim = 1, amrex_spacedim
          call amrex_multifab_destroy(fluxes(idim,ilev))
       end do
       call amrex_multifab_destroy(phiborder(ilev))
    end do

  end subroutine advance

end module evolve_module
