! Heat equation d(phi)/dt = D Laplacian(phi) on a periodic unit domain,
! written in Fortran with the AMReX Fortran interfaces.

program main

  use amrex_base_module
  use heat_module, only : init_phi, advance, max_error, pi

  implicit none

  integer :: n_cell, max_grid_size, plot_int, step
  real(amrex_real) :: diffusion_coef, cfl, stop_time
  real(amrex_real) :: dt, dt_step, time, decay, err, k2, tol
  real(amrex_real) :: space_err, time_err
  logical :: passed
  type(amrex_parmparse) :: pp
  type(amrex_box) :: domain
  type(amrex_boxarray) :: ba
  type(amrex_distromap) :: dm
  type(amrex_geometry) :: geom
  type(amrex_multifab) :: phi_old, phi_new

  ! Reads the inputs file and command-line parameters
  call amrex_init()

  n_cell = 32
  max_grid_size = 16
  diffusion_coef = 1.0_amrex_real
  cfl = 0.9_amrex_real
  stop_time = 0.01_amrex_real
  plot_int = -1 ! plotfile interval in steps; < 0 for no plotfiles

  call amrex_parmparse_build(pp)
  call pp%query("n_cell", n_cell)
  call pp%query("max_grid_size", max_grid_size)
  call pp%query("diffusion_coef", diffusion_coef)
  call pp%query("cfl", cfl)
  call pp%query("stop_time", stop_time)
  call pp%query("plot_int", plot_int)
  call amrex_parmparse_destroy(pp)

  ! amrex_geometry takes the physical domain and periodicity from the
  ! geometry.* parameters, so set them before building it.
  call amrex_geometry_set_prob_domain([0.0_amrex_real, 0.0_amrex_real, 0.0_amrex_real], &
                                      [1.0_amrex_real, 1.0_amrex_real, 1.0_amrex_real])
  call amrex_geometry_set_periodic([.true., .true., .true.])

  ! Indices are zero-based cell indices, as in C++.
  domain = amrex_box([0, 0, 0], [n_cell-1, n_cell-1, n_cell-1])
  call amrex_boxarray_build(ba, domain)
  call ba%maxSize(max_grid_size)
  call amrex_distromap_build(dm, ba)
  call amrex_geometry_build(geom, domain)

  ! One component and one ghost cell for the 2*dim+1 point stencil
  call amrex_multifab_build(phi_old, ba, dm, 1, 1)
  call amrex_multifab_build(phi_new, ba, dm, 1, 1)

  call init_phi(phi_new, geom)

  ! Explicit time step limit: dt <= dx^2 / (2 * dim * D)
  dt = cfl * geom%dx(1)**2 / (2.0_amrex_real * amrex_spacedim * diffusion_coef)

  time = 0.0_amrex_real
  step = 0
  if (plot_int > 0) call write_plotfile()

  do while (time < stop_time * (1.0_amrex_real - 1.e-6_amrex_real))
     ! Swap the data of the two MultiFabs without copying.
     call amrex_multifab_swap(phi_old, phi_new)
     dt_step = min(dt, stop_time - time)

     call advance(phi_old, phi_new, geom, dt_step, diffusion_coef)

     time = time + dt_step
     step = step + 1
     if (mod(step, 10) == 0 .and. amrex_parallel_ioprocessor()) then
        print '("Step ", i0, ", time = ", es16.10)', step, time
     end if
     if (plot_int > 0) then
        if (mod(step, plot_int) == 0) call write_plotfile()
     end if
  end do

  ! Exact solution: each sine mode decays as exp(-D |k|^2 t), |k|^2 = dim (2 pi)^2.
  k2 = amrex_spacedim * 4.0_amrex_real * pi**2
  decay = exp(-diffusion_coef * k2 * time)
  err = max_error(phi_new, geom, decay) / decay ! relative to the mode amplitude
  ! Allow twice the leading truncation errors of the Laplacian (dx^2)
  ! and of forward Euler (dt), plus round-off of about one epsilon per
  ! step on phi ~ 1, relative to the decaying mode amplitude.
  space_err = k2 / amrex_spacedim * geom%dx(1)**2 / 12.0_amrex_real
  time_err = diffusion_coef * k2 * dt / 2.0_amrex_real
  tol = 2.0_amrex_real * diffusion_coef * time * k2 * (space_err + time_err) &
       + step * epsilon(1.0_amrex_real) / decay
  passed = err <= tol

  if (amrex_parallel_ioprocessor()) then
     print '("Finished ", i0, " steps at time ", es16.10, ", mode amplitude = ", es16.10, &
           &", relative max error = ", es16.10)', step, time, decay, err
     if (passed) then
        print '(a)', "PASSED"
     else
        print '(a)', "FAILED"
     end if
  end if

  ! Objects in the main program are never finalized automatically, so
  ! destroy them before amrex_finalize.
  call amrex_multifab_destroy(phi_old)
  call amrex_multifab_destroy(phi_new)
  call amrex_geometry_destroy(geom)
  call amrex_distromap_destroy(dm)
  call amrex_boxarray_destroy(ba)

  call amrex_finalize()

  if (.not. passed) error stop 1

contains

  subroutine write_plotfile ()
    character(len=16) :: name
    type(amrex_string) :: varname(1)
    ! amrex_write_plotfile takes one MultiFab and Geometry per level.
    type(amrex_multifab) :: mf(1)
    type(amrex_geometry) :: gm(1)
    write(name, '("plt", i5.5)') step
    call amrex_string_build(varname(1), "phi")
    mf(1) = phi_new ! shallow copy
    gm(1) = geom
    call amrex_write_plotfile(trim(name), 1, mf, varname, gm, time, [step], [1])
  end subroutine write_plotfile

end program main
