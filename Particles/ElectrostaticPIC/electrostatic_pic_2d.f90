module electrostatic_pic_module

  use iso_c_binding
  use amrex_fort_module, only : amrex_real, amrex_particle_real

  implicit none

contains

  subroutine interpolate_cic_two_levels(particles, ns, np,      &
                                        Ex_p, Ey_p,             &
                                        Ex,   Ey,               &
                                        lo,   hi,   dx,         &
                                        cEx,  cEy,              &
                                        mask,                   &
                                        clo,  chi,  cdx,        &
                                        plo,  ng,   lev)        &
       bind(c,name='interpolate_cic_two_levels')
    integer, value,   intent(in)     :: ns, np
    real(amrex_particle_real), intent(in)     :: particles(ns,np)
    real(amrex_particle_real), intent(inout)  :: Ex_p(np), Ey_p(np)
    integer,          intent(in)     :: ng, lev
    integer,          intent(in)     :: lo(2), hi(2)
    integer,          intent(in)     :: clo(2), chi(2)
    real(amrex_real), intent(in)     :: Ex(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng)
    real(amrex_real), intent(in)     :: Ey(lo(1)-ng:hi(1)+ng, lo(2)-ng:hi(2)+ng)
    real(amrex_real), intent(in)     :: cEx(clo(1)-ng:chi(1)+ng, clo(2)-ng:chi(2)+ng)
    real(amrex_real), intent(in)     :: cEy(clo(1)-ng:chi(1)+ng, clo(2)-ng:chi(2)+ng)
    integer(c_int),   intent(in)     :: mask (lo(1):hi(1),lo(2):hi(2))
    real(amrex_real), intent(in)     :: plo(2)
    real(amrex_real), intent(in)     :: dx(2), cdx(2)

    integer i, j, n
    real(amrex_real) wx_lo, wy_lo, wx_hi, wy_hi
    real(amrex_real) lx, ly
    real(amrex_real) inv_dx(2), inv_cdx(2)
    inv_dx  = 1.0d0/dx
    inv_cdx = 1.0d0/cdx

    do n = 1, np

       lx = (particles(1, n) - plo(1))*inv_dx(1)
       ly = (particles(2, n) - plo(2))*inv_dx(2)

       i = floor(lx)
       j = floor(ly)

! use the coarse E if near the level boundary
       if (lev .eq. 1 .and. mask(i,j) .eq. 1) then

          lx = (particles(1, n) - plo(1))*inv_cdx(1)
          ly = (particles(2, n) - plo(2))*inv_cdx(2)

          i = floor(lx)
          j = floor(ly)

          wx_hi = lx - i
          wy_hi = ly - j

          wx_lo = 1.0d0 - wx_hi
          wy_lo = 1.0d0 - wy_hi

          Ex_p(n) = wx_lo*wy_lo*cEx(i,   j  ) + &
                    wx_lo*wy_hi*cEx(i,   j+1) + &
                    wx_hi*wy_lo*cEx(i+1, j  ) + &
                    wx_hi*wy_hi*cEx(i+1, j+1)

          Ey_p(n) = wx_lo*wy_lo*cEy(i,   j  ) + &
                    wx_lo*wy_hi*cEy(i,   j+1) + &
                    wx_hi*wy_lo*cEy(i+1, j  ) + &
                    wx_hi*wy_hi*cEy(i+1, j+1)

! otherwise use the fine
       else

          wx_hi = lx - i
          wy_hi = ly - j

          wx_lo = 1.0d0 - wx_hi
          wy_lo = 1.0d0 - wy_hi

          Ex_p(n) = wx_lo*wy_lo*Ex(i,   j  ) + &
                    wx_lo*wy_hi*Ex(i,   j+1) + &
                    wx_hi*wy_lo*Ex(i+1, j  ) + &
                    wx_hi*wy_hi*Ex(i+1, j+1)

          Ey_p(n) = wx_lo*wy_lo*Ey(i,   j  ) + &
                    wx_lo*wy_hi*Ey(i,   j+1) + &
                    wx_hi*wy_lo*Ey(i+1, j  ) + &
                    wx_hi*wy_hi*Ey(i+1, j+1)

       end if

    end do

  end subroutine interpolate_cic_two_levels

end module electrostatic_pic_module
