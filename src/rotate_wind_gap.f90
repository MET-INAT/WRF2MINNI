subroutine ROTATE_WIND_GAP
!------------------------------------------------------------------------------ 
! G. Calori, 3.1.2009
! Arianet 2009
! 
! LAST REV.:
! 01.08.2019   MA adapted for wrf2minni
! 10.11.2009   SF  modified to cope with gtpz conversion errors far from utm zone
! 03.05.2010   SF  removed bug in rotation (wrong initial x0 definition in grid loop)
!
! PURPOSE: Rotate wind components from geographic to target grid axes
!
!
! CALLS: GTPZ0
!-----------------------------------------------------------------------------
use general
use param_minni



implicit none

! coordinate conversion geographic parameters

integer :: out_sph = 12                  ! WGS84
integer :: ll_sys = 0                    ! proj code (Geographic)
integer :: out_sys = 1                   ! proj code (1=UTM; 4=Lambert conformal conic)
integer :: out_tpar=0.d0
integer ::  ll_zone = 99                  ! zone number
integer ::  out_zone = 32                 ! zone number
double precision :: ll_tpar(15)=0.d0     ! projection parameters
integer :: ll_unit = 4                   ! units (degrees)
integer :: out_unit = 2                  ! units (meters)
double precision :: crdin(2), crdio(2)
integer :: ipr, jpr, lemsg, lparm, ln27, ln83, length, iflg
character*128 :: fn27, fn83

real :: x0, y0, x1, y1
real :: angle, u, v, ue, ve

! locals

integer :: i, j, k

ipr = 1    ! do not print err msg
jpr = 1    ! do not print proj params


! look for wind components


do i = 1,nx 
  do j = 1,ny
    x0 = xstart *1000 + dx * real(i-1) * 1000
    y0 = ystart *1000 + dy * real(j-1) * 1000

    crdin(1) = dble(x0)
    crdin(2) = dble(y0)

!   compute grid point geographic coordinates
    call GTPZ0 &
         (crdin, out_sys, out_zone, out_tpar, out_unit, out_sph, ipr, jpr, &
          lemsg, lparm, crdio, ll_sys, ll_zone, ll_tpar, ll_unit, ln27, ln83, &
          fn27, fn83, length, iflg)

! re-compute x0, y0 from lat-lon coordinales

    crdin(1) = crdio(1)
    crdin(2) = crdio(2)


    call GTPZ0 &
         (crdin, ll_sys, ll_zone, ll_tpar, ll_unit, out_sph, ipr, jpr, &
          lemsg, lparm, crdio, out_sys, out_zone, out_tpar, out_unit, ln27, ln83, &
          fn27, fn83, length, iflg)

    x0 = sngl(crdio(1))
    y0 = sngl(crdio(2))
          
!   compute longitude incremented grid point output projection coordinates

    crdin(1) = crdin(1)+.1d0
    crdin(2) = crdin(2)

    call GTPZ0 &
         (crdin, ll_sys, ll_zone, ll_tpar, ll_unit, out_sph, ipr, jpr, &
          lemsg, lparm, crdio, out_sys, out_zone, out_tpar, out_unit, ln27, ln83, &
          fn27, fn83, length, iflg)

    x1 = sngl(crdio(1))
    y1 = sngl(crdio(2))

    angle = atan2(y1-y0,x1-x0)

    do k = 1,nz

      ue=ufarm(i,j,k)
      ve=vfarm(i,j,k)

      u = ue * cos(angle) - ve * sin(angle)
      v = ue * sin(angle) + ve * cos(angle)

      ufarm(i,j,k)=u
      vfarm(i,j,k)=v


    end do

  end do
end do

return

end subroutine ROTATE_WIND_GAP
