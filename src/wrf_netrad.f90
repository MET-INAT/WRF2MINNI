!Commento di Massimo D'Isidoro
! Da scambio Mail di Sandro a Gino il 3/3/2022:
!swdown=SWDOWN
!swup = SWDOWN*ALBEDO
!lwup = 5.67E-8 * TSK^4
!lwdown = GLW
!netrad=swdown - swup + lwdown - lwup

        subroutine wrf_netrad
        use param_wrf
        implicit none
        integer :: i,j
        integer :: LWUP
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(NRAD,HFX,LH,GRDFLX,south_north,west_east) &
!$OMP PRIVATE(i,j)
       do j=1,south_north
       do i=1,west_east
!        NRAD(i,j)=HFX(i,j)+LH(i,j)+GRDFLX(i,j)
        LWUP=EMISS(i,j)*5.67E-8*TSK(i,j)**4
        NRAD(i,j)=SWDOWN(i,j)*(1.-ALBEDO(i,j))+GLW(i,j)-LWUP
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

        end
