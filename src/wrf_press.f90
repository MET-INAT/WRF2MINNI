       subroutine wrf_press
!Compute total pressure
!input P perturbation pressure Pa
!input PB BASE STATE PRESSURE  Pa
!output press  total pressure  hPa

       use param_wrf

       IMPLICIT NONE
       integer :: i,j,k

!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(press,P,PB,bottom_top,south_north,west_east)     &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,south_north
       do i=1,west_east
       press(i,j,k)=(P(i,j,k)+PB(i,j,k))/100.0
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

       end subroutine
