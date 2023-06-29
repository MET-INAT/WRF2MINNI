       subroutine wrf_t
!Compute temperature from potential temperature
!input P perturbation pressure Pa
!input PB BASE STATE PRESSURE  Pa
!input T perturbation potential temperature (theta-t0) K
!output insituT temperature K 

       use param_wrf

       IMPLICIT NONE
       integer :: i,j,k

#ifdef debug
  print*,'debug: call compute T wrf_t.f90'
#endif


!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(insituT,T,P,PB,bottom_top,south_north,west_east) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,south_north
       do i=1,west_east
       insituT(i,j,k)=(T(i,j,k)+300.0)*                                 &
     &             ((P(i,j,k)+PB(i,j,k))*0.01/1000.0)**(2.0/7.0)
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

       end subroutine
