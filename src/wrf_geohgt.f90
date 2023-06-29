       subroutine wrf_geohgt
       use param_wrf
       implicit none
       integer :: i,j,k

#ifdef debug
  print*,'debug: call compute Z '
#endif


!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(geohgt_asl,phi,bottom_top,south_north,west_east) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,south_north
       do i=1,west_east
       geohgt_asl(i,j,k)=phi(i,j,k)/9.81 ! [m]
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER



       end
