       subroutine wrf_hgt_abovesl
    
       use param_wrf
      
       implicit none
       integer :: k,i,j
        
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(geohgt_asl,zwrf_asl,bottom_top,south_north,west_east) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,south_north
       do i=1,west_east
          zwrf_asl(i,j,k)=0.5*( geohgt_asl(i,j,k) + geohgt_asl(i,j,k+1))
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

       end subroutine
