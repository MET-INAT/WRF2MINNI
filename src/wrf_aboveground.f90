       subroutine wrf_aboveground(asl,agr,nk)
   
       use param_wrf
    
       implicit none
       integer :: k,nk,i,j
       real    :: asl(west_east,south_north,nk)
       real    :: agr(west_east,south_north,nk)
        

#ifdef debug
  print*,'debug: call wrf_aboveground'
#endif

!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(hgt,asl,agr,nk,south_north,west_east)    &
!$OMP PRIVATE(i,j,k)
       do k=1,nk
       do j=1,south_north
       do i=1,west_east
          agr(i,j,k)=asl(i,j,k) - hgt(i,j)
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

       end subroutine
