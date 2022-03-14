       subroutine wrf_totrad
       use param_wrf
       implicit none
       integer :: i,j
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(TRAD,SWDOWN,GLW,south_north,west_east) &
!$OMP PRIVATE(i,j)
       do j=1,south_north
       do i=1,west_east
       TRAD(i,j)=SWDOWN(i,j)+GLW(i,j)
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

       end
