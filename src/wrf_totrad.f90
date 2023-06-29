       subroutine wrf_totrad
       use param_wrf
       implicit none
       integer :: i,j


#ifdef debug
 print*,'debug: call totrad (actually SW)'
#endif

!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(TRAD,SWDOWN,GLW,south_north,west_east) &
!$OMP PRIVATE(i,j)
       do j=1,south_north
       do i=1,west_east
! Sembrerebbe non usata da surfpro, che invece prense  TOTRAD 
! che contiene solo short wave
!TOTRAD definita shortwave almeno nel python di sandro e dai riscontri di Mario
       TRAD(i,j)=SWDOWN(i,j)+GLW(i,j)
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

       end
