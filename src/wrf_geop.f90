       subroutine wrf_geop
!Compute total geopotential
!input PH perturbation geopotential m2 s-2
!input PHB  base-state geopotential m2 s-2
! output geop total geopotential m2 s-2

       use param_wrf

       IMPLICIT NONE
       integer :: i,j,k

!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(PH,PHB,phi,bottom_top,south_north,west_east)     &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,south_north
       do i=1,west_east
       phi(i,j,k)=PH(i,j,k)+PHB(i,j,k)
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

       end subroutine
