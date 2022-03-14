       subroutine check_prec(nx,ny,precfarm)
       implicit none
       integer :: i,j,nx,ny
       real    :: precfarm(nx,ny)
!$OMP SINGLE
       do j=1,ny
       do i=1,nx
          if (precfarm(i,j) .lt. 0.0 ) precfarm(i,j)=0
       enddo
       enddo
!$OMP END SINGLE


       end
