      subroutine mk_outgrid
      use general
      use param_minni
      implicit none

      integer :: i


! builds a regular grid
!$OMP SINGLE
      xfarm(1,:)=xstart
      yfarm(:,1)=ystart
      do i=2,nx
         xfarm(i,:)=xfarm(i-1,1)+dx
      enddo
      do i=2,ny 
         yfarm(:,i)=yfarm(1,i-1)+dy
      enddo
!$OMP END SINGLE


      end
