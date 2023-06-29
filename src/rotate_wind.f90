      subroutine rotate_wind
     
      use param_wrf
      implicit none
      real :: tmp
      integer :: k,i,j

#ifdef debug
    print*,'debug: call rotate wind'
#endif


!$OMP PARALLEL DO                         &
!$OMP  COLLAPSE(3)                        &
!$OMP DEFAULT(NONE)                       &
!$OMP SHARED(uint,cosalpha,vint,sinalpha,bottom_top,south_north,west_east) &                              
!$OMP FIRSTPRIVATE(tmp)                   &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,south_north
       do i=1,west_east
         tmp=uint(i,j,k)*cosalpha(i,j)-vint(i,j,k)*sinalpha(i,j)
         vint(i,j,k)=vint(i,j,k)*cosalpha(i,j)+uint(i,j,k)*sinalpha(i,j)
         uint(i,j,k)=tmp
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
!$OMP BARRIER


      end
