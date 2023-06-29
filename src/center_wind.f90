        subroutine center_wind
        use param_wrf
        IMPLICIT NONE
        integer :: i,j,k
        

#ifdef debug
  print*,'debug: call center_wind'
#endif


!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(U,V,vint,uint,bottom_top,south_north,west_east) &
!$OMP PRIVATE(i,j,k)
         do k=1,bottom_top
         do j=1,south_north
         do i=1,west_east
         uint(i,j,k)=                                                   &
     &      (U(i+1,j,k)+U(i,j,k))/2.0      
         vint(i,j,k)=                                                   &
     &      (V(i,j+1,k)+V(i,j,k))/2.0  
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

         end subroutine

