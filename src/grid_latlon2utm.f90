      subroutine grid_latlon2utm

      use general
      use param_wrf
      use param_minni

      implicit none
      integer :: i,j
      real    :: yeast,xnorth
      
      
      if ( utmcoordinate .eqv. .true. ) then
!$OMP PARALLEL         &
!$OMP DO COLLAPSE(2)   &
!$OMP DEFAULT (NONE)   &
!$OMP FIRSTPRIVATE(xnorth,yeast) &
!$OMP PRIVATE(i,j)    &
!$OMP SHARED(XWRF,YWRF,XLONG,XLAT,south_north,west_east) 
      do j=1,south_north
      do i=1,west_east
         call utm2geo(1,9.0,0.0,500000.0,XLAT(i,j),XLONG(i,j),xnorth,yeast,0)
         xwrf(i,j)=yeast/1000.0    ![km]
         ywrf(i,j)=xnorth/1000.0   ![km]
      enddo
      enddo
!$OMP END PARALLEL DO
!OMP BARRIER
      else
!$OMP PARALLEL         &
!$OMP DO COLLAPSE(2)   &
!$OMP DEFAULT (NONE)   &
!$OMP FIRSTPRIVATE(xnorth,yeast) &
!$OMP PRIVATE(i,j)    &
!$OMP SHARED(XWRF,YWRF,XLONG,XLAT,south_north,west_east) 
      do j=1,south_north
      do i=1,west_east
         xwrf(i,j)=XLONG(i,j)
         ywrf(i,j)=XLAT(i,j)
      enddo
      enddo
!$OMP END PARALLEL DO
!OMP BARRIER
      endif
      end
