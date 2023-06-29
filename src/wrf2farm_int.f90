      subroutine wrf2farm_int

      use general
      use param_wrf
      use param_minni

      IMPLICIT NONE

      integer    :: k,i,j
!wrf level height horizzontally interpolated on outgrid
      real         :: zwrf_agdint(nx,ny,bottom_top) 
      real,allocatable   :: dummy1(:,:,:) 
      real,allocatable   :: dummy2(:,:) 
      real,allocatable   :: dummy3(:,:,:) 
      real,allocatable   :: dummy4(:,:,:) 

#ifdef debug
  print*,'debug: wrf2farm:'
  print*,'debug: horizontal and vertical interpolation to farm grid'
#endif


!-------------OROGRAPHY---------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(hgtfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         hgtfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,HGT,                                     &
                     nx,ny,xfarm,yfarm,hgtfarm,.false.)
!-------------Sea Surface Temperature-------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(sstfarm,nx,ny)        &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         sstfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,SST,                                     &
                     nx,ny,xfarm,yfarm,sstfarm,.false.)
!-------------SNOW-------------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(snowfarm,nx,ny)        &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         snowfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,               &
     &               xwrf,ywrf,snowh,                                 &
                     nx,ny,xfarm,yfarm,snowfarm,.false.)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(snowfarm,nx,ny)       &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         if ( snowfarm(i,j) .lt. 0.001 ) then
             snowfarm(i,j)=0.0
         else
             snowfarm(i,j)=1.0
         endif
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!-------------Skin Sea Surface Temperature----------------------------------
       if ( fsstsk .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(sstskfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         sstskfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,SSTSK,                                   &
                     nx,ny,xfarm,yfarm,sstskfarm,.false.)
       endif
!-------------Rougthness Length ----------------------------------
       if ( fznt .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(z0farm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         z0farm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,         &
     &               xwrf,ywrf,znt,                            &
                     nx,ny,xfarm,yfarm,z0farm,.false.)
       endif

!-------------PRECIPITATION-----------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(precfarm,nx,ny)       &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         precfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,precflx,                                 &
                     nx,ny,xfarm,yfarm,precfarm,.false.)
         call check_prec(nx,ny,precfarm)
!-------------MONIN OBUKOV LENGTH-----------------------------------------------
      if ( flstar .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(lstarfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         lstarfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(RMOL,south_north,west_east) &
!$OMP PRIVATE(i,j)
       do j=1,south_north
       do i=1,west_east
         RMOL(i,j)=1.0/RMOL(i,j)
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,RMOL,                                    &
                     nx,ny,xfarm,yfarm,lstarfarm,.false.)
      endif
!-------------ALBEDO HEIGHT-----------------------------------------------
      if ( falbedo .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(albedofarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         albedofarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,ALBEDO,                                  &
                     nx,ny,xfarm,yfarm,albedofarm,.false.)
      endif
!-------------FRICTION VELOCITY-----------------------------------------------
      if ( fustar .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(ustarfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         ustarfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,UST,                                     &
                     nx,ny,xfarm,yfarm,ustarfarm,.false.)
      endif
!-------------PBL HEIGHT-----------------------------------------------
      if ( fpbl .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(pblfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         pblfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,PBLH,                                    &
                     nx,ny,xfarm,yfarm,pblfarm,.false.)
      endif
!------------NET HEAT FLUX----------------------------------------------
       if ( fnrad .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(netradfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         netradfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,nrad,                                    &
                     nx,ny,xfarm,yfarm,netradfarm,.false.)
       endif
!------------TOTAL HEAT FLUX----------------------------------------------
       if ( ftgrad .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(totradfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         totradfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,trad,                                    &
                     nx,ny,xfarm,yfarm,totradfarm,.false.)
       endif
!------------TOTAL SHORTWAVE HEAT FLUX----------------------------------------------
       if ( fsgrad .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(totsradfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         totsradfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,SWDOWN,                                  &
                     nx,ny,xfarm,yfarm,totsradfarm,.false.)
       endif
!------------TOTAL LONGWAVE HEAT FLUX----------------------------------------------
       if ( flgrad .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(totlradfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         totlradfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,GLW,                                     &
                     nx,ny,xfarm,yfarm,totlradfarm,.false.)
       endif
!------------SENSIBLE HEAT FLUX----------------------------------------------
       if ( fsh .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(shfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         shfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,HFX,                                     &
                     nx,ny,xfarm,yfarm,shfarm,.false.)
       endif
!------------LATENT HEAT FLUX----------------------------------------------
       if ( fsh .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(lhfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         lhfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,LH,                                     &
                     nx,ny,xfarm,yfarm,lhfarm,.false.)
       endif
!------------GROUND HEAT FLUX----------------------------------------------
       if ( fgh .eqv. .true. ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(ghfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         ghfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,GRDFLX,                                  &
                     nx,ny,xfarm,yfarm,ghfarm,.false.)
       endif
!-------------Total Cloud Cover-------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(tccfarm,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         tccfarm(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,tcc,                                     &
                     nx,ny,xfarm,yfarm,tccfarm,.false.)
!-------------LEVELS HEIGHT-------------------------------------------------
!Horizontal Interpoation from wrf to minni: height
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(zwrf_agdint,nx,ny,bottom_top ) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,ny
       do i=1,nx
         zwrf_agdint(i,j,k)=0.0
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_3dfields(west_east,south_north,bottom_top,       &
     &               xwrf,ywrf,zwrf_agd,                                &
                     nx,ny,xfarm,yfarm,zwrf_agdint,.false.)
!-------------TEMPERATURE---------------------------------------------------
!Horizontal Interpoation from wrf to minni: temperature
         if (.not. allocated(dummy1))                                   &
     &      allocate(dummy1(nx,ny,bottom_top))
         if (.not. allocated(dummy2))                                   &
     &      allocate(dummy2(nx,ny))
         if (.not. allocated(dummy3))                                   &
     &      allocate(dummy3(nx,ny,bottom_top+1))
         if (.not. allocated(dummy4))                                   &
     &      allocate(dummy4(nx,ny,bottom_top+1))
!if ( use_T2 ) then
!!$OMP PARALLEL DO            &
!!$OMP  COLLAPSE(2)           &
!!$OMP DEFAULT(NONE)          &
!!$OMP SHARED(dummy2,nx,ny) &
!!$OMP PRIVATE(i,j)
!   do j=1,ny
!     do i=1,nx
!       dummy2(i,j)=0.0
!     enddo
!   enddo
!!$OMP END PARALLEL DO
!!$OMP BARRIER
!!$OMP PARALLEL DO            &
!!$OMP  COLLAPSE(3)           &
!!$OMP DEFAULT(NONE)          &
!!$OMP SHARED(dummy3,dummy4,nx,ny,bottom_top_stag) &
!!$OMP PRIVATE(i,j,k)
!   do k=1,bottom_top_stag
!     do j=1,ny
!       do i=1,nx
!         dummy3(i,j,k)=0.0
!         dummy4(i,j,k)=0.0
!       enddo
!     enddo
!   enddo
!!$OMP END PARALLEL DO
!!$OMP BARRIER
!endif !if use T2, we have dummy2,3,4=0

!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,bottom_top,nx,ny) &
!$OMP PRIVATE(i,j,k)
do k=1,bottom_top
  do j=1,ny
    do i=1,nx
      dummy1(i,j,k)=0.0
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(tfarm,nx,ny,nz) &
!$OMP PRIVATE(i,j,k)
do k=1,nz
  do j=1,ny
    do i=1,nx
      tfarm(i,j,k)=0.0
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
! esce dummy1 che contiene la T wrf interpolata su griglia farm
call interp2d_3dfields(west_east,south_north,bottom_top,       &
                       xwrf,ywrf,insituT, nx,ny,xfarm,yfarm,dummy1,.false.)
!ho la mia dummy1 che e' la twrf interpolata su farm grid.
! se pero' voglio usare la T2, prosegui qua sotto 

!Horizontal Interpoation from wrf to minni: temperature 2 m
if ( use_T2 ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy2,nx,ny) &
!$OMP PRIVATE(i,j)
   do j=1,ny
     do i=1,nx
       dummy2(i,j)=0.0
     enddo
   enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,nx,ny,bottom_top_stag) &
!$OMP PRIVATE(i,j,k)
   do k=1,bottom_top_stag
     do j=1,ny
       do i=1,nx
         dummy3(i,j,k)=0.0
         dummy4(i,j,k)=0.0
       enddo
     enddo
   enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

! dummy2 e' la T2m di wrf su griglia farm (2d)
  call interp2d_2dfields(west_east,south_north, xwrf,ywrf,T2,           &
                       nx,ny,xfarm,yfarm,dummy2,.false.)
!Compose T2m and insituT on horizontal minni grid
!Compose wrf level height with 2m
! dummi3 diventa la t_wrf su grid farm
! dummy4 ci metto la z sul suolo
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,dummy1,zwrf_agdint,bottom_top_stag,ny,nx) &
!$OMP PRIVATE(i,j,k)
  do k=2,bottom_top_stag
    do j=1,ny
      do i=1,nx
        dummy3(i,j,k)=dummy1(i,j,k-1)
        dummy4(i,j,k)=zwrf_agdint(i,j,k-1)
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,dummy2,nx,ny) &
!$OMP PRIVATE(i,j)
! al livello piu' basso, metto in t_wrf (dummy3, la t2m)
! e la Z al livello piu' basso a 2m (dummy4)
  do j=1,ny
    do i=1,nx
      dummy3(i,j,1)=dummy2(i,j)
      dummy4(i,j,1)=2.0
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

!massimo
! goto 111
!Vertical intepolaton of insitu T taking into account T2m
  call wrf_interp_3d_z(nx,ny,bottom_top_stag,dummy3,dummy4,nz,zlev,tfarm)
else ! if use T2
!$OMP BARRIER
  call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,nz,zlev,tfarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(tfarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
  do j=1,ny
    do i=1,nx
      tfarm(i,j,1)=dummy1(i,j,1)
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
endif ! use t2

!-------------ZONAL VELOCITY---------------------------------------------------
!if ( use_W10 ) then
!!$OMP PARALLEL DO            &
!!$OMP  COLLAPSE(2)           &
!!$OMP DEFAULT(NONE)          &
!!$OMP SHARED(dummy2,nx,ny) &
!!$OMP PRIVATE(i,j)
!  do j=1,ny
!    do i=1,nx
!      dummy2(i,j)=0.0
!    enddo
!  enddo
!!$OMP END PARALLEL DO
!!$OMP BARRIER
!!$OMP PARALLEL DO            &
!!$OMP  COLLAPSE(3)           &
!!$OMP DEFAULT(NONE)          &
!!$OMP SHARED(dummy3,dummy4,nx,ny,bottom_top_stag) &
!!$OMP PRIVATE(i,j,k)
!  do k=1,bottom_top_stag
!    do j=1,ny
!      do i=1,nx
!        dummy3(i,j,k)=0.0
!        dummy4(i,j,k)=0.0
!      enddo
!    enddo
!  enddo
!!$OMP END PARALLEL DO
!!$OMP BARRIER
!endif
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,bottom_top,nx,ny) &
!$OMP PRIVATE(i,j,k)
do k=1,bottom_top
  do j=1,ny
    do i=1,nx
      dummy1(i,j,k)=0.0
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(ufarm,nx,ny,nz) &
!$OMP PRIVATE(i,j,k)
do k=1,nz
  do j=1,ny
    do i=1,nx
      ufarm(i,j,k)=0.0
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!Horizontal Interpoation from wrf to minni: zonal velocity
call interp2d_3dfields(west_east,south_north,bottom_top, xwrf,ywrf,uint,      &
                       nx,ny,xfarm,yfarm,dummy1,.false.)
!Horizontal Interpoation from wrf to minni: zonal velocity 10 m
if ( use_W10 ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy2,nx,ny) &
!$OMP PRIVATE(i,j)
  do j=1,ny
    do i=1,nx
      dummy2(i,j)=0.0
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,nx,ny,bottom_top_stag) &
!$OMP PRIVATE(i,j,k)
  do k=1,bottom_top_stag
    do j=1,ny
      do i=1,nx
        dummy3(i,j,k)=0.0
        dummy4(i,j,k)=0.0
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP BARRIER
  call interp2d_2dfields(west_east,south_north,xwrf,ywrf,U10,              &
                         nx,ny,xfarm,yfarm,dummy2,.false.)
!Compose U10 and uint on horizontal minni grid
!Compose wrf level height with 10m
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,dummy1,zwrf_agdint,bottom_top_stag,nx,ny) &
!$OMP PRIVATE(i,j,k)
  do k=2,bottom_top_stag
    do j=1,ny
      do i=1,nx
        dummy3(i,j,k)=dummy1(i,j,k-1)
        dummy4(i,j,k)=zwrf_agdint(i,j,k-1)
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,dummy2,nx,ny) &
!$OMP PRIVATE(i,j)
  do j=1,ny
    do i=1,nx
      dummy3(i,j,1)=dummy2(i,j)
      dummy4(i,j,1)=10.0
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!Vertical intepolaton of zonal velocity taking into account U10m
  call wrf_interp_3d_z(nx,ny,bottom_top_stag,dummy3,dummy4,nz,zlev,ufarm)
else !use W10
!$OMP BARRIER
  call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,nz,zlev,ufarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(ufarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
  do j=1,ny
    do i=1,nx
      ufarm(i,j,1)=dummy1(i,j,1)
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
endif ! use w10

!-----------MERIDIONAL VELOCITY-----------------------------------------------
!      if ( use_W10 ) then
!!$OMP PARALLEL DO            &
!!$OMP  COLLAPSE(2)           &
!!$OMP DEFAULT(NONE)          &
!!$OMP SHARED(dummy2,nx,ny) &
!!$OMP PRIVATE(i,j)
!       do j=1,ny
!       do i=1,nx
!         dummy2(i,j)=0.0
!       enddo
!       enddo
!!$OMP END PARALLEL DO
!!$OMP BARRIER
!!$OMP PARALLEL DO            &
!!$OMP  COLLAPSE(3)           &
!!$OMP DEFAULT(NONE)          &
!!$OMP SHARED(dummy3,dummy4,nx,ny,bottom_top_stag) &
!!$OMP PRIVATE(i,j,k)
!       do k=1,bottom_top_stag
!       do j=1,ny
!       do i=1,nx
!         dummy3(i,j,k)=0.0
!         dummy4(i,j,k)=0.0
!       enddo
!       enddo
!       enddo
!!$OMP END PARALLEL DO
!!$OMP BARRIER
!       endif
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
do k=1,bottom_top
  do j=1,ny
    do i=1,nx
      dummy1(i,j,k)=0.0
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(vfarm,nx,ny,nz) &
!$OMP PRIVATE(i,j,k)
!massimo debug
do k=1,nz
  do j=1,ny
    do i=1,nx
      vfarm(i,j,k)=0.0
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!Horizontal Interpoation from wrf to minni: meridional velocity
call interp2d_3dfields(west_east,south_north,bottom_top,xwrf,ywrf,vint,       &
                       nx,ny,xfarm,yfarm,dummy1,.false.)

if ( use_W10 ) then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy2,nx,ny) &
!$OMP PRIVATE(i,j)
  do j=1,ny
    do i=1,nx
      dummy2(i,j)=0.0
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,nx,ny,bottom_top_stag) &
!$OMP PRIVATE(i,j,k)
  do k=1,bottom_top_stag
    do j=1,ny
      do i=1,nx
        dummy3(i,j,k)=0.0
        dummy4(i,j,k)=0.0
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP BARRIER
!Horizontal Interpoation from wrf to minni: meridional velocity 10 m
  call interp2d_2dfields(west_east,south_north,xwrf,ywrf,V10,                 &
                        nx,ny,xfarm,yfarm,dummy2,.false.)
!Compose V10 and vint on horizontal minni grid
!Compose wrf level height with 10m
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,dummy1,zwrf_agdint,bottom_top_stag,nx,ny) &
!$OMP PRIVATE(i,j,k)
  do k=2,bottom_top_stag
    do j=1,ny
      do i=1,nx
        dummy3(i,j,k)=dummy1(i,j,k-1)
        dummy4(i,j,k)=zwrf_agdint(i,j,k-1)
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,dummy2,nx,ny) &
!$OMP PRIVATE(i,j)
  do j=1,ny
    do i=1,nx
      dummy3(i,j,1)=dummy2(i,j)
      dummy4(i,j,1)=10.0
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!Vertical intepolaton of meridional velocity taking into account  V10m
  call wrf_interp_3d_z(nx,ny,bottom_top_stag,dummy3,dummy4,nz,zlev,vfarm)
else
!$OMP BARRIER
  call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,nz,zlev,vfarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(vfarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
  do j=1,ny
    do i=1,nx
      vfarm(i,j,1)=dummy1(i,j,1)
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
endif ! use w10


!-------------RELATIVE HUMIDITY---------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
do k=1,bottom_top
  do j=1,ny
    do i=1,nx
      dummy1(i,j,k)=0.0
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(rhfarm,nx,ny,nz) &
!$OMP PRIVATE(i,j,k)
do k=1,nz
  do j=1,ny
    do i=1,nx
      rhfarm(i,j,k)=0.0
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!Horizontal Interpoation from wrf to minni: relative humidity
call interp2d_3dfields(west_east,south_north,bottom_top,xwrf,ywrf,rh,      &
                       nx,ny,xfarm,yfarm,dummy1,.false.)
call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,nz,zlev,rhfarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(rhfarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
do j=1,ny
  do i=1,nx
    rhfarm(i,j,1)=dummy1(i,j,1)
  enddo
enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

!---------------------------------------------------------------------------
!-------------HYDROMETEOR---------------------------------------------------
!---------------------------------------------------------------------------
if ( fhymet .eqv. .true. ) then
!-------------QCLOUD--------------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
  do k=1,bottom_top
    do j=1,ny
      do i=1,nx
        dummy1(i,j,k)=0.0
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(qcloudfarm,nx,ny,nz) &
!$OMP PRIVATE(i,j,k)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        qcloudfarm(i,j,k)=0.0
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!Horizontal Interpoation from wrf to minni: Cloud water mixing ratio
  call interp2d_3dfields(west_east,south_north,bottom_top,xwrf,ywrf,qcloud,   &
                         nx,ny,xfarm,yfarm,dummy1,.false.)
  call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,nz,zlev,qcloudfarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(qcloudfarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
  do j=1,ny
    do i=1,nx
      qcloudfarm(i,j,1)=dummy1(i,j,1)
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!-------------QICE--------------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
  do k=1,bottom_top
    do j=1,ny
      do i=1,nx
        dummy1(i,j,k)=0.0
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(qicefarm,nx,ny,nz) &
!$OMP PRIVATE(i,j,k)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        qicefarm(i,j,k)=0.0
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!Horizontal Interpoation from wrf to minni: Ice mixing ratio
  call interp2d_3dfields(west_east,south_north,bottom_top,xwrf,ywrf,qice,      &
                        nx,ny,xfarm,yfarm,dummy1,.false.)
  call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,nz,zlev,qicefarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(qicefarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
  do j=1,ny
    do i=1,nx
      qicefarm(i,j,1)=dummy1(i,j,1)
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!-------------QRAIN--------------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
  do k=1,bottom_top
    do j=1,ny
      do i=1,nx
        dummy1(i,j,k)=0.0
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(qrainfarm,nx,ny,nz) &
!$OMP PRIVATE(i,j,k)
  do k=1,nz
    do j=1,ny
      do i=1,nx
        qrainfarm(i,j,k)=0.0
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!Horizontal Interpoation from wrf to minni: Rain water mixing ratio
  call interp2d_3dfields(west_east,south_north,bottom_top,xwrf,ywrf,qrain,    &
     &               nx,ny,xfarm,yfarm,dummy1,.false.)
  call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,nz,zlev,qrainfarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(qrainfarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
  do j=1,ny
  do i=1,nx
    qrainfarm(i,j,1)=dummy1(i,j,1)
  enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!-------------QSNOW--------------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
  do k=1,bottom_top
  do j=1,ny
  do i=1,nx
    dummy1(i,j,k)=0.0
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(qsnowfarm,nx,ny,nz) &
!$OMP PRIVATE(i,j,k)
  do k=1,nz
  do j=1,ny
  do i=1,nx
    qsnowfarm(i,j,k)=0.0
  enddo
  enddo
  enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!Horizontal Interpoation from wrf to minni: Snow water mixing ratio
         call interp2d_3dfields(west_east,south_north,bottom_top,       &
     &               xwrf,ywrf,qsnow,                                   &
     &               nx,ny,xfarm,yfarm,dummy1,.false.)
         call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,      &
     &                        nz,zlev,qsnowfarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(qsnowfarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         qsnowfarm(i,j,1)=dummy1(i,j,1)
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!-------------QGRAUP--------------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,ny
       do i=1,nx
         dummy1(i,j,k)=0.0
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(qgraupfarm,nx,ny,nz) &
!$OMP PRIVATE(i,j,k)
       do k=1,nz
       do j=1,ny
       do i=1,nx
         qgraupfarm(i,j,k)=0.0
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!Horizontal Interpoation from wrf to minni: graupel mixing ratio
         call interp2d_3dfields(west_east,south_north,bottom_top,       &
     &               xwrf,ywrf,qgraup,                                   &
     &               nx,ny,xfarm,yfarm,dummy1,.false.)
         call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,      &
     &                        nz,zlev,qgraupfarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(qgraupfarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         qgraupfarm(i,j,1)=dummy1(i,j,1)
       enddo
       enddo
!$OMP END PARALLEL DO
!---------------------------------------------------------------------------------
!-----END---HYDROMETEOR-----------------------------------------------------------
!---------------------------------------------------------------------------------
      endif
!-------------VERTICAL VELOCITY---------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,ny
       do i=1,nx
         dummy1(i,j,k)=0.0
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,nx,ny,bottom_top_stag) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top_stag
       do j=1,ny
       do i=1,nx
         dummy3(i,j,k)=0.0
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_3dfields(west_east,south_north,bottom_top_stag,  &
     &               xwrf,ywrf,W,                                       &
                     nx,ny,xfarm,yfarm,dummy3,.false.)
!Vertical interpolation on centered grid
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,dummy3,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
         do k=1,bottom_top
         do j=1,ny
         do i=1,nx
            dummy1(i,j,k)=0.5*(dummy3(i,j,k)+dummy3(i,j,k+1))
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,      &
     &                        nz,zlev,wfarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(wfarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         wfarm(i,j,1)=dummy1(i,j,1)
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!------------------GEOPOTENTIAL---------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,ny
       do i=1,nx
         dummy1(i,j,k)=0.0
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,nx,ny,bottom_top_stag) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top_stag
       do j=1,ny
       do i=1,nx
         dummy3(i,j,k)=0.0
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_3dfields(west_east,south_north,bottom_top_stag,  &
     &               xwrf,ywrf,phi,                                     &
                     nx,ny,xfarm,yfarm,dummy3,.false.)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,dummy3,ny,nx,bottom_top) &
!$OMP PRIVATE(i,j,k)
         do k=1,bottom_top
         do j=1,ny
         do i=1,nx
            dummy1(i,j,k)=0.5*(dummy3(i,j,k)+dummy3(i,j,k+1))
         enddo
         enddo
         enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

         call wrf_interp_3d_z(nx,ny,bottom_top,dummy1,zwrf_agdint,      &
     &                        nz,zlev,phifarm)
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(phifarm,dummy1,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         phifarm(i,j,1)=dummy1(i,j,1)
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!----------------------PRESSURE---------------------------------------------------
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy2,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         dummy2(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,nx,ny,bottom_top_stag) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top_stag
       do j=1,ny
       do i=1,nx
         dummy3(i,j,k)=0.0
         dummy4(i,j,k)=0.0
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,bottom_top,nx,ny) &
!$OMP PRIVATE(i,j,k)
       do k=1,bottom_top
       do j=1,ny
       do i=1,nx
         dummy1(i,j,k)=0.0
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(PSFC,south_north,west_east) &
!$OMP PRIVATE(i,j,k)
       do j=1,south_north
       do i=1,west_east
         PSFC(i,j)=PSFC(i,j)/100.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call interp2d_3dfields(west_east,south_north,bottom_top_stag,  &
     &               xwrf,ywrf,press,                                   &
                     nx,ny,xfarm,yfarm,dummy3,.false.)
!Horizontal Interpoation from wrf to minni: surface pressure 
         call interp2d_2dfields(west_east,south_north,                  &
     &               xwrf,ywrf,PSFC,                                    &
                     nx,ny,xfarm,yfarm,dummy2,.false.)
!Vertical interpolation on centered grid
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy1,dummy3,nx,ny,bottom_top) &
!$OMP PRIVATE(i,j,k)
         do k=1,bottom_top
         do j=1,ny
         do i=1,nx
            dummy1(i,j,k)=0.5*(dummy3(i,j,k)+dummy3(i,j,k+1))
         enddo
         enddo
         enddo
!Compose surface pressure and pressure on horizontal minni grid
!Compose wrf level height with 0m
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,dummy1,zwrf_agdint,nx,ny,bottom_top_stag) &
!$OMP PRIVATE(i,j,k)
       do k=2,bottom_top_stag
       do j=1,ny
       do i=1,nx
         dummy3(i,j,k)=dummy1(i,j,k-1)
         dummy4(i,j,k)=zwrf_agdint(i,j,k-1)
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(dummy3,dummy4,dummy2,nx,ny) &
!$OMP PRIVATE(i,j)
       do j=1,ny
       do i=1,nx
         dummy3(i,j,1)=dummy2(i,j)
         dummy4(i,j,1)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
         call wrf_interp_3d_z(nx,ny,bottom_top_stag,dummy3,dummy4,      &
     &                        nz,zlev,pfarm)
!-------------------------------------------------------------------------------

         if (allocated(dummy1)) deallocate(dummy1)
         if (allocated(dummy2)) deallocate(dummy2)
         if (allocated(dummy3)) deallocate(dummy3)
         if (allocated(dummy4)) deallocate(dummy4)

      end subroutine wrf2farm_int
