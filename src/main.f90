      program main
! IMPORTANT bucket not implemented check if the value in namelist is -1 
! bucket_mm ecc....
      use general
      use param_wrf
      use param_minni
      use netcdf

      IMPLICIT NONE

      integer      :: ncid,idTime, idtmp
      integer      :: fnum,nfile,timestep,numtimestep
      character    :: cdummy
      character*256 :: filenl


!massimo: aggiunta namelist invece di parametri hardcoded
!namelist/general/dirin,dirout

!$OMP SINGLE
CALL GETARG(1,filenl)
!lettura namelits
open(99,file=filenl)
read(99,general_namelist)
read(99,output_namelist)
close(99)

PRINT general_namelist
PRINT  output_namelist

allocate(xfarm(nx,ny),yfarm(nx,ny), &
        hgtfarm(nx,ny), &
        sstfarm(nx,ny), &
        sstskfarm(nx,ny), &
        precfarm(nx,ny), &
        tccfarm(nx,ny), &
        pblfarm(nx,ny), &
        tfarm(nx,ny,nz), &
        ufarm(nx,ny,nz), &
        vfarm(nx,ny,nz), &
        wfarm(nx,ny,nz), &
        rhfarm(nx,ny,nz), &
        phifarm(nx,ny,nz), &
        pfarm(nx,ny,nz), &
        netradfarm(nx,ny), &
        totradfarm(nx,ny), &
        totlradfarm(nx,ny), &
        totsradfarm(nx,ny), &
        shfarm(nx,ny), &
        lhfarm(nx,ny), &
        ghfarm(nx,ny), &
        albedofarm(nx,ny), &
        ustarfarm(nx,ny), &
        qicefarm(nx,ny,nz),&
        qcloudfarm(nx,ny,nz),&
        qrainfarm(nx,ny,nz),&
        qsnowfarm(nx,ny,nz),&
        qgraupfarm(nx,ny,nz),&
        lstarfarm(nx,ny),    & 
        XA(nx,ny),    & 
        YA(nx,ny),    & 
        XB(nx,ny),    & 
        YB(nx,ny),    & 
        XC(nx,ny),    & 
        YC(nx,ny),    & 
        XD(nx,ny),    & 
        YD(nx,ny),    & 
        WA(nx,ny),    & 
        WB(nx,ny),    & 
        WC(nx,ny),    & 
        WD(nx,ny),    &
        SNOWFARM(nx,ny), &
        Z0FARM(nx,ny))

!mk output grid
      call mk_outgrid
!Get number of timestep in each file
      call check(nf90_open(trim(dirin)//trim(filein),nf90_nowrite, ncid))
      call check(nf90_inq_dimid(ncid, 'Time',idtime))              

      call check(nf90_inquire_dimension(ncid,idTime,cdummy,numtimestep))
      call check(nf90_inq_dimid(ncid, 'west_east',idtmp))              
      call check(nf90_inquire_dimension(ncid,idtmp,cdummy,west_east))
      call check(nf90_inq_dimid(ncid, 'south_north',idtmp))              
      call check(nf90_inquire_dimension(ncid,idtmp,cdummy,south_north))
      call check(nf90_inq_dimid(ncid, 'bottom_top',idtmp))              
      call check(nf90_inquire_dimension(ncid,idtmp,cdummy,bottom_top))
      call check(nf90_inq_dimid(ncid, 'bottom_top_stag',idtmp))              
      call check(nf90_inquire_dimension(ncid,idtmp,cdummy,bottom_top_stag))
!      call check(nf90_inq_dimid(ncid, 'soil_layers_stag',idtmp))              
!      call check(nf90_inquire_dimension(ncid,idtmp,cdummy,soil_layers_stag))
      call check(nf90_inq_dimid(ncid, 'west_east_stag',idtmp))              
      call check(nf90_inquire_dimension(ncid,idtmp,cdummy,west_east_stag))
      call check(nf90_inq_dimid(ncid, 'south_north_stag',idtmp))              
      call check(nf90_inquire_dimension(ncid,idtmp,cdummy,south_north_stag))

      call check(nf90_close(ncid))


allocate(XLAT(west_east, south_north),               &
         XLONG(west_east, south_north),              &
         XLAT_U(west_east_stag, south_north),   &
         XLONG_U(west_east_stag, south_north),   &
         XLAT_V(west_east, south_north_stag),   &
         XLONG_V(west_east, south_north_stag),   &
         U(west_east_stag,south_north,bottom_top),   &
         V(west_east,south_north_stag,bottom_top),   &
         W(west_east,south_north,bottom_top_stag),   &
         PH(west_east,south_north,bottom_top_stag),   &
         PHB(west_east,south_north,bottom_top_stag),   &
         T(west_east,south_north,bottom_top),   &
         PB(west_east,south_north,bottom_top),   &
         P(west_east,south_north,bottom_top),   &
         QVAPOR(west_east,south_north,bottom_top),   &
         RAINC(west_east,south_north),   &
         RAINNC(west_east,south_north),   &
         RAINSH(west_east,south_north),   &
         SNOWNC(west_east,south_north),   &
         HAILNC(west_east,south_north),   &
         GRAUPELNC(west_east,south_north),   &
         SST(west_east,south_north),   &
         SSTSK(west_east,south_north),   &
         HGT(west_east,south_north),   &
         COSALPHA(west_east,south_north),   &
         SINALPHA(west_east,south_north),   &
         T2(west_east,south_north),   &
         U10(west_east,south_north),   &
         V10(west_east,south_north),   &
         HFX(west_east,south_north),   &
         LH(west_east,south_north),   &
         GRDFLX(west_east,south_north),   &
         SWDOWN(west_east,south_north),   &
         GLW(west_east,south_north),   &
         PBLH(west_east,south_north),   &
         PSFC(west_east,south_north),   &
         ALBEDO(west_east,south_north),   &
         UST(west_east,south_north),   &
         RMOL(west_east,south_north), &
         insituT(west_east,south_north,bottom_top),  &
         press(west_east,south_north,bottom_top),  &
         rh(west_east,south_north,bottom_top),  &
         phi(west_east,south_north,bottom_top_stag),  &
         geohgt_asl(west_east,south_north,bottom_top_stag),  &
         geohgt_agd(west_east,south_north,bottom_top_stag),  &
         zwrf_asl(west_east,south_north,bottom_top),  &
         zwrf_agd(west_east,south_north,bottom_top),  &
         tcc(west_east,south_north),  &
         totprec(west_east,south_north),  &
         uint(west_east,south_north,bottom_top),  &
         vint(west_east,south_north,bottom_top),  &
         xwrf(west_east,south_north),  &
         ywrf(west_east,south_north),  &
         precflx(west_east,south_north),  &
         qcloud(west_east,south_north,bottom_top),  &
         qgraup(west_east,south_north,bottom_top),  &
         qice(west_east,south_north,bottom_top),  &
         qrain(west_east,south_north,bottom_top),  &
         qsnow(west_east,south_north,bottom_top),  &
         nrad(west_east,south_north),  &
         trad(west_east,south_north),  &
         snowh(west_east,south_north), &
         znt(west_east,south_north))


!$OMP END SINGLE
      do timestep=start_timestep,numtimestep
!read wrf data
!$OMP SINGLE
         call read_wrf(timestep)
!$OMP END SINGLE
!convert wrf lat lon grid to utm
         call grid_latlon2utm
!compute in situ temperature
         call wrf_t
!compute pressure
         call wrf_press
!compute geopotential
         call wrf_geop
!compute geopotential height
         call wrf_geohgt
!compute relative humidity
         call wrf_rh
!compute total cloud cover
         call wrf_tcc
!compute total precipitation
         call wrf_prec
         call wrf_precflux(timestep)
!compute level hight above sea level
         call wrf_hgt_abovesl
!compute level hight above ground
         call wrf_aboveground(zwrf_asl,zwrf_agd,bottom_top)
!compute geopotential hight above ground
         call wrf_aboveground(geohgt_asl,geohgt_agd,bottom_top_stag)
!Horizontal Interpolation from staggered to central grid 
         call center_wind
!Rotate wind         
      if (( rotate .eqv. .true. ) .and. ( utmcoordinate .eqv. .false. )) then
            call rotate_wind
      endif
      if ( fnrad .eqv. .true. ) then
         call wrf_netrad
      endif
      if ( ftgrad .eqv. .true. ) then
         call wrf_totrad
      endif
!Horizontal Interpolation from wrf to minni
         call wrf2farm_int
!OMP BARRIER         
!$OMP SINGLE
         if ( ( rotate .eqv. .true. ) .and. ( utmcoordinate .eqv. .true. )) then
            call rotate_wind_gap
         endif
         call wr_farm
!$OMP END SINGLE
      enddo !numtimestep
      end program main
!----------------------------------------------
      subroutine check(istatus)
      use netcdf
      IMPLICIT NONE
      integer,intent (in) :: istatus
      if (istatus /= nf90_noerr) then
         write(*,*) trim(adjustl(nf90_strerror(istatus)))
      endif
      end subroutine check
!----------------------------------------------
