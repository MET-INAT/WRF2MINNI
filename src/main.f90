program main
use general
use utils
use param_wrf
use param_minni
use netcdf

IMPLICIT NONE

integer      :: ncid,idTime, idtmp, idXTIME
integer      :: fnum,nfile,timestep,numtimestep
integer      :: year, month, day, hour, dt_step
integer      :: idummy,status
character    :: cdummy
character(len=256) :: filenl, time_units_dummy, csim_start, cfile_start
character(len=500) :: output_filename
character(len=4)   :: current_yyyy
character(len=2)   :: current_mm, current_dd, current_hh
character(len=10)   :: current_yyyymmdd
logical :: if_init, if_close
type(timedelta) :: dt_dummy


!$OMP SINGLE
CALL GETARG(1,filenl)
!lettura namelits
open(99,file=filenl)
read(99,general_namelist)
read(99,output_namelist)
close(99)

dirin=trim(dirin)//'/'

PRINT general_namelist
PRINT output_namelist


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

#ifdef debug
print*,'FILE IN INPUT: ',trim(dirin)//trim(filein)
#endif
!Get wrf grid dimensions (x,y,z,t) and close
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

call check(nf90_inq_dimid(ncid, 'west_east_stag',idtmp))              
call check(nf90_inquire_dimension(ncid,idtmp,cdummy,west_east_stag))

call check(nf90_inq_dimid(ncid, 'south_north_stag',idtmp))              
call check(nf90_inquire_dimension(ncid,idtmp,cdummy,south_north_stag))



allocate(XTIME(numtimestep),               &
         XLAT(west_east, south_north),               &
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
         SST(west_east,south_north),   &
         SSTSK(west_east,south_north),   &
         TSK(west_east,south_north),   &
         EMISS(west_east,south_north),   &
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

#ifdef debug
print*,'debug: get XTIME '
#endif
! ricavo il valore degli  steps temporali sottraendo il secondo col primo
call check(nf90_inq_varid(ncid,'XTIME',idXTIME))
call check(nf90_get_var(ncid,idXTIME, XTIME))
dt_step=XTIME(2)-XTIME(1) !time step in time units


! mi prendo anche le unità del tempo per conoscere gli steps a cosa si
! riferiscono (ore, minuti, secondi, giorni...)
call check(nf90_get_att(ncid, idXTIME, 'units', time_units_dummy))
time_units=string_trim(time_units_dummy,delimiter=' ')

#ifdef debug
print*,'debug: deltat units ', trim(time_units)
print*,'debug: deltat step ', dt_step
#endif

!check  bucket threshold
#ifdef debug
print*,'debug: get BUCKET Threshold '
#endif
status = nf90_inquire_attribute(ncid, nf90_global, "BUCKET_MM")
if (status == 0)then
  call check(nf90_get_att(ncid, NF90_GLOBAL, 'BUCKET_MM', bucket_threshold_mm))
  if ( bucket_threshold_mm >=  0 )then 
    if_bucket=.true. !this allows read of I_RAINC and I_RAINNC
    allocate( I_RAINC(west_east,south_north), I_RAINNC(west_east,south_north))
    print*,'BUCKET_THRESHOLD present (mm): ',bucket_threshold_mm
  else
    if_bucket=.false.
    print*,'BUCKET_THRESHOLD not present in wrfout'
  endif
endif

! get file start date
call check(nf90_get_att(ncid, NF90_GLOBAL, 'START_DATE', cfile_start))
! get start simulation date
call check(nf90_get_att(ncid, NF90_GLOBAL, 'SIMULATION_START_DATE', csim_start))

! set simulation start date (può non coincidere con la data iniziale del file ed
! essere precedente)
read(csim_start(1:4),'(I4)'),year
read(csim_start(6:7),'(I2)'),month
read(csim_start(9:10),'(I2)'),day
read(csim_start(12:13),'(I2)'),hour
simulation_start_date=datetime(year, month, day, hour)

! set actual file start date
read(cfile_start(1:4),'(I4)'),year
read(cfile_start(6:7),'(I2)'),month
read(cfile_start(9:10),'(I2)'),day
read(cfile_start(12:13),'(I2)'),hour
file_start_date=datetime(year, month, day, hour)


!close netcdf
call check(nf90_close(ncid))

! dalle info su step temporali e unità mi costruisco dei tipi datetime per
! gestire i tempi
!idummy è il numero di unità temporali (minuti, secondi, ore, ... ) contenute
!nel file e mi serve per calcolare quanti giorni ci sono nel file (dayxfile) che
!poi servirà per ricostruire il nome del file precedente per calcolare nel caso
!precflux
idummy=(numtimestep-1)*dt_step
select case (trim(time_units))
  case ('days')
    deltat_increment=timedelta(days=dt_step)
    deltat_since_file_start=timedelta(days=(start_timestep-1)*dt_step)
    dt_dummy=timedelta(days=idummy)
  case ('hours')
    deltat_increment=timedelta(hours=dt_step)
    deltat_since_file_start=timedelta(hours=(start_timestep-1)*dt_step)
    dt_dummy=timedelta(hours=idummy)
  case ('minutes')
    deltat_increment=timedelta(minutes=dt_step)
    deltat_since_file_start=timedelta(minutes=(start_timestep-1)*dt_step)
    dt_dummy=timedelta(minutes=idummy)
  case ('seconds')
    deltat_increment=timedelta(seconds=dt_step)
    deltat_since_file_start=timedelta(seconds=(start_timestep-1)*dt_step)
end select

dayxfile=dt_dummy%total_seconds()/86400
nhour=deltat_increment%total_seconds()/3600
current_date=file_start_date + deltat_since_file_start
print*,' '
print*,'NUMBER OF DAYS IN FILE INPUT: ',dayxfile
print*,'Time increment ',nhour,' hours'
print*,'FILE START DATE DATE: ',file_start_date%isoformat()
print*,'SIMULATION START DATE DATE: ',simulation_start_date%isoformat()
print*,'Start timestep ',start_timestep
print*,'End timestep   ',numtimestep
!$OMP END SINGLE
do timestep=start_timestep,numtimestep
!read wrf data
!$OMP SINGLE
  print*,' '
  print*,'Start timestep ',start_timestep
  print*,'End timestep   ',numtimestep
  print*,'Current timestep ',timestep
  print*,'CURRENT DATE: ',current_date%isoformat()
  year=current_date%GetYear()
  month=current_date%GetMonth()
  day=current_date%GetDay()
  hour=current_date%GetHour()
  write(current_yyyy,'(I4.4)'),year
  write(current_mm,'(I2.2)'),month
  write(current_dd,'(I2.2)'),day
  write(current_hh,'(I2.2)'),hour
  current_yyyymmdd=current_yyyy//current_mm//current_dd//current_hh
  output_filename=trim(dirout)//"/"//trim(out_prefix)//"_"//trim(current_yyyymmdd)//".nc"

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

  !compute total acc precipitation
  call wrf_prec

  !compute total precipitation
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

  ! Massimo: 
  ! 31 may 2022 Massimo: la rotate_wind_gap è concettualmente sbagliata secondo
  ! me e Gino. UTM conserva gli angoli, quindi va fatta la rotazione classica
  ! delle coordinate WRF alle coorfinate terrestri, chiamando rotate_wind, PRIMA
  ! dell'interpolazione orizzontale
  if ( ( rotate .eqv. .true. ) .and. ( utmcoordinate .eqv. .true. )) then
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

  ! scrive in output un file per ogni istante letto
  print*,'writing to output file ',trim(output_filename)
  print*,'---------------------------------------------------------------------'
  call wr_farm(output_filename)
!$OMP END SINGLE

  ! update the current_date
  current_date = current_date + deltat_increment
enddo !numtimestep
end program main

