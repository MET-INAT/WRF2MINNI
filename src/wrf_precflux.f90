      subroutine wrf_precflux(timestep) 
      use general
      use param_wrf
      use netcdf
      IMPLICIT NONE
      real         :: prec0
      real,dimension(west_east,south_north) :: dummy_rainc,dummy_irainc
      real,dimension(west_east,south_north) :: dummy_rainnc,dummy_irainnc
      real,dimension(west_east,south_north) :: dummy_rainsh
      integer      :: timestep,idtime,timetoload
      integer      :: i,j
      integer      :: year,month,day, hour, yearm1
      integer      :: ncid
      character    :: cdummy
      integer :: idRAINNC,idRAINC,idRAINSH,idi_rainnc, idi_rainc
      character*33 :: filetoread
      character*4 :: cyear
      character*2 :: cmonth
      character*2 :: cday
      character*2 :: chour
      logical :: file_exists
      type(datetime) :: old_date

#ifdef debug
  print*,'debug: call compute prec flux'
#endif


!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP
!$OMP SHARED(precflx,dummy_rainsh,dummy_irainc,dummy_rainc,dummy_irainnc,dummy_rainnc,west_east,south_north)       &
!$OMP PRIVATE(i,j)
       do j=1,south_north
       do i=1,west_east
         dummy_rainc(i,j)=0.0
         dummy_rainnc(i,j)=0.0
         dummy_rainsh(i,j)=0.0
         dummy_irainc(i,j)=0.0
         dummy_irainnc(i,j)=0.0
         precflx(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
       prec0=0.0
!OMP BARRIER   

!$OMP SINGLE
      if (( timestep .eq. 1 )  )  then
           if(dayxfile == 365 ) then
             year=current_date%getYear()
             yearm1=year-1
             if (isLeapYear(yearm1).eqv. .true.)then
               dayxfile=366
             else
               dayxfile=365
             endif
           endif
          old_date=current_date-timedelta(days=dayxfile) 
          year=old_date%getYear()
          month=old_date%getMonth()
          day=old_date%getDay()
          hour=old_date%getHour()
          write(cyear,'(I4)'),year
          write(cmonth,'(I2)'),month
          write(cday,'(I2)'),day
          write(chour,'(I2)'),hour
          if (month .le. 9) cmonth(1:1)="0"
          if (day .le. 9) cday(1:1)="0"
          if (hour .le. 9) chour(1:1)="0"
          !okkio, il file netcdf deve essere nominato assolutamente
          ! come un wrfout. Es: wrfout_d01_1980-01-01_00:00:00.nc
          filetoread=trim(filein(1:11))//cyear//"-"//cmonth//"-"//cday//"_"   &
     &         //chour//":00:00.nc"     
          inquire(file=filetoread, exist=file_exists)
          if(file_exists .eqv. .false.)then
            filetoread=trim(filein)
          endif
      else
         filetoread=trim(filein)
      endif ! first timestep
      !prec  is set to zero, but if current date is not 
      !start simulation date, then compute prec=prec(t)-prec(t-1)
      if ( current_date /= simulation_start_date ) then
         call check(nf90_open(trim(dirin)//trim(filetoread),               &
     &              nf90_nowrite, ncid))
         if (( timestep .eq. 1 ) )  then
            call check(nf90_inq_dimid(ncid, 'Time',idtime))
            call check(nf90_inquire_dimension(ncid,idTime,cdummy,       &
     &                 timetoload))
         else 
             timetoload = timestep-1
         endif
      print*,'COMPUTING precflux, from file ',trim(dirin)//trim(filetoread)
      print*,'timestep ',timetoload
         call check(nf90_inq_varid(ncid,'RAINC',idRAINC))
         call check(nf90_get_var(ncid,idRAINC,dummy_rainc,                   &
              start = (/ 1, 1, timetoload /),                           &
              count = (/ west_east, south_north, 1 /) ))


         call check(nf90_inq_varid(ncid,'RAINNC',idRAINNC))
         call check(nf90_get_var(ncid,idRAINNC,dummy_rainnc,                  &
              start = (/ 1, 1, timetoload /),                           &
              count = (/ west_east, south_north, 1 /) ))

 if (if_bucket .eqv. .true.)then
   call check(nf90_inq_varid(ncid,'I_RAINC',idi_rainc))
   call check(nf90_get_var(ncid,idi_rainc,dummy_irainc,                 &
              start = (/ 1, 1, timetoload /),                           &
              count = (/ west_east, south_north, 1 /) ))
   call check(nf90_inq_varid(ncid,'I_RAINNC',idi_rainnc))
   call check(nf90_get_var(ncid,idi_rainnc,dummy_irainnc,               &
              start = (/ 1, 1, timetoload /),                           &
              count = (/ west_east, south_north, 1 /) ))
 endif

!Massimo:
! per calcolare la precipitazione non occorrono le idrometeore
! da quanto ho capito nei vari forum di wrf, il totale e' in RAINC e
! RAINNC (al limite RAINSH che pero' e' in output solo
! se si attiva una particolare parametrizz della convez)
          call check(nf90_inq_varid(ncid,'RAINSH',idRAINSH))
          call check(nf90_get_var(ncid,idRAINSH,dummy_rainsh,                  &
               start = (/ 1, 1, timetoload /),                           &
               count = (/ west_east, south_north, 1 /) ))
endif
!$OMP END SINGLE
!compute precipitation flux
if (if_bucket .eqv. .true.)then
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(totprec,nhour,precflx,dummy_rainsh,dummy_irainc,dummy_rainc,dummy_irainnc,bucket_threshold_mm,dummy_rainnc,west_east,south_north)       &
!$OMP PRIVATE(i,j,prec0)
  do j=1,south_north
    do i=1,west_east
      dummy_rainnc(i,j)=dummy_rainnc(i,j)+bucket_threshold_mm*dummy_irainnc(i,j)
      dummy_rainc(i,j)=dummy_rainc(i,j)+bucket_threshold_mm*dummy_irainc(i,j)
      prec0=dummy_rainc(i,j)+dummy_rainnc(i,j) + dummy_rainsh(i,j)
      precflx(i,j)= ( totprec(i,j)-prec0 ) / nhour
    enddo
  enddo
!$OMP END PARALLEL DO
!OMP BARRIER         
else
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(totprec,nhour,precflx,dummy_rainsh,dummy_irainc,dummy_rainc,dummy_rainnc,west_east,south_north)       &
!$OMP PRIVATE(i,j,prec0)
  do j=1,south_north
    do i=1,west_east
      prec0=dummy_rainc(i,j)+dummy_rainnc(i,j) + dummy_rainsh(i,j)
      precflx(i,j)= ( totprec(i,j)-prec0 ) / nhour
    enddo
   enddo
!$OMP END PARALLEL DO
!OMP BARRIER         
endif


       end
