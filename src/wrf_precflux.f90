      subroutine wrf_precflux(timestep) 
      use general
      use param_wrf
      use netcdf
      IMPLICIT NONE
      real         :: prec0
      real         :: dummy1(west_east,south_north)
      real         :: dummy2(west_east,south_north)
      real         :: dummy3(west_east,south_north)
      real         :: dummy4(west_east,south_north)
      real         :: dummy5(west_east,south_north)
      real         :: dummy6(west_east,south_north)
      integer      :: timestep,idtime,timetoload
      integer      :: i,j
      integer      :: yr,mon,day,julday,jd
      integer      :: ncid
      character    :: cdummy
      integer :: idRAINNC,idRAINC,idRAINSH
      integer :: idGRAUPELNC,idHAILNC,idSNOWNC
      character*30 :: fileold
      character*4 :: cyr
      character*2 :: cmon
      character*2 :: cday
      character*2 :: chr

!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(precflx,dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,south_north,west_east) &
!$OMP PRIVATE(i,j)
       do j=1,south_north
       do i=1,west_east
         dummy1(i,j)=0.0
         dummy2(i,j)=0.0
         dummy3(i,j)=0.0
         dummy4(i,j)=0.0
         dummy5(i,j)=0.0
         dummy6(i,j)=0.0
         precflx(i,j)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
       prec0=0.0
!OMP BARRIER   

!$OMP SINGLE
      if (( timestep .eq. 1 )  )  then
          cyr=Times(1:4)
          cmon=Times(6:7)
          cday=Times(9:10)
          chr=Times(12:13)
          read(cyr,'(I4)'),yr
          read(cmon,'(I2)'),mon
          read(cday,'(I2)'),day
          jd=julday(mon,day,yr)     
          jd=jd - dayxfile
          call caldat(jd,mon,day,yr)
          write(cyr,'(I4)'),yr
          write(cmon,'(I2)'),mon
          write(cday,'(I2)'),day
          if (mon .le. 9) cmon(1:1)="0"
          if (day .le. 9) cday(1:1)="0"
          fileold=trim(filein(1:11))//cyr//"-"//cmon//"-"//cday//"_"   &
     &            //chr//":00:00"     
      else
         fileold=trim(filein)
      endif
      if ( Time_startsim /= Times ) then
         call check(nf90_open(trim(dirin)//trim(fileold),               &
     &              nf90_nowrite, ncid))
         if (( timestep .eq. 1 ) )  then
            call check(nf90_inq_dimid(ncid, 'Time',idtime))
            call check(nf90_inquire_dimension(ncid,idTime,cdummy,       &
     &                 timetoload))
         else 
             timetoload = timestep-1
         endif
         call check(nf90_inq_varid(ncid,'RAINC',idRAINC))
         call check(nf90_get_var(ncid,idRAINC,dummy1,                   &
              start = (/ 1, 1, timetoload /),                           &
              count = (/ west_east, south_north, 1 /) ))

         call check(nf90_inq_varid(ncid,'RAINNC',idRAINNC))
         call check(nf90_get_var(ncid,idRAINNC,dummy2,                  &
              start = (/ 1, 1, timetoload /),                           &
              count = (/ west_east, south_north, 1 /) ))

!Massimo:
! per calcolare la precipitazione non occorrono le idrometeore
! da quanto ho capito nei vari forum di wrf, il totale e' in RAINC e
! RAINNC (al limite RAINSH che pero' e' in output solo
! se si attiva una particolare parametrizz della convez)
          call check(nf90_inq_varid(ncid,'RAINSH',idRAINSH))
          call check(nf90_get_var(ncid,idRAINSH,dummy3,                  &
               start = (/ 1, 1, timetoload /),                           &
               count = (/ west_east, south_north, 1 /) ))
!
!         call check(nf90_inq_varid(ncid,'SNOWNC',idSNOWNC))
!         call check(nf90_get_var(ncid,idSNOWNC,dummy4,                  &
!              start = (/ 1, 1, timetoload /),                           &
!              count = (/ west_east, south_north, 1 /) ))
!
!         call check(nf90_inq_varid(ncid,'HAILNC',idHAILNC))
!         call check(nf90_get_var(ncid,idHAILNC,dummy5,                  &
!              start = (/ 1, 1, timetoload /),                           &
!              count = (/ west_east, south_north, 1 /) ))
!
!         call check(nf90_inq_varid(ncid,'GRAUPELNC',idGRAUPELNC))
!         call check(nf90_get_var(ncid,idGRAUPELNC,dummy6,               &
!              start = (/ 1, 1, timetoload /),                           &
!              count = (/ west_east, south_north, 1 /) ))
!         call check(nf90_close(ncid))
      endif
!$OMP END SINGLE
!compute precipitation flux
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(totprec,precflx,dummy1,dummy2,dummy3,nhour,south_north,west_east ) &
!$OMP PRIVATE(i,j,prec0)
!!!$OMP SHARED(totprec,precflx,dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,nhour,south_north,west_east ) &
       do j=1,south_north
       do i=1,west_east
         !prec0=dummy1(i,j)+dummy2(i,j)+dummy3(i,j)+dummy4(i,j)+         &
     &   !      dummy5(i,j)+dummy6(i,j)
         prec0=dummy1(i,j)+dummy2(i,j) + dummy3(i,j)
         precflx(i,j)= ( totprec(i,j)-prec0 ) / nhour
       enddo
       enddo
!$OMP END PARALLEL DO
!OMP BARRIER         


       end
