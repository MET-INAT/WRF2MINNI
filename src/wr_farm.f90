       subroutine wr_farm(fileout)
 
       use general
       use param_minni
       use param_wrf
       use netcdf

       implicit none
    
       type (datetime) :: datetime_farm, nc_ref_time
       type (timedelta) :: hours_since

       integer :: ncid
       integer :: idx,idy,idz,idtime
       integer :: idxvar,idyvar,idzvar,idtimevar
       integer :: idsst,idrel,idtcc,idprec,idw,idp,idphi
       integer :: idsstsk,idsnow
       integer :: idu,idv,idt,idrh
       integer :: idpbl,idnetrad,idtotrad,idstotrad,idltotrad
       integer :: idqcloud,idqice,idqrain,idqsnow,idqgraup
       integer :: idsh,idlh,idgh,idz0
       integer :: idalbedo,idustar,idlstar
       integer :: yr,mon,day,hr,farmtime
       integer,dimension(8) :: values
       character*4 :: cyr
       character*2 :: cmon
       character*2 :: cday
       character*2 :: chr
       character*2 :: cmin
       character*2 :: csec
       character*3 :: cmillisec
       character(len=*) ::  fileout

!----------------------------------------------------------------------
!$OMP SINGLE
! Open the model output file

      yr=current_date%getYear()
      mon=current_date%getMonth()
      day=current_date%getDay()
      hr=current_date%getHour()
      write(cyr,'(I4)'),yr
      write(cmon,'(I2)'),mon
      write(cday,'(I2)'),day
      write(chr,'(I2)'),hr

      hours_since=timedelta()

      datetime_farm=datetime(yr,mon,day,hr)
      nc_ref_time=datetime(1900,01,01,00)
      hours_since=datetime_farm-nc_ref_time 
      farmtime=hours_since%total_seconds()/3600.
#ifdef debug
print*,'debug: farmtime ',farmtime
#endif

      call check(nf90_create(trim(fileout),nf90_clobber, ncid))
!Define dimensions
#ifdef debug
print*,'debug: define dimensions '
#endif
      call check(nf90_def_dim(ncid,"x",nx,idx))
      call check(nf90_def_dim(ncid,"y",ny,idy))
      call check(nf90_def_dim(ncid,"z",nz,idz))
      call check(nf90_def_dim(ncid,"time",nf90_unlimited,idtime))
!Define variables
#ifdef debug
print*,'debug: define vars '
#endif
      call check(nf90_def_var(ncid,"x",nf90_float,                      &
     &  (/ idx  /),idxvar))
      call check(nf90_def_var(ncid,"y",nf90_float,                      &
     &  (/ idy  /),idyvar))
      call check(nf90_def_var(ncid,"z",nf90_float,                      &
     &  (/ idz  /),idzvar))
      call check(nf90_def_var(ncid,"time",nf90_double,                  &
     &  (/ idtime  /),idtimevar))
       if ( fznt .eqv. .true. ) then
      call check(nf90_def_var(ncid,"Z0",nf90_float,                  &
     &  (/ idx, idy, idtime /),idz0))
       endif
       if ( flstar .eqv. .true. ) then
      call check(nf90_def_var(ncid,"LSTAR",nf90_float,                  &
     &  (/ idx, idy, idtime /),idlstar))
       endif
       if ( fustar .eqv. .true. ) then
      call check(nf90_def_var(ncid,"USTAR",nf90_float,                  &
     &  (/ idx, idy, idtime /),idustar))
       endif
       if ( falbedo .eqv. .true. ) then
      call check(nf90_def_var(ncid,"ALBEDO",nf90_float,                 &
     &  (/ idx, idy, idtime /),idalbedo))
       endif
       if ( fsstsk .eqv. .true. ) then
      call check(nf90_def_var(ncid,"SSTSK",nf90_float,                  &
     &  (/ idx, idy, idtime /),idsstsk))
       endif
      call check(nf90_def_var(ncid,"SST",nf90_float,                    &
     &  (/ idx, idy, idtime /),idsst))
      call check(nf90_def_var(ncid,"SNOW",nf90_float,                   &
     &  (/ idx, idy, idtime /),idsnow))
      call check(nf90_def_var(ncid,"REL",nf90_float,                    &
     &  (/ idx, idy, idtime /),idrel))
      call check(nf90_def_var(ncid,"TCC",nf90_float,                    &
     &  (/ idx, idy, idtime /),idtcc))
      call check(nf90_def_var(ncid,"PREC",nf90_float,                   &
     &  (/ idx, idy, idtime /),idprec))
      if ( fpbl .eqv. .true. ) then
      call check(nf90_def_var(ncid,"HMIX",nf90_float,                   &
     &  (/ idx, idy, idtime /),idpbl))
      endif
       if ( fnrad .eqv. .true. ) then
      call check(nf90_def_var(ncid,"NETRAD",nf90_float,                 &
     &  (/ idx, idy, idtime /),idnetrad))
       endif
       if ( ftgrad .eqv. .true. ) then
      call check(nf90_def_var(ncid,"TRAD",nf90_float,                  &
     &  (/ idx, idy, idtime /),idtotrad))
       endif
       if ( fsgrad .eqv. .true. ) then
! TOTRAD è definito swdown in wrf2farm_int.f90
! ed e' quello che sembra si aspetti surfpro come TOTRAD (non al somma sw+lw)
! almeno il pytohn di sandro lo definisce "Total SOLAR radiazion" quindi SW
      call check(nf90_def_var(ncid,"TOTRAD",nf90_float,                 &
     &  (/ idx, idy, idtime /),idstotrad))
       endif
       if ( flgrad .eqv. .true. ) then
      call check(nf90_def_var(ncid,"LRAD",nf90_float,                  &
     &  (/ idx, idy, idtime /),idltotrad))
       endif
       if ( fsh .eqv. .true. ) then
      call check(nf90_def_var(ncid,"SHF",nf90_float,                  &
     &  (/ idx, idy, idtime /),idsh))
       endif
       if ( flh .eqv. .true. ) then
      call check(nf90_def_var(ncid,"LHF",nf90_float,                  &
     &  (/ idx, idy, idtime /),idlh))
       endif
       if ( fgh .eqv. .true. ) then
      call check(nf90_def_var(ncid,"GHF",nf90_float,                  &
     &  (/ idx, idy, idtime /),idgh))
       endif

      call check(nf90_def_var(ncid,"U",nf90_float,                      &
     &  (/ idx, idy, idz, idtime /),idu))
      call check(nf90_def_var(ncid,"V",nf90_float,                      &
     &  (/ idx, idy, idz, idtime /),idv))
      call check(nf90_def_var(ncid,"T",nf90_float,                      &
     &  (/ idx, idy, idz, idtime /),idt))
      call check(nf90_def_var(ncid,"RH",nf90_float,                     &
     &  (/ idx, idy, idz, idtime /),idrh))
      call check(nf90_def_var(ncid,"W",nf90_float,                      &
     &  (/ idx, idy, idz, idtime /),idw))
      call check(nf90_def_var(ncid,"P",nf90_float,                      &
     &  (/ idx, idy, idz, idtime /),idp))
      call check(nf90_def_var(ncid,"PHI",nf90_float,                    &
     &  (/ idx, idy, idz, idtime /),idphi))
      if ( fhymet .eqv. .true. ) then
      call check(nf90_def_var(ncid,"QCLOUD",nf90_float,                 &
     &  (/ idx, idy, idz, idtime /),idqcloud))
      call check(nf90_def_var(ncid,"QSNOW",nf90_float,                  &
     &  (/ idx, idy, idz, idtime /),idqsnow))
      call check(nf90_def_var(ncid,"QRAIN",nf90_float,                  &
     &  (/ idx, idy, idz, idtime /),idqrain))
      call check(nf90_def_var(ncid,"QICE",nf90_float,                   &
     &  (/ idx, idy, idz, idtime /),idqice))
      call check(nf90_def_var(ncid,"QGRAUP",nf90_float,                 &
     &  (/ idx, idy, idz, idtime /),idqgraup))
      endif

#ifdef debug
print*,'debug: define var atrributes '
#endif

       if ( flstar .eqv. .true. ) then
      call check(nf90_put_att(ncid,idlstar,"units","m")) 
       endif
       if ( fznt .eqv. .true. ) then
      call check(nf90_put_att(ncid,idz0,"units","m")) 
       endif
       if ( fustar .eqv. .true. ) then
      call check(nf90_put_att(ncid,idustar,"units","M/S")) 
       endif
       if ( falbedo .eqv. .true. ) then
      call check(nf90_put_att(ncid,idalbedo,"units","-")) 
       endif
       if ( fsstsk .eqv. .true. ) then
      call check(nf90_put_att(ncid,idsstsk,"units","C")) 
       endif
      call check(nf90_put_att(ncid,idsst,"units","C")) 
      call check(nf90_put_att(ncid,idsnow,"units","-")) 
       if ( fnrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idnetrad,"units","W/m2")) 
       endif
       if ( ftgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idtotrad,"units","W/m2")) 
       endif
       if ( fsgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idstotrad,"units","W/m2")) 
       endif
       if ( flgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idltotrad,"units","W/m2")) 
       endif
       if ( fsh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idsh,"units","W/m2")) 
       endif
       if ( flh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idlh,"units","W/m2")) 
       endif
       if ( fgh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idgh,"units","W/m2")) 
       endif
      if ( fpbl .eqv. .true. ) then
      call check(nf90_put_att(ncid,idpbl,"units","m")) 
      endif
      if ( fhymet .eqv. .true. ) then
      call check(nf90_put_att(ncid,idqcloud,"units","kg water/kg air")) 
      call check(nf90_put_att(ncid,idqice,"units","kg water/kg air")) 
      call check(nf90_put_att(ncid,idqrain,"units","kg water/kg air")) 
      call check(nf90_put_att(ncid,idqsnow,"units","kg water/kg air")) 
      call check(nf90_put_att(ncid,idqgraup,"units","kg water/kg air")) 
      endif
      call check(nf90_put_att(ncid,idrel,"units","m")) 
      call check(nf90_put_att(ncid,idtcc,"units","FREC")) 
      call check(nf90_put_att(ncid,idprec,"units","MM/H")) 
      call check(nf90_put_att(ncid,idphi,"units","m*g")) 
      call check(nf90_put_att(ncid,idu,"units","M/S")) 
      call check(nf90_put_att(ncid,idv,"units","M/S")) 
      call check(nf90_put_att(ncid,idw,"units","M/S")) 
      call check(nf90_put_att(ncid,idt,"units","DEG.K")) 
      call check(nf90_put_att(ncid,idp,"units","HPA.")) 
      call check(nf90_put_att(ncid,idrh,"units","PCT")) 
      if ( utmcoordinate .eqv. .true. ) then
      call check(nf90_put_att(ncid,idx,"units","km")) 
      call check(nf90_put_att(ncid,idy,"units","km")) 
      else
      call check(nf90_put_att(ncid,idx,"units","deg")) 
      call check(nf90_put_att(ncid,idy,"units","deg")) 
      endif
      call check(nf90_put_att(ncid,idz,"units","m")) 
      call check(nf90_put_att(ncid,idtime,"units",                      &
     &                        "hours since 1900-1-1 00:00:0.0")) 
      call check(nf90_put_att(ncid,idtime,"delta_t",                    &
     &                        "0000-00-00 01:00:00.00 +00:00")) 

       if ( flstar .eqv. .true. ) then
      call check(nf90_put_att(ncid,idlstar,"missing_value",real(-9.96921e+36))) 
       endif
       if ( fznt .eqv. .true. ) then
      call check(nf90_put_att(ncid,idz0,"missing_value",real(-9.96921e+36))) 
       endif
       if ( fustar .eqv. .true. ) then
      call check(nf90_put_att(ncid,idustar,"missing_value",real(-9.96921e+36))) 
       endif
       if ( falbedo .eqv. .true. ) then
      call check(nf90_put_att(ncid,idalbedo,"missing_value",real(-9.96921e+36))) 
       endif
       if ( fsstsk .eqv. .true. ) then
      call check(nf90_put_att(ncid,idsstsk,"missing_value",real(-9.96921e+36))) 
       endif
      call check(nf90_put_att(ncid,idsst,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idsnow,"missing_value",real(-9.96921e+36))) 
       if ( fnrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idnetrad,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idnetrad,"long_name","net radiation")) 
       endif
       if ( ftgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idtotrad,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idtotrad,"std_name","Total Radiation")) 
      call check(nf90_put_att(ncid,idtotrad,"long_name","Short and long wave radiation")) 
       endif
       if ( fsgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idstotrad,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idstotrad,"long_name","Total solar radiation")) 
      call check(nf90_put_att(ncid,idstotrad,"std_name","Global radiation")) 
       endif
       if ( flgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idltotrad,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idltotrad,"long_name","long wave radiation")) 
       endif
       if ( fsh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idsh,"missing_value",real(-9.96921e+36))) 
       endif
       if ( flh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idlh,"missing_value",real(-9.96921e+36))) 
       endif
       if ( fgh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idgh,"missing_value",real(-9.96921e+36))) 
       endif
      if ( fpbl .eqv. .true. ) then
      call check(nf90_put_att(ncid,idpbl,"missing_value",real(-9.96921e+36))) 
      endif
      call check(nf90_put_att(ncid,idrel,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idtcc,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idprec,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idphi,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idu,"missing_value",real(-9.96921e+36)))
      call check(nf90_put_att(ncid,idv,"missing_value",real(-9.96921e+36)))
      call check(nf90_put_att(ncid,idw,"missing_value",real(-9.96921e+36)))
      call check(nf90_put_att(ncid,idt,"missing_value",real(-9.96921e+36)))
      call check(nf90_put_att(ncid,idp,"missing_value",real(-9.96921e+36)))
      call check(nf90_put_att(ncid,idrh,"missing_value",real(-9.96921e+36))) 
      if ( fhymet .eqv. .true. ) then
      call check(nf90_put_att(ncid,idqcloud,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idqice,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idqrain,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idqsnow,"missing_value",real(-9.96921e+36))) 
      call check(nf90_put_att(ncid,idqgraup,"missing_value",real(-9.96921e+36))) 
      endif

       if ( flstar .eqv. .true. ) then
      call check(nf90_put_att(ncid,idlstar,"add_offset",real(0))) 
       endif
       if ( fznt .eqv. .true. ) then
      call check(nf90_put_att(ncid,idz0,"add_offset",real(0))) 
       endif
       if ( fustar .eqv. .true. ) then
      call check(nf90_put_att(ncid,idustar,"add_offset",real(0))) 
       endif
       if ( falbedo .eqv. .true. ) then
      call check(nf90_put_att(ncid,idalbedo,"add_offset",real(0))) 
       endif
       if ( fsstsk .eqv. .true. ) then
      call check(nf90_put_att(ncid,idsstsk,"add_offset",real(0))) 
       endif
      call check(nf90_put_att(ncid,idsst,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idsnow,"add_offset",real(0))) 
       if ( fnrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idnetrad,"add_offset",real(0))) 
       endif
       if ( ftgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idtotrad,"add_offset",real(0))) 
       endif
       if ( fsgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idstotrad,"add_offset",real(0))) 
       endif
       if ( flgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idltotrad,"add_offset",real(0))) 
       endif
       if ( fsh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idsh,"add_offset",real(0))) 
       endif
       if ( flh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idlh,"add_offset",real(0))) 
       endif
       if ( fgh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idgh,"add_offset",real(0))) 
       endif
      if ( fpbl .eqv. .true. ) then
      call check(nf90_put_att(ncid,idpbl,"add_offset",real(0))) 
      endif
      if ( fhymet .eqv. .true. ) then
      call check(nf90_put_att(ncid,idqcloud,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idqice,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idqrain,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idqsnow,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idqgraup,"add_offset",real(0))) 
      endif
      call check(nf90_put_att(ncid,idrel,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idtcc,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idprec,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idphi,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idu,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idv,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idw,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idt,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idp,"add_offset",real(0))) 
      call check(nf90_put_att(ncid,idrh,"add_offset",real(0))) 

       if ( flstar .eqv. .true. ) then
      call check(nf90_put_att(ncid,idlstar,"scale_factor",real(1))) 
       endif
       if ( fznt .eqv. .true. ) then
      call check(nf90_put_att(ncid,idz0,"scale_factor",real(1))) 
       endif
       if ( fustar .eqv. .true. ) then
      call check(nf90_put_att(ncid,idustar,"scale_factor",real(1))) 
       endif
       if ( falbedo .eqv. .true. ) then
      call check(nf90_put_att(ncid,idalbedo,"scale_factor",real(1))) 
       endif
       if ( fsstsk .eqv. .true. ) then
      call check(nf90_put_att(ncid,idsstsk,"scale_factor",real(1))) 
       endif
      call check(nf90_put_att(ncid,idsst,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idsnow,"scale_factor",real(1))) 
       if ( fnrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idnetrad,"scale_factor",real(1))) 
       endif
       if ( ftgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idtotrad,"scale_factor",real(1))) 
       endif
       if ( fsgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idstotrad,"scale_factor",real(1))) 
       endif
       if ( flgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idltotrad,"scale_factor",real(1))) 
       endif
       if ( fsh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idsh,"scale_factor",real(1))) 
       endif
       if ( flh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idlh,"scale_factor",real(1))) 
       endif
       if ( fgh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idgh,"scale_factor",real(1))) 
       endif
      if ( fpbl .eqv. .true. ) then
      call check(nf90_put_att(ncid,idpbl,"scale_factor",real(1))) 
      endif
      if ( fhymet .eqv. .true. ) then
      call check(nf90_put_att(ncid,idqcloud,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idqrain,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idqsnow,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idqice,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idqgraup,"scale_factor",real(1))) 
      endif
      call check(nf90_put_att(ncid,idrel,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idtcc,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idprec,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idphi,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idu,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idv,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idw,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idt,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idp,"scale_factor",real(1))) 
      call check(nf90_put_att(ncid,idrh,"scale_factor",real(1))) 

       if ( flstar .eqv. .true. ) then
      call check(nf90_put_att(ncid,idlstar,"actual_range",             &
     &     real((/minval(lstarfarm),maxval(lstarfarm)/)))) 
       endif
       if ( fznt .eqv. .true. ) then
      call check(nf90_put_att(ncid,idz0,"actual_range",             &
     &     real((/minval(z0farm),maxval(z0farm)/)))) 
       endif
       if ( fustar .eqv. .true. ) then
      call check(nf90_put_att(ncid,idustar,"actual_range",             &
     &     real((/minval(ustarfarm),maxval(ustarfarm)/)))) 
       endif
       if ( falbedo .eqv. .true. ) then
      call check(nf90_put_att(ncid,idalbedo,"actual_range",             &
     &     real((/minval(albedofarm),maxval(albedofarm)/)))) 
       endif
       if ( fsstsk .eqv. .true. ) then
      call check(nf90_put_att(ncid,idsst,"actual_range",                &
     &     real((/minval(sstskfarm),maxval(sstskfarm)/)))) 
       endif
      call check(nf90_put_att(ncid,idsst,"actual_range",                &
     &     real((/minval(sstfarm),maxval(sstfarm)/)))) 
      call check(nf90_put_att(ncid,idsnow,"actual_range",                &
     &     real((/minval(snowfarm),maxval(snowfarm)/)))) 
       if ( fnrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idnetrad,"actual_range",             &
     &     real((/minval(netradfarm),maxval(netradfarm)/)))) 
       endif
       if ( ftgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idtotrad,"actual_range",             &
     &     real((/minval(totradfarm),maxval(totradfarm)/)))) 
       endif
       if ( fsgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idstotrad,"actual_range",            &
     &     real((/minval(totsradfarm),maxval(totsradfarm)/)))) 
       endif
       if ( flgrad .eqv. .true. ) then
      call check(nf90_put_att(ncid,idltotrad,"actual_range",            &
     &     real((/minval(totlradfarm),maxval(totlradfarm)/)))) 
       endif
       if ( fsh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idsh,"actual_range",            &
     &     real((/minval(shfarm),maxval(shfarm)/)))) 
       endif
       if ( flh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idlh,"actual_range",            &
     &     real((/minval(lhfarm),maxval(lhfarm)/)))) 
       endif
       if ( fgh .eqv. .true. ) then
      call check(nf90_put_att(ncid,idgh,"actual_range",            &
     &     real((/minval(ghfarm),maxval(ghfarm)/)))) 
       endif
      if ( fpbl .eqv. .true. ) then
      call check(nf90_put_att(ncid,idpbl,"actual_range",                &
     &     real((/minval(pblfarm),maxval(pblfarm)/)))) 
      endif
      if ( fhymet .eqv. .true. ) then
      call check(nf90_put_att(ncid,idqcloud,"actual_range",             &
     &     real((/minval(qcloudfarm),maxval(qcloudfarm)/)))) 
      call check(nf90_put_att(ncid,idqice,"actual_range",               &
     &     real((/minval(qicefarm),maxval(qicefarm)/)))) 
      call check(nf90_put_att(ncid,idqrain,"actual_range",              &
     &     real((/minval(qrainfarm),maxval(qrainfarm)/)))) 
      call check(nf90_put_att(ncid,idqsnow,"actual_range",              &
     &     real((/minval(qsnowfarm),maxval(qsnowfarm)/)))) 
      call check(nf90_put_att(ncid,idqgraup,"actual_range",             &
     &     real((/minval(qgraupfarm),maxval(qgraupfarm)/)))) 
      endif
      call check(nf90_put_att(ncid,idrel,"actual_range",                &
     &     real((/minval(hgtfarm),maxval(hgtfarm)/)))) 
      call check(nf90_put_att(ncid,idtcc,"actual_range",                &
     &     real((/minval(tccfarm),maxval(tccfarm)/)))) 
      call check(nf90_put_att(ncid,idprec,"actual_range",               &
     &     real((/minval(precfarm),maxval(precfarm)/)))) 
      call check(nf90_put_att(ncid,idphi,"actual_range",                &
     &     real((/minval(phifarm),maxval(phifarm)/)))) 
      call check(nf90_put_att(ncid,idu,"actual_range",                  &
     &     real((/minval(ufarm),maxval(ufarm)/)))) 
      call check(nf90_put_att(ncid,idv,"actual_range",                  &
     &     real((/minval(vfarm),maxval(vfarm)/)))) 
      call check(nf90_put_att(ncid,idw,"actual_range",                  &
     &     real((/minval(wfarm),maxval(wfarm)/)))) 
      call check(nf90_put_att(ncid,idt,"actual_range",                  &
     &     real((/minval(tfarm),maxval(tfarm)/)))) 
      call check(nf90_put_att(ncid,idp,"actual_range",                  &
     &     real((/minval(pfarm),maxval(pfarm)/)))) 
      call check(nf90_put_att(ncid,idrh,"actual_range",                 &
     &     real((/minval(rhfarm),maxval(rhfarm)/)))) 

#ifdef debug
print*,'debug: define global attributes'
#endif
      call date_and_time(VALUES=values)
       write(cyr,'(I4)'),values(1)
       write(cmon,'(I2)'),values(2)
       write(cday,'(I2)'),values(3)
       write(chr,'(I2)'),values(5)
       write(cmin,'(I2)'),values(6)
       write(csec,'(I2)'),values(7)
       write(cmillisec,'(I3)'),values(8)
      call check(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","COARDS")) 
      call check(nf90_put_att(ncid,NF90_GLOBAL,"history",""))        
      call check(nf90_put_att(ncid,NF90_GLOBAL,"description",""))    
      call check(nf90_put_att(ncid,NF90_GLOBAL,"model","FARMZ"))    
      call check(nf90_put_att(ncid,NF90_GLOBAL,"lib_ver",""))    
      call check(nf90_put_att(ncid,NF90_GLOBAL,"creation_time",         &
     &cyr//" "//cmon//" "//cday//" H "//chr//"."//cmin//"."//csec//"."  &
     &//cmillisec))    


      call check(nf90_enddef(ncid))

#ifdef debug
print*,'debug: fill variables...'
#endif

!Fill in variable
      call check(nf90_put_var(ncid,idxvar,xfarm(:,1),                   &
     &                start = (/ 1 /),                                  &
                      count = (/ nx /)   ))
      call check(nf90_put_var(ncid,idyvar,yfarm(1,:),                   &
     &                start = (/ 1 /),                                  &
                      count = (/ ny /)   ))
      call check(nf90_put_var(ncid,idzvar,zlev,                         &
     &                start = (/ 1 /),                                  &
                      count = (/ nz /)   ))
      call check(nf90_put_var(ncid,idtimevar,dble(farmtime),            &
     &                start = (/ 1 /) ))                                
       if ( flstar .eqv. .true. ) then
      call check(nf90_put_var(ncid,idlstar,lstarfarm,                   &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
       if ( fznt .eqv. .true. ) then
      call check(nf90_put_var(ncid,idz0,z0farm,                   &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
       if ( fustar .eqv. .true. ) then
      call check(nf90_put_var(ncid,idustar,ustarfarm,                   &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
       if ( falbedo .eqv. .true. ) then
      call check(nf90_put_var(ncid,idalbedo,albedofarm,                 &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
       if ( fsstsk .eqv. .true. ) then
      call check(nf90_put_var(ncid,idsstsk,sstskfarm,                   &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
      call check(nf90_put_var(ncid,idsst,sstfarm,                       &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
      call check(nf90_put_var(ncid,idsnow,snowfarm,                       &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       if ( fnrad .eqv. .true. ) then
      call check(nf90_put_var(ncid,idnetrad,netradfarm,                 &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
       if ( ftgrad .eqv. .true. ) then
      call check(nf90_put_var(ncid,idtotrad,totradfarm,                 &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
       if ( fsgrad .eqv. .true. ) then
      call check(nf90_put_var(ncid,idstotrad,totsradfarm,               &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
       if ( flgrad .eqv. .true. ) then
      call check(nf90_put_var(ncid,idltotrad,totlradfarm,               &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
       if ( fsh .eqv. .true. ) then
      call check(nf90_put_var(ncid,idsh,shfarm,                         &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
       if ( flh .eqv. .true. ) then
      call check(nf90_put_var(ncid,idlh,lhfarm,                         &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
       if ( fgh .eqv. .true. ) then
      call check(nf90_put_var(ncid,idgh,ghfarm,                         &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
       endif
      call check(nf90_put_var(ncid,idrel,hgtfarm,                       &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
      call check(nf90_put_var(ncid,idprec,precfarm,                     &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
      call check(nf90_put_var(ncid,idtcc,tccfarm,                       &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
      if ( fpbl .eqv. .true. ) then
      call check(nf90_put_var(ncid,idpbl,pblfarm,                       &
     &                start = (/ 1, 1, 1 /),                     &
                      count = (/ nx, ny, 1 /)   ))
      endif
      call check(nf90_put_var(ncid,idu,ufarm,                           &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      call check(nf90_put_var(ncid,idv,vfarm,                           &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      call check(nf90_put_var(ncid,idt,tfarm,                           &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      call check(nf90_put_var(ncid,idrh,rhfarm,                         &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      call check(nf90_put_var(ncid,idw,wfarm,                           &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      call check(nf90_put_var(ncid,idp,pfarm,                           &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      call check(nf90_put_var(ncid,idphi,phifarm,                       &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      if ( fhymet .eqv. .true. ) then
      call check(nf90_put_var(ncid,idqcloud,qcloudfarm,                       &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      call check(nf90_put_var(ncid,idqice,qicefarm,                       &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      call check(nf90_put_var(ncid,idqrain,qrainfarm,                       &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      call check(nf90_put_var(ncid,idqsnow,qsnowfarm,                       &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      call check(nf90_put_var(ncid,idqgraup,qgraupfarm,                       &
     &                start = (/ 1, 1, 1,1 /),                   &
                      count = (/ nx, ny, nz, 1 /)   ))
      endif

#ifdef debug
print*,'debug: close file '
#endif
      call check(nf90_close(ncid))
!$OMP END SINGLE

      end
