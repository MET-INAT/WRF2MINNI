      subroutine read_wrf(timestep)

      use general
      use param_wrf
      use netcdf
       
      IMPLICIT NONE

      integer :: ncid
      integer :: idTimes,idXLAT,idXLONG
      integer :: idXLAT_V,idXLONG_V
      integer :: idXLAT_U,idXLONG_U
      integer :: idU,idV,idW,idPH,idPHB,idT,idP,idPB
      integer :: idQVAPOR
      integer :: idRAINNC,idRAINC,idRAINSH
      integer :: idGRAUPELNC,idHAILNC,idSNOWNC
      integer :: idSST,idSNOWH,idZNT
      integer :: idSSTSK
      integer :: idTSK, idEMISS
      integer :: idalbedo,idust
      integer :: idgrdflx,idhfx,idlh,idswdown,idglw
      integer :: idHGT,idT2,idU10,idV10,idPSFC,idPBLH
      integer :: idCOSALPHA,idSINALPHA
      integer :: idqcloud,idqice,idqrain,idqsnow,idqgraup
      character*500 :: fpd
      character :: cdummy
      integer :: timestep
      
!----------------------------------------------------------------------
! Open model output file
      fpd=trim(dirin)//trim(filein)
      call check(nf90_open(trim(fpd),nf90_nowrite, ncid))
!massimo debug
print*,"trim(dirin)//trim(filein)",trim(dirin)//trim(filein)
print*,'ncid ',ncid

print*,"debug: times"
      call check(nf90_inq_varid(ncid,'Times',idTimes))
      call check(nf90_get_var(ncid,idTimes, Times,                      &
           start = (/ 1, timestep /),                                   &
           count = (/ 19, 1 /) ))   
      print*,'Times, timestep ',Times,' ',timestep
      if ( timestep  .eq. start_timestep ) then

print*,'debug: xlat '
         call check(nf90_inq_varid(ncid,'XLAT',idXLAT))
         call check(nf90_get_var(ncid,idXLAT, XLAT,                     &
              start = (/ 1, 1, timestep /),                             &
              count = (/ west_east, south_north, 1 /) ))   

print*,'debug: xlong '
         call check(nf90_inq_varid(ncid,'XLONG',idXLONG))
         call check(nf90_get_var(ncid,idXLONG, XLONG,                   &
              start = (/ 1, 1, timestep /),                             &
              count = (/ west_east, south_north, 1 /) )) 

print*,'debug: lat_u '
         call check(nf90_inq_varid(ncid,'XLAT_U',idXLAT_U))
         call check(nf90_get_var(ncid,idXLAT_U, XLAT_U,                 &
              start = (/ 1, 1, timestep /),                             &
              count = (/ west_east_stag, south_north, 1 /) )) 

print*,'debug: long_u '
         call check(nf90_inq_varid(ncid,'XLONG_U',idXLONG_U))
         call check(nf90_get_var(ncid,idXLONG_U, XLONG_U,               &
              start = (/ 1, 1, timestep /),                             &
              count = (/ west_east_stag, south_north, 1 /) )) 

print*,'debug: lat_V '
         call check(nf90_inq_varid(ncid,'XLAT_V',idXLAT_V))
         call check(nf90_get_var(ncid,idXLAT_V, XLAT_V,                 &
              start = (/ 1, 1, timestep /),                             &
              count = (/ west_east, south_north_stag, 1 /) ))   

print*,'debug: long_v '
         call check(nf90_inq_varid(ncid,'XLONG_V',idXLONG_V))
         call check(nf90_get_var(ncid,idXLONG_V, XLONG_V,               &
              start = (/ 1, 1, timestep /),                             &
              count = (/ west_east, south_north_stag, 1 /) )) 
      endif
     
print*,'debug: u '
      call check(nf90_inq_varid(ncid,'U',idU))
      call check(nf90_get_var(ncid,idU, U,                              &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east_stag, south_north, bottom_top, 1 /) ))

print*,'debug: v '
      call check(nf90_inq_varid(ncid,'V',idV))
      call check(nf90_get_var(ncid,idV, V,                              &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north_stag, bottom_top, 1 /) ))

print*,'debug: w '
      call check(nf90_inq_varid(ncid,'W',idW))
      call check(nf90_get_var(ncid,idW, W,                              &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top_stag, 1 /) ))

print*,'debug: rh '
      call check(nf90_inq_varid(ncid,'PH',idPH))
      call check(nf90_get_var(ncid,idPH, PH,                            &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top_stag, 1 /) ))

print*,'debug: phb '
      call check(nf90_inq_varid(ncid,'PHB',idPHB))
      call check(nf90_get_var(ncid,idPHB, PHB,                          &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top_stag, 1 /) ))

print*,'debug: t '
      call check(nf90_inq_varid(ncid,'T',idT))
      call check(nf90_get_var(ncid,idT, T,                              &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top, 1 /) ))

print*,'debug: p '
      call check(nf90_inq_varid(ncid,'P',idP))
      call check(nf90_get_var(ncid,idP, P,                              &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top, 1 /) ))

print*,'debug: pb '
      call check(nf90_inq_varid(ncid,'PB',idPB))
      call check(nf90_get_var(ncid,idPB, PB,                            &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top, 1 /) ))

print*,'debug: qvapo '
      call check(nf90_inq_varid(ncid,'QVAPOR',idQVAPOR))
      call check(nf90_get_var(ncid,idQVAPOR,QVAPOR,                     &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top, 1 /) ))

print*,'debug: rainc '
      call check(nf90_inq_varid(ncid,'RAINC',idRAINC))
      call check(nf90_get_var(ncid,idRAINC,RAINC,                       &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))

print*,'debug: rinnc '
      call check(nf90_inq_varid(ncid,'RAINNC',idRAINNC))
      call check(nf90_get_var(ncid,idRAINNC,RAINNC,                     &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))

print*,'debug: rainsh '
      call check(nf90_inq_varid(ncid,'RAINSH',idRAINSH))
      call check(nf90_get_var(ncid,idRAINSH,RAINSH,                     &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))

print*,'debug: snowc '
      call check(nf90_inq_varid(ncid,'SNOWNC',idSNOWNC))
      call check(nf90_get_var(ncid,idSNOWNC,SNOWNC,                     &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))

print*,'debug: hail '
      call check(nf90_inq_varid(ncid,'HAILNC',idHAILNC))
      call check(nf90_get_var(ncid,idHAILNC,HAILNC,                     &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))

print*,'debug: groupel '
      call check(nf90_inq_varid(ncid,'GRAUPELNC',idGRAUPELNC))
      call check(nf90_get_var(ncid,idGRAUPELNC,GRAUPELNC,               &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))

print*,'debug: sst '
      call check(nf90_inq_varid(ncid,'SST',idSST))
      call check(nf90_get_var(ncid,idSST,SST,                           &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))

print*,'debug: snowh '
      call check(nf90_inq_varid(ncid,'SNOWH',idSNOWH))
      call check(nf90_get_var(ncid,idSNOWH,SNOWH,                           &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
print*,'debug: albed '
      if ( falbedo .eqv. .true. ) then
      call check(nf90_inq_varid(ncid,'ALBEDO',idalbedo))
      call check(nf90_get_var(ncid,idalbedo,ALBEDO,                     &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif

print*,'debug: znt '
      if ( fznt .eqv. .true. ) then
      call check(nf90_inq_varid(ncid,'ZNT',idZNT))
      call check(nf90_get_var(ncid,idZNT,ZNT,                           &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif
print*,'debug: sstsk '
      if ( fsstsk .eqv. .true. ) then
      call check(nf90_inq_varid(ncid,'SSTSK',idSSTSK))
      call check(nf90_get_var(ncid,idSSTSK,SSTSK,                       &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif
print*,'debug: hgt '
      call check(nf90_inq_varid(ncid,'HGT',idHGT))
      call check(nf90_get_var(ncid,idHGT,HGT,                           &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      if ( rotate .eqv. .true. ) then
print*,'debug: cosa '
          call check(nf90_inq_varid(ncid,'COSALPHA',idCOSALPHA))
          call check(nf90_get_var(ncid,idCOSALPHA,COSALPHA,             &
               start = (/ 1, 1, timestep /),                            &
               count = (/ west_east, south_north, 1 /) ))
print*,'debug: sina '
          call check(nf90_inq_varid(ncid,'SINALPHA',idSINALPHA))
          call check(nf90_get_var(ncid,idSINALPHA,SINALPHA,             &
               start = (/ 1, 1, timestep /),                            &
               count = (/ west_east, south_north, 1 /) ))
      endif

print*,'debug: t2 '
      call check(nf90_inq_varid(ncid,'T2',idT2))
      call check(nf90_get_var(ncid,idT2,T2,                             &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
print*,'debug: u10 '
      call check(nf90_inq_varid(ncid,'U10',idU10))
      call check(nf90_get_var(ncid,idU10,U10,                           &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
print*,'debug: v10 '
      call check(nf90_inq_varid(ncid,'V10',idV10))
      call check(nf90_get_var(ncid,idV10,V10,                           &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
print*,'debug: psfc '
      call check(nf90_inq_varid(ncid,'PSFC',idPSFC))
      call check(nf90_get_var(ncid,idPSFC,PSFC,                         &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      if ( fpbl .eqv. .true. ) then
print*,'debug: psfc '
      call check(nf90_inq_varid(ncid,'PBLH',idPBLH))
      call check(nf90_get_var(ncid,idPBLH,PBLH,                         &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif
      if ( ( fnrad .eqv. .true. ) .or. ( fgh .eqv. .true. ) ) then
print*,'debug: grdflx '
        call check(nf90_inq_varid(ncid,'GRDFLX',idgrdflx))
        call check(nf90_get_var(ncid,idgrdflx,GRDFLX,                     &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif
      if ( ( fnrad .eqv. .true. ) .or. ( fsh .eqv. .true. ) ) then
print*,'debug: hfx '
        call check(nf90_inq_varid(ncid,'HFX',idhfx))
        call check(nf90_get_var(ncid,idhfx,HFX,                           &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif
      if ( ( fnrad .eqv. .true. ) .or. ( flh .eqv. .true. ) ) then
print*,'debug: lh '
        call check(nf90_inq_varid(ncid,'LH',idlh))
        call check(nf90_get_var(ncid,idlh,LH,                             &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif

      if ( fnrad .eqv. .true. ) then
print*,'debug: tsk '
        call check(nf90_inq_varid(ncid,'TSK',idTSK))
        call check(nf90_get_var(ncid,idTSK,TSK,                       &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
print*,'debug: emiss '
        call check(nf90_inq_varid(ncid,'EMISS',idEMISS))
        call check(nf90_get_var(ncid,idEMISS,EMISS,                       &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif
      if ( ( ftgrad .eqv. .true. ) .or. ( flgrad .eqv. .true. ) ) then
print*,'debug: glw '
        call check(nf90_inq_varid(ncid,'GLW',idglw))
        call check(nf90_get_var(ncid,idglw,GLW,                           &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif
      if ( ( ftgrad .eqv. .true. ) .or. ( fsgrad .eqv. .true. ) ) then
print*,'debug: swdown '
        call check(nf90_inq_varid(ncid,'SWDOWN',idswdown))
        call check(nf90_get_var(ncid,idswdown,SWDOWN,                     &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif
      if ( fustar .eqv. .true. )  then
print*,'debug: ust '
        call check(nf90_inq_varid(ncid,'UST',idust))
        call check(nf90_get_var(ncid,idust,UST,                     &
           start = (/ 1, 1, timestep /),                                &
           count = (/ west_east, south_north, 1 /) ))
      endif
      if ( fhymet .eqv. .true. )  then
print*,'debug: qcloud '
      call check(nf90_inq_varid(ncid,'QCLOUD',idqcloud))
      call check(nf90_get_var(ncid,idqcloud, qcloud,                    &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top, 1 /) ))
print*,'debug: qice '
      call check(nf90_inq_varid(ncid,'QICE',idqice))
      call check(nf90_get_var(ncid,idqice, qice,                        &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top, 1 /) ))
print*,'debug: qrain '
      call check(nf90_inq_varid(ncid,'QRAIN',idqrain))   
      call check(nf90_get_var(ncid,idqrain, qrain,                      &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top, 1 /) ))
print*,'debug: qsnow '
      call check(nf90_inq_varid(ncid,'QSNOW',idqsnow))
      call check(nf90_get_var(ncid,idqsnow, qsnow,                      &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top, 1 /) ))
print*,'debug: qgraup '
      call check(nf90_inq_varid(ncid,'QGRAUP',idqgraup))
      call check(nf90_get_var(ncid,idqgraup, qgraup,                    &
           start = (/ 1, 1, 1, timestep /),                             &
           count = (/ west_east, south_north, bottom_top, 1 /) ))
      endif
! Close netcdf file

print*,'debug: CLOSE netcdf input '
      call check(nf90_close(ncid))
      end subroutine read_wrf
       
