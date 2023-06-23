       subroutine wrf_prec
!Compute total precipitation
!input RAINC convective precipitation             mm
!input RAINNC non convective precipitation        mm
!input RAINSH shallow convective precipitation    mm
!input GRAUPELNC non convective groupeln          mm
!input HAILNC non convective hail                 mm
!input SNOWNC non convective snow                 mm
!output totprec  totalprec                        mm 


!Massimo:
! per calcolare la precipitazione non occorrono le idrometeore
! da quanto ho capito nei vari forum di wrf, il totale e' in RAINC e
! RAINNC (al limite RAINSH che pero' e' in output solo
! se si attiva una particolare parametrizz della convez)

       use param_wrf

       IMPLICIT NONE
       integer :: i,j

!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(totprec,RAINC,RAINNC,RAINSH,south_north,west_east) &
!$OMP PRIVATE(i,j)
       do j=1,south_north
       do i=1,west_east
       totprec(i,j) = RAINC(i,j)+RAINNC(i,j)+RAINSH(i,j)
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

       end subroutine
