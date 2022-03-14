       subroutine wrf_rh
!Compute relative humidity
!input press total pressure hPa
!input insituT insitu temperature K
!input water vapor mixing ratio kg kg-1
!output relative humidity %

       use param_wrf

       IMPLICIT NONE
       real :: svp1,svp2,svp3,svpt0,ep_3
       real :: es,qvapors
       integer :: ii,jj,kk

       SVP1=0.6112
       SVP2=17.67
       SVP3=29.65
       SVPT0=273.15
       EP_3=0.622
     
!$OMP PARALLEL DO                                        &
!$OMP  COLLAPSE(3)                                       &
!$OMP DEFAULT(NONE)                                      &
!$OMP SHARED(rh,insituT,qvapor,press,bottom_top,south_north,west_east)                    &
!$OMP FIRSTPRIVATE(es,svp1,svp2,svpt0,svp3,qvapors,ep_3) &
!$OMP PRIVATE(ii,jj,kk)
       do kk=1,bottom_top
       do jj=1,south_north
       do ii=1,west_east
          es = 10.0*svp1*exp( svp2* (insituT(ii,jj,kk)-svpt0)/          &
     &                   (insituT(ii,jj,kk)-svp3 ) )
          qvapors = ep_3*es/(press(ii,jj,kk)-(1.0-ep_3)*es)
          rh(ii,jj,kk)=100.0*max(min(qvapor(ii,jj,kk)/qvapors,1.0),0.0)
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

       end subroutine
