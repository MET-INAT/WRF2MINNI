       subroutine wrf_tcc
!Compute total cloud cover 
!input press total pressure hPa
!input rh relative humidity %
!output tcc total cloud cover  [0-10]

       use param_wrf

       IMPLICIT NONE
      
       integer :: ii,jj,kk
       real    :: clfr_lo(west_east,south_north)
       real    :: clfr_mi(west_east,south_north)
       real    :: clfr_hi(west_east,south_north)
       real    :: pp

#ifdef debug
  print*,'debug: call compute tcc '
#endif

       
!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(tcc,clfr_lo,clfr_mi,clfr_hi,south_north,west_east) &
!$OMP PRIVATE(ii,jj)
       do jj=1,south_north
       do ii=1,west_east
       tcc(ii,jj)=0.0
       clfr_lo(ii,jj)=0.0
       clfr_mi(ii,jj)=0.0
       clfr_hi(ii,jj)=0.0
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER



!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(3)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(press,clfr_lo,rh,clfr_mi,clfr_hi,bottom_top,south_north,west_east ) &
!$OMP PRIVATE(ii,jj,kk)                        & 
!$OMP FIRSTPRIVATE(pp)
       do kk=1,bottom_top
       do jj=1,south_north
       do ii=1,west_east
          pp=press(ii,jj,kk)*100.0    ![hPa]->[Pa]
          if ( ( pp .ge. 80000.0 ) .and. ( pp .lt. 97000.0) ) then
             clfr_lo(ii,jj)=max(rh(ii,jj,kk),clfr_lo(ii,jj))
          elseif  ( ( pp .ge. 45000.0 ) .and. ( pp .lt. 80000.0) ) then
             clfr_mi(ii,jj)=max(rh(ii,jj,kk),clfr_mi(ii,jj))
          elseif  ( pp .gt. 45000.0 ) then
             clfr_hi(ii,jj)=max(rh(ii,jj,kk),clfr_hi(ii,jj))
          endif
       enddo
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER



!$OMP PARALLEL DO            &
!$OMP  COLLAPSE(2)           &
!$OMP DEFAULT(NONE)          &
!$OMP SHARED(clfr_lo,rh,clfr_mi,clfr_hi,tcc,south_north,west_east) &
!$OMP PRIVATE(ii,jj)
       do jj=1,south_north
       do ii=1,west_east
          clfr_lo(ii,jj)=4.0*clfr_lo(ii,jj)/100.0-3.0
          clfr_mi(ii,jj)=4.0*clfr_mi(ii,jj)/100.0-3.0
          clfr_hi(ii,jj)=2.5*clfr_hi(ii,jj)/100.0-1.5
          tcc(ii,jj)=max(clfr_lo(ii,jj),clfr_mi(ii,jj))         
          tcc(ii,jj)=max(tcc(ii,jj),clfr_hi(ii,jj))         
          tcc(ii,jj)=min(tcc(ii,jj),1.0)         
          tcc(ii,jj)=max(tcc(ii,jj),0.0)         
          tcc(ii,jj)=tcc(ii,jj)*10.0 !frac 0..10
          clfr_lo(ii,jj)=min(clfr_lo(ii,jj),1.0)
          clfr_lo(ii,jj)=max(clfr_lo(ii,jj),0.0)
          clfr_mi(ii,jj)=min(clfr_mi(ii,jj),1.0)
          clfr_mi(ii,jj)=max(clfr_mi(ii,jj),0.0)
          clfr_hi(ii,jj)=min(clfr_hi(ii,jj),1.0)
          clfr_hi(ii,jj)=max(clfr_hi(ii,jj),0.0)
       enddo
       enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

       end subroutine
