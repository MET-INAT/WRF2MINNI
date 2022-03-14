           subroutine wrf_interp_3d_z(nx,ny,nzi,V3D,Z,nzo,HEIGHT,V2D)
! V3D input field to be interpolated
! HEIGHT input height of the interpolated field
! V2D output interpolated field
      implicit none
      integer :: I,J,K
      integer :: nx,ny,nzi,nzo
      integer :: IP,IM,INTERP,KP
      real    :: V3D(nx,ny,nzi)
      real    :: V2D(nx,ny,nzo)
      real    :: HEIGHT(nzo)
      real    :: Z(nx,ny,nzi)
      real    :: W1,W2

! does vertical coordinate increase or decrease with increasing k?
! set offset appropriately

      IP = 0
      IM = 1
      if (Z(nx,ny,1) .gt. Z(nx,ny,nzi))  then
          IP = 1
          IM = 0
      endif 

!$OMP PARALLEL DO      &
!$OMP  COLLAPSE(3)     &
!$OMP DEFAULT (NONE)   &
!$OMP SHARED(IP,IM,nzo,ny,nx,nzi,HEIGHT,Z,V3D,V2D)  &
!$OMP FIRSTPRIVATE(INTERP,KP,W2,W1) &
!$OMP PRIVATE(I,J,K)
      do K = 1,nzo
      do J = 1,ny
      do I = 1,nx
         V2D(I,J,K) = 0.0
         INTERP = 0
         KP = nzi
         do  while ((INTERP .eq. 0) .and. (KP .ge. 2 ))
             if ((Z(I,J,KP-IM) .le. HEIGHT(K))  .and.                   &
     &            (Z(I,J,KP-IP) .gt. HEIGHT(K))) then
                   W2 = (HEIGHT(K)-Z(I,J,KP-IM))/                       &
     &                  (Z(I,J,KP-IP)-Z(I,J,KP-IM))
                   W1 = 1.0 - W2
                   V2D(I,J,K) = W1*V3D(I,J,KP-IM) + W2*V3D(I,J,KP-IP)
                   INTERP = 1
                  endif
             KP = KP - 1
         enddo
      enddo 
      enddo 
      enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

      end
