      subroutine interp2d_3dfields(nxi,nyi,nk,lonin,latin,fin, &
                                   nxo,nyo,lonout,latout,fout,patch)
     
      use param_minni, ONLY: XA,XB,XC,XD, YA,YB,YC,YD,WA,WB,WC,WD
      use general, ONLY : flag_interp
      implicit none

      integer :: nxo,nyo,nxi,nyi,i,j,nk,ii,jj,k
      integer :: sx,ex,sy,ey
      integer :: res(2),a(2),b(2),c(2),d(2)
      real :: lonin(nxi,nyi),latin(nxi,nyi),fin(nxi,nyi,nk)
      real :: lonout(nxo,nyo),latout(nxo,nyo),fout(nxo,nyo,nk)
      real :: x,y,dist,minimum,x1,y1
      real :: aox,aoy,box,boy,cox,coy,dox,doy,alpha,beta,t,s
      logical :: patch

      if ( patch ) then
         sx=2
         ex=nxo-1
         sy=2
         ey=nyo-1
      else
         sx=1
         ex=nxo
         sy=1
         ey=nyo
      endif
      if (flag_interp .eq. .true. ) then
!$OMP PARALLEL   DO    &
!$OMP COLLAPSE(2)      &
!$OMP DEFAULT (NONE)   &
!$OMP SHARED(fout,fin,latin,latout,lonin,lonout,nk,nxo,nyo,nxi,nyi,sx,ex,sy,ey)  &              
!$OMP SHARED(WA,WB,WC,WD,XA,XB,XC,XD,YA,YB,YC,YD)  &              
!$OMP FIRSTPRIVATE(x1,y1,x,y,dist,aox,aoy,box,boy,cox,coy,dox,doy,alpha,beta,t,s,res,a,b,c,d,minimum) &
!$OMP PRIVATE(i,j,ii,jj)
      do j=sy,ey
      do i=sx,ex
         minimum=9.0e10
         do jj=1,nyi 
         do ii=1,nxi 
         x=lonout(i,j)-lonin(ii,jj)
         y=latout(i,j)-latin(ii,jj)
         dist=x**2.0+y**2.0
         if (dist .lt. minimum ) then
            minimum=dist
            res(1)=ii
            res(2)=jj
            x1=x
            y1=y
         endif
         enddo
         enddo

         if (( x1 .ge. 0.0 ) .and.               &
     &       ( y1 .ge. 0.0 )) then 
!         c  x------------x d
!            |            |
!            | o          |
!    a=res-> x------------x b
!( (i,j) (i+1,j) (i,j+1) (i+1,j+1) )
             a=res
             b(1)=res(1)+1
             b(2)=res(2)
             c(1)=res(1)
             c(2)=res(2)+1
             d(1)=res(1)+1
             d(2)=res(2)+1
         else if (( x1 .le. 0.0 ) .and.           &
     &           (  y1 .le. 0.0 )) then 
!          c x------------x <-res=d
!            |          o |
!            |            |
!          a x------------x b
!( (i,j) (i-1,j) (i,j-1) (i-1,j-1) )
           a(1)=res(1)-1
           a(2)=res(2)-1
           b(1)=res(1)
           b(2)=res(2)-1
           c(1)=res(1)-1
           c(2)=res(2)
           d=res
         else if (( x1 .le. 0.0 ) .and.           &
     &           (  y1 .ge. 0.0 )) then 
!          c x------------x d
!            |            |
!            |          o |
!          a x------------x <-res=b
!( (i,j) (i-1,j) (i-1,j+1) (i,j+1)
         a(1)=res(1)-1
         a(2)=res(2)
         b=res
         c(1)=res(1)-1
         c(2)=res(2)+1
         d(1)=res(1)
         d(2)=res(2)+1
         else if (( x1 .ge. 0.0 ) .and.           &
     &           (  y1 .le. 0.0 )) then 
!    c=res-> x------------x d
!            | o          |
!            |            |
!         a  x------------x b
!( (i,j) (i,j-1) (i+1,j) (i+1,j-1) )
         a(1)=res(1)
         a(2)=res(2)-1
         b(1)=res(1)+1
         b(2)=res(2)-1
         c=res
         d(1)=res(1)+1
         d(2)=res(2)
         end if
!Now bilinear interpolation
        aox=abs(lonin(a(1),a(2))-lonout(i,j));
        aoy=abs(latin(a(1),a(2))-latout(i,j));
        box=abs(lonin(b(1),b(2))-lonout(i,j));
        boy=abs(latin(b(1),b(2))-latout(i,j));
        cox=abs(lonin(c(1),c(2))-lonout(i,j));
        coy=abs(latin(c(1),c(2))-latout(i,j));
        dox=abs(lonin(d(1),d(2))-lonout(i,j));
        doy=abs(latin(d(1),d(2))-latout(i,j));

        alpha=(aox+box+cox+dox)
        beta=(aoy+boy+coy+doy)
        s=(cox+aox)/alpha
        t=(coy+doy)/beta
        fout(i,j,:)=fin(a(1),a(2),:)*t*(1.0-s)+                         &
     &            fin(b(1),b(2),:)*t*s+                                 &
     &            fin(c(1),c(2),:)*(1.0-s)*(1.0-t)+                     &
                  fin(d(1),d(2),:)*(s)*(1.0-t)            
       XA(i,j)=a(1)
       YA(i,j)=a(2)
       WA(i,j)=t*(1.0-s)
       XB(i,j)=b(1)
       YB(i,j)=b(2)
       WB(i,j)=t*s
       XC(i,j)=c(1)
       YC(i,j)=c(2)
       WC(i,j)=(1.0-s)*(1.0-t)
       XD(i,j)=d(1)
       YD(i,j)=d(2)
       WD(i,j)=(s)*(1.0-t)
      enddo
      enddo
!$OMP END PARALLEL DO
!$OMP BARRIER
      flag_interp =.false.
      else

!$OMP PARALLEL DO      &
!$OMP  COLLAPSE(3)     &
!$OMP DEFAULT (NONE)   &
!$OMP SHARED(fout,fin,sy,sx,ey,ex,WA,WB,WC,WD,XA,XB,XC,XD,YA,YB,YC,YD) &
!$OMP SHARED(flag_interp,nk) &
!$OMP PRIVATE(i,j,k)
      do k=1,nk
      do j=sy,ey
      do i=sx,ex
        fout(i,j,k)=fin(XA(i,j),YA(i,j),k)*WA(i,j)+   &
                  fin(XB(i,j),YB(i,j),k)*WB(i,j)+   &
                  fin(XC(i,j),YC(i,j),k)*WC(i,j)+   &
                  fin(XD(i,j),YD(i,j),k)*WD(i,j)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
!$OMP BARRIER

      endif

      if ( patch) then
      fout(1,:,:)=fin(1,:,:)
      fout(nxo,:,:)=fin(nxi,:,:)
      fout(:,1,:)=fin(:,1,:)
      fout(:,nyo,:)=fin(:,nyi,:)
      endif
      end 
