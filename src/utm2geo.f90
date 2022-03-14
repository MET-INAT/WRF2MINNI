      SUBROUTINE utm2geo(num,cmeridian,fnorth,feast,lat,long,xnorth,yeast,dir)
!num : vector length
!cmeridian : 9 per UTM32
!fnorth : 0 north emisphere
!feast  : 500000 

      IMPLICIT none

!============================ INTERFACE ================================
!
! WRITTEN BY:        GEOTOOLS CORPORATION
!                    5808 BALCONES DRIVE, SUITE 202
!                    AUSTIN, TEXAS  78731
!
! DESCRIPTION: CONVERTS BETWEEN LAT/LONG AND X/Y.
!              LATITIUDE AND LONGITUDE ARE IN DECIMAL DEGREES.
!              XNORTH AND YEAST ARE IN METERS
!
!              TYPE=1 IS UTM, CLARKE 1866 SPHEROID.
!              DIR=0 MEANS CONVERT LATLONG TO UTM
!                 =1 MEANS CONVERT UTM TO LATLONG
!
!              BASLAT AND CMERIDIAN ARE THE BASE LATITUDE AND
!              CENTRAL MERIDIAN OF THE UTM ZONE CONTAINING THE 
!              POINT TO BE CONVERTED. THE INTERSECTION OF THE
!              BASE LATITUDE AND CENTRAL MERIDIAN IS THE ORIGIN
!              OF THE UTM X,Y COORDINATES WITHIN THE ZONE.
!              SEE USGS BULLETIN 1532, P.63-65.,AND USGS PP 1395,
!              P. 57-64.
!              
!
!              NOTE: FNORTH AND FEAST ARE PASSED IN, NORMALLY
!                    FNORTH IS 0. IN NORTHERN HEMISPHERE AND 10000 KM
!                    IN SOUTHERN HEMISPHERE. FEAST IS NORMALLY 500 KM.
!
!********************************************************************
!**********    NOTE: X AND Y IN THIS SUBROUTINE ARE REVERSED IN 
!**********          MEANING FROM THE UTM TRANSFORMATION EQUATIONS
!**********          GIVEN IN THE USGS PAPERS. THUS X IN THIS CODE
!**********          IS THE VERTICAL, OR LATITUDINAL, COORDINATE,
!**********          AND Y IS THE LONGITUDINAL COORDINATE.
!*********************************************************************
!
!   SEE: SNYDER,J.,1982,MAP PROJECTIONS USED BY THE USGS,
!        USGS BULLETIN 1532.
!        SNYDER,J.,1987,MAP PROJECTIONS - A WORKING MANUAL,
!        USGS PROFESSIONAL PAPER 1395.
!        NEWTON,G.D.,1985, COMPUTER PROGRAMS FOR COMMON MAP
!        PROJECTIONS, USGS BULLETIN 1642.

!============== PARAMETER AND LOCAL VARIABLE DECLARATIONS ==============
!
!     ---- PASS PARAMETERS
      INTEGER num,dir
      REAL cmeridian,fnorth,feast,xnorth(*),yeast(*),lat(*),long(*)
!
!     ---- LOCAL VARIABLES
      INTEGER i,j
      DOUBLE PRECISION latr,longr,blatr,cmeridr,dellong
      DOUBLE PRECISION n,t,c,a,m,m0,a2,a3,a4,a5,a6,t0
      DOUBLE PRECISION ars,ec,fa,fact,x3,x1
      DOUBLE PRECISION fp,fp0,sinfp,cosfp,term1,term2
      DOUBLE PRECISION p1,t1,t2,t4,t6,c1,c2,e2,e4,e6,e8
      DOUBLE PRECISION s1,s2,ecc,ee,eee,r1,an1
      DOUBLE PRECISION y1,yy,yyy,y2,y3,y4,y5,y6,y7,y8
      DOUBLE PRECISION x,y
      DOUBLE PRECISION ph1,ph2,ph3,ph4,f,ff,fff,el1,el2,el3,el4
!     
!     ---- DUMMY VARIABLES
      DOUBLE PRECISION phi,ec2,ec4,ec6,aa
!
!     ---- LOCAL CONSTANTS
      DOUBLE PRECISION esq,e4th,e6th,aradius,bradius,k0,ep2,degrad,cnvrt,pie
      REAL baslat
      PARAMETER (baslat=0.0)
!
!     ---- FUNCTIONS
!      LOGICAL REALNULL
      DOUBLE PRECISION trudist
!
!     ---- LOCAL STATEMENT FUNCTION TO COMPUTE M, THE TRUE DISTANCE
!          ALONG THE CENTRAL MERIDIAN FROM THE EQUATOR TO A GIVEN
!          LATITUDE PHI.  EQN. 3-21 IN (SNYDER, 1987, P.61). 
!
      TRUDIST(PHI,AA,EC2,EC4,EC6)=&
      AA*( (1.D0 -EC2/4.D0 -3.D0*EC4/64.D0 -5.D0*EC6/256.D0)*PHI&
         -(3.D0*EC2/8.D0 +3.D0*EC4/32.D0 +45.D0*EC6/1024.D0)&
          *DSIN(2.D0*PHI) &
         +(15.D0*EC4/256.D0 +45.D0*EC6/1024.D0)*DSIN(4.D0*PHI)&
         -(35.D0*EC6/3072.D0)*DSIN(6.D0*PHI)                    )
!
!
!     esq is the eccentricity, denoted as  "e squared" in the
!     transformation equations.  aradius and bradius are the equatorial
!     and polar radii, respectively. these are given in (snyder,1987),
!     pages 12 and 13, and are particular to the clarke 1866 ellipsoid.
!     to use a different ellipsoid with the utm projection, simply
!     assign the appropriate constant values to these three variables.
!     CLARKE 1866 ELLIPSOID      
!      esq=0.006768658d0
!      aradius=6.3782064d6
!      bradius=6.3565838d6

!     hayter ellipsoid
!      esq=0.00672267002233D0
!      aradius=6.378388D6
!      bradius=6.356912D6
!      feast=500000.0
!      fnorth=0.0

!        /*
!        do i=1,num
!         print*,long(i),lat(i)
!        end do
!        */

! WGS-84 ellipsoid
      esq=0.00669434D0
      aradius=6.378137D6
      bradius=6.356752D6
!      feast=500000.0
!      fnorth=10000000.0


      e4th=esq*esq
      e6th=e4th*esq
      ep2=esq/(1.0d0-esq)
      k0=9.996d-1
      pie=dacos(-1.0d0)
      degrad=pie/180.0d0
      cnvrt=180.0d0/pie

!
!=======================================================================
!
!        DIR=1
        if (dir.eq.0) then
!
!       ----convert lat-long to x-y
!
        cmeridr=cmeridian*degrad
        blatr=baslat*degrad

!       WRITE(*,*)'lAT AND lONG:', LAT(1), LONG(1)
!
        do 50 i=1,num
!
!          if (realnull(lat(i)) .or. realnull(long(i))) then
!
!           ----special handling for null value
!            xnorth(i)=rnull
!            yeast(i)=rnull
!
          if((lat(i).ge.90.0).or.(lat(i).le.-90.0)) then
!
!           ----special handling for either pole (snyder,1987,p.61)
            if (lat(i).ge.90.0) then
              latr=90.0*degrad
            else if (lat(i).le.-90.0) then
              latr=-90.0*degrad
            end if
            m=trudist(latr,aradius,esq,e4th,e6th)
            m0=0.        
            xnorth(i)=k0*(m-m0)+fnorth
            yeast(i)=0.+feast
!
          else
!
!           ----lat-long to x-y transformations, implementing eqns.
!               8-9,8-10,8-12 thru 8-15, and 4-20 in (snyder 87,p.61)
            latr=lat(i)*degrad
            longr=long(i)*degrad
            t0=dsin(latr)
            n=aradius/dsqrt(1.d0-(esq*t0*t0))
            t0=dtan(latr)
            t= t0*t0
            t0=dcos(latr)
            c=ep2*t0*t0
            dellong=longr-cmeridr
            a=dcos(latr)*dellong
            a2=a*a
            a4=a2*a2
            a6=a4*a2
            a3=a2*a
            a5=a3*a2
            m0=0.
            m =trudist(latr,aradius,esq,e4th,e6th) 
            y=k0*n*(a+(1.d0-t+c)*a3/6.d0+(5.d0-18.d0*t+t*t+ &
              72.d0*c-58.d0*ep2)*a5/120.d0)
            x=k0*(m -m0 +n*dtan(latr)*( a2/2.d0+(5.d0-t+9.d0*c+ &
              4.d0*c*c)*a4/24.d0 +(61.d0-58.d0*t+t*t+600.d0*c- &
              330.d0*ep2)*a6/720.d0)  )
!
!           ----add 10000 km in southern hemisphere to keep northing >0
            xnorth(i)=x+fnorth
!
!           ----add 500 km to keep easting positive
            yeast(i)=y+feast

!           WRITE(*,*)'nORTH:',FNORTH, X, XNORTH(I)
!           WRITE(*,*)'eAST:',FEAST,Y,YEAST(I)
          end if
 50     continue
!
      else if (dir.eq.1) then
!
!       ----convert x-y to lat-long
        aradius=aradius/1000.d0
        bradius=bradius/1000.d0
        ars=k0*aradius
        ec=1.d0-(bradius*bradius)/(aradius*aradius)
        fa=1.d0-ec
        fact=ec/fa
        do 250 i=1,num
!          if (realnull(xnorth(i)) .or. realnull(yeast(i))) then
!
!           ----special handling for null value
!            lat(i)=rnull
!            long(i)=rnull
!          else
!
!         ----calculation of footpoint latitude is solution of elliptic
!               integral expanded to 4th power in eccentricity
            x3=(xnorth(i)-fnorth)/1000.d0
            x1=dabs(x3)
            x1=x1/(ars*fa)
            fp=x1
            do 210 j=1,4
              fp0=fp
              sinfp=dsin(fp)
              cosfp=dcos(fp)
              term1=0.5d0*(fp-sinfp*cosfp)
              term2=0.75d0*term1-0.25d0*(sinfp**3.d0)*cosfp
              fp=x1-1.5d0*ec*term1-1.875d0*ec*ec*term2
  210       continue
!
            p1=fp
            t1=dsin(p1)/dcos(p1)
            t2=t1*t1
            t4=t2*t2
            t6=t4*t2
            c1=dcos(p1)
            c2=c1*c1
            e2=fact*c2
            e4=e2*e2
            e6=e4*e2
            e8=e6*e2
            s1=dsin(p1)
            s2=s1*s1
            ecc=ec*s2
            ee=1.d0-ecc
            ee=dsqrt(ee)
            eee=ee*ee*ee
            r1=(ars*fa)/eee
            an1=ars/ee
!
            y1=(yeast(i)-feast)/1000.d0
            yy=y1/r1
            yyy=y1/an1
            y2=yy*yyy
            y3=yyy*yyy*yyy
            y4=yy*y3
            y5=yyy*yyy*y3
            y6=yy*y5
            y7=yyy*yyy*y5
            y8=yy*y7
!
            ph1=-y2/2.d0
            ph2=y4/24.d0
            ph2=ph2*(5.d0+3.d0*t2+e2-4.d0*e4-9.d0*e2*t2)
            ph3=-y6/720.d0
            f=  61.d0+90.d0*t2+46.d0*e2+45.d0*t4-252.d0*t2*e2
            ff= -3.d0*e4+100.d0*e6-66.d0*t2*e4-90.d0*t4*e2
            fff=88.d0*e8+225.d0*t4*e4+84.d0*t2*e6-192.d0*t2*e8
            ph3=ph3*(f+ff+fff)
            ph4=y8/40320.d0
            ph4=ph4*(1385.d0+3633.d0*t2+4095.d0*t4+1575.d0*t6)
            lat(i)=cnvrt*(p1+t1*(ph1+ph2+ph3+ph4))
            if (x3 .lt. 0.d0) lat(i)=-lat(i)
!
            el1=yyy
            el2=-y3/6.d0
            el2=el2*(1.d0+2.d0*t2+e2)
            el3=y5/120.d0
            f=5.d0+6.d0*e2+28.d0*t2-3.d0*e4+8.d0*t2*e2
            ff=24.d0*t4-4.d0*e6+4.d0*t2*e4+24.d0*t2*e6
            el3=el3*(f+ff)
            el4=-y7/5040.d0
            el4=el4*(61.d0+662.d0*t2+1320.d0*t4+720.d0*t6)
            long(i)=cmeridian+cnvrt*(el1+el2+el3+el4)/c1
!          end if
  250   continue
      end if

      return
      end

