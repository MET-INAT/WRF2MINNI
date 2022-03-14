!----------------------------------------------------------------------
! MGM, 26.01.2007 - CVS check
!----------------------------------------------------------------------
!     GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE - VERSION 2.0.2      
!     FORTRAN 77 LANGUAGE FOR IBM, AMDAHL, ENCORE, VAX, CONCURRENT, AND
!     DATA GENERAL COMPUTERS                                           
!                   ADJLZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION ADJLZ0 (LON)                           
!                                                                      
! FUNCTION TO ADJUST LONGITUDE ANGLE TO MODULO 180 DEGREES.            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      DATA TWO,PI /2.0D0,3.14159265358979323846D0/                     
!                                                                      
  020 ADJLZ0 = LON                                                     
      IF (DABS(LON) .LE. PI) RETURN                                    
      TWOPI = TWO * PI                                                 
      LON = LON - DSIGN (TWOPI,LON)                                    
      GO TO 020                                                        
!                                                                      
      END                                                              
!                   ASINZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION ASINZ0 (CON)                           
!                                                                      
! THIS FUNCTION ADJUSTS FOR ROUND-OFF ERRORS IN COMPUTING ARCSINE      
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      DATA ONE /1.0D0/                                                 
!                                                                      
      IF (DABS(CON) .GT. ONE) THEN                                     
       CON = DSIGN (ONE,CON)                                           
       ENDIF                                                           
      ASINZ0 = DASIN (CON)                                             
      RETURN                                                           
!                                                                      
      END                                                              
!                   DMSPZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION DMSPZ0 (SGNA,DEGS,MINS,SECS)           
!                                                                      
! SUBROUTINE TO CONVERT UNPACKED DMS TO PACKED DMS ANGLE               
! SGNA : SIGN OF ANGLE                                                 
! DEGS : DEGREES PORTION OF ANGLE                                      
! MINS : MINUTES PORTION OF ANGLE                                      
! SECS : SECONDS PORTION OF ANGLE                                      
!                                                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      REAL*4 SECS                                                      
      INTEGER*4 DEGS,MINS                                              
      CHARACTER*1 SGNA,NEG                                             
      DATA CON1,CON2 /1000000.0D0,1000.0D0/                            
      DATA NEG /'-'/                                                   
!                                                                      
      CON = DBLE (DEGS) * CON1 + DBLE (MINS) * CON2 + DBLE (SECS)      
      IF (SGNA .EQ. NEG) CON = - CON                                   
      DMSPZ0 = CON                                                     
      RETURN                                                           
!                                                                      
      END                                                              
!                   E0FNZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION E0FNZ0 (ECCNTS)                        
!                                                                      
! FUNCTION TO COMPUTE CONSTANT (E0).                                   
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      DATA QUART,ONE,ONEQ,THREE,SIXT /0.25D0,1.0D0,1.25D0,3.0D0,16.0D0/
!                                                                      
      E0FNZ0 = ONE - QUART * ECCNTS * (ONE + ECCNTS / SIXT *            &
     &         (THREE + ONEQ * ECCNTS))                                
!                                                                      
      RETURN                                                           
      END                                                              
!                   E1FNZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION E1FNZ0 (ECCNTS)                        
!                                                                      
! FUNCTION TO COMPUTE CONSTANT (E1).                                   
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      DATA CON1,CON2,CON3 /0.375D0,0.25D0,0.46875D0/                   
      DATA ONE /1.0D0/                                                 
!                                                                      
      E1FNZ0 = CON1 * ECCNTS * (ONE + CON2 * ECCNTS *                   &
     &         (ONE + CON3 * ECCNTS))                                  
!                                                                      
      RETURN                                                           
      END                                                              
!                   E2FNZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION E2FNZ0 (ECCNTS)                        
!                                                                      
! FUNCTION TO COMPUTE CONSTANT (E2).                                   
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      DATA CON1,CON2 /0.05859375D0,0.75D0/                             
      DATA ONE /1.0D0/                                                 
!                                                                      
      E2FNZ0 = CON1 * ECCNTS * ECCNTS * (ONE + CON2 * ECCNTS)          
!                                                                      
      RETURN                                                           
      END                                                              
!                   E3FNZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION E3FNZ0 (ECCNTS)                        
!                                                                      
! FUNCTION TO COMPUTE CONSTANT (E3).                                   
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
!                                                                      
      E3FNZ0 = ECCNTS*ECCNTS*ECCNTS*(35.D0/3072.D0)                    
!                                                                      
      RETURN                                                           
      END                                                              
!                   E4FNZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION E4FNZ0 (ECCENT)                        
!                                                                      
! FUNCTION TO COMPUTE CONSTANT (E4).                                   
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      DATA ONE /1.0D0/                                                 
!                                                                      
      CON = ONE + ECCENT                                               
      COM = ONE - ECCENT                                               
      E4FNZ0 = DSQRT ((CON ** CON) * (COM ** COM))                     
!                                                                      
      RETURN                                                           
      END                                                              
!                   GTPZ0                                              
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      SUBROUTINE GTPZ0(CRDIN,INSYS,INZONE,TPARIN,INUNIT,INSPH,IPR,JPR,  &
     &     LEMSG,LPARM,CRDIO,IOSYS,IOZONE,TPARIO,IOUNIT,LN27,LN83,      &
     &     FN27,FN83,LENGTH,IFLG)                                      
!                                                                      
! *********************************************************************
!                                                                      
! INPUT ***************************************************************
! CRDIN  : COORDINATES IN INPUT SYSTEM (2 DP WORDS ARRAY).             
! INSYS  : CODE NUMBER OF INPUT COORDINATE SYSTEM (INTEGER).           
!            =  0 , GEOGRAPHIC                                         
!            =  1 , U T M                                              
!            =  2 , STATE PLANE                                        
!            =  3 , ALBERS CONICAL EQUAL-AREA                          
!            =  4 , LAMBERT CONFORMAL CONIC                            
!            =  5 , MERCATOR                                           
!            =  6 , POLAR STEREOGRAPHIC                                
!            =  7 , POLYCONIC                                          
!            =  8 , EQUIDISTANT CONIC                                  
!            =  9 , TRANSVERSE MERCATOR                                
!            = 10 , STEREOGRAPHIC                                      
!            = 11 , LAMBERT AZIMUTHAL EQUAL-AREA                       
!            = 12 , AZIMUTHAL EQUIDISTANT                              
!            = 13 , GNOMONIC                                           
!            = 14 , ORTHOGRAPHIC                                       
!            = 15 , GENERAL VERTICAL NEAR-SIDE PERSPECTIVE             
!            = 16 , SINUSOIDAL                                         
!            = 17 , EQUIRECTANGULAR (PLATE CARREE)                     
!            = 18 , MILLER CYLINDRICAL                                 
!            = 19 , VAN DER GRINTEN I                                  
!            = 20 , OBLIQUE MERCATOR (HOTINE)                          
!            = 21 , ROBINSON                                           
!            = 22 , SPACE OBLIQUE MERCATOR                             
!            = 23 , MODIFIED-STEREOGRAPHIC CONFORMAL (ALASKA)          
! INZONE : CODE NUMBER OF INPUT COORDINATE ZONE (INTEGER).             
! TPARIN : PARAMETERS OF INPUT REFERENCE SYSTEM (15 DP WORDS ARRAY).   
! INUNIT : CODE NUMBER OF UNITS OF MEASURE FOR INPUT COORDINATES (I*4) 
!            = 0 , RADIANS.                                            
!            = 1 , U.S. FEET.                                          
!            = 2 , METERS.                                             
!            = 3 , SECONDS OF ARC.                                     
!            = 4 , DEGREES OF ARC.                                     
!            = 5 , INTERNATIONAL FEET.                                 
!            = 6 , USE LEGISLATED DISTANCE UNITS FROM NADUT TABLE      
! INSPH  : INPUT SPHEROID CODE.  SEE SPHDZ0 FOR PROPER CODES.          
! IPR    : PRINTOUT FLAG FOR ERROR MESSAGES. 0=YES, 1=NO               
! JPR    : PRINTOUT FLAG FOR PROJECTION PARAMETERS 0=YES, 1=NO         
! LEMSG  : LOGICAL UNIT FOR LISTING ERROR MESSAGES IF IPR = 0          
! LPARM  : LOGICAL UNIT FOR LISTING PROJECTION PARAMETERS IF JPR = 0   
! LN27   : LOGICAL UNIT FOR NAD 1927 SPCS PARAMETER FILE               
! FN27   : FILE NAME OF NAD 1927 SPCS PARAMETERS                       
! LN83   : LOGICAL UNIT FOR NAD 1983 SPCS PARAMETER FILE               
! FN83   : FILE NAME OF NAD 1983 SPCS PARAMETERS                       
! LENGTH : RECORD LENGTH OF NAD1927 AND NAD1983 PARAMETER FILES        
! OUTPUT ***                                                       ****
! IOSYS  : CODE NUMBER OF OUTPUT COORDINATE SYSTEM (INTEGER).          
! IOZONE : CODE NUMBER OF OUTPUT COORDINATE ZONE (INTEGER).            
! TPARIO : PARAMETERS OF OUTPUT REFERENCE SYSTEM (15 DP WORDS ARRAY).  
! IOUNIT : CODE NUMBER OF UNITS OF MEASURE FOR OUTPUT COORDINATES (I*4)
! CRDIO  : COORDINATES IN OUTPUT REFERENCE SYSTEM (2 DP WORDS ARRAY).  
! IFLG   : RETURN FLAG (INTEGER).                                      
!            = 0 , SUCCESSFUL TRANSFORMATION.                          
!            = 1 , ILLEGAL INPUT SYSTEM CODE.                          
!            = 2 , ILLEGAL OUTPUT SYSTEM CODE.                         
!            = 3 , ILLEGAL INPUT UNIT CODE.                            
!            = 4 , ILLEGAL OUTPUT UNIT CODE.                           
!            = 5 , INCONSISTENT UNIT AND SYSTEM CODES FOR INPUT.       
!            = 6 , INCONSISTENT UNIT AND SYSTEM CODES FOR OUTPUT.      
!            = 7 , ILLEGAL INPUT ZONE CODE.                            
!            = 8 , ILLEGAL OUTPUT ZONE CODE.                           
!      OTHERWISE , ERROR CODE FROM PROJECTION COMPUTATIONAL MODULE.    
!                                                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      INTEGER*4 NAD27(134), NAD83(134), NADUT(54), SPTYPE(134)         
      INTEGER*4 SYSUNT(24), SWITCH(23), ITER                           
      INTEGER*2 INMOD, IOMOD, FWD, INV                                 
      CHARACTER*128 FN27, FN83, FILE27, FILE83                         
      DIMENSION CRDIN(2),CRDIO(2),TPARIN(15),TPARIO(15),COORD(2)       
      DIMENSION DUMMY(15), PDIN(15), PDIO(15)                          
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z,E4Z                    
      COMMON /PROJZ0/ IPROJ                                            
      COMMON /SPCS/ ISPHER,LU27,LU83,LEN,MSYS,FILE27,FILE83            
      COMMON /TOGGLE/ SWITCH                                           
!                                                                      
      PARAMETER (MAXUNT=6, MAXSYS=23)                                  
      PARAMETER (FWD=0, INV=1)                                         
      DATA SYSUNT / 0 , 23*2 /                                         
      DATA PDIN/15*0.0D0/, PDIO/15*0.0D0/                              
      DATA INSP/999/, INPJ/999/, INZN/99999/                           
      DATA IOSP/999/, IOPJ/999/, IOZN/99999/                           
      DATA ITER /0/                                                    
      DATA JFLAG/0/                                                    
!                                                                      
      DATA NAD27/0101,0102,5010,5300,0201,0202,0203,0301,0302,0401,0402, &
     &           0403,0404,0405,0406,0407,0501,0502,0503,0600,0700,0901, &
     &           0902,0903,1001,1002,5101,5102,5103,5104,5105,1101,1102, &
     &           1103,1201,1202,1301,1302,1401,1402,1501,1502,1601,1602, &
     &           1701,1702,1703,1801,1802,1900,2001,2002,2101,2102,2103, &
     &           2111,2112,2113,2201,2202,2203,2301,2302,2401,2402,2403, &
     &           2501,2502,2503,2601,2602,2701,2702,2703,2800,2900,3001, &
     &           3002,3003,3101,3102,3103,3104,3200,3301,3302,3401,3402, &
     &           3501,3502,3601,3602,3701,3702,3800,3901,3902,4001,4002, &
     &           4100,4201,4202,4203,4204,4205,4301,4302,4303,4400,4501, &
     &           4502,4601,4602,4701,4702,4801,4802,4803,4901,4902,4903, &
     &           4904,5001,5002,5003,5004,5005,5006,5007,5008,5009,5201, &
     &           5202,5400/                                             
!                                                                      
      DATA NAD83/0101,0102,5010,5300,0201,0202,0203,0301,0302,0401,0402, &
     &           0403,0404,0405,0406,0000,0501,0502,0503,0600,0700,0901, &
     &           0902,0903,1001,1002,5101,5102,5103,5104,5105,1101,1102, &
     &           1103,1201,1202,1301,1302,1401,1402,1501,1502,1601,1602, &
     &           1701,1702,1703,1801,1802,1900,2001,2002,2101,2102,2103, &
     &           2111,2112,2113,2201,2202,2203,2301,2302,2401,2402,2403, &
     &           2500,0000,0000,2600,0000,2701,2702,2703,2800,2900,3001, &
     &           3002,3003,3101,3102,3103,3104,3200,3301,3302,3401,3402, &
     &           3501,3502,3601,3602,3701,3702,3800,3900,0000,4001,4002, &
     &           4100,4201,4202,4203,4204,4205,4301,4302,4303,4400,4501, &
     &           4502,4601,4602,4701,4702,4801,4802,4803,4901,4902,4903, &
     &           4904,5001,5002,5003,5004,5005,5006,5007,5008,5009,5200, &
     &           0000,5400/                                            
!                                                                      
!     TABLE OF UNIT CODES AS SPECIFIED BY STATE LAWS AS OF 2/1/92      
!     FOR NAD 1983 SPCS - 1 = U.S. SURVEY FEET, 2 = METERS,            
!                         5 = INTERNATIONAL FEET                       
!                                                                      
!     NADUT - UNIT CODES FOR THE STATES ARRANGED IN STATE NUMBER ORDER 
!              (FIRST TWO DIGITS OF ZONE NUMBER)                       
!                                                                      
      DATA NADUT /1, 5, 1, 1, 5, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2, 2, &  
     &            1, 1, 5, 2, 1, 2, 5, 1, 2, 2, 2, 1, 1, 1, 5, 2, 1, 5, &
     &            2, 2, 5, 2, 1, 1, 5, 2, 2, 1, 2, 1, 2, 2, 1, 2, 2, 2/ 
!                                                                      
!     TABLE OF STATE PLANE ZONE TYPES:  4 = LAMBERT, 7 = POLYCONIC,    
!     9 = TRANSVERSE MERCATOR, AND 20 = OBLIQUE MERCATOR               
!                                                                      
      DATA SPTYPE / 9, 9, 4, 4, 9, 9, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,   & 
     &              4, 4, 4, 9, 9, 9, 4, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, &
     &              9, 9, 9, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 9, 9, 4, 4, &
     &              4, 9, 9, 9, 4, 4, 4, 4, 4, 4, 9, 9, 9, 9, 9, 4, 4, &
     &              4, 4, 4, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 4, 4, 4, &
     &              4, 4, 4, 4, 4, 4, 4, 4, 4, 9, 4, 4, 4, 4, 4, 4, 4, &
     &              4, 4, 4, 4, 4, 4, 9, 4, 4, 4, 4, 4, 4, 4, 4, 4, 9, &
     &              9, 9, 9,20, 9, 9, 9, 9, 9, 9, 9, 9, 4, 4, 7/       
!                                                                      
!     SETUP                                                            
!                                                                      
      IOSPH = INSPH                                                    
      IPEMSG = IPR                                                     
      IPPARM = JPR                                                     
      IPELUN = LEMSG                                                   
      IPPLUN = LPARM                                                   
      IPROJ = INSYS                                                    
      LU27 = LN27                                                      
      FILE27 = FN27                                                    
      LU83 = LN83                                                      
      FILE83 = FN83                                                    
      LEN = LENGTH                                                     
!                                                                      
!     INITIALIZE SWITCH FOR EACH PROJECTION TO ZERO                    
!                                                                      
      ITER = ITER + 1                                                  
      IF (ITER .LE. 1) THEN                                            
         DO 5 I=1,15                                                   
            DUMMY(I) = 0.0D0                                           
    5    CONTINUE                                                      
         MSYS = 2                                                      
      END IF                                                           
      INSPCS = 2                                                       
      IOSPCS = 2                                                       
      IF (JFLAG.NE.0) GO TO 10                                         
      EZ = 0.0D0                                                       
      ESZ = 0.0D0                                                      
      CALL SPHDZ0(0,DUMMY)                                             
      JFLAG = 1                                                        
!                                                                      
! CHECK VALIDITY OF CODES FOR REFERENCE SYSTEMS.                       
!                                                                      
   10 IF (INSYS.LT.0 .OR. INSYS.GT.MAXSYS) THEN                        
         IF (IPEMSG .NE. 0) WRITE (IPELUN,2000) INSYS                  
 2000    FORMAT (' ILLEGAL SOURCE REFERENCE SYSTEM CODE = ',I6)        
         IFLG = 1                                                      
         RETURN                                                        
      END IF                                                           
!                                                                      
      IF (IOSYS.LT.0 .OR. IOSYS.GT.MAXSYS) THEN                        
         IF (IPEMSG .NE. 0) WRITE (IPELUN,2010) IOSYS                  
 2010    FORMAT (' ILLEGAL TARGET REFERENCE SYSTEM CODE = ',I6)        
         IFLG = 2                                                      
         RETURN                                                        
      END IF                                                           
!                                                                      
!     FORCE INITIALIZATION OF PROJECTIONS IF SPHEROID OR PROJECTION    
!     HAS CHANGED FROM PREVIOUS INPUT - OUTPUT SET                     
!                                                                      
      IF (INSPH .NE. INSP) THEN                                        
         DO 11 I = 1,MAXSYS                                            
            SWITCH(I) = 0                                              
   11    CONTINUE                                                      
      END IF                                                           
!                                                                      
      IF (INSYS .GT. 0) THEN                                           
         IF (INSYS .NE. INPJ .AND. INSYS .NE. IOPJ) SWITCH(INSYS) = 0  
         IF (SWITCH(INSYS) .NE. INZONE .AND. SWITCH(INSYS) .NE. IOZONE) & 
     &      SWITCH(INSYS) = 0                                          
      END IF                                                           
!                                                                      
      IF (IOSYS .GT. 0) THEN                                           
         IF (IOSYS .NE. INPJ .AND. IOSYS .NE. IOPJ) SWITCH(IOSYS) = 0  
         IF (SWITCH(IOSYS) .NE. INZONE .AND. SWITCH(IOSYS) .NE. IOZONE) & 
     &      SWITCH(IOSYS) = 0                                          
      END IF                                                           
!                                                                      
!     CHECK FOR REPEAT OF INPUT SYSTEM                                 
!                                                                      
      INMOD = 1                                                        
      IF (INSYS .EQ. 2) THEN                                           
         IF (INZONE .GT. 0) THEN                                       
            ID = 0                                                     
            IF (INSPH .EQ. 0) THEN                                     
               DO 12 I = 1,134                                         
                  IF (INZONE .EQ. NAD27(I)) ID = I                     
   12          CONTINUE                                                
            END IF                                                     
            IF (INSPH .EQ. 8) THEN                                     
               DO 13 I = 1,134                                         
                  IF (INZONE .EQ. NAD83(I)) ID = I                     
   13          CONTINUE                                                
            END IF                                                     
            IF (ID .NE. 0) INSPCS = SPTYPE(ID)                         
            IF (INZONE .NE. SWITCH(INSPCS)) GO TO 15                   
         END IF                                                        
      END IF                                                           
      IF (INSP .NE. INSPH) GO TO 15                                    
      IF (INPJ .NE. INSYS) GO TO 15                                    
      IF (INZN .NE. INZONE) GO TO 15                                   
      IF (INSYS .GE. 3) THEN                                           
         DO 14 I=1,15                                                  
            IF (TPARIN(I) .NE. PDIN(I)) GO TO 15                       
   14    CONTINUE                                                      
      END IF                                                           
      INMOD = 0                                                        
      GO TO 30                                                         
!                                                                      
!     SAVE INPUT SYSTEM PARAMETERS                                     
!                                                                      
   15 INSP = INSPH                                                     
      INPJ = INSYS                                                     
      INZN = INZONE                                                    
      DO 16 I=1,15                                                     
   16 PDIN(I) = TPARIN(I)                                              
!                                                                      
! CHECK CONSISTENCY BETWEEN UNITS OF MEASURE                           
!                                                                      
      IF (INUNIT.LT.0 .OR. INUNIT.GT.MAXUNT) THEN                      
         IF (IPEMSG .NE. 0) WRITE (IPELUN,2020) INUNIT                 
 2020    FORMAT (' ILLEGAL SOURCE UNIT CODE = ',I6)                    
         IFLG = 3                                                      
         RETURN                                                        
      END IF                                                           
!                                                                      
!     CHECK FOR REPEAT OF OUTPUT SYSTEM                                
!                                                                      
   30 IOMOD = 1                                                        
      IF (IOSYS .EQ. 2) THEN                                           
         IF (IOZONE .GT. 0) THEN                                       
            ID = 0                                                     
            IF (IOSPH .EQ. 0) THEN                                     
               DO 32 I = 1,134                                         
                  IF (IOZONE .EQ. NAD27(I)) ID = I                     
   32          CONTINUE                                                
            END IF                                                     
            IF (IOSPH .EQ. 8) THEN                                     
               DO 33 I = 1,134                                         
                  IF (IOZONE .EQ. NAD83(I)) ID = I                     
   33          CONTINUE                                                
            END IF                                                     
            IF (ID .NE. 0) IOSPCS = SPTYPE(ID)                         
            IF (IOZONE .NE. SWITCH(INSPCS)) GO TO 35                   
         END IF                                                        
      END IF                                                           
      IF (IOSP .NE. INSPH) GO TO 35                                    
      IF (IOSP .NE. IOSPH) GO TO 35                                    
      IF (IOPJ .NE. IOSYS) GO TO 35                                    
      IF (IOZN .NE. IOZONE) GO TO 35                                   
      IF (IOSYS .GE. 3) THEN                                           
         DO 34 I=1,15                                                  
            IF (TPARIO(I) .NE. PDIO(I)) GO TO 35                       
   34    CONTINUE                                                      
      END IF                                                           
      IOMOD = 0                                                        
      GO TO 80                                                         
!                                                                      
!     SAVE OUTPUT SYSTEM PARAMETERS                                    
!                                                                      
   35 IOSP = INSPH                                                     
      IOPJ = IOSYS                                                     
      IOZN = IOZONE                                                    
      DO 36 I=1,15                                                     
   36 PDIO(I) = TPARIO(I)                                              
!                                                                      
! CHECK CONSISTENCY BETWEEN UNITS OF MEASURE                           
!                                                                      
      IF (IOUNIT.LT.0 .OR. IOUNIT.GT.MAXUNT) THEN                      
         IF (IPEMSG .NE. 0) WRITE (IPELUN,2030) IOUNIT                 
 2030    FORMAT (' ILLEGAL TARGET UNIT CODE = ',I6)                    
         IFLG = 4                                                      
         RETURN                                                        
      END IF                                                           
!                                                                      
   80 IUNIT = SYSUNT(INSYS + 1)                                        
!                                                                      
!     CHANGE UNITS TO LEGISLATED UNITS USING TABLE                     
!                                                                      
      IF (INSPH .EQ. 0 .AND. INSYS .EQ. 2 .AND. INUNIT .EQ. 6) INUNIT=1
      IF (INSPH .EQ. 8 .AND. INSYS .EQ. 2 .AND. INUNIT .EQ. 6) THEN    
         IND = 0                                                       
         DO 90 I = 1,134                                               
            IF (INZONE .EQ. NAD83(I)) IND = I                          
   90    CONTINUE                                                      
         IF (IND .NE. 0) INUNIT = NADUT( INT(INZONE/100))              
      END IF                                                           
      CALL UNTFZ0 (INUNIT,IUNIT,FACTOR,IFLG)                           
      IF (IFLG .EQ. 0) GO TO 100                                       
      IFLG = 5                                                         
      RETURN                                                           
  100 COORD(1) = FACTOR * CRDIN(1)                                     
      COORD(2) = FACTOR * CRDIN(2)                                     
      IUNIT = SYSUNT(IOSYS + 1)                                        
!                                                                      
!     CHANGE UNITS TO LEGISLATED UNITS USING TABLE                     
!                                                                      
      IF (INSPH .EQ. 0 .AND. IOSYS .EQ. 2 .AND. IOUNIT .EQ. 6) IOUNIT=1
      IF (INSPH .EQ. 8 .AND. IOSYS .EQ. 2 .AND. IOUNIT .EQ. 6) THEN    
         IND = 0                                                       
         DO 110 I = 1,134                                              
            IF (IOZONE .EQ. NAD83(I)) IND = I                          
  110    CONTINUE                                                      
         IF (IND .NE. 0) IOUNIT = NADUT( INT(IOZONE/100))              
      END IF                                                           
      CALL UNTFZ0 (IUNIT,IOUNIT,FACTOR,IFLG)                           
      IF (IFLG .EQ. 0) GO TO 120                                       
      IFLG = 6                                                         
      RETURN                                                           
  120 IF (INSYS.NE.IOSYS.OR.INZONE.NE.IOZONE.OR.INZONE.LE.0) GO TO 140 
      CRDIO(1) = FACTOR * COORD(1)                                     
      CRDIO(2) = FACTOR * COORD(2)                                     
      RETURN                                                           
!                                                                      
! COMPUTE TRANSFORMED COORDINATES AND ADJUST THEIR UNITS.              
!                                                                      
  140 IF (INSYS .EQ. 0) GO TO 520                                      
      IF (INZONE.GT.60 .OR. INSYS.EQ.1) GO TO 200                      
      IF (IPEMSG .NE. 0) WRITE (IPELUN,2040) INZONE                    
 2040 FORMAT (' ILLEGAL SOURCE ZONE NUMBER = ',I6)                     
      IFLG = 7                                                         
      RETURN                                                           
!                                                                      
! INVERSE TRANSFORMATION.                                              
!                                                                      
  200 IPROJ=INSYS                                                      
      ISPHER = INSPH                                                   
      IF (INSYS.GE.3) CALL SPHDZ0(INSPH,TPARIN)                        
!                                                                      
!     CHECK FOR CHANGE IN ZONE FROM LAST USE OF THE INPUT PROJECTION   
!                                                                      
      IF (INSYS .EQ. 1 .AND. INZONE .NE. SWITCH(9)) THEN               
         SWITCH(1) = 0                                                 
         INMOD = 1                                                     
      END IF                                                           
      IF (INSYS .EQ. 2 .AND. INZONE .NE. SWITCH(INSPCS)) THEN          
         SWITCH(2) = 0                                                 
         INMOD = 1                                                     
      END IF                                                           
      IF (INZONE .NE. SWITCH(INSYS)) THEN                              
         SWITCH(INSYS) = 0                                             
         INMOD = 1                                                     
      END IF                                                           
!                                                                      
      IF (INSYS .EQ. 1) THEN                                           
         IF (INZONE.EQ.0.AND.TPARIN(1).NE.0.0D0) GO TO 211             
         TPARIN(1) = 1.0D6*DBLE(6*INZONE-183)                          
         TPARIN(2) = DSIGN(4.0D7,DBLE(INZONE))                         
  211    CALL SPHDZ0(INSPH,DUMMY)                                      
         TPARIN(14) = DUMMY(1)                                         
         TPARIN(15) = DUMMY(2)                                         
         IF (INMOD .NE. 0) THEN                                        
            CALL PJINIT (INSYS,INZONE,TPARIN)                          
            IF (IERROR .NE. 0) INZN = 99999                            
            IF (IERROR .NE. 0) GO TO 500                               
         END IF                                                        
         CALL PJ01Z0 (COORD,CRDIO,INV)                                 
      END IF                                                           
!                                                                      
      IF (INSYS .GT. 1) THEN                                           
         IF (INMOD .NE. 0) THEN                                        
            MSYS = INSPCS                                              
            CALL PJINIT (INSYS,INZONE,TPARIN)                          
            IF (IERROR .NE. 0) INZN = 99999                            
            IF (IERROR .NE. 0) GO TO 500                               
         END IF                                                        
         IF (INSYS .EQ. 2) CALL PJ02Z0 (COORD,CRDIO,INV)               
         IF (INSYS .EQ. 3) CALL PJ03Z0 (COORD,CRDIO,INV)               
         IF (INSYS .EQ. 4) CALL PJ04Z0 (COORD,CRDIO,INV)               
         IF (INSYS .EQ. 5) CALL PJ05Z0 (COORD,CRDIO,INV)               
         IF (INSYS .EQ. 6) CALL PJ06Z0 (COORD,CRDIO,INV)               
         IF (INSYS .EQ. 7) CALL PJ07Z0 (COORD,CRDIO,INV)               
         IF (INSYS .EQ. 8) CALL PJ08Z0 (COORD,CRDIO,INV)               
         IF (INSYS .EQ. 9) CALL PJ09Z0 (COORD,CRDIO,INV)               
         IF (INSYS .EQ. 10) CALL PJ10Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 11) CALL PJ11Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 12) CALL PJ12Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 13) CALL PJ13Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 14) CALL PJ14Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 15) CALL PJ15Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 16) CALL PJ16Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 17) CALL PJ17Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 18) CALL PJ18Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 19) CALL PJ19Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 20) CALL PJ20Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 21) CALL PJ21Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 22) CALL PJ22Z0 (COORD,CRDIO,INV)              
         IF (INSYS .EQ. 23) CALL PJ23Z0 (COORD,CRDIO,INV)              
      END IF                                                           
!                                                                      
  500 IFLG = IERROR                                                    
      DO 510 I = 1,15                                                  
  510 TPARIN(I) = PDIN(I)                                              
      IF (IFLG .NE. 0) RETURN                                          
      CRDIO(1) = ADJLZ0(CRDIO(1))                                      
      IF (IOSYS .EQ. 0) GO TO 920                                      
      COORD(1) = CRDIO(1)                                              
      COORD(2) = CRDIO(2)                                              
  520 IF (INSYS .EQ. 0 .AND. IOSYS .EQ. 0) THEN                        
         CRDIO(1) = COORD(1)                                           
         CRDIO(2) = COORD(2)                                           
         GO TO 920                                                     
      END IF                                                           
      IF (IOZONE.GT.60 .OR. IOSYS.EQ.1) GO TO 540                      
      IF (IPEMSG .NE. 0) WRITE (IPELUN,2050) IOSYS                     
 2050 FORMAT (' ILLEGAL TARGET ZONE NUMBER = ',I6)                     
      IFLG = 8                                                         
      RETURN                                                           
!                                                                      
! FORWARD TRANSFORMATION.                                              
!                                                                      
  540 IPROJ=IOSYS                                                      
      ISPHER = INSPH                                                   
      IF (IOSYS.GE.3) CALL SPHDZ0(INSPH,TPARIO)                        
!                                                                      
!     CHECK FOR CHANGE IN ZONE FROM LAST USE OF THE OUTPUT PROJECTION  
!                                                                      
      IF (IOSYS .EQ. 1 .AND. IOZONE .NE. SWITCH(9)) THEN               
         SWITCH(1) = 0                                                 
         IOMOD = 1                                                     
      END IF                                                           
      IF (IOSYS .EQ. 2 .AND. IOZONE .NE. SWITCH(IOSPCS)) THEN          
         SWITCH(2) = 0                                                 
         IOMOD = 1                                                     
      END IF                                                           
      IF (IOZONE .NE. SWITCH(IOSYS)) THEN                              
         SWITCH(IOSYS) = 0                                             
         IOMOD = 1                                                     
      END IF                                                           
!                                                                      
      IF (IOSYS .EQ. 1) THEN                                           
         TPARIO(1) = COORD(1)                                          
         TPARIO(2) = COORD(2)                                          
         CALL SPHDZ0(INSPH,DUMMY)                                      
         TPARIO(14) = DUMMY(1)                                         
         TPARIO(15) = DUMMY(2)                                         
         IF (IOMOD .NE. 0) THEN                                        
            CALL PJINIT (IOSYS,IOZONE,TPARIO)                          
            IF (IERROR .NE. 0) IOZN = 99999                            
            IF (IERROR .NE. 0) GO TO 900                               
         END IF                                                        
         CALL PJ01Z0 (COORD,CRDIO,FWD)                                 
      END IF                                                           
!                                                                      
      IF (IOSYS .GT. 1) THEN                                           
         IF (IOMOD .NE. 0) THEN                                        
            MSYS = IOSPCS                                              
            CALL PJINIT (IOSYS,IOZONE,TPARIO)                          
            IF (IERROR .NE. 0) IOZN = 99999                            
            IF (IERROR .NE. 0) GO TO 900                               
         END IF                                                        
         IF (IOSYS .EQ. 2) CALL PJ02Z0 (COORD,CRDIO,FWD)               
         IF (IOSYS .EQ. 3) CALL PJ03Z0 (COORD,CRDIO,FWD)               
         IF (IOSYS .EQ. 4) CALL PJ04Z0 (COORD,CRDIO,FWD)               
         IF (IOSYS .EQ. 5) CALL PJ05Z0 (COORD,CRDIO,FWD)               
         IF (IOSYS .EQ. 6) CALL PJ06Z0 (COORD,CRDIO,FWD)               
         IF (IOSYS .EQ. 7) CALL PJ07Z0 (COORD,CRDIO,FWD)               
         IF (IOSYS .EQ. 8) CALL PJ08Z0 (COORD,CRDIO,FWD)               
         IF (IOSYS .EQ. 9) CALL PJ09Z0 (COORD,CRDIO,FWD)               
         IF (IOSYS .EQ. 10) CALL PJ10Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 11) CALL PJ11Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 12) CALL PJ12Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 13) CALL PJ13Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 14) CALL PJ14Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 15) CALL PJ15Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 16) CALL PJ16Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 17) CALL PJ17Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 18) CALL PJ18Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 19) CALL PJ19Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 20) CALL PJ20Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 21) CALL PJ21Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 22) CALL PJ22Z0 (COORD,CRDIO,FWD)              
         IF (IOSYS .EQ. 23) CALL PJ23Z0 (COORD,CRDIO,FWD)              
      END IF                                                           
!                                                                      
  900 IFLG = IERROR                                                    
      DO 910 I = 1,15                                                  
  910 TPARIO(I) = PDIO(I)                                              
  920 CRDIO(1) = FACTOR * CRDIO(1)                                     
      CRDIO(2) = FACTOR * CRDIO(2)                                     
      RETURN                                                           
!                                                                      
      END                                                              
!                   MLFNZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION MLFNZ0 (E0,E1,E2,E3,PHI)               
!                                                                      
! FUNCTION TO COMPUTE CONSTANT (M).                                    
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      DATA TWO,FOUR,SIX /2.0D0,4.0D0,6.0D0/                            
!                                                                      
      MLFNZ0 = E0 * PHI - E1 * DSIN (TWO * PHI) + E2 * DSIN (FOUR * PHI) & 
     & - E3 * DSIN (SIX * PHI)                                         
!                                                                      
      RETURN                                                           
      END                                                              
!                   MSFNZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION MSFNZ0 (ECCENT,SINPHI,COSPHI)          
!                                                                      
! FUNCTION TO COMPUTE CONSTANT (SMALL M).                              
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      DATA ONE /1.0D0/                                                 
!                                                                      
      CON = ECCENT * SINPHI                                            
      MSFNZ0 = COSPHI / DSQRT (ONE - CON * CON)                        
!                                                                      
      RETURN                                                           
      END                                                              
!                   PAKCZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION PAKCZ0 (PAK)                           
!                                                                      
! SUBROUTINE TO CONVERT 2 DIGIT PACKED DMS TO 3 DIGIT PACKED DMS ANGLE.
!                                                                      
! SGNA : SIGN OF ANGLE                                                 
! DEGS : DEGREES PORTION OF ANGLE                                      
! MINS : MINUTES PORTION OF ANGLE                                      
! SECS : SECONDS PORTION OF ANGLE                                      
!                                                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      INTEGER*4 DEGS,MINS                                              
      CHARACTER*1 SGNA,IBLANK,NEG                                      
      DATA ZERO,CON1,CON2 /0.0D0,10000.0D0,100.0D0/                    
      DATA CON3,CON4 /1000000.0D0,1000.0D0/                            
      DATA TOL /1.0D-3/                                                
      DATA IBLANK,NEG /' ','-'/                                        
!                                                                      
      SGNA = IBLANK                                                    
      IF (PAK .LT. ZERO) SGNA = NEG                                    
      CON = DABS (PAK)                                                 
      DEGS = IDINT ((CON / CON1) + TOL)                                
      CON = DMOD ( CON , CON1)                                         
      MINS = IDINT ((CON / CON2) + TOL)                                
      SECS = DMOD (CON , CON2)                                         
!                                                                      
      CON = DBLE (DEGS) * CON3 + DBLE (MINS) * CON4 + SECS             
      IF (SGNA .EQ. NEG) CON = - CON                                   
      PAKCZ0 = CON                                                     
      RETURN                                                           
!                                                                      
      END                                                              
!                   PAKDZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      SUBROUTINE PAKDZ0 (PAK,SGNA,DEGS,MINS,SECS)                      
!                                                                      
! SUBROUTINE TO CONVERT PACKED DMS TO UNPACKED DMS ANGLE.              
!                                                                      
! SGNA : SIGN OF ANGLE                                                 
! DEGS : DEGREES PORTION OF ANGLE                                      
! MINS : MINUTES PORTION OF ANGLE                                      
! SECS : SECONDS PORTION OF ANGLE                                      
!                                                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      REAL*4 SECS                                                      
      INTEGER*4 DEGS,MINS                                              
      CHARACTER*1 SGNA,IBLANK,NEG                                      
      DATA ZERO,CON1,CON2 /0.0D0,1000000.0D0,1000.0D0/                 
      DATA TOL /1.0D-4/                                                
      DATA IBLANK,NEG /' ','-'/                                        
!                                                                      
      SGNA = IBLANK                                                    
      IF (PAK .LT. ZERO) SGNA = NEG                                    
      CON = DABS (PAK)                                                 
      DEGS = IDINT ((CON / CON1) + TOL)                                
      CON = DMOD ( CON , CON1)                                         
      MINS = IDINT ((CON / CON2) + TOL)                                
      SECS = SNGL ( DMOD (CON , CON2))                                 
      RETURN                                                           
!                                                                      
      END                                                              
!                   PAKRZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION PAKRZ0 (ANG)                           
!                                                                      
! FUNCTION TO CONVERT DMS PACKED ANGLE INTO RADIANS.                   
!                                                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      DATA SECRAD /0.4848136811095359D-5/                              
!                                                                      
! CONVERT ANGLE TO SECONDS OF ARC                                      
!                                                                      
      SEC = PAKSZ0 (ANG)                                               
!                                                                      
! CONVERT ANGLE TO RADIANS.                                            
!                                                                      
      PAKRZ0 = SEC * SECRAD                                            
!                                                                      
      RETURN                                                           
      END                                                              
!                   PAKSZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION PAKSZ0 (ANG)                           
!                                                                      
! FUNCTION TO CONVERT DMS PACKED ANGLE INTO SECONDS OF ARC.            
!                                                                      
      IMPLICIT REAL*8 (A-H,M-Z)                                        
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      DIMENSION CODE(2)                                                
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      DATA CODE /1000000.0D0,1000.0D0/                                 
      DATA ZERO,ONE /0.0D0,1.0D0/                                      
      DATA C1,C2 /3600.0D0,60.0D0/                                     
      DATA TOL /1.0D-4/                                                
!                                                                      
! SEPARATE DEGREE FIELD.                                               
!                                                                      
      FACTOR = ONE                                                     
      IF (ANG .LT. ZERO) FACTOR = - ONE                                
      SEC = DABS(ANG)                                                  
      TMP = CODE(1)                                                    
      I = IDINT ((SEC / TMP) + TOL)                                    
      IF (I .GT. 360) GO TO 020                                        
      DEG = DBLE (I)                                                   
!                                                                      
! SEPARATE MINUTES FIELD.                                              
!                                                                      
      SEC = SEC - DEG * TMP                                            
      TMP = CODE(2)                                                    
      I = IDINT ((SEC / TMP) + TOL)                                    
      IF (I .GT. 60) GO TO 020                                         
      MIN = DBLE (I)                                                   
!                                                                      
! SEPARATE SECONDS FIELD.                                              
!                                                                      
      SEC = SEC - MIN * TMP                                            
      IF (SEC .GT. C2) GO TO 020                                       
      SEC = FACTOR * (DEG * C1 + MIN * C2 + SEC)                       
      GO TO 040                                                        
!                                                                      
! ERROR DETECTED IN DMS FORM.                                          
!                                                                      
  020 WRITE (IPELUN,2000) ANG                                          
 2000 FORMAT ('0ERROR PAKSZ0'/                                     &
     &        ' ILLEGAL DMS FIELD =',F15.3)                            
      STOP 16                                                          
!                                                                      
  040 PAKSZ0 = SEC                                                     
!                                                                      
      RETURN                                                           
      END                                                              
!                   PHI1Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION PHI1Z0 (ECCENT,QS)                     
!                                                                      
! FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-1).                          
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 II,NIT                                                 
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      DATA HALF,ONE /0.5D0,1.0D0/                                      
      DATA EPSLN,TOL,NIT /1.0D-7,1.0D-10,15/                           
!                                                                      
      PHI1Z0 = ASINZ0 (HALF * QS)                                      
      IF (ECCENT .LT. EPSLN) RETURN                                    
!                                                                      
      ECCNTS = ECCENT * ECCENT                                         
      PHI = PHI1Z0                                                     
      DO 020 II = 1,NIT                                                
      SINPI = DSIN (PHI)                                               
      COSPI = DCOS (PHI)                                               
      CON = ECCENT * SINPI                                             
      COM = ONE - CON * CON                                            
      DPHI = HALF * COM * COM / COSPI * (QS / (ONE - ECCNTS) -       &  
     &       SINPI / COM + HALF / ECCENT * DLOG ((ONE - CON) /        & 
     &       (ONE + CON)))                                             
      PHI = PHI + DPHI                                                 
      IF (DABS(DPHI) .GT. TOL) GO TO 020                               
      PHI1Z0 = PHI                                                     
      RETURN                                                           
  020 CONTINUE                                                         
!                                                                      
      IF (IPEMSG .EQ. 0) WRITE (IPELUN,2000) NIT,ECCENT,QS             
 2000 FORMAT ('0ERROR PHI1Z0' /                                       & 
     &        ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/  & 
     &        ' ECCENTRICITY =',D25.16,'   QS =',D25.16)               
      IERROR = 001                                                     
      RETURN                                                           
!                                                                      
      END                                                              
!                   PHI2Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION PHI2Z0 (ECCENT,TS)                     
!                                                                      
! FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-2).                          
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 II,NIT                                                 
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      DATA HALF,ONE,TWO /0.5D0,1.0D0,2.0D0/                            
      DATA TOL,NIT /1.0D-10,15/                                        
      DATA HALFPI /1.5707963267948966D0/                               
!                                                                      
      ECCNTH = HALF * ECCENT                                           
      PHI = HALFPI - TWO * DATAN (TS)                                  
      DO 020 II = 1,NIT                                                
      SINPI = DSIN (PHI)                                               
      CON = ECCENT * SINPI                                             
      DPHI = HALFPI - TWO * DATAN (TS * ((ONE - CON) /   &              
     &       (ONE + CON)) ** ECCNTH) - PHI                             
      PHI = PHI + DPHI                                                 
      IF (DABS(DPHI) .GT. TOL) GO TO 020                               
      PHI2Z0 = PHI                                                     
      RETURN                                                           
  020 CONTINUE                                                         
!                                                                      
      IF (IPEMSG .EQ. 0) WRITE (IPELUN,2000) NIT,ECCENT,TS             
 2000 FORMAT ('0ERROR PHI2Z0' /                                       & 
     &        ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/  & 
     &        ' ECCENTRICITY =',D25.16,'   TS =',D25.16)               
      IERROR = 002                                                     
      RETURN                                                           
!                                                                      
      END                                                              
!                   PHI3Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION PHI3Z0 (ML,E0,E1,E2,E3)                
!                                                                      
! FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-3).                          
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 II,NIT                                                 
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      DATA TWO,FOUR,SIX /2.0D0,4.0D0,6.0D0/                            
      DATA TOL,NIT /1.0D-10,15/                                        
!                                                                      
      PHI = ML                                                         
      DO 020 II = 1,NIT                                                
      DPHI = (ML + E1 * DSIN (TWO * PHI) - E2 * DSIN (FOUR * PHI)  &  
     &       + E3 * DSIN (SIX * PHI)) / E0 - PHI                       
      PHI = PHI + DPHI                                                 
      IF (DABS(DPHI) .GT. TOL) GO TO 020                               
      PHI3Z0 = PHI                                                     
      RETURN                                                           
  020 CONTINUE                                                         
!                                                                      
      IF (IPEMSG .EQ. 0) WRITE (IPELUN,2000) NIT,ML,E0,E1,E2,E3        
 2000 FORMAT ('0ERROR PHI3Z0' /                                      &  
     &        ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/  & 
     &        ' ML =',D25.16,'   E0 =',D25.16/                        & 
     &        ' E1 =',D25.16,'   E2 =',D25.16,'   E3=',D25.16)         
      IERROR = 003                                                     
      RETURN                                                           
!                                                                      
      END                                                              
!                   PHI4Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      SUBROUTINE PHI4Z0 (ECCNTS,E0,E1,E2,E3,A,B,C,PHI)                 
!                                                                      
! FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-4).                          
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 II,NIT                                                 
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      DATA ONE,TWO,FOUR,SIX /1.0D0,2.0D0,4.0D0,6.0D0/                  
      DATA TOL,NIT /1.0D-10,15/                                        
!                                                                      
      PHI = A                                                          
      DO 020 II = 1,NIT                                                
      SINPHI = DSIN (PHI)                                              
      TANPHI = DTAN (PHI)                                              
      C = TANPHI * DSQRT (ONE - ECCNTS * SINPHI * SINPHI)              
      SIN2PH = DSIN (TWO * PHI)                                        
      ML = E0 * PHI - E1 * SIN2PH + E2 * DSIN (FOUR * PHI)             &
     &      - E3 * DSIN (SIX * PHI)                                    
      MLP = E0 - TWO * E1 * DCOS (TWO * PHI) + FOUR * E2 *             &
     &      DCOS (FOUR * PHI) - SIX * E3 * DCOS (SIX * PHI)            
      CON1 = TWO * ML + C * (ML * ML + B) - TWO * A *                  &
     &       (C * ML + ONE)                                            
      CON2 = ECCNTS * SIN2PH * (ML * ML + B - TWO * A * ML) / (TWO * C)
      CON3 = TWO * (A - ML) * (C * MLP - TWO / SIN2PH) - TWO * MLP     
      DPHI = CON1 / (CON2 + CON3)                                      
      PHI = PHI + DPHI                                                 
      IF (DABS(DPHI) .GT. TOL) GO TO 020                               
      RETURN                                                           
  020 CONTINUE                                                         
!                                                                      
      IF (IPEMSG .EQ. 0) WRITE (IPELUN,2000) NIT,E0,E1,E2,E3,A,B,C,   & 
     & ECCNTS                                                          
 2000 FORMAT ('0ERROR PHI4Z0' /                                        &
     &        ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/   &
     &        ' E0 =',D25.16,'   E1 =',D25.16/                         &
     &        ' E2 =',D25.16,'   E3 =',D25.16/                         &
     &        ' A  =',D25.16,'   B  =',D25.16/                         &
     &        ' C  =',D25.16/                                          &
     &        ' ECCENTRICITY SQUARE =',D25.16)                         
      IERROR = 004                                                     
      RETURN                                                           
!                                                                      
      END                                                              
!                   PJINIT                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      SUBROUTINE PJINIT (ISYS,ZONE,DATA)                               
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      REAL*4 SECS(5)                                                   
      real*8 save9
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN,ITEMP               
      INTEGER*4 LAND, PATH, LIMIT, IND02, IND06, IND09, ISYS, KEEPZN   
      INTEGER*4 SWITCH(23),I,ZONE,DEGS(5),MINS(5)                      
      INTEGER*4 ID, IND, ITEM, ITYPE, MODE, N, MSYS                    
      INTEGER*4 ISPHER, LUNIT, LU27, LU83, LEN, NAD27(134), NAD83(134) 
      CHARACTER*128 DATUM, FILE27, FILE83                              
      CHARACTER*32 PNAME                                               
      CHARACTER*1 SGNA(5)                                              
!                                                                      
      DIMENSION DATA(15),BUFFL(15)                                     
      DIMENSION TABLE(9)                                               
      DIMENSION PR(20),XLR(20)                                         
      DIMENSION ACOEF(6),BCOEF(6)                                      
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z,E4Z                    
      COMMON /SPHRZ0/ AZZ                                              
      COMMON /NORM/ Q,T,U,W,ES22,P22,SA,CA,XJ                          
      COMMON /SPCS/ ISPHER,LU27,LU83,LEN,MSYS,FILE27,FILE83            
      COMMON /PJ02/ ITYPE                                              
      COMMON /PJ03/ A03,LON003,X003,Y003,C,E03,ES03,NS03,RH003         
      COMMON /PJ04/ A04,LON004,X004,Y004,E04,F04,NS04,RH004            
      COMMON /PJ05/ A05,LON005,X005,Y005,E05,M1                        
      COMMON /PJ06/ A06,LON006,X006,Y006,E06,E4,FAC,MCS,TCS,IND06      
      COMMON /PJ07/ A07,LON007,X007,Y007,E07,E007,E107,E207,E307,ES07, &
     &              ML007                                              
      COMMON /PJ08/ A08,LON008,X008,Y008,E008,E108,E208,E308,GL,NS08,  &
     &              RH008                                              
      COMMON /PJ09/ A09,LON009,X009,Y009,ES09,ESP,E009,E109,E209,E309, &
     &              KS009,LAT009,ML009,IND09                           
      COMMON /PJ10/ A10,LON010,X010,Y010,COSP10,LAT010,SINP10          
      COMMON /PJ11/ A11,LON011,X011,Y011,COSP11,LAT011,SINP11          
      COMMON /PJ12/ A12,LON012,X012,Y012,COSP12,LAT012,SINP12          
      COMMON /PJ13/ A13,LON013,X013,Y013,COSP13,LAT013,SINP13          
      COMMON /PJ14/ A14,LON014,X014,Y014,COSP14,LAT014,SINP14          
      COMMON /PJ15/ A15,LON015,X015,Y015,COSP15,LAT015,P,SINP15        
      COMMON /PJ16/ A16,LON016,X016,Y016                               
      COMMON /PJ17/ A17,LON017,X017,Y017,LAT1                          
      COMMON /PJ18/ A18,LON018,X018,Y018                               
      COMMON /PJ19/ A19,LON019,X019,Y019                               
      COMMON /PJ20/ LON020,X020,Y020,AL,BL,COSALF,COSGAM,E20,EL,SINALF,  & 
     &              SINGAM,U0                                          
      COMMON /PJ21/ A21,LON021,X021,Y021,PR,XLR                        
      COMMON /PJ22/ A22,X022,Y022,A2,A4,B,C1,C3,LAND,PATH              
      COMMON /PJ23/ A23,LON023,X023,Y023,ACOEF,BCOEF,EC,LAT023,    &    
     &              CCHIO,SCHIO,N                                      
      COMMON /TOGGLE/ SWITCH                                           
!                                                                      
      data save9 /1.0D128/
      DATA PI /3.14159265358979323846D0/                               
      DATA HALFPI /1.5707963267948966D0/                               
      DATA ZERO,HALF,ONE,TWO /0.0D0,0.5D0,1.0D0,2.0D0/                 
      DATA EPSLN /1.0D-10/                                             
      DATA TOL /1.0D-7/                                                
      DATA TOL09 /1.0D-5/                                              
      DATA NINTYD /90000000.0D0/                                       
      DATA DG1 /0.01745329252D0/                                       
!                                                                      
      DATA NAD27/0101,0102,5010,5300,0201,0202,0203,0301,0302,0401,0402, & 
     &           0403,0404,0405,0406,0407,0501,0502,0503,0600,0700,0901, & 
     &           0902,0903,1001,1002,5101,5102,5103,5104,5105,1101,1102, &  
     &           1103,1201,1202,1301,1302,1401,1402,1501,1502,1601,1602, & 
     &           1701,1702,1703,1801,1802,1900,2001,2002,2101,2102,2103, & 
     &           2111,2112,2113,2201,2202,2203,2301,2302,2401,2402,2403, &  
     &           2501,2502,2503,2601,2602,2701,2702,2703,2800,2900,3001, &  
     &           3002,3003,3101,3102,3103,3104,3200,3301,3302,3401,3402, & 
     &           3501,3502,3601,3602,3701,3702,3800,3901,3902,4001,4002, & 
     &           4100,4201,4202,4203,4204,4205,4301,4302,4303,4400,4501, & 
     &           4502,4601,4602,4701,4702,4801,4802,4803,4901,4902,4903, & 
     &           4904,5001,5002,5003,5004,5005,5006,5007,5008,5009,5201, & 
     &           5202,5400/                                            
!                                                                      
      DATA NAD83/0101,0102,5010,5300,0201,0202,0203,0301,0302,0401,0402, & 
     &           0403,0404,0405,0406,0000,0501,0502,0503,0600,0700,0901, & 
     &           0902,0903,1001,1002,5101,5102,5103,5104,5105,1101,1102, & 
     &           1103,1201,1202,1301,1302,1401,1402,1501,1502,1601,1602, & 
     &           1701,1702,1703,1801,1802,1900,2001,2002,2101,2102,2103, & 
     &           2111,2112,2113,2201,2202,2203,2301,2302,2401,2402,2403, & 
     &           2500,0000,0000,2600,0000,2701,2702,2703,2800,2900,3001, & 
     &           3002,3003,3101,3102,3103,3104,3200,3301,3302,3401,3402, & 
     &           3501,3502,3601,3602,3701,3702,3800,3900,0000,4001,4002, & 
     &           4100,4201,4202,4203,4204,4205,4301,4302,4303,4400,4501, & 
     &           4502,4601,4602,4701,4702,4801,4802,4803,4901,4902,4903, & 
     &           4904,5001,5002,5003,5004,5005,5006,5007,5008,5009,5200, & 
     &           0000,5400/                                            
! .................................................................... 
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                              .  U T M  .                             
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 1) THEN                                            
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(1).NE.0 .AND. SWITCH(1).EQ.ZONE) RETURN            
         SWITCH(1) = ZONE                                              
         IF (SWITCH(9).NE.0.AND.SWITCH(9).EQ.ZONE.AND.DATA(14).EQ.SAVE) &  
     &   RETURN                                                        
         KEEPZN = ZONE                                                 
         ZONE = IABS(ZONE)                                             
         SAVE = DATA(1)                                                
         IF (ZONE .EQ. 0) THEN                                         
            ZONE = IDINT( ( (DATA(1) * 180.0D0 / PI)                    &
     &             + (TOL09 / 3600.D0) ) / 6.D0 )                      
            IND = 1                                                    
            IF (DATA(1) .LT. ZERO) IND = 0                             
            ZONE = MOD ((ZONE + 30), 60) + IND                         
            KEEPZN = ZONE                                              
            IF (DATA(2) .LT. ZERO) KEEPZN = -ZONE                      
         END IF                                                        
         IF (ZONE.LT.1 .OR. ZONE.GT.60) THEN                           
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,140) KEEPZN               
  140       FORMAT ('0ERROR PJ01Z0'/                           &        
     &              ' ILLEGAL ZONE NO. : ',I10)                        
            IERROR = 011                                               
            RETURN                                                     
         END IF                                                        
         BUFFL(1) = DATA(14)                                           
         BUFFL(2) = DATA(15)                                           
         BUFFL(3) = 0.9996D0                                           
         BUFFL(4) = ZERO                                               
         BUFFL(5) = DBLE (6 * ZONE - 183) * 1.0D6                      
         BUFFL(6) = ZERO                                               
         BUFFL(7) = 500000.0D0                                         
         BUFFL(8) = ZERO                                               
         IF (DATA(2) .LT. ZERO) BUFFL(8) = 10000000.0D0                
         IF (KEEPZN .LT. 0) BUFFL(8) = 10000000.0D0                    
         IF (BUFFL(1).NE.0.0D0.AND.BUFFL(1).NE.SAVE9) SWITCH(9) = 0    
         SAVE9 = BUFFL(1)                                              
         ITEMP = IPPARM                                                
         IPPARM = 1                                                    
         DO 145 I=1,8                                                  
            DATA(I) = BUFFL(I)                                         
  145    CONTINUE                                                      
         AZ = DATA(14)                                                 
         EZ = DATA(15)                                                 
         SWITCH(9) = 0                                                 
         GO TO 900                                                     
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                           .  STATE PLANE  .                          
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 2) THEN                                            
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(2).NE.0 .AND. SWITCH(2).EQ.ZONE) RETURN            
         IF (ISPHER .NE. 0 .AND. ISPHER .NE. 8) THEN                   
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,205) ISPHER               
  205       FORMAT('0ERROR PJ02Z0'/                                   &
     &             ' SPHEROID NO. ',I4,' IS INVALID FOR STATE PLANE', &
     &             ' TRANSFORMATIONS')                                 
            IERROR = 020                                               
            RETURN                                                     
         END IF                                                        
         IF (ZONE .GT. 0) THEN                                         
            IND02 = 0                                                  
            IF (ISPHER .EQ. 0) THEN                                    
               DO 210 I = 1,134                                        
                  IF (ZONE .EQ. NAD27(I)) IND02 = I                    
  210          CONTINUE                                                
            END IF                                                     
            IF (ISPHER .EQ. 8) THEN                                    
               DO 220 I = 1,134                                        
                  IF (ZONE .EQ. NAD83(I)) IND02 = I                    
  220          CONTINUE                                                
            END IF                                                     
            IF (IND02 .EQ. 0) THEN                                     
               IF (IPEMSG .EQ. 0) WRITE (IPELUN,240) ZONE, ISPHER      
               IERROR = 021                                            
               RETURN                                                  
            END IF                                                     
         ELSE                                                          
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,240) ZONE, ISPHER         
            IERROR = 021                                               
            RETURN                                                     
         END IF                                                        
         IF (ISPHER .EQ. 0) THEN                                       
            LUNIT = LU27                                               
            DATUM = FILE27                                             
         END IF                                                        
         IF (ISPHER .EQ. 8) THEN                                       
            LUNIT = LU83                                               
            DATUM = FILE83                                             
         END IF                                                        
         OPEN (UNIT=LUNIT,FILE=DATUM,STATUS='OLD',ACCESS='DIRECT',    &
     &   RECL=LEN)                                                     
         READ (UNIT=LUNIT,REC=IND02) PNAME,ID,TABLE                    
         CLOSE (UNIT=LUNIT,STATUS='KEEP')                              
         IF (ID .LE. 0) THEN                                           
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,240) ZONE, ISPHER         
  240       FORMAT('0ERROR PJ02Z0'/                                   &
     &             ' ILLEGAL ZONE NO. : ',I8,' FOR SPHEROID NO. : ',I4)
            IERROR = 021                                               
            RETURN                                                     
         END IF                                                        
         ITYPE = ID                                                    
         AZ = TABLE(1)                                                 
         ES = TABLE(2)                                                 
         ESZ = ES                                                      
         EZ  = DSQRT(ES)                                               
         E0Z = E0FNZ0(ES)                                              
         E1Z = E1FNZ0(ES)                                              
         E2Z = E2FNZ0(ES)                                              
         E3Z = E3FNZ0(ES)                                              
         E4Z = E4FNZ0(EZ)                                              
         ITEMP = IPPARM                                                
         IPPARM = 1                                                    
!                                                                      
!     TRANSVERSE MERCATOR PROJECTION                                   
!                                                                      
         IF (ITYPE .EQ. 1) THEN                                        
            DATA(3) = TABLE(4)                                         
            DATA(5) = PAKCZ0(TABLE(3))                                 
            DATA(6) = PAKCZ0(TABLE(7))                                 
            DATA(7) = TABLE(8)                                         
            DATA(8) = TABLE(9)                                         
            MSYS = 9                                                   
            SWITCH(MSYS) = 0                                           
            GO TO 900                                                  
         END IF                                                        
!                                                                      
!     LAMBERT CONFORMAL PROJECTION                                     
!                                                                      
         IF (ITYPE .EQ. 2) THEN                                        
            DATA(3) = PAKCZ0(TABLE(6))                                 
            DATA(4) = PAKCZ0(TABLE(5))                                 
            DATA(5) = PAKCZ0(TABLE(3))                                 
            DATA(6) = PAKCZ0(TABLE(7))                                 
            DATA(7) = TABLE(8)                                         
            DATA(8) = TABLE(9)                                         
            MSYS = 4                                                   
            SWITCH(MSYS) = 0                                           
            GO TO 400                                                  
         END IF                                                        
!                                                                      
!     POLYCONIC PROJECTION                                             
!                                                                      
         IF (ITYPE .EQ. 3) THEN                                        
            DATA(5) = PAKCZ0(TABLE(3))                                 
            DATA(6) = PAKCZ0(TABLE(4))                                 
            DATA(7) = TABLE(5)                                         
            DATA(8) = TABLE(6)                                         
            MSYS = 7                                                   
            SWITCH(MSYS) = 0                                           
            GO TO 700                                                  
         END IF                                                        
!                                                                      
!     OBLIQUE MERCATOR PROJECTION                                      
!                                                                      
         IF (ITYPE .EQ. 4) THEN                                        
            DATA(3) = TABLE(4)                                         
            DATA(4) = PAKCZ0(TABLE(6))                                 
            DATA(5) = PAKCZ0(TABLE(3))                                 
            DATA(6) = PAKCZ0(TABLE(7))                                 
            DATA(7) = TABLE(8)                                         
            DATA(8) = TABLE(9)                                         
            DATA(13) = ONE                                             
            MSYS = 20                                                  
            SWITCH(MSYS) = 0                                           
            GO TO 2000                                                 
         END IF                                                        
!                                                                      
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                    .  ALBERS CONICAL EQUAL AREA  .                   
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 3) THEN                                            
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(3).NE.0 .AND. SWITCH(3).EQ.ZONE) RETURN            
         SWITCH(3) = 0                                                 
         A03 = AZ                                                      
         E03 = EZ                                                      
         ES03 = ESZ                                                    
         LAT1 = PAKRZ0 (DATA(3))                                       
         LAT2 = PAKRZ0 (DATA(4))                                       
         IF (DABS(LAT1+LAT2) .LT. EPSLN) THEN                          
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,340)                      
  340       FORMAT ('0ERROR PJ03Z0'/                                  &
     &              ' EQUAL LATITUDES FOR ST. PARALLELS ON OPPOSITE', &
     &              ' SIDES OF EQUATOR')                               
            IERROR = 031                                               
            RETURN                                                     
         END IF                                                        
         LON003 = PAKRZ0 (DATA(5))                                     
         LAT003 = PAKRZ0 (DATA(6))                                     
         X003 = DATA(7)                                                
         Y003 = DATA(8)                                                
         SINP03 = DSIN (LAT1)                                          
         CON = SINP03                                                  
         COSP03 = DCOS (LAT1)                                          
         MS1 = MSFNZ0 (E03,SINP03,COSP03)                              
         QS1 = QSFNZ0 (E03,SINP03,COSP03)                              
         SINP03 = DSIN (LAT2)                                          
         COSP03 = DCOS (LAT2)                                          
         MS2 = MSFNZ0 (E03,SINP03,COSP03)                              
         QS2 = QSFNZ0 (E03,SINP03,COSP03)                              
         SINP03 = DSIN (LAT003)                                        
         COSP03 = DCOS (LAT003)                                        
         QS0 = QSFNZ0 (E03,SINP03,COSP03)                              
         IF (DABS(LAT1-LAT2) .GE. EPSLN) THEN                          
            NS03 = (MS1 * MS1 - MS2 * MS2) / (QS2 - QS1)               
         ELSE                                                          
            NS03 = CON                                                 
         END IF                                                        
         C = MS1 * MS1 + NS03 * QS1                                    
         RH003 = A03 * DSQRT (C - NS03 * QS0) / NS03                   
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LAT1,SGNA(1),DEGS(1),MINS(1),SECS(1))            
         CALL RADDZ0 (LAT2,SGNA(2),DEGS(2),MINS(2),SECS(2))            
         CALL RADDZ0 (LON003,SGNA(3),DEGS(3),MINS(3),SECS(3))          
         CALL RADDZ0 (LAT003,SGNA(4),DEGS(4),MINS(4),SECS(4))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,350) A03,ES03,              &
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,4),            &
     &            X003,Y003                                            
  350   FORMAT ('0INITIALIZATION PARAMETERS (ALBERS CONICAL EQUAL-AREA'  & 
     &           ' PROJECTION)'/                                       &
     &           ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/    &
     &           ' ECCENTRICITY SQUARED         =',F12.9/              &
     &           ' LATITUDE OF 1ST ST. PARALLEL = ',A1,2I3,F7.3/       &
     &           ' LATITUDE OF 2ND ST. PARALLEL = ',A1,2I3,F7.3/       &
     &           ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/       &
     &           ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/       &
     &           ' FALSE EASTING                =',F12.2,' METERS'/    &
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A03                                                 
         DATA(2) = ES03                                                
         SWITCH(3) = ZONE                                              
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                     .  LAMBERT CONFORMAL CONIC  .                    
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 4) THEN                                            
!                                                                      
  400    IERROR = 0                                                    
         IF (SWITCH(4).NE.0 .AND. SWITCH(4).EQ.ZONE) RETURN            
         SWITCH(4) = 0                                                 
         A04 = AZ                                                      
         E04 = EZ                                                      
         ES = ESZ                                                      
         LAT1 = PAKRZ0 (DATA(3))                                       
         LAT2 = PAKRZ0 (DATA(4))                                       
         IF (DABS(LAT1+LAT2) .LT. EPSLN) THEN                          
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,440)                      
  440       FORMAT ('0ERROR PJ04Z0'/                                  & 
     &              ' EQUAL LATITUDES FOR ST. PARALLELS ON OPPOSITE', & 
     &              ' SIDES OF EQUATOR')                               
            IERROR = 041                                               
            RETURN                                                     
         END IF                                                        
         LON004 = PAKRZ0 (DATA(5))                                     
         LAT004 = PAKRZ0 (DATA(6))                                     
         X004 = DATA(7)                                                
         Y004 = DATA(8)                                                
         SINP04 = DSIN (LAT1)                                          
         CON = SINP04                                                  
         COSP04 = DCOS (LAT1)                                          
         MS1 = MSFNZ0 (E04,SINP04,COSP04)                              
         TS1 = TSFNZ0 (E04,LAT1,SINP04)                                
         SINP04 = DSIN (LAT2)                                          
         COSP04 = DCOS (LAT2)                                          
         MS2 = MSFNZ0 (E04,SINP04,COSP04)                              
         TS2 = TSFNZ0 (E04,LAT2,SINP04)                                
         SINP04 = DSIN (LAT004)                                        
         TS0 = TSFNZ0 (E04,LAT004,SINP04)                              
         IF (DABS(LAT1-LAT2) .GE. EPSLN) THEN                          
            NS04 = DLOG (MS1 / MS2) / DLOG (TS1 / TS2)                 
         ELSE                                                          
            NS04 = CON                                                 
         END IF                                                        
         F04 = MS1 / (NS04 * TS1 ** NS04)                              
         RH004 = A04 * F04 * TS0 ** NS04                               
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LAT1,SGNA(1),DEGS(1),MINS(1),SECS(1))            
         CALL RADDZ0 (LAT2,SGNA(2),DEGS(2),MINS(2),SECS(2))            
         CALL RADDZ0 (LON004,SGNA(3),DEGS(3),MINS(3),SECS(3))          
         CALL RADDZ0 (LAT004,SGNA(4),DEGS(4),MINS(4),SECS(4))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,450) A04,ES,         &       
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,4),     &       
     &            X004,Y004                                            
  450    FORMAT ('0INITIALIZATION PARAMETERS (LAMBERT CONFORMAL CONIC', & 
     &           ' PROJECTION)'/                                      &
     &           ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/   &
     &           ' ECCENTRICITY SQUARED         =',F12.9/             &
     &           ' LATITUDE OF 1ST ST. PARALLEL = ',A1,2I3,F7.3/      &
     &           ' LATITUDE OF 2ND ST. PARALLEL = ',A1,2I3,F7.3/      &
     &           ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/      &
     &           ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/      &
     &           ' FALSE EASTING                =',F12.2,' METERS'/   &
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A04                                                 
         DATA(2) = ES                                                  
         SWITCH(4) = ZONE                                              
!                                                                      
!     LIST STATE PLANE INITIALIZATION PARAMETERS IF NECESSARY          
!                                                                      
         IF (ISYS .EQ. 2) THEN                                         
            IPPARM = ITEMP                                             
            IF (IERROR .NE. 0) RETURN                                  
            IF (IPPARM .EQ. 0) WRITE (IPPLUN,470) ZONE, PNAME          
  470    FORMAT (' INITIALIZATION PARAMETERS (STATE PLANE PROJECTION)'/ & 
     &           ' ZONE NUMBER = ',I4,5X,' ZONE NAME = ',A32)          
            SWITCH(2) = ZONE                                           
            RETURN                                                     
         END IF                                                        
!                                                                      
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                            .  MERCATOR  .                            
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 5) THEN                                            
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(5).NE.0 .AND. SWITCH(5).EQ.ZONE) RETURN            
         SWITCH(5) = 0                                                 
         A05 = AZ                                                      
         E05 = EZ                                                      
         ES = ESZ                                                      
         LON005 = PAKRZ0 (DATA(5))                                     
         LAT1 = PAKRZ0 (DATA(6))                                       
         M1 = DCOS(LAT1) / (DSQRT( ONE - ES * DSIN(LAT1) **2))         
         X005 = DATA(7)                                                
         Y005 = DATA(8)                                                
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LAT1,SGNA(1),DEGS(1),MINS(1),SECS(1))            
         CALL RADDZ0 (LON005,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,550) A05,ES,                  & 
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),              & 
     &            X005,Y005                                            
  550    FORMAT ('0INITIALIZATION PARAMETERS (MERCATOR',              & 
     &           ' PROJECTION)'/                                       &
     &           ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/    &
     &           ' ECCENTRICITY SQUARED         =',F12.9/              &
     &           ' LATITUDE OF TRUE SCALE       = ',A1,2I3,F7.3/       &
     &           ' CENTRAL LONGITUDE            = ',A1,2I3,F7.3/       &
     &           ' FALSE EASTING                =',F12.2,' METERS'/    &
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A05                                                 
         DATA(2) = ES                                                  
         SWITCH(5) = ZONE                                              
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                       .  POLAR STEREOGRAPHIC  .                      
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 6) THEN                                            
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(6).NE.0 .AND. SWITCH(6).EQ.ZONE) RETURN            
         SWITCH(6) = 0                                                 
         A06 = AZ                                                      
         E06 = EZ                                                      
         ES = ESZ                                                      
         E4 = E4Z                                                      
         LON006 = PAKRZ0 (DATA(5))                                     
         SAVE = DATA(6)                                                
         LATC = PAKRZ0 (SAVE)                                          
         X006 = DATA(7)                                                
         Y006 = DATA(8)                                                
         FAC = ONE                                                     
         IF (SAVE .LT. ZERO) FAC =-ONE                                 
         IND06 = 0                                                     
         IF (DABS(SAVE) .NE. NINTYD) THEN                              
            IND06 = 1                                                  
            CON1 = FAC * LATC                                          
            SINPHI = DSIN (CON1)                                       
            COSPHI = DCOS (CON1)                                       
            MCS = MSFNZ0 (E06,SINPHI,COSPHI)                           
            TCS = TSFNZ0 (E06,CON1,SINPHI)                             
         END IF                                                        
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON006,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         CALL RADDZ0 (LATC,SGNA(2),DEGS(2),MINS(2),SECS(2))            
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,650) A06,ES,               &  
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),            & 
     &            X006,Y006                                            
  650    FORMAT ('0INITIALIZATION PARAMETERS (POLAR STEREOGRAPHIC',   & 
     &           ' PROJECTION)'/                                       &
     &           ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/   & 
     &           ' ECCENTRICITY SQUARED         =',F12.9/             & 
     &           ' LONGITUDE OF Y-AXIS          = ',A1,2I3,F7.3/      & 
     &           ' LATITUDE OF TRUE SCALE       = ',A1,2I3,F7.3/      & 
     &           ' FALSE EASTING                =',F12.2,' METERS'/   & 
     &           ' FALSE NORTHING               =',F12.2,' METERS')   
         DATA(1) = A06                                                 
         DATA(2) = ES                                                  
         SWITCH(6) = ZONE                                              
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                            .  POLYCONIC  .                           
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 7) THEN                                            
!                                                                      
  700    IERROR = 0                                                    
         IF (SWITCH(7).NE.0 .AND. SWITCH(7).EQ.ZONE) RETURN            
         SWITCH(7) = 0                                                 
         A07 = AZ                                                      
         E07 = EZ                                                      
         ES07 = ESZ                                                    
         E007 = E0Z                                                    
         E107 = E1Z                                                    
         E207 = E2Z                                                    
         E307 = E3Z                                                    
         LON007 = PAKRZ0 (DATA(5))                                     
         LAT007 = PAKRZ0 (DATA(6))                                     
         X007 = DATA(7)                                                
         Y007 = DATA(8)                                                
         ML007 = MLFNZ0 (E007,E107,E207,E307,LAT007)                   
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON007,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         CALL RADDZ0 (LAT007,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,750) A07,ES07,             &  
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),            & 
     &            X007,Y007                                            
  750    FORMAT ('0INITIALIZATION PARAMETERS (POLYCONIC',            &  
     &           ' PROJECTION)'/                                     &  
     &           ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/  &  
     &           ' ECCENTRICITY SQUARED         =',F12.9/            &  
     &           ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/     &  
     &           ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/     &  
     &           ' FALSE EASTING                =',F12.2,' METERS'/  &  
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A07                                                 
         DATA(2) = ES07                                                
         SWITCH(7) = ZONE                                              
!                                                                      
!     LIST STATE PLANE INITIALIZATION PARAMETERS IF NECESSARY          
!                                                                      
         IF (ISYS .EQ. 2) THEN                                         
            IPPARM = ITEMP                                             
            IF (IERROR .NE. 0) RETURN                                  
            IF (IPPARM .EQ. 0) WRITE (IPPLUN,470) ZONE, PNAME          
            SWITCH(2) = ZONE                                           
            RETURN                                                     
         END IF                                                        
!                                                                      
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                        .  EQUIDISTANT CONIC  .                       
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 8) THEN                                            
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(8).NE.0 .AND. SWITCH(8).EQ.ZONE) RETURN            
         SWITCH(8) = 0                                                 
         A08 = AZ                                                      
         E = EZ                                                        
         ES = ESZ                                                      
         E008 = E0Z                                                    
         E108 = E1Z                                                    
         E208 = E2Z                                                    
         E308 = E3Z                                                    
         LAT1 = PAKRZ0 (DATA(3))                                       
         LAT2 = PAKRZ0 (DATA(4))                                       
         IF (DABS(LAT1+LAT2) .LT. EPSLN) THEN                          
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,840)                      
  840       FORMAT ('0ERROR PJ08Z0'/                                  & 
     &              ' EQUAL LATITUDES FOR ST. PARALLELS ON OPPOSITE', & 
     &              ' SIDES OF EQUATOR')                               
            IERROR = 081                                               
            RETURN                                                     
         END IF                                                        
         LON008 = PAKRZ0 (DATA(5))                                     
         LAT0 = PAKRZ0 (DATA(6))                                       
         X008 = DATA(7)                                                
         Y008 = DATA(8)                                                
         SINPHI = DSIN (LAT1)                                          
         COSPHI = DCOS (LAT1)                                          
         MS1 = MSFNZ0 (E,SINPHI,COSPHI)                                
         ML1 = MLFNZ0 (E008,E108,E208,E308,LAT1)                       
         IND = 0                                                       
         IF (DATA(9) .NE. ZERO) THEN                                   
            IND = 1                                                    
            SINPHI = DSIN (LAT2)                                       
            COSPHI = DCOS (LAT2)                                       
            MS2 = MSFNZ0 (E,SINPHI,COSPHI)                             
            ML2 = MLFNZ0 (E008,E108,E208,E308,LAT2)                    
            IF (DABS(LAT1-LAT2) .GE. EPSLN) THEN                       
               NS08 = (MS1 - MS2) / (ML2 - ML1)                        
            ELSE                                                       
               NS08 = SINPHI                                           
            END IF                                                     
         ELSE                                                          
            NS08 = SINPHI                                              
         END IF                                                        
         GL = ML1 + MS1 / NS08                                         
         ML0 = MLFNZ0 (E008,E108,E208,E308,LAT0)                       
         RH008 = A08 * (GL - ML0)                                      
!                                                                      
!    LIST RESULTS OF PARAMETER INITIALIZATION.                         
!                                                                      
         CALL RADDZ0 (LAT1,SGNA(1),DEGS(1),MINS(1),SECS(1))            
         CALL RADDZ0 (LAT2,SGNA(2),DEGS(2),MINS(2),SECS(2))            
         CALL RADDZ0 (LON008,SGNA(3),DEGS(3),MINS(3),SECS(3))          
         CALL RADDZ0 (LAT0,SGNA(4),DEGS(4),MINS(4),SECS(4))            
         IF (IND .NE. 0) THEN                                          
            IF (IPPARM .EQ. 0) WRITE (IPPLUN,850) A08,ES,     &         
     &               (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,4),  &        
     &               X008,Y008                                         
  850       FORMAT ('0INITIALIZATION PARAMETERS (EQUIDISTANT CONIC',   &
     &              ' PROJECTION)'/                                    &
     &              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/ &
     &              ' ECCENTRICITY SQUARED         =',F12.9/           &
     &              ' LATITUDE OF 1ST ST. PARALLEL = ',A1,2I3,F7.3/    &
     &              ' LATITUDE OF 2ND ST. PARALLEL = ',A1,2I3,F7.3/    &
     &              ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/    &
     &              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/    &
     &              ' FALSE EASTING                =',F12.2,' METERS'/ &
     &              ' FALSE NORTHING               =',F12.2,' METERS') 
         ELSE                                                          
            IF (IPPARM .EQ. 0) WRITE (IPPLUN,860) A08,ES,           &   
     &                SGNA(1),DEGS(1),MINS(1),SECS(1),              &   
     &               (SGNA(I),DEGS(I),MINS(I),SECS(I),I=3,4),       &   
     &               X008,Y008                                         
  860       FORMAT ('0INITIALIZATION PARAMETERS (EQUIDISTANT CONIC',   & 
     &              ' PROJECTION)'/                                    & 
     &              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/ & 
     &              ' ECCENTRICITY SQUARED         =',F12.9/           & 
     &              ' LATITUDE OF ST. PARALLEL     = ',A1,2I3,F7.3/    & 
     &              ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/    & 
     &              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/    & 
     &              ' FALSE EASTING                =',F12.2,' METERS'/ & 
     &              ' FALSE NORTHING               =',F12.2,' METERS') 
         END IF                                                        
         DATA(1) = A08                                                 
         DATA(2) = ES                                                  
         SWITCH(8) = ZONE                                              
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                       .  TRANSVERSE MERCATOR  .                      
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 9) THEN                                            
!                                                                      
  900    IERROR = 0                                                    
         IF (DATA(1).NE.0.0D0.AND.DATA(1).NE.SAVE) SWITCH(9) = 0       
         IF (SWITCH(9).NE.0 .AND. SWITCH(9).EQ.ZONE) RETURN            
         SWITCH(9) = 0                                                 
         SAVE = DATA(1)                                                
         A09 = AZ                                                      
         E09 = EZ                                                      
         ES09 = ESZ                                                    
         E009 = E0Z                                                    
         E109 = E1Z                                                    
         E209 = E2Z                                                    
         E309 = E3Z                                                    
         KS009 = DATA(3)                                               
         LON009 = PAKRZ0 (DATA(5))                                     
         LAT009 = PAKRZ0 (DATA(6))                                     
         X009 = DATA(7)                                                
         Y009 = DATA(8)                                                
         ML009 = A09 * MLFNZ0 (E009,E109,E209,E309,LAT009)             
         IND09 = 1                                                     
         ESP = ES09                                                    
         IF (E09 .GE. TOL09) THEN                                      
            IND09 = 0                                                  
            ESP = ES09 / (ONE - ES09)                                  
         END IF                                                        
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON009,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         CALL RADDZ0 (LAT009,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,950) A09,ES09,KS009, &        
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),     &       
     &            X009,Y009                                            
  950    FORMAT ('0INITIALIZATION PARAMETERS (TRANSVERSE MERCATOR', &   
     &           ' PROJECTION)'/                                    &   
     &           ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/ &   
     &           ' ECCENTRICITY SQUARED         =',F12.9/           &   
     &           ' SCALE FACTOR AT C. MERIDIAN  =',F9.6/            &   
     &           ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/    &   
     &           ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/    &   
     &           ' FALSE EASTING                =',F12.2,' METERS'/ &   
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A09                                                 
         DATA(2) = ES09                                                
         SWITCH(9) = ZONE                                              
!                                                                      
!     LIST UTM PROJECTION INITIALIZATION PARAMETERS IF NECESSARY       
!                                                                      
         IF (ISYS .EQ. 1) THEN                                         
            IPPARM = ITEMP                                             
            BUFFL(1) = A09                                             
            BUFFL(2) = ES09                                            
            ZONE = KEEPZN                                              
            SWITCH(9) = ZONE                                           
            IF (IERROR .NE. 0) RETURN                                  
            IF (IPPARM .EQ. 0) WRITE (IPPLUN,960) ZONE,BUFFL(1),       &
     &            BUFFL(2),BUFFL(3),                                   &
     &            SGNA(1),DEGS(1),MINS(1),SECS(1),                     & 
     &            BUFFL(7),BUFFL(8)                                    
  960          FORMAT ('0INITIALIZATION PARAMETERS (U T M PROJECTION)'/ & 
     &            ' ZONE = ',I3/                                       &
     &            ' SEMI-MAJOR AXIS OF ELLIPSOID = ',F12.2,' METERS'/  &
     &            ' ECCENTRICITY SQUARED         = ',F18.15/           &
     &            ' SCALE FACTOR AT C. MERIDIAN  = ',F9.6/             &
     &            ' LONGITUDE OF CENTRAL MERIDIAN= ',A1,2I3,F7.3/      &
     &            ' FALSE EASTING                = ',F12.2,' METERS'/  &
     &            ' FALSE NORTHING               = ',F12.2,' METERS')  
            SWITCH(1) = ZONE                                           
            RETURN                                                     
         END IF                                                        
!                                                                      
!     LIST STATE PLANE INITIALIZATION PARAMETERS IF NECESSARY          
!                                                                      
         IF (ISYS .EQ. 2) THEN                                         
            IPPARM = ITEMP                                             
            IF (IERROR .NE. 0) RETURN                                  
            IF (IPPARM .EQ. 0) WRITE (IPPLUN,470) ZONE, PNAME          
            SWITCH(2) = ZONE                                           
            RETURN                                                     
         END IF                                                        
!                                                                      
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                          .  STEREOGRAPHIC  .                         
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 10) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(10).NE.0 .AND. SWITCH(10).EQ.ZONE) RETURN          
         SWITCH(10) = 0                                                
         A10 = AZZ                                                     
         LON010 = PAKRZ0 (DATA(5))                                     
         LAT010 = PAKRZ0 (DATA(6))                                     
         X010 = DATA(7)                                                
         Y010 = DATA(8)                                                
         SINP10 = DSIN (LAT010)                                        
         COSP10 = DCOS (LAT010)                                        
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON010,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         CALL RADDZ0 (LAT010,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,1050) A10,          &         
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),    &         
     &            X010,Y010                                            
 1050    FORMAT ('0INITIALIZATION PARAMETERS (STEREOGRAPHIC',         & 
     &           ' PROJECTION)'/                                      & 
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/   & 
     &           ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/      & 
     &           ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/      & 
     &           ' FALSE EASTING                =',F12.2,' METERS'/   & 
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A10                                                 
         SWITCH(10) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                   .  LAMBERT AZIMUTHAL EQUAL-AREA  .                 
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 11) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(11).NE.0 .AND. SWITCH(11).EQ.ZONE) RETURN          
         SWITCH(11) = 0                                                
         A11 = AZZ                                                     
         LON011 = PAKRZ0 (DATA(5))                                     
         LAT011 = PAKRZ0 (DATA(6))                                     
         X011 = DATA(7)                                                
         Y011 = DATA(8)                                                
         SINP11 = DSIN (LAT011)                                        
         COSP11 = DCOS (LAT011)                                        
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON011,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         CALL RADDZ0 (LAT011,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,1150) A11,                &
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),       &   
     &            X011,Y011                                            
 1150 FORMAT ('0INITIALIZATION PARAMETERS (LAMBERT AZIMUTHAL EQUAL-AREA'  & 
     &          ,' PROJECTION)'/                                     &  
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/ &   
     &           ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/    &   
     &           ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/    &   
     &           ' FALSE EASTING                =',F12.2,' METERS'/ &   
     &           ' FALSE NORTHING               =',F12.2,' METERS') 
         DATA(1) = A11                                                 
         SWITCH(11) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                      .  AZIMUTHAL EQUIDISTANT  .                     
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 12) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(12).NE.0 .AND. SWITCH(12).EQ.ZONE) RETURN          
         SWITCH(12) = 0                                                
         A12 = AZZ                                                     
         LON012 = PAKRZ0 (DATA(5))                                     
         LAT012 = PAKRZ0 (DATA(6))                                     
         X012 = DATA(7)                                                
         Y012 = DATA(8)                                                
         SINP12 = DSIN (LAT012)                                        
         COSP12 = DCOS (LAT012)                                        
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON012,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         CALL RADDZ0 (LAT012,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,1250) A12,                  & 
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),            & 
     &            X012,Y012                                            
 1250    FORMAT ('0INITIALIZATION PARAMETERS (AZIMUTHAL EQUIDISTANT', & 
     &           ' PROJECTION)'/                                      & 
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/   & 
     &           ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/      & 
     &           ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/      & 
     &           ' FALSE EASTING                =',F12.2,' METERS'/   & 
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A12                                                 
         SWITCH(12) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                            .  GNOMONIC  .                            
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 13) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(13).NE.0 .AND. SWITCH(13).EQ.ZONE) RETURN          
         SWITCH(13) = 0                                                
         A13 = AZZ                                                     
         LON013 = PAKRZ0 (DATA(5))                                     
         LAT013 = PAKRZ0 (DATA(6))                                     
         X013 = DATA(7)                                                
         Y013 = DATA(8)                                                
         SINP13 = DSIN (LAT013)                                        
         COSP13 = DCOS (LAT013)                                        
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON013,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         CALL RADDZ0 (LAT013,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,1350) A13,                &   
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),          &   
     &            X013,Y013                                            
 1350    FORMAT ('0INITIALIZATION PARAMETERS (GNOMONIC',             &  
     &           ' PROJECTION)'/                                     &  
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/  &  
     &           ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/     &  
     &           ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/     &  
     &           ' FALSE EASTING                =',F12.2,' METERS'/  &  
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A13                                                 
         SWITCH(13) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                          .  ORTHOGRAPHIC  .                          
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 14) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(14).NE.0 .AND. SWITCH(14).EQ.ZONE) RETURN          
         SWITCH(14) = 0                                                
         A14 = AZZ                                                     
         LON014 = PAKRZ0 (DATA(5))                                     
         LAT014 = PAKRZ0 (DATA(6))                                     
         X014 = DATA(7)                                                
         Y014 = DATA(8)                                                
         SINP14 = DSIN (LAT014)                                        
         COSP14 = DCOS (LAT014)                                        
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON014,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         CALL RADDZ0 (LAT014,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,1450) A14,                 &  
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),           &  
     &            X014,Y014                                            
 1450    FORMAT ('0INITIALIZATION PARAMETERS (ORTHOGRAPHIC',         &  
     &           ' PROJECTION)'/                                     &  
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/  &  
     &           ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/     &  
     &           ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/     &  
     &           ' FALSE EASTING                =',F12.2,' METERS'/  &  
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A14                                                 
         SWITCH(14) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!             .  GENERAL VERTICAL NEAR-SIDE PERSPECTIVE  .             
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 15) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(15).NE.0 .AND. SWITCH(15).EQ.ZONE) RETURN          
         SWITCH(15) = 0                                                
         A15 = AZZ                                                     
         P = ONE + DATA(3) / A15                                       
         LON015 = PAKRZ0 (DATA(5))                                     
         LAT015 = PAKRZ0 (DATA(6))                                     
         X015 = DATA(7)                                                
         Y015 = DATA(8)                                                
         SINP15 = DSIN (LAT015)                                        
         COSP15 = DCOS (LAT015)                                        
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON015,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         CALL RADDZ0 (LAT015,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,1550) A15,DATA(3),         &  
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),           &  
     &            X015,Y015                                            
 1550 FORMAT ('0INITIALIZATION PARAMETERS (GENERAL VERTICAL NEAR-SIDE',&
     &           ' PERSPECTIVE PROJECTION)'/                           &
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/    &
     &           ' HEIGHT OF PERSPECTIVE POINT'/                       &
     &           ' ABOVE SPHERE                 =',F12.2,' METERS'/    &
     &           ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/       &
     &           ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/       &
     &           ' FALSE EASTING                =',F12.2,' METERS'/    &
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A15                                                 
         SWITCH(15) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                           .  SINUSOIDAL  .                           
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 16) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(16).NE.0 .AND. SWITCH(16).EQ.ZONE) RETURN          
         SWITCH(16) = 0                                                
         A16 = AZZ                                                     
         LON016 = PAKRZ0 (DATA(5))                                     
         X016 = DATA(7)                                                
         Y016 = DATA(8)                                                
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON016,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,1650) A16,                 &  
     &            SGNA(1),DEGS(1),MINS(1),SECS(1),                   &  
     &            X016,Y016                                            
 1650    FORMAT ('0INITIALIZATION PARAMETERS (SINUSOIDAL',           &  
     &           ' PROJECTION)'/                                     &  
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/  &  
     &           ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/     &  
     &           ' FALSE EASTING                =',F12.2,' METERS'/  &  
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A16                                                 
         SWITCH(16) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                         .  EQUIRECTANGULAR  .                        
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 17) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(17).NE.0 .AND. SWITCH(17).EQ.ZONE) RETURN          
         SWITCH(17) = 0                                                
         A17 = AZZ                                                     
         LAT1 = PAKRZ0 (DATA(6))                                       
         LON017 = PAKRZ0 (DATA(5))                                     
         X017 = DATA(7)                                                
         Y017 = DATA(8)                                                
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LAT1,SGNA(1),DEGS(1),MINS(1),SECS(1))            
         CALL RADDZ0 (LON017,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,1750) A17,                &   
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),          &   
     &            X017,Y017                                            
 1750 FORMAT ('0INITIALIZATION PARAMETERS (EQUIRECTANGULAR PROJECTION)'&
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/   & 
     &           ' LATITUDE OF TRUE SCALE       = ',A1,2I2,F7.3/      & 
     &           ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/      & 
     &           ' FALSE EASTING                =',F12.2,' METERS'/   & 
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A17                                                 
         SWITCH(17) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                       .  MILLER CYLINDRICAL  .                       
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 18) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(18).NE.0 .AND. SWITCH(18).EQ.ZONE) RETURN          
         SWITCH(18) = 0                                                
         A18 = AZZ                                                     
         LON018 = PAKRZ0 (DATA(5))                                     
         X018 = DATA(7)                                                
         Y018 = DATA(8)                                                
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON018,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,1850) A18,               &    
     &             SGNA(1),DEGS(1),MINS(1),SECS(1),                &    
     &             X018,Y018                                           
 1850    FORMAT ('0INITIALIZATION PARAMETERS (MILLER CYLINDRICAL',  &   
     &           ' PROJECTION)'/                                    &   
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/  &  
     &           ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/     &  
     &           ' FALSE EASTING                =',F12.2,' METERS'/  &  
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A18                                                 
         SWITCH(18) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                        .  VAN DER GRINTEN I  .                       
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 19) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(19).NE.0 .AND. SWITCH(19).EQ.ZONE) RETURN          
         SWITCH(19) = 0                                                
         A19 = AZZ                                                     
         LON019 = PAKRZ0 (DATA(5))                                     
         X019 = DATA(7)                                                
         Y019 = DATA(8)                                                
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON019,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,1950) A19,                 &  
     &             SGNA(1),DEGS(1),MINS(1),SECS(1),                  &  
     &             X019,Y019                                           
 1950    FORMAT ('0INITIALIZATION PARAMETERS (VAN DER GRINTEN I',    &  
     &           ' PROJECTION)'/                                     &  
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/  &  
     &           ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/     &  
     &           ' FALSE EASTING                =',F12.2,' METERS'/  &  
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A19                                                 
         SWITCH(19) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                    .  OBLIQUE MERCATOR (HOTINE)  .                   
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 20) THEN                                           
!                                                                      
 2000    IERROR = 0                                                    
         IF (SWITCH(20).NE.0 .AND. SWITCH(20).EQ.ZONE) RETURN          
         SWITCH(20) = 0                                                
         MODE = 0                                                      
         IF (DATA(13) .NE. ZERO) MODE = 1                              
         A = AZ                                                        
         E20 = EZ                                                      
         ES = ESZ                                                      
         KS0 = DATA(3)                                                 
         LAT0 = PAKRZ0 (DATA(6))                                       
         X020 = DATA(7)                                                
         Y020 = DATA(8)                                                
         SINPH0 = DSIN (LAT0)                                          
         COSPH0 = DCOS (LAT0)                                          
         CON = ONE - ES * SINPH0 * SINPH0                              
         COM = DSQRT (ONE - ES)                                        
         BL = DSQRT (ONE + ES * COSPH0 ** 4 / (ONE - ES))              
         AL = A * BL * KS0 * COM / CON                                 
         IF (DABS(LAT0).LT.EPSLN) TS0 = 1.0D0                          
         IF (DABS(LAT0).LT.EPSLN) D=1.0D0                              
         IF (DABS(LAT0).LT.EPSLN) EL=1.0D0                             
         IF (DABS(LAT0).GE.EPSLN) THEN                                 
            TS0 = TSFNZ0 (E20,LAT0,SINPH0)                             
            CON = DSQRT (CON)                                          
            D = BL * COM / (COSPH0 * CON)                              
            F = D + DSIGN (DSQRT (DMAX1 ((D * D - ONE), 0.0D0)) , LAT0)
            EL = F * TS0 ** BL                                         
         END IF                                                        
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,2050) A,ES,KS0               
 2050 FORMAT ('0INITIALIZATION PARAMETERS (OBLIQUE MERCATOR ''HOTINE'''&
     &           ' PROJECTION)'/                                      & 
     &           ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/   & 
     &           ' ECCENTRICITY SQUARED         =',F12.9/             & 
     &           ' SCALE AT CENTER              =',F12.9)              
         IF (MODE .NE. 0) THEN                                         
            ALPHA = PAKRZ0 (DATA(4))                                   
            LONC = PAKRZ0 (DATA(5))                                    
            G = HALF * (F - ONE / F)                                   
            GAMMA = ASINZ0 (DSIN (ALPHA) / D)                          
            LON020 = LONC - ASINZ0 (G * DTAN (GAMMA)) / BL             
!                                                                      
!     LIST INITIALIZATION PARAMETERS (CASE B).                         
!                                                                      
            CALL RADDZ0 (ALPHA,SGNA(1),DEGS(1),MINS(1),SECS(1))        
            CALL RADDZ0 (LONC,SGNA(2),DEGS(2),MINS(2),SECS(2))         
            CALL RADDZ0 (LAT0,SGNA(3),DEGS(3),MINS(3),SECS(3))         
            IF (IPPARM .EQ. 0) WRITE (IPPLUN,2060)               &       
     &               (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,3)           
 2060       FORMAT (' AZIMUTH OF CENTRAL LINE      = ',A1,2I3,F7.3/   &  
     &              ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/    & 
     &              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3)    
            CON = DABS (LAT0)                                          
            IF (CON.GT.EPSLN .AND. DABS(CON - HALFPI).GT.EPSLN) THEN   
               SINGAM = DSIN (GAMMA)                                   
               COSGAM = DCOS (GAMMA)                                   
               SINALF = DSIN (ALPHA)                                   
               COSALF = DCOS (ALPHA)                                   
               U0 = DSIGN((AL/BL)*DATAN(DSQRT(D*D-ONE)/COSALF),LAT0)   
               IF (IPPARM .EQ. 0) WRITE (IPPLUN,2080) X020,Y020        
               DATA(1) = A                                             
               DATA(2) = ES                                            
               SWITCH(20) = ZONE                                       
!                                                                      
!     LIST STATE PLANE INITIALIZATION PARAMETERS IF NECESSARY          
!                                                                      
               IF (ISYS .EQ. 2) THEN                                   
                  IPPARM = ITEMP                                       
                  IF (IERROR .NE. 0) RETURN                            
                  IF (IPPARM .EQ. 0) WRITE (IPPLUN,470) ZONE, PNAME    
                  SWITCH(2) = ZONE                                     
                  RETURN                                               
               END IF                                                  
!                                                                      
               RETURN                                                  
            ELSE                                                       
               IF (IPEMSG .EQ. 0) WRITE (IPELUN,2040)                  
 2040          FORMAT ('0ERROR PJ20Z0'/                  &               
     &                 ' INPUT DATA ERROR')                            
               IERROR = 201                                            
               RETURN                                                  
            END IF                                                     
         END IF                                                        
         LON1 = PAKRZ0 (DATA(9))                                       
         LAT1 = PAKRZ0 (DATA(10))                                      
         LON2 = PAKRZ0 (DATA(11))                                      
         LAT2 = PAKRZ0 (DATA(12))                                      
         SINPHI = DSIN (LAT1)                                          
         TS1 = TSFNZ0 (E20,LAT1,SINPHI)                                
         SINPHI = DSIN (LAT2)                                          
         TS2 = TSFNZ0 (E20,LAT2,SINPHI)                                
         H = TS1 ** BL                                                 
         L = TS2 ** BL                                                 
         F = EL / H                                                    
         G = HALF * (F - ONE / F)                                      
         J = (EL * EL - L * H) / (EL * EL + L * H)                     
         P = (L - H) / (L + H)                                         
         CALL RADDZ0 (LON2,SGNA(3),DEGS(3),MINS(3),SECS(3))            
         DLON = LON1 - LON2                                            
         IF (DLON .LT. -PI) LON2 = LON2 - 2.D0 * PI                    
         IF (DLON .GT.  PI) LON2 = LON2 + 2.D0 * PI                    
         DLON = LON1 - LON2                                            
         LON020 = HALF * (LON1 + LON2) - DATAN (J * DTAN (HALF * BL * & 
     &          DLON) / P) / BL                                        
         DLON = ADJLZ0 (LON1 - LON020)                                 
         GAMMA = DATAN (DSIN (BL * DLON) / G)                          
         ALPHA = ASINZ0 (D * DSIN (GAMMA))                             
         CALL RADDZ0 (LON1,SGNA(1),DEGS(1),MINS(1),SECS(1))            
         CALL RADDZ0 (LAT1,SGNA(2),DEGS(2),MINS(2),SECS(2))            
         CALL RADDZ0 (LAT2,SGNA(4),DEGS(4),MINS(4),SECS(4))            
         CALL RADDZ0 (LAT0,SGNA(5),DEGS(5),MINS(5),SECS(5))            
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,2070)       &                 
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,5)              
 2070    FORMAT (' LONGITUDE OF 1ST POINT       = ',A1,2I3,F7.3/  &  
     &           ' LATITUDE OF 1ST POINT        = ',A1,2I3,F7.3/  &    
     &           ' LONGITUDE OF 2ND POINT       = ',A1,2I3,F7.3/  &     
     &           ' LATITUDE OF 2ND POINT        = ',A1,2I3,F7.3/  &     
     &           ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3)       
         IF (DABS(LAT1 - LAT2) .LE. EPSLN) THEN                        
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,2040)                     
            IERROR = 202                                               
            RETURN                                                     
        ELSE                                                           
            CON = DABS (LAT1)                                          
         END IF                                                        
         IF (CON.LE.EPSLN .OR. DABS(CON - HALFPI).LE.EPSLN) THEN       
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,2040)                     
            IERROR = 202                                               
            RETURN                                                     
         ELSE                                                          
            IF (DABS(DABS(LAT0) - HALFPI) .LE. EPSLN) THEN             
               IF (IPEMSG .EQ. 0) WRITE (IPELUN,2040)                  
               IERROR = 202                                            
               RETURN                                                  
            END IF                                                     
         END IF                                                        
         SINGAM = DSIN (GAMMA)                                         
         COSGAM = DCOS (GAMMA)                                         
         SINALF = DSIN (ALPHA)                                         
         COSALF = DCOS (ALPHA)                                         
         U0 = DSIGN((AL/BL)*DATAN(DSQRT(D*D-ONE)/COSALF),LAT0)         
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,2080) X020,Y020              
 2080    FORMAT (' FALSE EASTING                =',F12.2,' METERS'/  &    
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A                                                   
         DATA(2) = ES                                                  
         SWITCH(20) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                       .       ROBINSON       .                       
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 21) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(21).NE.0 .AND. SWITCH(21).EQ.ZONE) RETURN          
         SWITCH(21) = 0                                                
         A21 = AZZ                                                     
         LON021 = PAKRZ0 (DATA(5))                                     
         X021 = DATA(7)                                                
         Y021 = DATA(8)                                                
         PR(1)=-0.062D0                                                
         XLR(1)=0.9986D0                                               
         PR(2)=0.D0                                                    
         XLR(2)=1.D0                                                   
         PR(3)=0.062D0                                                 
         XLR(3)=0.9986D0                                               
         PR(4)=0.124D0                                                 
         XLR(4)=0.9954D0                                               
         PR(5)=0.186D0                                                 
         XLR(5)=0.99D0                                                 
         PR(6)=0.248D0                                                 
         XLR(6)=0.9822D0                                               
         PR(7)=0.31D0                                                  
         XLR(7)=0.973D0                                                
         PR(8)=0.372D0                                                 
         XLR(8)=0.96D0                                                 
         PR(9)=0.434D0                                                 
         XLR(9)=0.9427D0                                               
         PR(10)=0.4958D0                                               
         XLR(10)=0.9216D0                                              
         PR(11)=0.5571D0                                               
         XLR(11)=0.8962D0                                              
         PR(12)=0.6176D0                                               
         XLR(12)=0.8679D0                                              
         PR(13)=0.6769D0                                               
         XLR(13)=0.835D0                                               
         PR(14)=0.7346D0                                               
         XLR(14)=0.7986D0                                              
         PR(15)=0.7903D0                                               
         XLR(15)=0.7597D0                                              
         PR(16)=0.8435D0                                               
         XLR(16)=0.7186D0                                              
         PR(17)=0.8936D0                                               
         XLR(17)=0.6732D0                                              
         PR(18)=0.9394D0                                               
         XLR(18)=0.6213D0                                              
         PR(19)=0.9761D0                                               
         XLR(19)=0.5722D0                                              
         PR(20)=1.0D0                                                  
         XLR(20)=0.5322D0                                              
         DO 2110 I=1,20                                                
 2110    XLR(I)=XLR(I) * 0.9858D0                                      
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON021,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,2150) A21,        &           
     &             SGNA(1),DEGS(1),MINS(1),SECS(1),          &          
     &             X021,Y021                                           
 2150    FORMAT ('0INITIALIZATION PARAMETERS (ROBINSON',            &   
     &           ' PROJECTION)'/                                     &  
     &           ' RADIUS OF SPHERE             =',F12.2,' METERS'/   & 
     &           ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/       &
     &           ' FALSE EASTING                =',F12.2,' METERS'/    &
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A21                                                 
         SWITCH(21) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!                      .  SPACE OBLIQUE MERCATOR  .                    
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 22) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(22).NE.0 .AND. SWITCH(22).EQ.ZONE) RETURN          
         SWITCH(22) = 0                                                
         A22 = AZ                                                      
         E = EZ                                                        
         ES22 = ESZ                                                    
         X022 = DATA(7)                                                
         Y022 = DATA(8)                                                
         LAND = IDINT(DATA(3)+TOL)                                     
         PATH = IDINT(DATA(4)+TOL)                                     
!                                                                      
!        CHECK IF LANDSAT NUMBER IS WITHIN RANGE 1 - 5                 
!                                                                      
         IF (LAND .GT. 0 .AND. LAND .LE. 5) THEN                       
            IF (LAND .LE. 3) LIMIT = 251                               
            IF (LAND .GE. 4) LIMIT = 233                               
         ELSE                                                          
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,2240) LAND, PATH          
            IERROR = 221                                               
            RETURN                                                     
         END IF                                                        
!                                                                      
!        CHECK IF PATH NUMBER IS WITHIN RANGE 1 - 251 FOR LANDSATS 1 - 
!        OR RANGE 1 - 233 FOR LANDSATS 4 - 5                           
!                                                                      
         IF (PATH .LE. 0 .OR. PATH .GT. LIMIT) THEN                    
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,2240) LAND, PATH          
 2240       FORMAT ('0ERROR PJ22Z0'/                                  & 
     &              ' LANDSAT NUMBER ',I2,' AND / OR PATH NUMBER ',I4, &
     &              ' ARE OUT OF RANGE')                               
            IERROR = 221                                               
            RETURN                                                     
         END IF                                                        
         P1=1440.0D0                                                   
         IF (LAND.LE.3) THEN                                           
            P2=103.2669323D0                                           
            ALF=99.092D0*DG1                                           
         ELSE                                                          
            P2=98.8841202D0                                            
            ALF=98.20D0*DG1                                            
         END IF                                                        
         SA=DSIN(ALF)                                                  
         CA=DCOS(ALF)                                                  
         IF (DABS(CA).LT.1.D-9) CA=1.D-9                               
         ESC=ES22*CA*CA                                                
         ESS=ES22*SA*SA                                                
         W=((ONE-ESC)/(ONE-ES22))**TWO-ONE                             
         Q=ESS/(ONE-ES22)                                              
         T=(ESS*(TWO-ES22))/(ONE-ES22)**TWO                            
         U=ESC/(ONE-ES22)                                              
         XJ=(ONE-ES22)**3                                              
         P22=P2/P1                                                     
!                                                                      
!        COMPUTE FOURIER COEFFICIENTS.  LAM IS CURRENT VALUE OF        
!        LAMBDA DOUBLE-PRIME.                                          
!                                                                      
         LAM=0                                                         
         CALL SERAZ0 (FB,FA2,FA4,FC1,FC3,LAM)                          
         SUMA2=FA2                                                     
         SUMA4=FA4                                                     
         SUMB=FB                                                       
         SUMC1=FC1                                                     
         SUMC3=FC3                                                     
         DO 2210 I=9,81,18                                             
         LAM=DBLE(I)                                                   
         CALL SERAZ0 (FB,FA2,FA4,FC1,FC3,LAM)                          
         SUMA2=SUMA2+4.0D0*FA2                                         
         SUMA4=SUMA4+4.0D0*FA4                                         
         SUMB=SUMB+4.0D0*FB                                            
         SUMC1=SUMC1+4.0D0*FC1                                         
         SUMC3=SUMC3+4.0D0*FC3                                         
 2210    CONTINUE                                                      
         DO 2220 I=18,72,18                                            
         LAM=DBLE(I)                                                   
         CALL SERAZ0 (FB,FA2,FA4,FC1,FC3,LAM)                          
         SUMA2=SUMA2+TWO*FA2                                           
         SUMA4=SUMA4+TWO*FA4                                           
         SUMB=SUMB+TWO*FB                                              
         SUMC1=SUMC1+TWO*FC1                                           
         SUMC3=SUMC3+TWO*FC3                                           
 2220    CONTINUE                                                      
         LAM=90.0D0                                                    
         CALL SERAZ0 (FB,FA2,FA4,FC1,FC3,LAM)                          
         SUMA2=SUMA2+FA2                                               
         SUMA4=SUMA4+FA4                                               
         SUMB=SUMB+FB                                                  
         SUMC1=SUMC1+FC1                                               
         SUMC3=SUMC3+FC3                                               
!                                                                      
!        THESE ARE THE VALUES OF FOURIER CONSTANTS.                    
!                                                                      
         A2=SUMA2/30.D0                                                
         A4=SUMA4/60.D0                                                
         B=SUMB/30.D0                                                  
         C1=SUMC1/15.D0                                                
         C3=SUMC3/45.D0                                                
!                                                                      
!        LIST RESULTS OF PARAMETER INITIALIZATION.                     
!                                                                      
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,2250) A22,ES22,LAND,PATH,   & 
     &                                          X022,Y022              
 2250    FORMAT ('0INITIALIZATION PARAMETERS (SPACE OBL. MERCATOR', &   
     &           ' PROJECTION)'/                                    &   
     &           ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/ &   
     &           ' ECCENTRICITY SQUARED         =',F12.9/           &   
     &           ' LANDSAT NO.                  = ',I3/             &   
     &           ' PATH                         = ',I5/             &   
     &           ' FALSE EASTING                =',F12.2,' METERS'/ &   
     &           ' FALSE NORTHING               =',F12.2,' METERS'/)   
         DATA(1) = A22                                                 
         DATA(2) = ES22                                                
         SWITCH(22) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!             .  INITIALIZATION OF PROJECTION PARAMETERS  .            
!                                                                      
!          .  MODIFIED-STEREOGRAPHIC CONFORMAL (FOR ALASKA)  .         
! .....................................................................
!                                                                      
      IF (ISYS .EQ. 23) THEN                                           
!                                                                      
         IERROR = 0                                                    
         IF (SWITCH(23).NE.0 .AND. SWITCH(23).EQ.ZONE) RETURN          
         SWITCH(23) = 0                                                
         A23 = AZ                                                      
         EC2 = 0.6768657997291094D-02                                  
         EC  = DSQRT (EC2)                                             
         N=6                                                           
         LON023 = -152.0D0*DG1                                         
         LAT023 = 64.0D0*DG1                                           
         X023 = DATA(7)                                                
         Y023 = DATA(8)                                                
         ACOEF(1)=0.9945303D0                                          
         ACOEF(2)=0.0052083D0                                          
         ACOEF(3)=0.0072721D0                                          
         ACOEF(4)=-0.0151089D0                                         
         ACOEF(5)=0.0642675D0                                          
         ACOEF(6)=0.3582802D0                                          
         BCOEF(1)=0.0D0                                                
         BCOEF(2)=-.0027404D0                                          
         BCOEF(3)=0.0048181D0                                          
         BCOEF(4)=-0.1932526D0                                         
         BCOEF(5)=-0.1381226D0                                         
         BCOEF(6)=-0.2884586D0                                         
         ESPHI=EC*DSIN(LAT023)                                         
         CHIO=TWO*DATAN(DTAN((HALFPI+LAT023)/TWO)*((ONE-ESPHI)/    &    
     &       (ONE+ESPHI))**(EC/TWO)) - HALFPI                          
         SCHIO=DSIN(CHIO)                                              
         CCHIO=DCOS(CHIO)                                              
!                                                                      
!     LIST RESULTS OF PARAMETER INITIALIZATION.                        
!                                                                      
         CALL RADDZ0 (LON023,SGNA(1),DEGS(1),MINS(1),SECS(1))          
         CALL RADDZ0 (LAT023,SGNA(2),DEGS(2),MINS(2),SECS(2))          
         IF (IPPARM .EQ. 0) WRITE (IPPLUN,2350) A23,EC2,              & 
     &            (SGNA(I),DEGS(I),MINS(I),SECS(I),I=1,2),            & 
     &            X023,Y023                                            
 2350    FORMAT ('0INITIALIZATION PARAMETERS (MOD. STEREOGRAPHIC',    & 
     &           ' CONFORMAL PROJECTION, ALASKA)'/                    & 
     &           ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/   & 
     &           ' ECCENTRICITY SQUARED         =',F12.9/             & 
     &           ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/      & 
     &           ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/      & 
     &           ' FALSE EASTING                =',F12.2,' METERS'/   & 
     &           ' FALSE NORTHING               =',F12.2,' METERS')    
         DATA(1) = A23                                                 
         SWITCH(23) = ZONE                                             
         RETURN                                                        
      END IF                                                           
!                                                                      
!     INITIALIZATION OF PROJECTION COMPLETED                           
!                                                                      
      END                                                              
!                   PJ01Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                              *  U T M  *                             
! *********************************************************************
!                                                                      
      SUBROUTINE PJ01Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC, FWD, INV                                        
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /TOGGLE/ SWITCH                                           
      PARAMETER (FWD=0, INV=1)                                         
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(1) .NE. 0) GO TO 140                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT ('0ERROR PJ01Z0'/                     &                 
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 013                                                  
         RETURN                                                        
  140    CALL PJ09Z0 (GEOG,PROJ,FWD)                                   
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(1) .NE. 0) GO TO 160                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
         IERROR = 014                                                  
         RETURN                                                        
  160    CALL PJ09Z0 (PROJ,GEOG,INV)                                   
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ02Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                           *  STATE PLANE  *                          
! *********************************************************************
!                                                                      
      SUBROUTINE PJ02Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23), ITYPE                                      
      INTEGER*2 INDIC, FWD, INV                                        
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ02/ ITYPE                                              
      COMMON /TOGGLE/ SWITCH                                           
!                                                                      
      PARAMETER (FWD=0, INV=1)                                         
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(2) .EQ. 0) THEN                                    
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,250)                      
  250       FORMAT ('0ERROR PJ02Z0'/                    &               
     &              ' PROJECTION WAS NOT INITIALIZED')                 
            IERROR = 023                                               
            RETURN                                                     
         END IF                                                        
!                                                                      
!     TRANSVERSE MERCATOR PROJECTION                                   
!                                                                      
         IF (ITYPE .EQ. 1) THEN                                        
            CALL PJ09Z0 (GEOG,PROJ,FWD)                                
         END IF                                                        
!                                                                      
!     LAMBERT CONFORMAL PROJECTION                                     
!                                                                      
         IF (ITYPE .EQ. 2) THEN                                        
            CALL PJ04Z0 (GEOG,PROJ,FWD)                                
         END IF                                                        
!                                                                      
!     POLYCONIC PROJECTION                                             
!                                                                      
         IF (ITYPE .EQ. 3) THEN                                        
            CALL PJ07Z0 (GEOG,PROJ,FWD)                                
         END IF                                                        
!                                                                      
!     OBLIQUE MERCATOR PROJECTION                                      
!                                                                      
         IF (ITYPE .EQ. 4) THEN                                        
            CALL PJ20Z0 (GEOG,PROJ,FWD)                                
         END IF                                                        
!                                                                      
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(2) .EQ. 0) THEN                                    
            IF (IPEMSG .EQ. 0) WRITE (IPELUN,250)                      
            IERROR = 025                                               
            RETURN                                                     
         END IF                                                        
!                                                                      
!     TRANSVERSE MERCATOR PROJECTION                                   
!                                                                      
         IF (ITYPE .EQ. 1) THEN                                        
            CALL PJ09Z0 (PROJ,GEOG,INV)                                
         END IF                                                        
!                                                                      
!     LAMBERT CONFORMAL PROJECTION                                     
!                                                                      
         IF (ITYPE .EQ. 2) THEN                                        
            CALL PJ04Z0 (PROJ,GEOG,INV)                                
         END IF                                                        
!                                                                      
!     POLYCONIC PROJECTION                                             
!                                                                      
         IF (ITYPE .EQ. 3) THEN                                        
            CALL PJ07Z0 (PROJ,GEOG,INV)                                
         END IF                                                        
!                                                                      
!     OBLIQUE MERCATOR PROJECTION                                      
!                                                                      
         IF (ITYPE .EQ. 4) THEN                                        
            CALL PJ20Z0 (PROJ,GEOG,INV)                                
         END IF                                                        
!                                                                      
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ03Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                    *  ALBERS CONICAL EQUAL AREA  *                   
! *********************************************************************
!                                                                      
      SUBROUTINE PJ03Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,E,ES,LAT1,LAT2,LON0,LAT0,X0,Y0,NS,C,RH0 ******
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ03/ A,LON0,X0,Y0,C,E,ES,NS,RH0                         
      COMMON /TOGGLE/ SWITCH                                           
      DATA TOL /1.0D-7/                                                
      DATA HALFPI /1.5707963267948966D0/                               
      DATA ZERO,HALF,ONE /0.0D0,0.5D0,1.0D0/                           
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(3) .NE. 0) GO TO 220                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT ('0ERROR PJ03Z0'/                                      & 
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 033                                                  
         RETURN                                                        
  220    SINPHI = DSIN (GEOG(2))                                       
         COSPHI = DCOS (GEOG(2))                                       
         QS = QSFNZ0 (E,SINPHI,COSPHI)                                 
         RH = A * DSQRT (C - NS * QS) / NS                             
         THETA = NS * ADJLZ0 (GEOG(1) - LON0)                          
         PROJ(1) = X0 + RH * DSIN (THETA)                              
         PROJ(2) = Y0 + RH0 - RH * DCOS (THETA)                        
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(3) .NE. 0) GO TO 240                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
         IERROR = 034                                                  
         RETURN                                                        
  240    X = PROJ(1) - X0                                              
         Y = RH0 - PROJ(2) + Y0                                        
         RH = DSIGN (DSQRT (X * X + Y * Y) , NS)                       
         THETA = ZERO                                                  
         CON = DSIGN (ONE , NS)                                        
         IF (RH .NE. ZERO) THETA = DATAN2 (CON * X , CON * Y)          
         CON = RH * NS / A                                             
         QS = (C - CON * CON) / NS                                     
         IF (E .LT. TOL) GO TO 260                                     
         CON = ONE - HALF * (ONE - ES) * DLOG ((ONE - E) /            & 
     &         (ONE + E)) / E                                          
         IF ((DABS(CON) - DABS(QS)) .GT. TOL) GO TO 260                
         GEOG(2) = DSIGN (HALFPI , QS)                                 
         GO TO 280                                                     
  260    GEOG(2) = PHI1Z0 (E,QS)                                       
         IF (IERROR .EQ. 0) GO TO 280                                  
         IERROR = 035                                                  
         RETURN                                                        
  280    GEOG(1) = ADJLZ0 (THETA / NS + LON0)                          
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ04Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                     *  LAMBERT CONFORMAL CONIC  *                    
! *********************************************************************
!                                                                      
      SUBROUTINE PJ04Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,E,ES,LAT1,LAT2,LON0,LAT0,X0,Y0,NS,F,RH0 ******
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ04/ A,LON0,X0,Y0,E,F,NS,RH0                            
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN /1.0D-10/                                             
      DATA ZERO,ONE /0.0D0,1.0D0/                                      
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(4) .NE. 0) GO TO 200                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT ('0ERROR PJ04Z0'/                                       & 
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 043                                                  
         RETURN                                                        
  200    CON = DABS (DABS (GEOG(2)) - HALFPI)                          
         IF (CON .GT. EPSLN) GO TO 220                                 
         CON = GEOG(2) * NS                                            
         IF (CON .GT. ZERO) GO TO 210                                  
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2030)                        
 2030    FORMAT ('0ERROR PJ04Z0'/                                       &
     &           ' POINT CANNOT BE PROJECTED')                         
         IERROR = 044                                                  
         RETURN                                                        
  210    RH = ZERO                                                     
         GO TO 230                                                     
  220    SINPHI = DSIN (GEOG(2))                                       
         TS = TSFNZ0 (E,GEOG(2),SINPHI)                                
         RH = A * F * TS ** NS                                         
  230    THETA = NS * ADJLZ0 (GEOG(1) - LON0)                          
         PROJ(1) = X0 + RH * DSIN (THETA)                              
         PROJ(2) = Y0 + RH0 - RH * DCOS (THETA)                        
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(4) .NE. 0) GO TO 240                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
         IERROR = 045                                                  
         RETURN                                                        
  240    X = PROJ(1) - X0                                              
         Y = RH0 - PROJ(2) + Y0                                        
         RH = DSIGN (DSQRT (X*X + Y*Y) , NS)                           
         THETA = ZERO                                                  
         CON = DSIGN (ONE , NS)                                        
         IF (RH .NE. ZERO) THETA = DATAN2 (CON * X , CON * Y)          
         IF (RH.NE.ZERO .OR. NS.GT.ZERO) GO TO 250                     
         GEOG(2) = - HALFPI                                            
         GO TO 260                                                     
  250    CON = ONE / NS                                                
         TS = (RH / (A * F)) ** CON                                    
         GEOG(2) = PHI2Z0 (E,TS)                                       
         IF (IERROR .EQ. 0) GO TO 260                                  
         IERROR = 046                                                  
         RETURN                                                        
  260    GEOG(1) = ADJLZ0 (THETA / NS + LON0)                          
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ05Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                            *  MERCATOR  *                            
! *********************************************************************
!                                                                      
      SUBROUTINE PJ05Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,E,ES,LON0,X0,Y0,NS,F,RH0,LAT1,M1 *************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ05/ A,LON0,X0,Y0,E,M1                                  
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN /1.0D-10/                                             
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(5) .NE. 0) GO TO 220                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ05Z0'/                                     & 
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 052                                                  
         RETURN                                                        
  220    IF (DABS(DABS(GEOG(2)) - HALFPI) .GT. EPSLN) GO TO 240        
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT ('0ERROR PJ05Z0'/                                     & 
     &           ' TRANSFORMATION CANNOT BE COMPUTED AT THE POLES')    
         IERROR = 053                                                  
         RETURN                                                        
  240    SINPHI = DSIN (GEOG(2))                                       
         TS = TSFNZ0 (E,GEOG(2),SINPHI)                                
         PROJ(1) = X0 + A * M1 * ADJLZ0 (GEOG(1) - LON0)               
         PROJ(2) = Y0 - A * M1 * DLOG (TS)                             
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(5) .NE. 0) GO TO 260                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 054                                                  
         RETURN                                                        
  260    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         TS = DEXP (- Y / (A * M1))                                    
         GEOG(2) = PHI2Z0 (E,TS)                                       
         IF (IERROR .EQ. 0) GO TO 280                                  
         IERROR = 055                                                  
         RETURN                                                        
  280    GEOG(1) = ADJLZ0 (LON0 + X / (A * M1))                        
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ06Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                       *  POLAR STEREOGRAPHIC  *                      
! *********************************************************************
!                                                                      
      SUBROUTINE PJ06Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23),IND                                         
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,E,ES,LON0,LATC,X0,Y0,E4,MCS,TCS,FAC,IND ******
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ06/ A,LON0,X0,Y0,E,E4,FAC,MCS,TCS,IND                  
      COMMON /TOGGLE/ SWITCH                                           
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                            
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(6) .NE. 0) GO TO 220                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ06Z0'/                                 &     
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 062                                                  
         RETURN                                                        
  220    CON1 = FAC * ADJLZ0 (GEOG(1) - LON0)                          
         CON2 = FAC * GEOG(2)                                          
         SINPHI = DSIN (CON2)                                          
         TS = TSFNZ0 (E,CON2,SINPHI)                                   
         IF (IND .EQ. 0) GO TO 240                                     
         RH = A * MCS * TS / TCS                                       
         GO TO 260                                                     
  240    RH = TWO * A * TS / E4                                        
  260    PROJ(1) = X0 + FAC * RH * DSIN (CON1)                         
         PROJ(2) = Y0 - FAC * RH * DCOS (CON1)                         
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(6) .NE. 0) GO TO 320                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 063                                                  
         RETURN                                                        
  320    X = FAC * (PROJ(1) - X0)                                      
         Y = FAC * (PROJ(2) - Y0)                                      
         RH = DSQRT (X * X + Y * Y)                                    
         IF (IND .EQ. 0) GO TO 340                                     
         TS = RH * TCS / (A * MCS)                                     
         GO TO 360                                                     
  340    TS = RH * E4 / (TWO * A)                                      
  360    GEOG(2) = FAC * PHI2Z0 (E,TS)                                 
         IF (IERROR .EQ. 0) GO TO 380                                  
         IERROR = 064                                                  
         RETURN                                                        
  380    IF (RH .NE. ZERO) GO TO 400                                   
         GEOG(1) = FAC * LON0                                          
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  400    GEOG(1) = ADJLZ0 (FAC * DATAN2 (X , -Y) + LON0)               
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ07Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                            *  POLYCONIC  *                           
! *********************************************************************
!                                                                      
      SUBROUTINE PJ07Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,E,ES,LON0,LAT0,X0,Y0,E0,E1,E2,ML0 ************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ07/ A,LON0,X0,Y0,E,E0,E1,E2,E3,ES,ML0                  
      COMMON /TOGGLE/ SWITCH                                           
      DATA TOL /1.0D-7/                                                
      DATA ZERO,ONE /0.0D0,1.0D0/                                      
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(7) .NE. 0) GO TO 220                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ07Z0'/                                    &  
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 072                                                  
         RETURN                                                        
  220    CON = ADJLZ0 (GEOG(1) - LON0)                                 
         IF (DABS(GEOG(2)) .GT. TOL) GO TO 240                         
         PROJ(1) = X0 + A * CON                                        
         PROJ(2) = Y0 - A * ML0                                        
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
  240    SINPHI = DSIN (GEOG(2))                                       
         COSPHI = DCOS (GEOG(2))                                       
         ML = MLFNZ0 (E0,E1,E2,E3,GEOG(2))                             
         MS = MSFNZ0 (E,SINPHI,COSPHI)                                 
         CON = CON * SINPHI                                            
         PROJ(1) = X0 + A * MS * DSIN (CON) / SINPHI                   
         PROJ(2) = Y0 + A * (ML - ML0 + MS * (ONE - DCOS(CON)) / SINPHI)
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(7) .NE. 0) GO TO 320                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 073                                                  
         RETURN                                                        
  320    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         AL = ML0 + Y / A                                              
         IF (DABS (AL) .GT. TOL) GO TO 340                             
         GEOG(1) = X / A + LON0                                        
         GEOG(2) = ZERO                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  340    B = AL * AL + (X / A) ** 2                                    
         CALL PHI4Z0 (ES,E0,E1,E2,E3,AL,B,C,GEOG(2))                   
         IF (IERROR .EQ. 0) GO TO 360                                  
         IERROR = 074                                                  
         RETURN                                                        
  360    GEOG(1) = ADJLZ0 (ASINZ0 (X * C / A) / DSIN (GEOG(2)) + LON0) 
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ08Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                        *  EQUIDISTANT CONIC  *                       
! *********************************************************************
!                                                                      
      SUBROUTINE PJ08Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! ** PARAMETERS * A,E,ES,LAT1,LAT2,LON0,LAT0,X0,Y0,E0,E1,E2,E3,NS,GL,RH
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ08/ A,LON0,X0,Y0,E0,E1,E2,E3,GL,NS,RH0                 
      COMMON /TOGGLE/ SWITCH                                           
      DATA ZERO,ONE /0.0D0,1.0D0/                                      
      DATA EPSLN /1.0D-10/                                             
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(8) .NE. 0) GO TO 300                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2030)                        
 2030    FORMAT ('0ERROR PJ08Z0'/                                    &  
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 083                                                  
         RETURN                                                        
  300    ML = MLFNZ0 (E0,E1,E2,E3,GEOG(2))                             
         RH = A * (GL - ML)                                            
         THETA = NS * ADJLZ0 (GEOG(1) - LON0)                          
         PROJ(1) = X0 + RH * DSIN (THETA)                              
         PROJ(2) = Y0 + RH0 - RH * DCOS (THETA)                        
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(8) .NE. 0) GO TO 320                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2030)                        
         IERROR = 084                                                  
         RETURN                                                        
  320    X = PROJ(1) - X0                                              
         Y = RH0 - PROJ(2) + Y0                                        
         RH = DSIGN (DSQRT (X * X + Y * Y) , NS)                       
         THETA = ZERO                                                  
         CON = DSIGN (ONE , NS)                                        
         IF (RH .NE. ZERO) THETA = DATAN2 (CON * X , CON * Y)          
         ML = GL - RH / A                                              
         GEOG(2) = PHI3Z0 (ML,E0,E1,E2,E3)                             
         IF (IERROR .EQ. 0) GO TO 340                                  
         IERROR = 085                                                  
         RETURN                                                        
  340    GEOG(1) = ADJLZ0 (LON0 + THETA / NS)                          
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ09Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                       *  TRANSVERSE MERCATOR  *                      
! *********************************************************************
!                                                                      
      SUBROUTINE PJ09Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23),I,IND,NIT                                   
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS ** A,E,ES,KS0,LON0,LAT0,X0,Y0,E0,E1,E2,E3,ESP,ML0,IND
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ09/ A,LON0,X0,Y0,ES,ESP,E0,E1,E2,E3,KS0,LAT0,ML0,IND   
      COMMON /TOGGLE/ SWITCH                                           
      DATA ZERO,HALF,ONE,TWO,THREE /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/     
      DATA FOUR,FIVE,SIX,EIGHT,NINE /4.0D0,5.0D0,6.0D0,8.0D0,9.0D0/    
      DATA HALFPI /1.5707963267948966D0/                               
      DATA TEN /10.0D0/                                                
      DATA EPSLN,NIT /1.0D-10,6/                                       
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(9) .NE. 0) GO TO 220                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ09Z0'/                                    &  
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 092                                                  
         RETURN                                                        
  220    DLON = ADJLZ0 (GEOG(1) - LON0)                                
         LAT = GEOG(2)                                                 
         IF (IND .EQ. 0) GO TO 240                                     
         COSPHI = DCOS (LAT)                                           
         B = COSPHI * DSIN (DLON)                                      
         IF (DABS(DABS(B) - ONE) .GT. EPSLN) GO TO 230                 
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT ('0ERROR PJ09Z0'/                                     & 
     &           ' POINT PROJECTS INTO INFINITY')                      
         IERROR = 093                                                  
         RETURN                                                        
  230    PROJ(1) = HALF * A * KS0 * DLOG ((ONE + B) / (ONE - B)) + X0  
         CON = DACOS (COSPHI * DCOS (DLON) / DSQRT (ONE - B * B))      
         IF (LAT .LT. ZERO) CON =-CON                                  
         PROJ(2) = A * KS0 * (CON - LAT0) + Y0                         
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
!                                                                      
  240    SINPHI = DSIN (LAT)                                           
         COSPHI = DCOS (LAT)                                           
         AL = COSPHI * DLON                                            
         ALS = AL * AL                                                 
         C = ESP * COSPHI * COSPHI                                     
         TQ = DTAN (LAT)                                               
         T = TQ * TQ                                                   
         N = A / DSQRT (ONE - ES * SINPHI * SINPHI)                    
         ML = A * MLFNZ0 (E0,E1,E2,E3,LAT)                             
         PROJ(1) = KS0 * N * AL * (ONE + ALS / SIX * (ONE - T + C +     &
     &             ALS / 20.0D0 * (FIVE - 18.0D0 * T + T * T + 72.0D0 * &
     &             C - 58.0D0 * ESP))) + X0                            
         PROJ(2) = KS0 *(ML - ML0 + N * TQ *(ALS *(HALF + ALS / 24.0D0*  &
     &             (FIVE - T + NINE * C + FOUR * C * C + ALS / 30.0D0 *  & 
     &             (61.0D0 - 58.0D0 * T + T * T + 600.0D0 * C -           &    
     &             330.0D0 * ESP))))) + Y0                             
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(9) .NE. 0) GO TO 320                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 094                                                  
         RETURN                                                        
  320    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         IF (IND .EQ. 0) GO TO 340                                     
         F = DEXP (X / (A * KS0))                                      
         G = HALF * (F - ONE / F)                                      
         TEMP = LAT0 + Y / (A * KS0)                                   
         H = DCOS (TEMP)                                               
         CON = DSQRT ((ONE - H * H) / (ONE + G * G))                   
         GEOG(2) = ASINZ0 (CON)                                        
         IF (TEMP .LT. ZERO) GEOG(2) =-GEOG(2)                         
         IF (G.NE.ZERO .OR. H.NE.ZERO) GO TO 330                       
         GEOG(1) = LON0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  330    GEOG(1) = ADJLZ0 (DATAN2 (G,H) + LON0)                        
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
!                                                                      
  340    CON = (ML0 + Y / KS0) / A                                     
         PHI = CON                                                     
         DO 360 I = 1,NIT                                              
         DPHI = ((CON + E1 * DSIN (TWO * PHI) - E2 * DSIN (FOUR * PHI)  &
     &          + E3 * DSIN (SIX * PHI)) / E0) - PHI                   
         PHI = PHI + DPHI                                              
         IF (DABS(DPHI) .LE. EPSLN) GO TO 380                          
  360    CONTINUE                                                      
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2030) NIT                    
 2030    FORMAT ('0ERROR PI09Z0' /                                      & 
     &           ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS')
         IERROR = 095                                                  
         RETURN                                                        
  380    IF (DABS(PHI) .LT. HALFPI) GO TO 400                          
         GEOG(2) = DSIGN (HALFPI , Y)                                  
         GEOG(1) = LON0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  400    SINPHI = DSIN (PHI)                                           
         COSPHI = DCOS (PHI)                                           
         TANPHI = DTAN (PHI)                                           
         C = ESP * COSPHI * COSPHI                                     
         CS = C * C                                                    
         T = TANPHI * TANPHI                                           
         TS = T * T                                                    
         CON = ONE - ES * SINPHI * SINPHI                              
         N = A / DSQRT (CON)                                           
         R = N * (ONE - ES) / CON                                      
         D = X / (N * KS0)                                             
         DS = D * D                                                    
         GEOG(2) = PHI - (N * TANPHI * DS / R) * (HALF - DS / 24.0D0 *   &
     &             (FIVE + THREE * T + TEN * C - FOUR * CS - NINE * ESP  &
     &             - DS / 30.0D0 * (61.0D0 + 90.0D0 * T + 298.0D0 * C +  &
     &             45.0D0 * TS - 252.0D0 * ESP - THREE * CS)))         
         GEOG(1) = ADJLZ0 (LON0 + (D * (ONE - DS / SIX * (ONE + TWO *    &
     &             T + C - DS / 20.0D0 * (FIVE - TWO * C + 28.0D0 * T -  &
     &             THREE * CS + EIGHT * ESP + 24.0D0 * TS))) / COSPHI))
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ10Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                          *  STEREOGRAPHIC  *                         
! *********************************************************************
!                                                                      
      SUBROUTINE PJ10Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 ****************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ10/ A,LON0,X0,Y0,COSPH0,LAT0,SINPH0                    
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN /1.0D-10/                                             
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                            
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(10) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ10Z0'/                              &          
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 102                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         SINPHI = DSIN (GEOG(2))                                       
         COSPHI = DCOS (GEOG(2))                                       
         COSLON = DCOS (LON)                                           
         G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                
         IF (DABS(G + ONE) .GT. EPSLN) GO TO 140                       
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT ('0ERROR PJ10Z0'/                                      &  
     &           ' POINT PROJECTS INTO INFINITY')                      
         IERROR = 103                                                  
         RETURN                                                        
  140    KSP = TWO / (ONE + G)                                         
         PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                  
         PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI *  &
     &             COSLON)                                             
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(10) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 104                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         RH = DSQRT (X * X + Y * Y)                                    
         Z = TWO * DATAN (RH / (TWO * A))                              
         SINZ = DSIN (Z)                                               
         COSZ = DCOS (Z)                                               
         GEOG(1) = LON0                                                
         IF (DABS(RH) .GT. EPSLN) GO TO 240                            
         GEOG(2) = LAT0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  240    GEOG(2) = ASINZ0 (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)     
         CON = DABS (LAT0) - HALFPI                                    
         IF (DABS (CON) .GT. EPSLN) GO TO 260                          
         IF (LAT0 .LT. ZERO) GO TO 250                                 
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  250    GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  260    CON = COSZ - SINPH0 * DSIN (GEOG(2))                          
         IF (DABS(CON).LT.EPSLN.AND.DABS(X).LT.EPSLN) RETURN           
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH))) 
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ11Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                   *  LAMBERT AZIMUTHAL EQUAL-AREA  *                 
! *********************************************************************
!                                                                      
      SUBROUTINE PJ11Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 ****************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ11/ A,LON0,X0,Y0,COSPH0,LAT0,SINPH0                    
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN /1.0D-10/                                             
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                            
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(11) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ11Z0'/                              &       
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 112                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         SINPHI = DSIN (GEOG(2))                                       
         COSPHI = DCOS (GEOG(2))                                       
         COSLON = DCOS (LON)                                           
         G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                
         IF (G .NE. -ONE) GO TO 140                                    
         CON = TWO * A                                                 
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020) CON                    
 2020    FORMAT (' POINT PROJECTS INTO A CIRCLE OF RADIUS =',F12.2,   & 
     &           ' METERS')                                            
         IERROR = 113                                                  
         RETURN                                                        
  140    KSP = DSQRT (TWO / (ONE + G))                                 
         PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                  
         PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI * &
     &             COSLON)                                             
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(11) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 114                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         RH = DSQRT (X * X + Y * Y)                                    
         CON = RH / (TWO * A)                                          
         IF (CON .LE. ONE) GO TO 230                                   
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2030)                        
 2030    FORMAT ('0ERROR PJ11Z0'/                                 &     
     &           ' INPUT DATA ERROR')                                  
         IERROR = 115                                                  
         RETURN                                                        
  230    Z = TWO * ASINZ0 (CON)                                        
         SINZ = DSIN (Z)                                               
         COSZ = DCOS (Z)                                               
         GEOG(1) = LON0                                                
         IF (DABS(RH) .GT. EPSLN) GO TO 240                            
         GEOG(2) = LAT0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  240    GEOG(2) = ASINZ0 (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)     
         CON = DABS (LAT0) - HALFPI                                    
         IF (DABS (CON) .GT. EPSLN) GO TO 260                          
         IF (LAT0 .LT. ZERO) GO TO 250                                 
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  250    GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  260    CON = COSZ - SINPH0 * DSIN (GEOG(2))                          
         IF (CON .EQ. ZERO) RETURN                                     
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH))) 
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ12Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                      *  AZIMUTHAL EQUIDISTANT  *                     
! *********************************************************************
!                                                                      
      SUBROUTINE PJ12Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 ****************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ12/ A,LON0,X0,Y0,COSPH0,LAT0,SINPH0                    
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN /1.0D-10/                                             
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                            
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(12) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ12Z0'/                                     & 
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 122                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         SINPHI = DSIN (GEOG(2))                                       
         COSPHI = DCOS (GEOG(2))                                       
         COSLON = DCOS (LON)                                           
         G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                
         IF (DABS(DABS(G) - ONE) .GE. EPSLN) GO TO 140                 
         KSP = ONE                                                     
         IF (G .GE. ZERO) GO TO 160                                    
         CON = TWO * HALFPI * A                                        
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020) CON                    
 2020    FORMAT (' POINT PROJECTS INTO CIRCLE OF RADIUS =',F12.2,     & 
     &           ' METERS')                                            
         IERROR = 123                                                  
         RETURN                                                        
  140    Z = DACOS (G)                                                 
         KSP = Z / DSIN (Z)                                            
  160    PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                  
         PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI *  &
     &             COSLON)                                             
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(12) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 124                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         RH = DSQRT (X * X + Y * Y)                                    
         IF (RH .LE. (TWO * HALFPI * A)) GO TO 230                     
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2030)                        
 2030    FORMAT ('0ERROR PJ12Z0'/                    &                   
     &           ' INPUT DATA ERROR')                                  
         IERROR = 125                                                  
         RETURN                                                        
  230    Z = RH / A                                                    
         SINZ = DSIN (Z)                                               
         COSZ = DCOS (Z)                                               
         GEOG(1) = LON0                                                
         IF (DABS(RH) .GT. EPSLN) GO TO 240                            
         GEOG(2) = LAT0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  240    GEOG(2) = ASINZ0 (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)     
         CON = DABS (LAT0) - HALFPI                                    
         IF (DABS (CON) .GT. EPSLN) GO TO 260                          
         IF (LAT0 .LT. ZERO) GO TO 250                                 
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  250    GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  260    CON = COSZ - SINPH0 * DSIN (GEOG(2))                          
         IF (DABS(CON).LT.EPSLN.AND.DABS(X).LT.EPSLN) RETURN           
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH))) 
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ13Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                            *  GNOMONIC  *                            
! *********************************************************************
!                                                                      
      SUBROUTINE PJ13Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 ****************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ13/ A,LON0,X0,Y0,COSPH0,LAT0,SINPH0                    
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN /1.0D-10/                                             
      DATA ZERO,ONE /0.0D0,1.0D0/                                      
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(13) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ13Z0'/                                     & 
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 132                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         SINPHI = DSIN (GEOG(2))                                       
         COSPHI = DCOS (GEOG(2))                                       
         COSLON = DCOS (LON)                                           
         G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                
         IF (G .GT. ZERO) GO TO 140                                    
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT (' POINT PROJECTS INTO INFINITY')                      
         IERROR = 133                                                  
         RETURN                                                        
  140    KSP = ONE / G                                                 
         PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                  
         PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI * &  
     &             COSLON)                                             
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(13) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 134                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         RH = DSQRT (X * X + Y * Y)                                    
         Z = DATAN (RH / A)                                            
         SINZ = DSIN (Z)                                               
         COSZ = DCOS (Z)                                               
         GEOG(1) = LON0                                                
         IF (DABS(RH) .GT. EPSLN) GO TO 240                            
         GEOG(2) = LAT0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  240    GEOG(2) = ASINZ0 (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)     
         CON = DABS (LAT0) - HALFPI                                    
         IF (DABS (CON) .GT. EPSLN) GO TO 260                          
         IF (LAT0 .LT. ZERO) GO TO 250                                 
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  250    GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  260    CON = COSZ - SINPH0 * DSIN (GEOG(2))                          
         IF (DABS(CON).LT.EPSLN.AND.DABS(X).LT.EPSLN) RETURN           
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH))) 
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ14Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                          *  ORTHOGRAPHIC  *                          
! *********************************************************************
!                                                                      
      SUBROUTINE PJ14Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 ****************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ14/ A,LON0,X0,Y0,COSPH0,LAT0,SINPH0                    
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN /1.0D-10/                                             
      DATA ZERO,ONE /0.0D0,1.0D0/                                      
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(14) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ14Z0'/                                 &      
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 142                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         SINPHI = DSIN (GEOG(2))                                       
         COSPHI = DCOS (GEOG(2))                                       
         COSLON = DCOS (LON)                                           
         G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                
         KSP = ONE                                                     
         IF (G.GT.ZERO .OR. DABS(G).LE.EPSLN) GO TO 140                
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT (' POINT CANNOT BE PROJECTED')                         
         IERROR = 143                                                  
         RETURN                                                        
  140    PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                  
         PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI *  &
     &             COSLON)                                             
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(14) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 144                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         RH = DSQRT (X * X + Y * Y)                                    
         IF (RH .LE. A) GO TO 230                                      
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2030)                        
 2030    FORMAT ('0ERROR PJ14Z0'/                               &        
     &           ' INPUT DATA ERROR')                                  
         IERROR = 145                                                  
         RETURN                                                        
  230    Z = ASINZ0 (RH / A)                                           
         SINZ = DSIN (Z)                                               
         COSZ = DCOS (Z)                                               
         GEOG(1) = LON0                                                
         IF (DABS(RH) .GT. EPSLN) GO TO 240                            
         GEOG(2) = LAT0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  240    GEOG(2) = ASINZ0 (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)     
         CON = DABS (LAT0) - HALFPI                                    
         IF (DABS (CON) .GT. EPSLN) GO TO 260                          
         IF (LAT0 .LT. ZERO) GO TO 250                                 
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  250    GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  260    CON = COSZ - SINPH0 * DSIN (GEOG(2))                          
         IF (DABS(CON).LT.EPSLN.AND.DABS(X).LT.EPSLN) RETURN           
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH))) 
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ15Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!              *  GENERAL VERTICAL NEAR-SIDE PERSPECTIVE  *            
! *********************************************************************
!                                                                      
      SUBROUTINE PJ15Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,P,LON0,LAT0,X0,Y0,SINPH0,COSPH0 **************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ15/ A,LON0,X0,Y0,COSPH0,LAT0,P,SINPH0                  
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN /1.0D-10/                                             
      DATA ZERO,ONE /0.0D0,1.0D0/                                      
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(15) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ15Z0'/                                     & 
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 152                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         SINPHI = DSIN (GEOG(2))                                       
         COSPHI = DCOS (GEOG(2))                                       
         COSLON = DCOS (LON)                                           
         G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                
         IF (G .GE. (ONE / P)) GO TO 140                               
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT (' POINT CANNOT BE PROJECTED')                         
         IERROR = 153                                                  
         RETURN                                                        
  140    KSP = (P - ONE) / (P - G)                                     
         PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                  
         PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI * &
     &             COSLON)                                             
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(15) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 154                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         RH = DSQRT (X * X + Y * Y)                                    
         R = RH / A                                                    
         CON = P - ONE                                                 
         COM = P + ONE                                                 
         IF (R .LE. DSQRT (CON / COM)) GO TO 230                       
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2030)                        
 2030    FORMAT ('0ERROR PJ15Z0'/                                    &  
     &           ' INPUT DATA ERROR')                                  
         IERROR = 155                                                  
         RETURN                                                        
  230    SINZ = (P - DSQRT (ONE - R * R * COM / CON)) /              &  
     &          (CON / R + R / CON)                                    
         Z = ASINZ0 (SINZ)                                             
         SINZ = DSIN (Z)                                               
         COSZ = DCOS (Z)                                               
         GEOG(1) = LON0                                                
         IF (DABS(RH) .GT. EPSLN) GO TO 240                            
         GEOG(2) = LAT0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  240    GEOG(2) = ASINZ0 (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)     
         CON = DABS (LAT0) - HALFPI                                    
         IF (DABS (CON) .GT. EPSLN) GO TO 260                          
         IF (LAT0 .LT. ZERO) GO TO 250                                 
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  250    GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  260    CON = COSZ - SINPH0 * DSIN (GEOG(2))                          
         IF (DABS(CON).LT.EPSLN.AND.DABS(X).LT.EPSLN) RETURN           
         GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH))) 
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ16Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                           *  SINUSOIDAL  *                           
! *********************************************************************
!                                                                      
      SUBROUTINE PJ16Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,LON0,X0,Y0 ***********************************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ16/ A,LON0,X0,Y0                                       
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN /1.0D-10/                                             
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(16) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ16Z0'/                                 &     
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 162                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         PROJ(1) = X0 + A * LON * DCOS (GEOG(2))                       
         PROJ(2) = Y0 + A * GEOG(2)                                    
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(16) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 163                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         GEOG(2) = Y / A                                               
         IF (DABS(GEOG(2)) .LE. HALFPI) GO TO 230                      
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT ('0ERROR PJ16Z0'/                                  &    
     &           ' INPUT DATA ERROR')                                  
         IERROR = 164                                                  
         RETURN                                                        
  230    CON = DABS (GEOG(2)) - HALFPI                                 
         IF (DABS (CON) .GT. EPSLN) GO TO 240                          
         GEOG(1) = LON0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  240    GEOG(1) = ADJLZ0 (LON0 + X / (A * DCOS (GEOG(2))))            
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ17Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                  *  EQUIRECTANGULAR   *                              
! *********************************************************************
!                                                                      
      SUBROUTINE PJ17Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,LON0,X0,Y0,LAT1 ******************************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ17/ A,LON0,X0,Y0,LAT1                                  
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(17) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ17Z0'/                               &       
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 172                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         PROJ(1) = X0 + A * LON * DCOS(LAT1)                           
         PROJ(2) = Y0 + A * GEOG(2)                                    
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(17) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 173                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         GEOG(2) = Y / A                                               
         IF (DABS(GEOG(2)) .LE. HALFPI) GO TO 240                      
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2020)                        
 2020    FORMAT ('0ERROR PJ17Z0'/                                    &  
     &           ' INPUT DATA ERROR')                                  
         IERROR = 174                                                  
         RETURN                                                        
  240    GEOG(1) = ADJLZ0 (LON0 + X / (A * DCOS(LAT1) ))               
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ18Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                       *  MILLER CYLINDRICAL  *                       
! *********************************************************************
!                                                                      
      SUBROUTINE PJ18Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,LON0,X0,Y0 ***********************************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ18/ A,LON0,X0,Y0                                       
      COMMON /TOGGLE/ SWITCH                                           
      DATA FORTPI /0.78539816339744833D0/                              
      DATA ZERO,ONEQ,TWOH /0.0D0,1.25D0,2.5D0/                         
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(18) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ18Z0'/                       &               
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 182                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         PROJ(1) = X0 + A * LON                                        
         PROJ(2) = Y0 + A * DLOG (DTAN (FORTPI + GEOG(2) / TWOH)) * ONE
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(18) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 183                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         GEOG(1) = ADJLZ0 (LON0 + X / A)                               
         GEOG(2) = TWOH * DATAN (DEXP (Y / A / ONEQ)) - FORTPI * TWOH  
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ19Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                        *  VAN DER GRINTEN I  *                       
! *********************************************************************
!                                                                      
      SUBROUTINE PJ19Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,LON0,X0,Y0 ***********************************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ19/ A,LON0,X0,Y0                                       
      COMMON /TOGGLE/ SWITCH                                           
      DATA PI /3.14159265358979323846D0/                               
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN/1.0D-10/                                              
      DATA ZERO,HALF,ONE,TWO,THREE/0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/      
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(19) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ19Z0'/                                     & 
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 192                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         LAT = GEOG(2)                                                 
         IF (DABS(LAT) .GT. EPSLN) GO TO 140                           
         PROJ(1) = X0 + A * LON                                        
         PROJ(2) = Y0                                                  
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
  140    THETA = ASINZ0 (DMIN1(DABS (LAT /HALFPI),ONE))                
         IF (DABS(LON).GT.EPSLN.AND.DABS(DABS(LAT)-HALFPI).GT.EPSLN)  & 
     &       GO TO 160                                                 
         PROJ(1) = X0                                                  
         PROJ(2) = Y0 + PI * A * DSIGN( DTAN (HALF * THETA), LAT)      
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
  160    AL = HALF * DABS (PI / LON - LON / PI)                        
         ASQ = AL * AL                                                 
         SINTHT = DSIN (THETA)                                         
         COSTHT = DCOS (THETA)                                         
         G = COSTHT / (SINTHT + COSTHT - ONE)                          
         GSQ = G * G                                                   
         M = G * (TWO / SINTHT - ONE)                                  
         MSQ = M * M                                                   
         CON = PI * A * (AL * (G - MSQ) + DSQRT (ASQ * (G - MSQ)**2 -  &
     &         (MSQ + ASQ) * (GSQ - MSQ))) / (MSQ + ASQ)               
         CON = DSIGN (CON , LON)                                       
         PROJ(1) = X0 + CON                                            
         CON = DABS (CON / (PI * A))                                   
         PROJ(2) = Y0 + DSIGN (PI * A * DSQRT (ONE - CON * CON -       &
     &             TWO * AL * CON) , LAT)                              
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
!     ALGORITHM DEVELOPED BY D.P. RUBINCAM, THE AMERICAN CARTOGRAPHER, 
!                1981, V. 8, NO. 2, P. 177-180.                        
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(19) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 193                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         CON = PI * A                                                  
         XX = X / CON                                                  
         YY = Y / CON                                                  
         XYS = XX * XX + YY * YY                                       
         C1 = -DABS(YY) * (ONE + XYS)                                  
         C2 = C1 - TWO * YY * YY + XX * XX                             
         C3 = -TWO * C1 + ONE + TWO * YY * YY + XYS*XYS                
         D = YY * YY / C3 + (TWO * C2 * C2 * C2/ C3/ C3/ C3 - 9.0D0 * C &
     &       * C2/ C3/ C3) / 27.0D0                                    
         A1 = (C1 - C2 * C2/ THREE/ C3)/ C3                            
         M1 = TWO * DSQRT(-A1/ THREE)                                  
         CON = ((THREE * D) / A1) / M1                                 
         IF (DABS(CON).GT.ONE) CON = DSIGN(ONE,CON)                    
         TH1 = DACOS(CON)/THREE                                        
         GEOG(2) = (-M1 * DCOS(TH1 + PI/ THREE) - C2/ THREE/ C3)        &
     &   * DSIGN(PI,Y)                                                 
         IF (DABS(XX).GE.EPSLN) GO TO 230                              
         GEOG(1) = LON0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  230    CONTINUE                                                      
         GEOG(1) = LON0 + PI * (XYS - ONE + DSQRT(ONE + TWO * (XX * XX  &
     &      - YY * YY) + XYS * XYS))/ TWO/ XX                          
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ20Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                    *  OBLIQUE MERCATOR (HOTINE)  *                   
! *********************************************************************
!                                                                      
      SUBROUTINE PJ20Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN                     
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,E,ES,KS0,ALPHA,LONC,LON1,LAT1,LON2,LAT2,LAT0 *
! ********************** X0,Y0,GAMMA,LON0,AL,BL,EL ********************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ20/ LON0,X0,Y0,AL,BL,COSALF,COSGAM,E,EL,SINALF,SINGAM,U
      COMMON /TOGGLE/ SWITCH                                           
      DATA PI /3.14159265358979323846D0/                               
      DATA HALFPI /1.5707963267948966D0/                               
      DATA TOL,EPSLN /1.0D-7,1.0D-10/                                  
      DATA ZERO,HALF,ONE /0.0D0,0.5D0,1.0D0/                           
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(20) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2050)                        
 2050    FORMAT ('0ERROR PJ20Z0'/                                &       
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 204                                                  
         RETURN                                                        
  220    SINPHI = DSIN (GEOG(2))                                       
         DLON = ADJLZ0 (GEOG(1) - LON0)                                
         VL = DSIN (BL * DLON)                                         
         IF (DABS(DABS(GEOG(2)) - HALFPI) .GT. EPSLN) GO TO 230        
         UL = SINGAM * DSIGN (ONE , GEOG(2))                           
         US = AL * GEOG(2) / BL                                        
         GO TO 250                                                     
  230    TS = TSFNZ0 (E,GEOG(2),SINPHI)                                
         Q = EL / TS ** BL                                             
         S = HALF * (Q - ONE / Q)                                      
         T = HALF * (Q + ONE / Q)                                      
         UL = (S * SINGAM - VL * COSGAM) / T                           
         CON = DCOS (BL * DLON)                                        
         IF (DABS(CON) .LT. TOL) GO TO 240                             
         US = AL * DATAN ((S * COSGAM + VL * SINGAM) / CON) / BL       
         IF (CON .LT. ZERO) US = US + PI * AL / BL                     
         GO TO 250                                                     
  240    US = AL * BL * DLON                                           
  250    IF (DABS(DABS(UL) - ONE) .GT. EPSLN) GO TO 255                
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2060)                        
 2060    FORMAT ('0ERROR PJ20Z0'/                                    &  
     &           ' POINT PROJECTS INTO INFINITY')                      
         IERROR = 205                                                  
         RETURN                                                        
  255    VS = HALF * AL * DLOG ((ONE - UL) / (ONE + UL)) / BL          
         US = US - U0                                                  
         PROJ(1) = X0 + VS * COSALF + US * SINALF                      
         PROJ(2) = Y0 + US * COSALF - VS * SINALF                      
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(20) .NE. 0) GO TO 280                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2050)                        
         IERROR = 206                                                  
         RETURN                                                        
  280    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         VS = X * COSALF - Y * SINALF                                  
         US = Y * COSALF + X * SINALF                                  
         US = US + U0                                                  
         Q = DEXP (- BL * VS / AL)                                     
         S = HALF * (Q - ONE / Q)                                      
         T = HALF * (Q + ONE / Q)                                      
         VL = DSIN (BL * US / AL)                                      
         UL = (VL * COSGAM + S * SINGAM) / T                           
         IF (DABS (DABS (UL) - ONE) .GE. EPSLN) GO TO 300              
         GEOG(1) = LON0                                                
         GEOG(2) = DSIGN (HALFPI , UL)                                 
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  300    CON = ONE / BL                                                
         TS = (EL / DSQRT ((ONE + UL) / (ONE - UL))) ** CON            
         GEOG(2) = PHI2Z0 (E,TS)                                       
         CON = DCOS (BL * US / AL)                                     
         LON = LON0 - DATAN2 ((S * COSGAM - VL * SINGAM) , CON) / BL   
         GEOG(1) = ADJLZ0 (LON)                                        
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ21Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                       *       ROBINSON       *                       
! *********************************************************************
!                                                                      
      SUBROUTINE PJ21Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN,IP1,NN              
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2),                    & 
     & PR(20),XLR(20)                                                  
! **** PARAMETERS **** A,LON0,X0,Y0 ***********************************
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ21/ A,LON0,X0,Y0,PR,XLR                                
      COMMON /TOGGLE/ SWITCH                                           
      DATA DG1 /0.01745329252D0/                                       
      DATA PI /3.14159265358979323846D0/                               
      DATA EPSLN /1.0D-10/                                             
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(21) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ21Z0'/                               &       
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 212                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
         P2=DABS(GEOG(2)/5.0D0/DG1)                                    
         IP1=IDINT(P2-EPSLN)                                           
!                                                                      
!        STIRLING'S INTERPOLATION FORMULA (USING 2ND DIFF.)            
!            USED WITH LOOKUP TABLE TO COMPUTE RECTANGULAR COORDINATES 
!            FROM LAT/LONG.                                            
!                                                                      
         P2=P2-DBLE(IP1)                                               
         X=A*(XLR(IP1+2)+P2*(XLR(IP1+3)-XLR(IP1+1))/2.0D0              &
     &     +P2*P2*(XLR(IP1+3)-2.0D0*XLR(IP1+2)+XLR(IP1+1))/2.0D0)*LON  
         Y=A*(PR(IP1+2)+P2*(PR(IP1+3)-PR(IP1+1))/2.0D0                 &
     &     +P2*P2*(PR(IP1+3)-2.0D0*PR(IP1+2)+PR(IP1+1))/2.0D0)*PI/2.0D0 & 
     &     *DSIGN(1.0D0,GEOG(2))                                       
         PROJ(1) = X0 + X                                              
         PROJ(2) = Y0 + Y                                              
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(21) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 213                                                  
         RETURN                                                        
  220    X = PROJ(1) - X0                                              
         Y = PROJ(2) - Y0                                              
         YY = 2.0D0 * Y / PI / A                                       
         PHID = YY * 90.0D0                                            
         P2 = DABS(PHID / 5.0D0)                                       
         IP1 = IDINT(P2 - EPSLN)                                       
         IF (IP1.EQ.0) IP1 = 1                                         
         NN = 0                                                        
!                                                                      
!     STIRLING'S INTERPOLATION FORMULA AS USED IN FORWARD TRANSFORMATIO
!     IS REVERSED FOR FIRST ESTIMATION OF LAT. FROM RECTANGULAR        
!     COORDINATES.  LAT. IS THEN ADJUSTED BY ITERATION UNTIL USE OF    
!     FORWARD SERIES PROVIDES CORRECT VALUE OF Y WITHIN TOLERANCE.     
!                                                                      
  230    U = PR(IP1 + 3) - PR(IP1 + 1)                                 
         V = PR(IP1 + 3) - 2.0D0 * PR(IP1 + 2) + PR(IP1 + 1)           
         T = 2.0D0 * (DABS(YY) - PR(IP1 + 2))/ U                       
         C = V / U                                                     
         P2 = T * (1.0D0 - C * T * (1.0D0 - 2.0D0 * C * T))            
         IF (P2.LT.0.0D0.AND.IP1.NE.1) GO TO 240                       
         PHID = DSIGN((P2 + DBLE(IP1)) * 5.0D0, Y)                     
  235    P2 = DABS(PHID / 5.0D0)                                       
         IP1 = IDINT(P2 - EPSLN)                                       
         P2 = P2 - DBLE(IP1)                                           
         Y1=A*(PR(IP1+2)+P2*(PR(IP1+3)-PR(IP1+1))/2.0D0                 & 
     &     +P2*P2*(PR(IP1+3)-2.0D0*PR(IP1+2)+PR(IP1+1))/2.0D0)*PI/2.0D0 & 
     &     * DSIGN(1.0D0,Y)                                            
         PHID = PHID - 180.0D0* (Y1 - Y) / PI / A                      
         NN = NN + 1                                                   
         IF (NN.LE.20) GO TO 237                                       
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,245)                         
         IERROR = 214                                                  
         RETURN                                                        
  237    IF (DABS(Y1 - Y).GT.0.00001D0) GO TO 235                      
         GO TO 250                                                     
  240    IP1 = IP1 - 1                                                 
         GO TO 230                                                     
  245    FORMAT ('0ERROR PJ21Z0'/                                     & 
     &           ' TOO MANY ITERATIONS FOR INVERSE ROBINSON')          
  250    GEOG(2) = PHID * DG1                                          
!                                                                      
!        CALCULATE LONG. USING FINAL LAT. WITH TRANSPOSED FORWARD      
!        STIRLING'S INTERPOLATION FORMULA.                             
!                                                                      
         GEOG(1)=LON0+X/A/(XLR(IP1+2)+P2*(XLR(IP1+3)-XLR(IP1+1))/2.0D0 & 
     &     +P2*P2*(XLR(IP1+3)-2.0D0*XLR(IP1+2)+XLR(IP1+1))/2.0D0)      
         GEOG(1) = ADJLZ0(GEOG(1))                                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ22Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!                      *  SPACE OBLIQUE MERCATOR  *                    
! *********************************************************************
!                                                                      
      SUBROUTINE PJ22Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN,PATH,LAND,NN,L      
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2)                      
! **** PARAMETERS **** A,E,ES,LON0,LATC,X0,Y0,MCS,TCS,FAC,IND *********
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /NORM/ Q,T,U,W,ES,P22,SA,CA,XJ                            
      COMMON /PJ22/ A,X0,Y0,A2,A4,B,C1,C3,LAND,PATH                    
      COMMON /TOGGLE/ SWITCH                                           
      DATA TOL /1.0D-7/                                                
      DATA DG1 /0.01745329252D0/                                       
      DATA PI /3.14159265358979323846D0/                               
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                            
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(22) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ22Z0'/                              &        
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 222                                                  
         RETURN                                                        
  220    IF (LAND.GE.4) GO TO 225                                      
         LON=GEOG(1)-128.87D0*DG1+PI*TWO/251.D0*DBLE(PATH)             
         GO TO 230                                                     
  225    LON=GEOG(1)-129.30D0*DG1+PI*TWO/233.D0*DBLE(PATH)             
  230    LAT=GEOG(2)                                                   
!                                                                      
!        TEST FOR LAT. AND LONG. APPROACHING 90 DEGREES.               
!                                                                      
         IF (LAT.GT.1.570796D0) LAT=1.570796D0                         
         IF (LAT.LT.-1.570796D0) LAT =-1.570796D0                      
         IF (LAT.GE.0) LAMPP=PI/TWO                                    
         IF (LAT.LT.0) LAMPP=1.5D0*PI                                  
         NN=0                                                          
  231    SAV=LAMPP                                                     
         L=0                                                           
         LAMTP=LON+P22*LAMPP                                           
         CL=DCOS(LAMTP)                                                
         IF (DABS(CL).LT.TOL) LAMTP=LAMTP-TOL                          
         FAC=LAMPP-(DSIGN(ONE,CL))*DSIN(LAMPP)*PI/TWO                  
  232    LAMT=LON+P22*SAV                                              
         C=DCOS(LAMT)                                                  
         IF (DABS(C).LT.TOL) THEN                                      
            LAMDP = SAV                                                
            GO TO 233                                                  
         END IF                                                        
         XLAM=((ONE-ES)*DTAN(LAT)*SA+DSIN(LAMT)*CA)/C                  
         LAMDP=DATAN(XLAM)                                             
         LAMDP=LAMDP+FAC                                               
         DIF=DABS(SAV)-DABS(LAMDP)                                     
         IF (DABS(DIF).LT.TOL) GO TO 233                               
         SAV=LAMDP                                                     
         L=L+1                                                         
         IF (L.GT.50) GO TO 234                                        
         GO TO 232                                                     
!                                                                      
!        ADJUST FOR LANDSAT ORIGIN.                                    
!                                                                      
  233    RLM=PI*(16.D0/31.D0+ONE/248.D0)                               
         RLM2=RLM+TWO*PI                                               
         NN=NN+1                                                       
         IF (NN.GE.3) GO TO 236                                        
         IF (LAMDP.GT.RLM.AND.LAMDP.LT.RLM2) GO TO 236                 
         IF (LAMDP.LE.RLM) LAMPP=2.5D0*PI                              
         IF (LAMDP.GE.RLM2) LAMPP=PI/TWO                               
         GO TO 231                                                     
  234    IF (IPEMSG .EQ. 0) WRITE (IPELUN,235)                         
  235    FORMAT ('0ERROR PJ22Z0'/                                     & 
     &           ' 50 ITERATIONS WITHOUT CONVERGENCE.')                
         IERROR = 223                                                  
  236    CONTINUE                                                      
!                                                                      
!        LAMDP COMPUTED.  NOW COMPUTE PHIDP.                           
!                                                                      
         SP=DSIN(LAT)                                                  
         PHIDP=ASINZ0(((ONE-ES)*CA*SP-SA*DCOS(LAT)*DSIN(LAMT))/DSQRT(ON & 
     &     -ES*SP*SP))                                                 
!                                                                      
!        COMPUTE X AND Y                                               
!                                                                      
         TANPH=DLOG(DTAN(PI/4.0D0+PHIDP/TWO))                          
         SD=DSIN(LAMDP)                                                
         SDSQ=SD*SD                                                    
         S=P22*SA*DCOS(LAMDP)*DSQRT((ONE+T*SDSQ)/((ONE+W*SDSQ)*(ONE    & 
     &     +Q*SDSQ)))                                                  
         D=DSQRT(XJ*XJ+S*S)                                            
         X=B*LAMDP+A2*DSIN(TWO*LAMDP)+A4*DSIN(4.0D0*LAMDP)-TANPH*S/D   
         X=A*X                                                         
         Y=C1*SD+C3*DSIN(3.0D0*LAMDP)+TANPH*XJ/D                       
         Y=A*Y                                                         
         PROJ(1)=X+X0                                                  
         PROJ(2)=Y+Y0                                                  
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(22) .NE. 0) GO TO 320                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 224                                                  
         RETURN                                                        
  320    X = PROJ(1) -X0                                               
         Y = PROJ(2) -Y0                                               
!                                                                      
!        COMPUTE TRANSFORMED LAT/LONG AND GEODETIC LAT/LONG, GIVEN X,Y.
!                                                                      
!        BEGIN INVERSE COMPUTATION WITH APPROXIMATION FOR LAMDP.  SOLVE
!        FOR TRANSFORMED LONG.                                         
!                                                                      
         LAMDP=X/A/B                                                   
         NN=0                                                          
  325    SAV=LAMDP                                                     
         SD=DSIN(LAMDP)                                                
         SDSQ=SD*SD                                                    
         S=P22*SA*DCOS(LAMDP)*DSQRT((ONE+T*SDSQ)/((ONE+W*SDSQ)*(ONE+Q  & 
     &     *SDSQ)))                                                    
         LAMDP=X/A+Y/A*S/XJ-A2*DSIN(TWO*LAMDP)-A4*DSIN(4.0D0*LAMDP)   & 
     &     -(S/XJ)*(C1*DSIN(LAMDP)+C3*DSIN(3.0D0*LAMDP))               
         LAMDP=LAMDP/B                                                 
         DIF=LAMDP-SAV                                                 
         IF (DABS(DIF).LT.TOL) GO TO 330                               
         NN=NN+1                                                       
         IF (NN.EQ.50) GO TO 330                                       
         GO TO 325                                                     
!                                                                      
!        COMPUTE TRANSFORMED LAT.                                      
!                                                                      
  330    SL=DSIN(LAMDP)                                                
         FAC=DEXP(DSQRT(ONE+S*S/XJ/XJ)*(Y/A-C1*SL-C3*DSIN(3.0D0*LAMDP)))
         ACTAN=DATAN(FAC)                                              
         PHIDP=TWO*(ACTAN-PI/4.0D0)                                    
!                                                                      
!        COMPUTE GEODETIC LATITUDE.                                    
!                                                                      
         DD=SL*SL                                                      
         IF (DABS(DCOS(LAMDP)).LT.TOL) LAMDP=LAMDP-TOL                 
         SPP=DSIN(PHIDP)                                               
         SPPSQ=SPP*SPP                                                 
         LAMT=DATAN(((ONE-SPPSQ/(ONE-ES))*DTAN(LAMDP)*CA-SPP*SA*DSQRT(( &
     &   ONE+Q*DD)*(ONE-SPPSQ)-SPPSQ*U)/DCOS(LAMDP))/(ONE-SPPSQ*(ONE+U)) &
     &   )                                                             
!                                                                      
!        CORRECT INVERSE QUADRANT.                                     
!                                                                      
         IF (LAMT.GE.0) SL=ONE                                         
         IF (LAMT.LT.0) SL=-ONE                                        
         IF (DCOS(LAMDP).GE.0) SCL=ONE                                 
         IF (DCOS(LAMDP).LT.0) SCL=-ONE                                
         LAMT=LAMT-PI/TWO*(ONE-SCL)*SL                                 
         LON=LAMT-P22*LAMDP                                            
!                                                                      
!        COMPUTE GEODETIC LATITUDE.                                    
!                                                                      
         IF (DABS(SA).LT.TOL) LAT=ASINZ0(SPP/DSQRT((ONE-ES)*(ONE-ES)   & 
     &      +ES*SPPSQ))                                                
         IF (DABS(SA).LT.TOL) GO TO 335                                
         LAT=DATAN((DTAN(LAMDP)*DCOS(LAMT)-CA*DSIN(LAMT))/((ONE-ES)*SA))
  335    CONTINUE                                                      
         IF (LAND.GE.4) GO TO 370                                      
         GEOG(1)=LON+128.87D0*DG1-PI*TWO/251.D0*DBLE(PATH)             
         GO TO 380                                                     
  370    GEOG(1)=LON+129.30D0*DG1-PI*TWO/233.D0*DBLE(PATH)             
  380    GEOG(2)=LAT                                                   
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   PJ23Z0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! **              MATHEMATICAL ANALYSIS BY JOHN SNYDER                *
! *********************************************************************
!            * MODIFIED-STEREOGRAPHIC CONFORMAL (FOR ALASKA) *         
! *********************************************************************
!                                                                      
      SUBROUTINE PJ23Z0 (COORD,CRDIO,INDIC)                            
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      INTEGER*4 IERROR,IPEMSG,IPELUN,IPPARM,IPPLUN,N,J,NN              
      INTEGER*4 SWITCH(23)                                             
      INTEGER*2 INDIC                                                  
      DIMENSION GEOG(2),PROJ(2),COORD(2),CRDIO(2),                    & 
     & ACOEF(6),BCOEF(6)                                               
! **** PARAMETERS **** A,E,ES,LON0,LAT0,X0,Y0,SINPH0,COSPH0 ***********
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PJ23/ A,LON0,X0,Y0,ACOEF,BCOEF,EC,LAT0,CCHIO,SCHIO,N     
      COMMON /TOGGLE/ SWITCH                                           
      DATA HALFPI /1.5707963267948966D0/                               
      DATA EPSLN /1.0D-10/                                             
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                            
!                                                                      
! .....................................................................
!                      .  FORWARD TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 0) THEN                                           
!                                                                      
         GEOG(1) = COORD(1)                                            
         GEOG(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(23) .NE. 0) GO TO 120                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
 2010    FORMAT ('0ERROR PJ23Z0'/                                   &   
     &           ' PROJECTION WAS NOT INITIALIZED')                    
         IERROR = 232                                                  
         RETURN                                                        
  120    LON = ADJLZ0 (GEOG(1) - LON0)                                 
!                                                                      
!     CALCULATE X-PRIME AND Y-PRIME FOR OBLIQUE STEREOGRAPHIC PROJ.    
!          FROM LAT/LONG.                                              
!                                                                      
         SINLON = DSIN (LON)                                           
         COSLON = DCOS (LON)                                           
         ESPHI = EC *DSIN(GEOG(2))                                     
         CHI=TWO*DATAN(DTAN((HALFPI+GEOG(2))/TWO)*((ONE-ESPHI)/(ONE   & 
     &      +ESPHI))**(EC/TWO)) - HALFPI                               
         SCHI=DSIN(CHI)                                                
         CCHI=DCOS(CHI)                                                
         G=SCHIO*SCHI+CCHIO*CCHI*COSLON                                
         S=TWO/(ONE+G)                                                 
         XP=S*CCHI*SINLON                                              
         YP=S*(CCHIO*SCHI-SCHIO*CCHI*COSLON)                           
!                                                                      
!     USE KNUTH ALGORITHM FOR SUMMING COMPLEX TERMS, TO CONVERT        
!     OBLIQUE STEREOGRAPHIC TO MODIFIED-STEREOGRAPHIC COORD.           
!                                                                      
         R=XP+XP                                                       
         S=XP*XP+YP*YP                                                 
         AR=ACOEF(N)                                                   
         AI=BCOEF(N)                                                   
         BR=ACOEF(N-1)                                                 
         BI=BCOEF(N-1)                                                 
         DO 140 J=2,N                                                  
         ARN=BR+R*AR                                                   
         AIN=BI+R*AI                                                   
         IF (J.EQ.N) GO TO 140                                         
         BR=ACOEF(N-J)-S*AR                                            
         BI=BCOEF(N-J)-S*AI                                            
         AR=ARN                                                        
         AI=AIN                                                        
  140    CONTINUE                                                      
         BR=-S*AR                                                      
         BI=-S*AI                                                      
         AR=ARN                                                        
         AI=AIN                                                        
         X=XP*AR-YP*AI+BR                                              
         Y=YP*AR+XP*AI+BI                                              
         PROJ(1)=X*A+X0                                                
         PROJ(2)=Y*A+Y0                                                
         CRDIO(1) = PROJ(1)                                            
         CRDIO(2) = PROJ(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
! .....................................................................
!                      .  INVERSE TRANSFORMATION  .                    
! .....................................................................
!                                                                      
      IF (INDIC .EQ. 1) THEN                                           
!                                                                      
         PROJ(1) = COORD(1)                                            
         PROJ(2) = COORD(2)                                            
         IERROR = 0                                                    
         IF (SWITCH(23) .NE. 0) GO TO 220                              
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,2010)                        
         IERROR = 234                                                  
         RETURN                                                        
  220    X = (PROJ(1) - X0)/A                                          
         Y = (PROJ(2) - Y0)/A                                          
         XP=X                                                          
         YP=Y                                                          
         NN=0                                                          
!                                                                      
!     USE KNUTH ALGORITHM FOR SUMMING COMPLEX TERMS, TO CONVERT        
!     MODIFIED-STEREOGRAPHIC CONFORMAL TO OBLIQUE STEREOGRAPHIC        
!     COORDINATES (XP,YP).                                             
!                                                                      
  225    R=XP+XP                                                       
         S=XP*XP+YP*YP                                                 
         AR=ACOEF(N)                                                   
         AI=BCOEF(N)                                                   
         BR=ACOEF(N-1)                                                 
         BI=BCOEF(N-1)                                                 
         CR=DBLE(N)*AR                                                 
         CI=DBLE(N)*AI                                                 
         DR=(DBLE(N-1))*BR                                             
         DI=(DBLE(N-1))*BI                                             
         DO 230 J=2,N                                                  
         ARN=BR+R*AR                                                   
         AIN=BI+R*AI                                                   
         IF (J.EQ.N) GO TO 230                                         
         BR=ACOEF(N-J)-S*AR                                            
         BI=BCOEF(N-J)-S*AI                                            
         AR=ARN                                                        
         AI=AIN                                                        
         CRN=DR+R*CR                                                   
         CIN=DI+R*CI                                                   
         DR=DBLE(N-J)*ACOEF(N-J)-S*CR                                  
         DI=DBLE(N-J)*BCOEF(N-J)-S*CI                                  
         CR=CRN                                                        
         CI=CIN                                                        
  230    CONTINUE                                                      
         BR=-S*AR                                                      
         BI=-S*AI                                                      
         AR=ARN                                                        
         AI=AIN                                                        
         FXYR=XP*AR-YP*AI+BR-X                                         
         FXYI=YP*AR+XP*AI+BI-Y                                         
         FPXYR=XP*CR-YP*CI+DR                                          
         FPXYI=YP*CR+XP*CI+DI                                          
         DEN=FPXYR*FPXYR+FPXYI*FPXYI                                   
         DXP=-(FXYR*FPXYR+FXYI*FPXYI)/DEN                              
         DYP=-(FXYI*FPXYR-FXYR*FPXYI)/DEN                              
         XP=XP+DXP                                                     
         YP=YP+DYP                                                     
         DS=DABS(DXP)+DABS(DYP)                                        
         NN=NN+1                                                       
         IF (NN.LE.20) GO TO 237                                       
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,235)                         
  235    FORMAT ('0ERROR PJ23Z0'/                                     & 
     &           ' TOO MANY ITERATIONS IN ITERATING INVERSE')          
         IERROR = 235                                                  
         GO TO 238                                                     
  237    IF (DS.GT.EPSLN) GO TO 225                                    
!                                                                      
!     CONVERT OBLIQUE STEREOGRAPHIC COORDINATES TO LAT/LONG.           
!                                                                      
  238    RH = DSQRT (XP * XP + YP * YP)                                
         Z = TWO * DATAN (RH / TWO)                                    
         SINZ = DSIN (Z)                                               
         COSZ = DCOS (Z)                                               
         GEOG(1) = LON0                                                
         IF (DABS(RH) .GT. EPSLN) GO TO 240                            
         GEOG(2) = LAT0                                                
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
  240    CHI = ASINZ0 (COSZ * SCHIO + YP *SINZ * CCHIO / RH)           
         NN=0                                                          
         PHI=CHI                                                       
  250    ESPHI=EC*DSIN(PHI)                                            
         DPHI=TWO*DATAN(DTAN((HALFPI+CHI)/TWO)*((ONE+ESPHI)/(ONE-ESPHI)) &
     &      **(EC/TWO)) - HALFPI - PHI                                 
         PHI = PHI + DPHI                                              
         NN = NN + 1                                                   
         IF (NN.LE.20) GO TO 257                                       
         IF (IPEMSG .EQ. 0) WRITE (IPELUN,255)                         
  255    FORMAT ('0ERROR PJ23Z0'/                                     & 
     &           ' TOO MANY ITERATIONS IN CALCULATING PHI FROM CHI')   
         IERROR = 236                                                  
         GO TO 260                                                     
  257    IF (DABS(DPHI).GT.EPSLN) GO TO 250                            
  260    GEOG(2)=PHI                                                   
         GEOG(1) = ADJLZ0 (LON0 + DATAN2(XP*SINZ, RH*CCHIO*COSZ-YP*SCHI &
     &     *SINZ))                                                     
         CRDIO(1) = GEOG(1)                                            
         CRDIO(2) = GEOG(2)                                            
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
!                   QSFNZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION QSFNZ0 (ECCENT,SINPHI,COSPHI)          
!                                                                      
! FUNCTION TO COMPUTE CONSTANT (SMALL Q).                              
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      DATA HALF,ONE,TWO /0.5D0,1.0D0,2.0D0/                            
      DATA EPSLN /1.0D-7/                                              
!                                                                      
      IF (ECCENT .LT. EPSLN) GO TO 020                                 
      CON = ECCENT * SINPHI                                            
      QSFNZ0 = (ONE - ECCENT * ECCENT) * (SINPHI / (ONE - CON * CON) - &
     &         (HALF / ECCENT) * DLOG ((ONE - CON) / (ONE + CON)))     
      RETURN                                                           
!                                                                      
  020 QSFNZ0 = TWO * SINPHI                                            
      RETURN                                                           
      END                                                              
!                   RADDZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      SUBROUTINE RADDZ0 (RAD,SGNA,DEGS,MINS,SECS)                      
!                                                                      
! SUBROUTINE TO CONVERT ANGLE FROM RADIANS TO SIGNED DMS               
! SGNA : SIGN OF ANGLE                                                 
! DEGS : DEGREES PORTION OF ANGLE                                      
! MINS : MINUTES PORTION OF ANGLE                                      
! SECS : SECONDS PORTION OF ANGLE                                      
!                                                                      
      REAL*8 RAD,CON,RADSEC,ZERO,TOL                                   
      REAL*4 SECS                                                      
      INTEGER*4 DEGS,MINS                                              
      CHARACTER*1 SGNA,BLANK,NEG                                       
      DATA RADSEC /206264.806247D0/                                    
      DATA ZERO,TOL /0.0D0,1.0D-4/                                     
      DATA BLANK,NEG /' ','-'/                                         
!                                                                      
! CONVERT THE ANGLE TO SECONDS.                                        
!                                                                      
      CON = DABS(RAD) * RADSEC                                         
      ISEC = IDINT(CON + TOL)                                          
!                                                                      
! DETERMINE THE SIGN OF THE ANGLE.                                     
!                                                                      
      SGNA = BLANK                                                     
      IF (RAD .LT. ZERO .AND. CON .GE. 0.00005D0) SGNA = NEG           
      IF (CON .LT. 0.00005D0) CON = ZERO                               
!                                                                      
! COMPUTE DEGREES PART OF THE ANGLE.                                   
!                                                                      
      INTG = ISEC / 3600                                               
      DEGS = INTG                                                      
      ISEC = INTG * 3600                                               
      CON = CON - DBLE(ISEC)                                           
      ISEC = IDINT(CON + TOL)                                          
!                                                                      
! COMPUTE MINUTES PART OF THE ANGLE.                                   
!                                                                      
      MINS = ISEC / 60                                                 
      ISEC = MINS * 60                                                 
      CON = CON - DBLE(ISEC)                                           
!                                                                      
! COMPUTE SECONDS PART OF THE ANGLE.                                   
!                                                                      
      SECS = SNGL(CON)                                                 
!                                                                      
!     INCREASE MINS IF SECS CLOSE TO 60.000                            
!                                                                      
      IF(SECS .LT. 59.9995D0) RETURN                                   
      MINS = MINS + 1                                                  
      SECS = 0.0                                                       
!                                                                      
!     INCREASE DEGS IF MINS EQUAL 60                                   
!                                                                      
      IF(MINS .LE. 59) RETURN                                          
      MINS = 0                                                         
      DEGS = DEGS + 1                                                  
!                                                                      
      RETURN                                                           
      END                                                              
!                   SERAZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      SUBROUTINE SERAZ0 (FB,FA2,FA4,FC1,FC3,LAM)                       
!                                                                      
!     COMPUTES INTEGRAL FUNCTION OF TRANSFORMED LONG. FOR FOURIER      
!     CONSTANTS A2, A4, B, C1, AND C3.                                 
!     LAM IS INTEGRAL VALUE OF TRANSFORMED LONG.                       
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      COMMON /NORM/ Q,T,U,W,ES,P22,SA,CA,XJ                            
      DATA DG1 /0.01745329252D0/                                       
      DATA ONE,TWO /1.0D0,2.0D0/                                       
      LAM=LAM*DG1                                                      
      SD=DSIN(LAM)                                                     
      SDSQ=SD*SD                                                       
      S=P22*SA*DCOS(LAM)*DSQRT((ONE+T*SDSQ)/((ONE+W*SDSQ)             & 
     &  *(ONE+Q*SDSQ)))                                                
      H=DSQRT((ONE+Q*SDSQ)/(ONE+W*SDSQ))*(((ONE+W*SDSQ)/              & 
     &   ((ONE+Q*SDSQ)**TWO))-P22*CA)                                  
      SQ=DSQRT(XJ*XJ+S*S)                                              
      FB=(H*XJ-S*S)/SQ                                                 
      FA2=FB*DCOS(TWO*LAM)                                             
      FA4=FB*DCOS(4.0D0*LAM)                                           
      FC=S*(H+XJ)/SQ                                                   
      FC1=FC*DCOS(LAM)                                                 
      FC3=FC*DCOS(3.0D0*LAM)                                           
      RETURN                                                           
      END                                                              
!                   SPHDZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      SUBROUTINE SPHDZ0(ISPH,PARM)                                     
!                                                                      
!     SUBROUTINE TO COMPUTE SPHEROID PARAMETERS                        
!                                                                      
!     ISPH IS THE SPHEROID CODE FROM THE FOLLOWING LIST:               
!     0 = CLARKE 1866           1 = CLARKE 1880                        
!     2 = BESSEL                3 = NEW INTERNATIONAL 1967             
!     4 = INTERNATIONAL 1909    5 = WGS 72                             
!     6 = EVEREST               7 = WGS 66                             
!     8 = GRS 1980              9 = AIRY                               
!    10 = MODIFIED EVEREST     11 = MODIFIED AIRY                      
!    12 = WGS 84               13 = SOUTHEAST ASIA                     
!    14 = AUSTRALIAN NATIONAL  15 = KRASSOVSKY                         
!    16 = HOUGH                17 = MERCURY 1960                       
!    18 = MODIFIED MERC 1968   19 = SPHERE OF RADIUS 6370997 M         
!                                                                      
!    PARM IS ARRAY OF PROJECTION PARAMETERS:                           
!       PARM(1) IS THE SEMI-MAJOR AXIS                                 
!       PARM(2) IS THE ECCENTRICITY SQUARED                            
!                                                                      
!     IF ISPH IS NEGATIVE, USER SPECIFIED PROJECTION PARAMETERS ARE TO 
!     DEFINE THE RADIUS OF SPHERE OR ELLIPSOID CONSTANTS AS APPROPRIATE
!                                                                      
!     IF ISPH = 0 , THE DEFAULT IS RESET TO CLARKE 1866                
!                                                                      
! ****                                                             ****
!                                                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      DIMENSION PARM(15),AXIS(20),BXIS(20)                             
!                                                                      
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z,E4Z                    
      COMMON /SPHRZ0/ AZZ                                              
      COMMON /ERRMZ0/ IERROR                                           
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      COMMON /PROJZ0/ IPROJ                                            
!                                                                      
      DATA ZERO,ONE /0.0D0,1.0D0/                                      
!                                                                      
      DATA AXIS/6378206.4D0,6378249.145D0,6377397.155D0,6378157.5D0,      &
     & 6378388.0D0,6378135.0D0,6377276.3452D0,6378145.0D0,6378137.0D0,    &
     & 6377563.396D0,6377304.063D0,6377340.189D0,6378137.0D0,6378155.D0,  &
     & 6378160.0D0,6378245.0D0,6378270.0D0,6378166.0D0,6378150.0D0,       &
     & 6370997.0D0/                                                    
!                                                                      
      DATA BXIS/6356583.8D0,6356514.86955D0,6356078.96284D0,             & 
     & 6356772.2D0,6356911.94613D0,6356750.519915D0,6356075.4133D0,      & 
     & 6356759.769356D0,6356752.314140D0,6356256.91D0,6356103.039D0,     & 
     & 6356034.448D0,6356752.314245D0,6356773.3205D0,6356774.719D0,      & 
     & 6356863.0188D0,6356794.343479D0,6356784.283666D0,6356768.337303D0, &  
     & 6370997.0D0/                                                   
!                                                                      
      IF (ISPH.GE.0) GO TO 5                                           
!                                                                      
!     INITIALIZE USER SPECIFIED SPHERE AND ELLIPSOID PARAMETERS        
!                                                                      
      AZZ = ZERO                                                       
      AZ = ZERO                                                        
      EZ = ZERO                                                        
      ESZ = ZERO                                                       
      E0Z = ZERO                                                       
      E1Z = ZERO                                                       
      E2Z = ZERO                                                       
      E3Z = ZERO                                                       
      E4Z = ZERO                                                       
!                                                                      
!     FETCH FIRST TWO USER SPECIFIED PROJECTION PARAMETERS             
!                                                                      
      A = DABS(PARM(1))                                                
      B = DABS(PARM(2))                                                
      IF (A .GT. ZERO .AND. B .GT. ZERO) GO TO 13                      
      IF (A .GT. ZERO .AND. B .LE. ZERO) GO TO 12                      
      IF (A .LE. ZERO .AND. B .GT. ZERO) GO TO 11                      
!                                                                      
!     DEFAULT NORMAL SPHERE AND CLARKE 1866 ELLIPSOID                  
!                                                                      
      JSPH = 1                                                         
      GO TO 10                                                         
!                                                                      
!     DEFAULT CLARKE 1866 ELLIPSOID                                    
!                                                                      
   11 A = AXIS(1)                                                      
      B = BXIS(1)                                                      
      GO TO 14                                                         
!                                                                      
!     USER SPECIFIED RADIUS OF SPHERE                                  
!                                                                      
   12 AZZ = A                                                          
      GO TO 15                                                         
!                                                                      
!     USER SPECIFIED SEMI-MAJOR AND SEMI-MINOR AXES OF ELLIPSOID       
!                                                                      
   13 IF (B .LE. ONE) GO TO 15                                         
   14 ES = ONE - (B / A)**2                                            
      GO TO 16                                                         
!                                                                      
!     USER SPECIFIED SEMI-MAJOR AXIS AND ECCENTRICITY SQUARED          
!                                                                      
   15 ES = B                                                           
   16 AZ = A                                                           
      ESZ = ES                                                         
      EZ  = DSQRT(ES)                                                  
      E0Z = E0FNZ0(ES)                                                 
      E1Z = E1FNZ0(ES)                                                 
      E2Z = E2FNZ0(ES)                                                 
      E3Z = E3FNZ0(ES)                                                 
      E4Z = E4FNZ0(EZ)                                                 
      PARM(1) = A                                                      
      PARM(2) = ES                                                     
      RETURN                                                           
!                                                                      
!     CHECK FOR VALID SPHEROID SELECTION                               
!                                                                      
    5 IF (PARM(1).NE.ZERO.AND.IPROJ.NE.1) RETURN                       
      JSPH = IABS(ISPH) + 1                                            
      IF (JSPH.LE.20) GO TO 10                                         
      IERROR = 999                                                     
      IF (IPEMSG .EQ. 0) WRITE (IPELUN,1) ISPH                         
    1 FORMAT('0ERROR SPHDZ0:  SPHEROID CODE OF ',I5,' RESET TO 0')     
      ISPH = 0                                                         
      JSPH = 1                                                         
!                                                                      
!     RETRIEVE A AND B AXES FOR SELECTED SPHEROID                      
!                                                                      
   10 A = AXIS(JSPH)                                                   
      B = BXIS(JSPH)                                                   
      ES = ONE - (B / A)**2                                            
!                                                                      
!     SET COMMON BLOCK PARAMETERS FOR SELECTED SPHEROID                
!                                                                      
      AZZ = 6370997.0D0                                                
      EZ  = DSQRT(ES)                                                  
      E0Z = E0FNZ0(ES)                                                 
      E1Z = E1FNZ0(ES)                                                 
      E2Z = E2FNZ0(ES)                                                 
      E3Z = E3FNZ0(ES)                                                 
      E4Z = E4FNZ0(EZ)                                                 
      AZ  = A                                                          
      ESZ = ES                                                         
      IF (ES.EQ.ZERO) AZZ=A                                            
!                                                                      
      PARM(1) = A                                                      
      PARM(2) = ES                                                     
      RETURN                                                           
      END                                                              
!                   TSFNZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      DOUBLE PRECISION FUNCTION TSFNZ0 (ECCENT,PHI,SINPHI)             
!                                                                      
! FUNCTION TO COMPUTE CONSTANT (SMALL T).                              
!                                                                      
      IMPLICIT REAL*8 (A-Z)                                            
      DATA HALF,ONE /0.5D0,1.0D0/                                      
      DATA HALFPI /1.5707963267948966D0/                               
!                                                                      
      CON = ECCENT * SINPHI                                            
      COM = HALF * ECCENT                                              
      CON = ((ONE - CON) / (ONE + CON)) ** COM                         
      TSFNZ0 = DTAN (HALF * (HALFPI - PHI)) / CON                      
!                                                                      
      RETURN                                                           
      END                                                              
!                   UNTFZ0                                             
! *********************************************************************
! ** GENERAL CARTOGRAPHIC TRANSFORMATION PACKAGE (GCTP) VERSION 2.0.2 *
! ** U. S. GEOLOGICAL SURVEY - SNYDER, ELASSAL, AND LINCK    06/08/94 *
! *********************************************************************
!                                                                      
      SUBROUTINE UNTFZ0 (INUNIT,IOUNIT,FACTOR,IFLG)                    
!                                                                      
! SUBROUTINE TO DETERMINE CONVERGENCE FACTOR BETWEEN TWO LINEAL UNITS  
!                                                                      
! * INPUT ........                                                     
! * INUNIT * UNIT CODE OF SOURCE.                                      
! * IOUNIT * UNIT CODE OF TARGET.                                      
!                                                                      
! * OUTPUT .......                                                     
! * FACTOR * CONVERGENCE FACTOR FROM SOURCE TO TARGET.                 
! * IFLG   * RETURN FLAG .EQ. 0 , NORMAL RETURN.                       
!            RETURN FLAG .NE. 0 , ABNORMAL RETURN.                     
!                                                                      
      IMPLICIT REAL*8 (A-H,O-Z)                                        
      DIMENSION FACTRS(6,6)                                            
      COMMON /PRINZ0/ IPEMSG,IPELUN,IPPARM,IPPLUN                      
      PARAMETER (ZERO = 0.0D0, MAXUNT = 6)                             
      DATA FACTRS /0.1000000000000000D01 , 0.0000000000000000D00 , &    
     &             0.0000000000000000D00 , 0.2062648062470963D06 , &    
     &             0.5729577951308231D02 , 0.0000000000000000D00 , &    
     &             0.0000000000000000D00 , 0.1000000000000000D01 , &    
     &             0.3048006096012192D00 , 0.0000000000000000D00 , &    
     &             0.0000000000000000D00 , 0.1000002000004000D01 , &    
     &             0.0000000000000000D00 , 0.3280833333333333D01 , &    
     &             0.1000000000000000D01 , 0.0000000000000000D00 , &    
     &             0.0000000000000000D00 , 0.3280839895013124D01 , &    
     &             0.4848136811095360D-5 , 0.0000000000000000D00 , &    
     &             0.0000000000000000D00 , 0.1000000000000000D01 , &    
     &             0.2777777777777778D-3 , 0.0000000000000000D00 , &    
     &             0.1745329251994330D-1 , 0.0000000000000000D00 , &    
     &             0.0000000000000000D00 , 0.3600000000000000D04 , &    
     &             0.1000000000000000D01 , 0.0000000000000000D00 , &    
     &             0.0000000000000000D00 , 0.9999980000000000D00 , &    
     &             0.3048000000000000D00 , 0.0000000000000000D00 , &    
     &             0.0000000000000000D00 , 0.1000000000000000D01 /     
!                                                                      
      IF (INUNIT .GE. 0 .AND. INUNIT .LT. MAXUNT .AND.                & 
     &    IOUNIT .GE. 0 .AND. IOUNIT .LT. MAXUNT) THEN                 
         FACTOR = FACTRS(IOUNIT+1 , INUNIT+1)                          
         IF (FACTOR .NE. ZERO) THEN                                    
            IFLG = 0                                                   
            RETURN                                                     
         ELSE                                                          
            IF (IPEMSG .NE. 0) WRITE (IPELUN,2000) INUNIT,IOUNIT       
 2000       FORMAT (' INCONSISTENT UNIT CODES = ',I6,' / ',I6)         
            IFLG = 12                                                  
            RETURN                                                     
         END IF                                                        
      ELSE                                                             
         IF (INUNIT.LT.0 .OR. INUNIT.GE.MAXUNT) THEN                   
            IF (IPEMSG .NE. 0) WRITE (IPELUN,2010) INUNIT,IOUNIT       
 2010       FORMAT (' ILLEGAL SOURCE OR TARGET UNIT CODE = ',I6,' / ', &
     &              I6)                                                
         END IF                                                        
         IF (IOUNIT.LT.0 .OR. IOUNIT.GE.MAXUNT) THEN                   
            IF (IPEMSG .NE. 0) WRITE (IPELUN,2010) IOUNIT,IOUNIT       
         END IF                                                        
         IFLG = 11                                                     
         RETURN                                                        
      END IF                                                           
!                                                                      
      END                                                              
