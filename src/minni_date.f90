      INTEGER FUNCTION  minni_date (YEAR,MONTH,DAY,HR)
!
!---COMPUTES THE JULIAN DATE (JD) GIVEN A GREGORIAN CALENDAR
!   DATE (YEAR,MONTH,DAY).
!
      INTEGER YEAR,MONTH,DAY,I,J,K,HR
!
      I= YEAR
      J= MONTH
      K= DAY
!
      minni_date = (K-32075+1461*(I+4800+(J-14)/12)/4+367               &
     &    *(J-2-(J-14)/12*12)                                           &
     &    /12-3*((I+4900+(J-14)/12)/100)/4-1721424)*24+HR
!
      RETURN
      end function minni_date


