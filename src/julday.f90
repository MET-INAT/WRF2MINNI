      function julday(mm,id,iyyy)
      integer julday,id,iyyy,mm,IGREG
      parameter (IGREG=15+31*(10+12*1582))
      integer ja,jm,jy
      jy=iyyy
! massimo: per ovviare pause che e' obsoleto:
      !if (jy.eq.0) pause 'julday: there is no year zero'
      if (jy.eq.0) write(*,*) 'julday: there is no year zero'
      if (jy.eq.0) read(*,*)
! massimo: fine

      if (jy.lt.0) jy=jy+1
      if (mm.gt.2) then
          jm=mm+1
      else
          jy=jy-1
	  jm=mm+13
      endif
      julday=365*jy+int(0.25d0*jy+2000.d0)+int(30.6001d0*jm)+id+1718995
      if (id+31*(mm+12*iyyy).ge.IGREG) then
          ja=int(0.01d0*jy)
	  julday=julday+2-ja+int(0.25d0*ja)
      endif
      return
      end
