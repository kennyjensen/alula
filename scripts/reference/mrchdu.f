C     Original XFOIL MRCHDU on supplied surface BL states, without a flow solve.
C     Both independent surface arrays receive the same prescribed profile.
C     No wakes or early trips: XIFSET uses the exact terminal station.
      PROGRAM MRCHDUREF
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      READ(*,*) NCASES
      DO IC=1,NCASES
        CALL BLPINI
        READ(*,*) REINF,MINF,AMCRIT,NS,ITOLD
        IF(NS+1.GT.IVX) STOP 1
        GAMBL=1.4
        GM1BL=0.4
        QINFBL=1.0
        HVRAT=0.35
        IDAMPV=0
        BULE=1.0
        TR=1.0+0.5*GM1BL*MINF**2
        HSTINV=GM1BL*MINF**2/TR
        HSTINV_MS=GM1BL/TR**2
        RSTBL=TR**(1.0/GM1BL)
        RSTBL_MS=0.5*RSTBL/TR
        TKBL=0.0
        TKBL_MS=0.0
        HERAT=1.0-0.5*HSTINV
        HERAT_MS=-0.5*HSTINV_MS
        RFAC=SQRT(HERAT**3)*(1.0+HVRAT)/(HERAT+HVRAT)
        REYBL=REINF*RFAC
        REYBL_RE=RFAC
        REYBL_MS=REYBL*(1.5/HERAT-1.0/(HERAT+HVRAT))*HERAT_MS
        ANTE=0.0
        DO IS=1,2
          ACRIT(IS)=AMCRIT
          XSTRIP(IS)=1.0
          IBLTE(IS)=NS+1
          NBL(IS)=NS+1
          ITRAN(IS)=ITOLD+2
          DO IBL=1,NS+1
            CTAU(IBL,IS)=0.0
            THET(IBL,IS)=0.0
            DSTR(IBL,IS)=0.0
          ENDDO
        ENDDO
        DO IBL=2,NS+1
          READ(*,*) XSSI(IBL,1),UEDG(IBL,1),CTAU(IBL,1),
     &              THET(IBL,1),DSTR(IBL,1)
          XSSI(IBL,2)=XSSI(IBL,1)
          UEDG(IBL,2)=UEDG(IBL,1)
          CTAU(IBL,2)=CTAU(IBL,1)
          THET(IBL,2)=THET(IBL,1)
          DSTR(IBL,2)=DSTR(IBL,1)
          DO IS=1,2
            MASS(IBL,IS)=DSTR(IBL,IS)*UEDG(IBL,IS)
          ENDDO
        ENDDO
        WRITE(*,*) 'CASE',IC
        CALL MRCHDU
        DO IS=1,2
          WRITE(*,*) 'TRANSITION',IC,IS,ITRAN(IS)-2,
     &               XSSITR(IS),TFORCE(IS)
          DO IBL=2,NS+1
            WRITE(*,*) 'STATION',IC,IS,IBL-2,XSSI(IBL,IS),
     &        UEDG(IBL,IS),CTAU(IBL,IS),THET(IBL,IS),DSTR(IBL,IS)
          ENDDO
        ENDDO
      ENDDO
      END
