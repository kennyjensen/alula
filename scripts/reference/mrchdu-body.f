C     Unmodified MRCHDU on differing upper/lower states and their wake.
C     Sharp TE, natural transition or transition at TE; no early trips.
C     IBLPAN: upper NBL ends at TE; lower NBL includes the merged wake.
C     XICALC: first wake arclength equals the lower TE arclength.
      PROGRAM MRCHDUBODYREF
      INCLUDE 'XFOIL.INC'
      INCLUDE 'XBL.INC'
      INTEGER NS(2),ITOLD(2)
      READ(*,*) NCASES
      DO IC=1,NCASES
        CALL BLPINI
        READ(*,*) REINF,MINF,AMCRIT,NS(1),NS(2),NW,
     &            ITOLD(1),ITOLD(2)
        IF(MAX(NS(1),NS(2))+NW+1.GT.IVX) STOP 1
        IF(NW.GT.IWX) STOP 2
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
        DO IW=1,NW
          WGAP(IW)=0.0
        ENDDO
        DO IS=1,2
          ACRIT(IS)=AMCRIT
          XSTRIP(IS)=1.0
          IBLTE(IS)=NS(IS)+1
          NBL(IS)=IBLTE(IS)
          IF(IS.EQ.2) NBL(IS)=IBLTE(IS)+NW
          ITRAN(IS)=ITOLD(IS)+2
          XSSI(1,IS)=0.0
          CTAU(1,IS)=0.0
          THET(1,IS)=0.0
          DSTR(1,IS)=0.0
          DO IBL=2,NBL(IS)
            READ(*,*) XSSI(IBL,IS),UEDG(IBL,IS),CTAU(IBL,IS),
     &                THET(IBL,IS),DSTR(IBL,IS)
            MASS(IBL,IS)=DSTR(IBL,IS)*UEDG(IBL,IS)
          ENDDO
        ENDDO
        WRITE(*,*) 'CASE',IC
        CALL MRCHDU
        DO IS=1,2
          WRITE(*,*) 'TRANSITION',IC,IS,ITRAN(IS)-2,
     &               XSSITR(IS),TFORCE(IS)
          DO IBL=2,NBL(IS)
            WRITE(*,*) 'STATION',IC,IS,IBL-2,XSSI(IBL,IS),
     &        UEDG(IBL,IS),CTAU(IBL,IS),THET(IBL,IS),DSTR(IBL,IS)
          ENDDO
        ENDDO
      ENDDO
      END
