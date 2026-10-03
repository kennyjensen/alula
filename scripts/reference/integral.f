C     Native BLSYS/TRDIF/TESYS oracle for arbitrary outer-flow edge states.
      PROGRAM INTEGRAL
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
      READ(*,*) NCASES
      DO IC=1,NCASES
        CALL BLPINI
        READ(*,*) ITYP,REINF,MINF,IMODE,AMCRIT,XIFORC,BULE
        GAMBL=1.4
        GM1BL=0.4
        QINFBL=1.0
        HVRAT=0.35
        IDAMPV=0
        TR=1.0+0.5*GM1BL*MINF**2
        HSTINV=GM1BL*MINF**2/TR
        HSTINV_MS=GM1BL/TR**2
        RSTBL=TR**(1.0/GM1BL)
        RSTBL_MS=0.5*RSTBL/TR
        TKBL=0.0
        TKBL_MS=0.0
        IF(IMODE.EQ.1) THEN
          BETA=SQRT(1.0-MINF**2)
          TKBL=(1.0-BETA)/(1.0+BETA)
          TKBL_MS=1.0/(BETA*(1.0+BETA)**2)
        ENDIF
        HERAT=1.0-0.5*HSTINV
        HERAT_MS=-0.5*HSTINV_MS
        RFAC=SQRT(HERAT**3)*(1.0+HVRAT)/(HERAT+HVRAT)
        REYBL=REINF*RFAC
        REYBL_RE=RFAC
        REYBL_MS=REYBL*(1.5/HERAT-1.0/(HERAT+HVRAT))*HERAT_MS
        SIMI=ITYP.EQ.0
        WAKE=ITYP.EQ.3
        TURB=ITYP.EQ.2.OR.ITYP.EQ.3
        TRAN=ITYP.EQ.4
        TRFORC=.FALSE.
        TRFREE=.FALSE.
        XT=0.0
        KTYPE=1
        IF(ITYP.EQ.2) KTYPE=2
        IF(ITYP.EQ.3) KTYPE=3
        DO J=1,2
          READ(*,*) XSI,AMI,CTI,THI,DSI,DSWI,UEI
          CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWI,UEI)
          CALL BLKIN
          IF(J.EQ.1) THEN
            CALL BLVAR(KTYPE)
            DO K=1,NCOM
              COM1(K)=COM2(K)
            ENDDO
          ENDIF
        ENDDO
        READ(*,*) CTE,TTE,DTE
        IF(ITYP.EQ.4) CALL TRCHEK
        IF(ITYP.EQ.5) THEN
          CALL TESYS(CTE,TTE,DTE)
        ELSE
          CALL BLSYS
        ENDIF
        WRITE(*,*) 'META',IC,XT,TRFORC,TRFREE
        DO K=1,3
          WRITE(*,*) 'R',IC,K,-VSREZ(K),VSR(K),VSM(K),VSX(K)
          DO L=1,5
            WRITE(*,*) 'J',IC,K,L,VS1(K,L),VS2(K,L)
          ENDDO
        ENDDO
      ENDDO
      END
