C     Original XFOIL TRCHEK and BLVAR values for auxiliary-state transfer.
      PROGRAM TRTRANSFER
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
      READ(*,*) NCASES
      DO IC=1,NCASES
        CALL BLPINI
        READ(*,*) REINF,MINF,AMCRIT,XIFORC
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
        SIMI=.FALSE.
        WAKE=.FALSE.
        TURB=.FALSE.
        TRAN=.FALSE.
        DO J=1,2
          READ(*,*) XSI,AMI,THI,DSI,UEI
          IF(J.EQ.1) AMUP=AMI
          CALL BLPRV(XSI,AMUP,0.03,THI,DSI,0.0,UEI)
          CALL BLKIN
          IF(J.EQ.1) THEN
            CALL BLVAR(1)
            DO K=1,NCOM
              COM1(K)=COM2(K)
            ENDDO
          ENDIF
        ENDDO
        CALL TRCHEK
        WRITE(*,*) 'AMPL',IC,AMPL2,TRAN,TRFORC
        CALL BLPRV(XSI,0.0,0.03,THI,DSI,0.0,UEI)
        CALL BLKIN
        TURB=.TRUE.
        TRAN=.FALSE.
        CALL BLVAR(2)
        CTINIT=CTRCON*EXP(-CTRCEX/(HK2-1.0))*CQ2
        WRITE(*,*) 'SHEAR',IC,CTINIT
      ENDDO
      END
