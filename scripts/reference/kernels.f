C     Local residual and analytic Jacobian oracle from original XFOIL.
      PROGRAM KERNELS
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
      READ(*,*) NCASES
      DO IC=1,NCASES
        CALL BLPINI
        READ(*,*) ITYP, REYBL
        QINFBL=1.0
        TKBL=0.0
        TKBL_MS=0.25
        RSTBL=1.0
        RSTBL_MS=0.5
        HSTINV=0.0
        HSTINV_MS=0.4
        REYBL_MS=0.0
        REYBL_RE=1.0
        GAMBL=1.4
        GM1BL=0.4
        HVRAT=0.35
        AMCRIT=9.0
        BULE=1.0
        IDAMPV=0
        DO J=1,2
          READ(*,*) XSI,AMI,CTI,THI,DSI,DSWI,UEI
          CALL BLPRV(XSI,AMI,CTI,THI,DSI,DSWI,UEI)
          CALL BLKIN
          CALL BLVAR(MAX(1,ITYP))
          IF(J.EQ.1) THEN
            DO K=1,NCOM
              COM1(K)=COM2(K)
            ENDDO
          ENDIF
        ENDDO
        CALL BLMID(MAX(1,ITYP))
        CALL BLDIF(ITYP)
        DO K=1,4
          WRITE(*,*) 'R',IC,K,VSREZ(K),VSR(K),VSM(K),VSX(K)
          DO L=1,5
            WRITE(*,*) 'J',IC,K,L,VS1(K,L),VS2(K,L)
          ENDDO
        ENDDO
      ENDDO
      READ(*,*) NCASES
      DO IC=1,NCASES
        READ(*,*) HK,RT,MSQ
        CALL CFT(HK,RT,MSQ,CF,CF_HK,CF_RT,CF_MSQ)
        WRITE(*,*) 'CFT',IC,CF,CF_HK,CF_RT,CF_MSQ
      ENDDO
      END
