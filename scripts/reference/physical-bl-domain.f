C     Original BLPRV/BLKIN/HKIN, before BLVAR's closure safeguards.
      PROGRAM BLDOMAIN
      IMPLICIT REAL(M)
      INCLUDE 'XBL.INC'
      READ(*,*) NCASES
      DO IC=1,NCASES
        READ(*,*) MINF,GAMBL,THI,DSI,DWI,UEI
        GM1BL=GAMBL-1.0
        QINFBL=1.0
        TKBL=0.0
        TKBL_MS=0.0
        HVRAT=0.35
        TR=1.0+0.5*GM1BL*MINF**2
        HSTINV=GM1BL*MINF**2/TR
        HSTINV_MS=GM1BL/TR**2
        RSTBL=TR**(1.0/GM1BL)
        RSTBL_MS=0.5*RSTBL/TR
        REYBL=1.0E6
        REYBL_RE=1.0
        REYBL_MS=0.0
        CALL BLPRV(0.5,0.0,0.03,THI,DSI,DWI,UEI)
        CALL BLKIN
        WRITE(*,*) IC,M2,M2_U2,HK2,HK2_T2,HK2_D2,HK2_U2
      ENDDO
      END
