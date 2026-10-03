C     Native DILW values, native partials, and independent differences of
C     the unchanged native value (the printed partial is not an oracle
C     for the mathematical derivative).
      PROGRAM WAKEDISSIPATION
      REAL HK,RT,DI,DHK,DRT,H,DP2,DP1,DM1,DM2,DUMMY1,DUMMY2
      REAL FD(3)
      READ(*,*) NCASES
      DO IC=1,NCASES
        READ(*,*) HK,RT
        CALL DILW(HK,RT,DI,DHK,DRT)
        DO K=1,3
          H=HK*10.0**(-K-2)
          CALL DILW(HK+2*H,RT,DP2,DUMMY1,DUMMY2)
          CALL DILW(HK+H,RT,DP1,DUMMY1,DUMMY2)
          CALL DILW(HK-H,RT,DM1,DUMMY1,DUMMY2)
          CALL DILW(HK-2*H,RT,DM2,DUMMY1,DUMMY2)
          FD(K)=(DM2-8*DM1+8*DP1-DP2)/(12*H)
        ENDDO
        WRITE(*,*) IC,DI,DHK,DRT,FD
      ENDDO
      END
