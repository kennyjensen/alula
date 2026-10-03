C     Finite-Mach comparison; all numerical routines are original XFOIL.
C     Separate driver preserves the Mach-zero fixture provenance.
      PROGRAM REFERENCE
      INCLUDE 'XFOIL.INC'
      CALL INIT
      READ(*,*) N, ALFA, REINF1, ACRIT(1), ACRIT(2),
     &          XSTRIP(1), XSTRIP(2), ITMAX, MINF1
      DO I=1,N
        READ(*,*) X(I),Y(I)
      ENDDO
      CALL SCALC(X,Y,S,N)
      CALL SEGSPL(X,XP,S,N)
      CALL SEGSPL(Y,YP,S,N)
      CALL NCALC(X,Y,S,N,NX,NY)
      CALL LEFIND(SLE,X,XP,Y,YP,S,N)
      XLE=SEVAL(SLE,X,XP,S,N)
      YLE=SEVAL(SLE,Y,YP,S,N)
      XTE=0.5*(X(1)+X(N))
      YTE=0.5*(Y(1)+Y(N))
      CHORD=SQRT((XTE-XLE)**2+(YTE-YLE)**2)
      CALL TECALC
      CALL APCALC
      ALFA=ALFA*DTOR
      LALFA=.TRUE.
      LVISC=.TRUE.
      VACCEL=0.01
      CALL SPECAL
      CALL VISCAL(ITMAX)
      WRITE(*,*) 'RESULT',LVCONV,CL,CM,CD,CDF,CDP,RMSBL,
     &           XOCTR(1),XOCTR(2)
      DO I=1,N+NW
        WRITE(*,*) 'CP',I,X(I),Y(I),CPI(I),CPV(I),QVIS(I)
      ENDDO
      DO IS=1,2
        DO IBL=2,NBL(IS)
          WRITE(*,*) 'BL',IS,IBL,IPAN(IBL,IS),
     &    XSSI(IBL,IS),UEDG(IBL,IS),THET(IBL,IS),DSTR(IBL,IS),
     &    CTAU(IBL,IS),MASS(IBL,IS),TAU(IBL,IS),DIS(IBL,IS)
        ENDDO
      ENDDO
      END
C     Only INIT's graphics setup is replaced. No physics stubs.
      SUBROUTINE PLINITIALIZE
      END
      SUBROUTINE PLPORT
      END
      SUBROUTINE COLORSPECTRUMHUES(N,SCHEME)
      CHARACTER*(*) SCHEME
      END
