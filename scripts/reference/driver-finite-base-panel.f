C     Original XFOIL inviscid GGCALC/TECALC/PSILIN, no viscous solve.
      PROGRAM REFERENCE
      INCLUDE 'XFOIL.INC'
      INTEGER NQUERY, IQ
      REAL PXQ, PYQ, PSI1, UQ, VQ
      CALL INIT
      READ(*,*) N, ALFA
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
      LVISC=.FALSE.
      MINF1=0.
      CALL SPECAL
      CALL TECALC
      WRITE(*,*) 'RESULT',CL,CM,CDP,CHORD,ANTE,ASTE,DSTE,
     & SIGTE,GAMTE,SHARP
      DO I=1,N
        WRITE(*,*) 'NODE',I,GAM(I),XP(I),YP(I)
      ENDDO
      READ(*,*) NQUERY
      IF(NQUERY.LT.0.OR.NQUERY.GT.100) STOP 2
      DO IQ=1,NQUERY
        READ(*,*) PXQ,PYQ
        CALL PSILIN(0,PXQ,PYQ,0.0,1.0,PSI1,UQ,.FALSE.,.FALSE.)
        CALL PSILIN(0,PXQ,PYQ,-1.0,0.0,PSI1,VQ,.FALSE.,.FALSE.)
        WRITE(*,*) 'FIELD',IQ,PXQ,PYQ,PSI1,UQ,VQ
      ENDDO
      END
C     Only INIT's graphics setup is replaced; all physics is original.
      SUBROUTINE PLINITIALIZE
      END
      SUBROUTINE PLPORT
      END
      SUBROUTINE COLORSPECTRUMHUES(N,SCHEME)
      CHARACTER*(*) SCHEME
      END
