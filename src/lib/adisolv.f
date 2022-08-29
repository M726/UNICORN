      SUBROUTINE ADISOLV(IEQN,ISOR,RELX,TOLR)                   
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LIJ=LI+LJ,
     1     LPD=55*LE-15*LE-4*LIJ-4)
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/DUMMY/ UTMP(LJ,LI),AP(LJ,LI),AEE(LJ,LI),AE(LJ,LI),
     1 AW(LJ,LI),AWW(LJ,LI),ANN(LJ,LI),AN(LJ,LI),AS(LJ,LI),ASS(LJ,LI),
     2 RHS1(LJ,LI),RHS(LJ,LI),FPA(LJ+2,LI+2),AU(LJ,LI),BU(LJ,LI),
     3 AWM(LJ),ASM(LI),APM(LIJ),DUM(LPD)
      ENORM=DFLOAT(LE)
C                      RELAX BY ADI SCHEME
      SOR1=RELX
      SOR2=(1.0-RELX)
      JBIG=2
      IBIG=2
      IF(IEQN.EQ.1) IBIG=3
      IF(IEQN.EQ.2) JBIG=3
      DO 100 ISORP=1,ISOR
      RSD=0.0
C                      SWEEP IN X - DIRECTION
      DO 20 I=IBIG,LI-1
      II=I+1
      DO J=JBIG,LJ-1
      JJ=J+1
      RHS(J,I)=SOR1*(RHS1(J,I)-(ANN(J,I)*FPA(JJ+2,II)
     1 +AN(J,I)*FPA(JJ+1,II)+AS(J,I)*FPA(JJ-1,II)
     2 +ASS(J,I)*FPA(JJ-2,II)))+SOR2*(AP(J,I)*FPA(JJ,II)
     3 +AEE(J,I)*FPA(JJ,II+2)+AE(J,I)*FPA(JJ,II+1)
     4 +AW(J,I)*FPA(JJ,II-1)+AWW(J,I)*FPA(JJ,II-2))
      ENDDO
   20 CONTINUE
      DO 22 J=JBIG,LJ-1
      JJ=J+1
      A1=AP(J,IBIG)
      AU(J,IBIG)=AE(J,IBIG)/A1
      BU(J,IBIG)=AEE(J,IBIG)/A1
      RHS(J,IBIG)=(RHS(J,IBIG)-FPA(JJ,IBIG)*AW(J,IBIG)
     1         -FPA(JJ,IBIG-1)*AWW(J,IBIG))/A1
      RHS(J,IBIG+1)=(RHS(J,IBIG+1)-FPA(JJ,IBIG)*AWW(J,IBIG+1))
      RHS(J,LI-2)=RHS(J,LI-2)-AEE(J,LI-2)*FPA(JJ,LI+1)
      AWM(J)=AW(J,IBIG+1)
      APM(J)=AP(J,IBIG+1)
   22 CONTINUE
      DO 24 I=IBIG+1,LI-2
      II=I+1
      DO J=JBIG,LJ-1
      JJ=J+1
      A1=APM(J)-AU(J,I-1)*AWM(J)
      AU(J,I)=(AE(J,I)-BU(J,I-1)*AWM(J))/A1
      BU(J,I)=AEE(J,I)/A1
      RHS(J,I)=(RHS(J,I)-RHS(J,I-1)*AWM(J))/A1
      AWM(J)=AW(J,I+1)-AU(J,I-1)*AWW(J,I+1)
      APM(J)=AP(J,I+1)-BU(J,I-1)*AWW(J,I+1)
      RHS(J,I+1)=RHS(J,I+1)-RHS(J,I-1)*AWW(J,I+1)
      ENDDO
   24 CONTINUE
      DO 26 J=JBIG,LJ-1
      JJ=J+1
      A1=APM(J)-AU(J,LI-2)*AWM(J)
      RHS(J,LI-1)=(RHS(J,LI-1)-AEE(J,LI-1)*FPA(JJ,LI+2)
     1         -AE(J,LI-1)*FPA(JJ,LI+1)-RHS(J,LI-2)*AWM(J))/A1
      RHS(J,LI-2)=RHS(J,LI-2)-RHS(J,LI-1)*AU(J,LI-2)
      FPA(JJ,LI)=RHS(J,LI-1)
      FPA(JJ,LI-1)=RHS(J,LI-2)
   26 CONTINUE
      DO 28 I=LI-3,IBIG,-1
      DO J=JBIG,LJ-1
      RHS(J,I)=RHS(J,I)-RHS(J,I+1)*AU(J,I)-RHS(J,I+2)*BU(J,I)
      FPA(J+1,I+1)=RHS(J,I)
      ENDDO
   28 CONTINUE
C                      SWEEP IN Y - DIRECTION
      DO 30 J=JBIG,LJ-1
      JJ=J+1
      DO I=IBIG,LI-1
      II=I+1
      RHS(J,I)=SOR1*(RHS1(J,I)-(AEE(J,I)*FPA(JJ,II+2)
     1 +AE(J,I)*FPA(JJ,II+1)+AW(J,I)*FPA(JJ,II-1)
     2 +AWW(J,I)*FPA(JJ,II-2)))+SOR2*(AP(J,I)*FPA(JJ,II)
     3 +ANN(J,I)*FPA(JJ+2,II)+AN(J,I)*FPA(JJ+1,II)
     4 +AS(J,I)*FPA(JJ-1,II)+ASS(J,I)*FPA(JJ-2,II))
      ENDDO
   30 CONTINUE
      DO 32 I=IBIG,LI-1
      II=I+1
      A1=AP(JBIG,I)
      AU(JBIG,I)=AN(JBIG,I)/A1
      BU(JBIG,I)=ANN(JBIG,I)/A1
      RHS(JBIG,I)=(RHS(JBIG,I)-FPA(JBIG,II)*AS(JBIG,I)
     1         -FPA(JBIG-1,II)*ASS(JBIG,I))/A1
      RHS(JBIG+1,I)=(RHS(JBIG+1,I)-FPA(JBIG,II)*ASS(JBIG+1,I))
      RHS(LJ-2,I)=RHS(LJ-2,I)-ANN(LJ-2,I)*FPA(LJ+1,II)
      ASM(I)=AS(JBIG+1,I)
      APM(I)=AP(JBIG+1,I)
   32 CONTINUE
      DO 34 J=JBIG+1,LJ-2
      JJ=J+1
      DO I=IBIG,LI-1
      II=I+1
      A1=APM(I)-AU(J-1,I)*ASM(I)
      AU(J,I)=(AN(J,I)-BU(J-1,I)*ASM(I))/A1
      BU(J,I)=ANN(J,I)/A1
      RHS(J,I)=(RHS(J,I)-RHS(J-1,I)*ASM(I))/A1
      ASM(I)=AS(J+1,I)-AU(J-1,I)*ASS(J+1,I)
      APM(I)=AP(J+1,I)-BU(J-1,I)*ASS(J+1,I)
      RHS(J+1,I)=RHS(J+1,I)-RHS(J-1,I)*ASS(J+1,I)
      ENDDO
   34 CONTINUE
      DO 36 I=IBIG,LI-1
      II=I+1
      A1=APM(I)-AU(LJ-2,I)*ASM(I)
      RHS(LJ-1,I)=(RHS(LJ-1,I)-ANN(LJ-1,I)*FPA(LJ+2,II)
     1         -AN(LJ-1,I)*FPA(LJ+1,II)-RHS(LJ-2,I)*ASM(I))/A1
      RHS(LJ-2,I)=RHS(LJ-2,I)-RHS(LJ-1,I)*AU(LJ-2,I)
   36 CONTINUE
      DO 38 J=LJ-3,JBIG,-1
      DO I=IBIG,LI-1
      RHS(J,I)=RHS(J,I)-RHS(J+1,I)*AU(J,I)-RHS(J+2,I)*BU(J,I)
      ENDDO
   38 CONTINUE
      DO 40 J=JBIG,LJ-1 
      DO I=IBIG,LI-1
      RSD=RSD+DABS(FPA(J+1,I+1)-RHS(J,I))
      FPA(J+1,I+1)=RHS(J,I)
      ENDDO
   40 CONTINUE
      IF((RSD/ENORM).LE.TOLR) GO TO 101
  100 CONTINUE
      ISOR=ISORP-1
  101 CONTINUE
      ISOR=ISORP
      DO 200 J=1,LJ
      DO I=1,LI
      RHS(J,I)=FPA(J+1,I+1)
      ENDDO
  200 CONTINUE
      RETURN
      END