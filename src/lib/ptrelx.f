      SUBROUTINE PTRELX(IEQN,ISOR,RELX,TOLR)                   
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LIJ=LI+LJ,
     1   LPD=55*LE-13*LE-2*LIJ-4) 
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/DUMMY/ UTMP(LJ,LI),AP(LJ,LI),AEE(LJ,LI),AE(LJ,LI),
     1 AW(LJ,LI),AWW(LJ,LI),ANN(LJ,LI),AN(LJ,LI),AS(LJ,LI),ASS(LJ,LI),
     2 RHS1(LJ,LI),FP(LJ,LI),FPA(LJ+2,LI+2),DUM(LPD)
      ENORM=DFLOAT(LE)
C                      RELAX AT POINT BY POINT
      SOR1=DABS(RELX)
      SOR2=(1.0-RELX)
      JBIG=2
      IBIG=2
      IF(IEQN.EQ.1) IBIG=3
      IF(IEQN.EQ.2) JBIG=3
      DO 100 ISORP=1,ISOR
      RSD=0.0
      DO 20 J=JBIG,LJ-1
      DO I=IBIG,LI-1
      JJ=J+1
      II=I+1
      FPAO=FPA(JJ,II)
      FPA(JJ,II)=SOR2*FPAO+(RHS1(J,I)
     1  -ANN(J,I)*FPA(JJ+2,II)-ASS(J,I)*FPA(JJ-2,II)
     2  -AN(J,I)*FPA(JJ+1,II)-AS(J,I)*FPA(JJ-1,II)
     3  -AEE(J,I)*FPA(JJ,II+2)-AWW(J,I)*FPA(JJ,II-2)
     4  -AE(J,I)*FPA(JJ,II+1)-AW(J,I)*FPA(JJ,II-1))*SOR1/AP(J,I)
      RSD=RSD+DABS(FPA(JJ,II)-FPAO)
      ENDDO
   20 CONTINUE
      IF((RSD/ENORM).LE.TOLR) GO TO 101
  100 CONTINUE
      ISOR=ISORP-1
  101 CONTINUE
      ISOR=ISORP
      DO 200 J=JBIG,LJ-1 
      DO I=IBIG,LI-1 
      FP(J,I)=FPA(J+1,I+1)
      ENDDO
  200 CONTINUE
      RETURN
      END