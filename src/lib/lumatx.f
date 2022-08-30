      SUBROUTINE LUMATX(KBOUND)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,
     1 LB=LJ-2,LL=(LI-2)*LB,LL1=LB*LL-(LB+1)*LB/2)
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/LUMAT/ PSL(LL1),PSD(LL),IBAND(LB)
C-------------  KBOUND=0 --- FOR DIRICHLET CONDITIONS--------
C-------------  KBOUND=1 --- FOR NEUMANN CONDITIONS--------
      LEND=LL-KBOUND
      SYM=DFLOAT(ISYM)
      DO 601 I=1,LB
      IA=(I-1)*I
      IBAND(I)=(I-1)*LL-IA/2
  601 CONTINUE
      IF(KBOUND.EQ.1) THEN
C------------------------  NEUMANN CONDITION  --------------------------
      DO 608 I=2,LI-1
      DO 611 J=2,LJ-1
      IF(I.EQ.LI-1.AND.J.EQ.LJ-1) GO TO 608
      AE=0.0
      AW=0.0
      AN=0.0
      AS=0.0
      II=(I-2)*(LJ-2)+J-1
      IF(I.GT.2) THEN
                 AW=Y(J,I)*YYC(J)/XXS(I)
                 IF(ISYM.EQ.0) AW=YYC(J)/XXS(I)
                 IA=IBAND(LB)+II-LB
                 PSL(IA)=AW
                 END IF
      IF(I.LT.LI-1) THEN 
                 AE=Y(J,I)*YYC(J)/XXS(I+1)
                 IF(ISYM.EQ.0) AE=YYC(J)/XXS(I+1)
                 END IF
      IF(J.GT.2) THEN
                 AS=0.5*(Y(J-1,I)+Y(J,I))*XXC(I)/YYS(J)
                 IF(ISYM.EQ.0) AS=XXC(I)/YYS(J)
                 IA=II-1
                 PSL(IA)=AS
                 END IF
      IF(J.LT.LJ-1) THEN
                 AN=0.5*(Y(J,I)+Y(J+1,I))*XXC(I)/YYS(J+1)
                 IF(ISYM.EQ.0) AN=XXC(I)/YYS(J+1)
                 END IF
      PSD(II)=-(AE+AW+AN+AS)
  611 CONTINUE
  608 CONTINUE
      END IF
      IF(KBOUND.EQ.0) THEN
C------------------------  DIRICHLET CONDITION  ------------------------
      DO 609 I=2,LI-1
      DO 613 J=2,LJ-1
      AE=1.0/XXS(I+1)/XXC(I)
      AW=1.0/XXS(I)/XXC(I)
      AN=(1.0/YYC(J)+0.5*SYM/Y(J,I))/YYS(J+1)
      AS=(1.0/YYC(J)-0.5*SYM/Y(J,I))/YYS(J)
      II=(I-2)*(LJ-2)+J-1
      IF(I.GT.2) THEN
                 IA=IBAND(LB)+II-LB
                 PSL(IA)=AW
                 END IF
      IF(J.GT.2) THEN
                 IA=II-1
                 PSL(IA)=AS
                 END IF
      PSD(II)=-(AE+AW+AN+AS)
  613 CONTINUE
  609 CONTINUE
      END IF
C------------------------     ELEMINATION     --------------------------
      DO 610 I=1,LEND
      DIAG=PSD(I)
      IJEND=I+LB
      IF(IJEND.GT.LEND) IJEND=LEND
      DO 612 J=I+1,IJEND
      ID=IBAND(J-I)+I
      FACT=PSL(ID)/DIAG
      DO 614 K=I+1,J-1
      IC=IBAND(K-I)+I
      IA=IBAND(J-K)+K
      PSL(IA)=PSL(IA)-PSL(ID)*PSL(IC)
  614 CONTINUE
      PSD(J)=PSD(J)-FACT*PSL(ID)
      PSL(ID)=FACT
  612 CONTINUE
  610 CONTINUE 
      RETURN                                                            
      END