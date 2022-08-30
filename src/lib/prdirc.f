      SUBROUTINE PRDIRC(KBOUND,RESDM)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LB=LJ-2,LL=(LI-2)*LB,
     1    LL1=LB*LL-(LB+1)*LB/2,LPD=55*LE-4*LE-LL)                 
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),
     1     W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB05/ RHONP(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/LUMAT/ PSL(LL1),PSD(LL),IBAND(LB)
      COMMON/CB09/ IBOT(LI),ITOP(LI),JLFT(LJ),JRGT(LJ),
     1       FBOT(8+LSP,LI),FTOP(8+LSP,LI),FLFT(8+LSP,LJ),FRGT(8+LSP,LJ)
      COMMON/CB10/ FBXP(4,LJ),FBXM(4,LJ),FBYP(4,LI),FBYM(4,LI)
      COMMON/DUMMY/ UOLD(LJ,LI),AA(LJ,LI),BB(LJ,LI),RHS(LL),
     1              YSWCH(LJ,LI),DUM(LPD)
C-------------  KBOUND=0 --- FOR DIRICHLET CONDITIONS--------
C-------------  KBOUND=1 --- FOR NEUMANN CONDITIONS--------
      LEND=LL-KBOUND
      SYM=DFLOAT(ISYM)
      DO 8 I=1,LL
      RHS(I)=0.0
    8 CONTINUE
      DO 101 I=1,LI
      DO J=1,LJ
      YSWCH(J,I)=Y(J,I)
      IF(ISYM.EQ.0) YSWCH(J,I)=1.0
      ENDDO
  101 CONTINUE
      DO 102 I=2,LI
      DO J=2,LJ
      AA(J,I)=0.5*(RHONP(J,I)+RHONP(J,I-1))
      BB(J,I)=0.5*(RHONP(J,I)+RHONP(J-1,I))
      ENDDO
  102 CONTINUE
      DO 20 I=2,LI-1
      DO J=2,LJ-1
      AMDOT=YSWCH(J,I)*XXC(I)*YYC(J)*(RHONP(J,I)-RHO(J,I))/DT
     1     +YYC(J)*YSWCH(J,I)*(AA(J,I+1)*U(J,I+1)-AA(J,I)*U(J,I))
     2     +XXC(I)*0.5*((YSWCH(J+1,I)+YSWCH(J,I))*BB(J+1,I)*V(J+1,I)
     3     -(YSWCH(J,I)+YSWCH(J-1,I))*BB(J,I)*V(J,I))
      II=(I-2)*(LJ-2)+J-1
      IF(ISKIP(J,I).GE.3.AND.ISKIP(J,I).NE.11) AMDOT=0.0
      RHS(II)=AMDOT/DT
      ENDDO
   20 CONTINUE
      IF(KBOUND.EQ.1) RHS(LL)=0.0
C------------------------------ SOLVE BY LU ----------------------------
      DO 620 I=1,LEND
      JEND=I+LB
      IF(JEND.GT.LEND) JEND=LEND
      FACT=RHS(I)
      DO 622 J=I+1,JEND
      IA=IBAND(J-I)+I
      RHS(J)=RHS(J)-PSL(IA)*FACT
  622 CONTINUE
      RHS(I)=RHS(I)/PSD(I)
  620 CONTINUE
C
      DO 630 I=LEND,1,-1
      JBIG=I-LB
      IF(JBIG.LT.1) JBIG=1
      FACT=RHS(I)
      DO 632 J=JBIG,I-1
      IB=IBAND(I-J)+J
      RHS(J)=RHS(J)-PSL(IB)*FACT
  632 CONTINUE
  630 CONTINUE
C--------------------------------- OVER --------------------------------
      DO 310 I=2,LI-1
      DO J=2,LJ-1
      II=(I-2)*(LJ-2)+J-1
      P(J,I)=RHS(II)
      ENDDO
  310 CONTINUE
      IF(KBOUND.EQ.1) P(LJ-1,LI-1)=0.0
      DO 312 I=1,LI
      P(1,I)=P(2,I)
      P(LJ,I)=P(LJ-1,I)
  312 CONTINUE
      DO 314 J=1,LJ
      P(J,1)=P(J,2)
      P(J,LI)=P(J,LI-1)
  314 CONTINUE
      DO 320 I=3,LI-1
      DO J=2,LJ-1
      U(J,I)=U(J,I)-DT*(P(J,I)-P(J,I-1))/XXS(I)/AA(J,I)
      ENDDO
  320 CONTINUE
      DO 322 I=2,LI-1
      DO J=3,LJ-1
      V(J,I)=V(J,I)-DT*(P(J,I)-P(J-1,I))/YYS(J)/BB(J,I)
      ENDDO
  322 CONTINUE
      DO 400 I=2,LI-1
      DO 401 J=2,LJ-1
      IF(I.EQ.LI-1.AND.J.EQ.LJ-1) GO TO 401
      AMDOT=YSWCH(J,I)*XXC(I)*YYC(J)*(RHONP(J,I)-RHO(J,I))/DT
     1     +YYC(J)*YSWCH(J,I)*(AA(J,I+1)*U(J,I+1)-AA(J,I)*U(J,I))
     2     +XXC(I)*0.5*((YSWCH(J+1,I)+YSWCH(J,I))*BB(J+1,I)*V(J+1,I)
     3     -(YSWCH(J,I)+YSWCH(J-1,I))*BB(J,I)*V(J,I))
      IF(ISKIP(J,I).GE.3.AND.ISKIP(J,I).NE.11) AMDOT=0.0
      RESDM=RESDM+DABS(AMDOT)
  401 CONTINUE
  400 CONTINUE
      RETURN
      END