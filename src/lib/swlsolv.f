      SUBROUTINE SWLSOLV(ISOR,RELXW,TOLRW,SIGW,RESDW)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LIJ=LI+LJ,
     1     LPD=55*LE-8*LE) 
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LE,LSP),U(LJ,LI),V(LJ,LI),
     1  W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB05/ RHONP(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDFSP(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/CB09/ IBOT(LI),ITOP(LI),JLFT(LJ),JRGT(LJ),
     1       FBOT(8+LSP,LI),FTOP(8+LSP,LI),FLFT(8+LSP,LJ),FRGT(8+LSP,LJ)
      COMMON/CB10/ FBXP(4,LJ),FBXM(4,LJ),FBYP(4,LI),FBYM(4,LI)
      COMMON/BODY/IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10),NBODY
      COMMON/FINJ/IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30),NFINJ
      COMMON/KEPS/ AKREF,EPREF,C1KE,C2KE,CMKE,CMKEQ,EE,AVON,PCON
      COMMON/DUMMY/ AE(LJ,LI),AW(LJ,LI),AN(LJ,LI),AS(LJ,LI),
     1 AP(LJ,LI),FPA(LJ,LI),EMU(LJ,LI),RHS1(LJ,LI),DUM(LPD)
      SYM=DFLOAT(ISYM)
      ENORM=DFLOAT(LE)
      DO 10 I=1,LI
      DO J=1,LJ
      EMU(J,I)=REI*(VMU(J,I)+TMU(J,I)/SIGW)
      FPA(J,I)=W(J,I)
      ENDDO
   10 CONTINUE
C----------------------------COEFFICIENTS-------------------------------
      DO 20 I=2,LI-1
      DO J=2,LJ-1
      DXE=XXS(I+1)
      DXW=XXS(I)
      DXCE=XXC(I+1)
      DXCP=XXC(I)
      DXCW=XXC(I-1)
C
      DYN=YYS(J+1)
      DYS=YYS(J)
      DYCN=YYC(J+1)
      DYCP=YYC(J)
      DYCS=YYC(J-1)
C
      DTBRX=DT/DXCP
      DTBRY=DT/DYCP
      RHOP=RHO(J,I)
      RHOUE=0.5*(RHO(J,I+1)+RHO(J,I))*U(J,I+1)
      RHOUW=0.5*(RHO(J,I)+RHO(J,I-1))*U(J,I)
      RHOVN=0.5*(RHO(J+1,I)+RHO(J,I))*V(J+1,I)
      RHOVS=0.5*(RHO(J,I)+RHO(J-1,I))*V(J,I)
      SYMT=SYM*DT*0.5*EMU(J,I)/Y(J,I)
C
      ACONV=ABS(RHOUE)
      ADIFU=(EMU(J,I+1)+EMU(J,I))/DXE
      IF(ADIFU.GE.ACONV) ACONV=ADIFU
      AE(J,I)=DTBRX*0.5*(RHOUE-ACONV)
C
      ACONV=ABS(RHOUW)
      ADIFU=(EMU(J,I)+EMU(J,I-1))/DXW
      IF(ADIFU.GE.ACONV) ACONV=ADIFU
      AW(J,I)=DTBRX*0.5*(-RHOUW-ACONV)
C
      ACONV=ABS(RHOVN)
      ADIFU=(EMU(J+1,I)+EMU(J,I))/DYN
      IF(ADIFU.GE.ACONV) ACONV=ADIFU
      AN(J,I)=DTBRY*0.5*(RHOVN-ACONV)-SYMT/DYN
C
      ACONV=ABS(RHOVS)
      ADIFU=(EMU(J,I)+EMU(J-1,I))/DYS
      IF(ADIFU.GE.ACONV) ACONV=ADIFU
      AS(J,I)=DTBRY*0.5*(-RHOVS-ACONV)+SYMT/DYS
C
      AP(J,I)=RHO(J,I)-(AE(J,I)+AW(J,I)+AN(J,I)+AS(J,I))
C--------------------SOURCE TERMS IN SWIRL EQUATION---------------------
      SW=-(0.5*RHO(J,I)*(V(J,I)+V(J+1,I))
     1   +0.5*((EMU(J+1,I)-EMU(J,I))/DYN+(EMU(J,I)-EMU(J-1,I))/DYS)
     2   +EMU(J,I)/Y(J,I))*W(J,I)/Y(J,I)
      RHS1(J,I)=RHO(J,I)*W(J,I)+DT*SW
      ENDDO
   20 CONTINUE
C----------------------------- INSERT BODY -----------------------------
      IF(NBODY.GE.1) THEN
          DO 28 N=1,NBODY
          IBODYM=IBDM(N)
          IBODYP=IBDP(N)
          JBODYM=JBDM(N)
          JBODYP=JBDP(N)
          DO 24 I=IBODYM,IBODYP
          IF(ISKIP(JBODYP,I).GE.10) THEN
          JA=JBODYP+1
          RHS1(JA,I)=RHS1(JA,I)+AS(JA,I)*FPA(JA,I)
          AS(JA,I)=0.0
          END IF
          IF(ISKIP(JBODYM,I).GE.10) THEN
          JA=JBODYM-1
          RHS1(JA,I)=RHS1(JA,I)-AN(JA,I)*FPA(JA,I)
          AN(JA,I)=0.0
          END IF
   24     CONTINUE
          DO 26 J=JBODYM,JBODYP
          IF(ISKIP(J,IBODYP).GE.10) THEN
          IA=IBODYP+1
          RHS1(J,IA)=RHS1(J,IA)-AW(J,IA)*FPA(J,IA)
          AW(J,IA)=0.0
          END IF
          IF(ISKIP(J,IBODYM).GE.10) THEN
          IA=IBODYM-1
          RHS1(J,IA)=RHS1(J,IA)-AE(J,IA)*FPA(J,IA)
          AE(J,IA)=0.0
          END IF
   26     CONTINUE
   28     CONTINUE
          END IF
C-----------------------------------------------------------------------
      SOR1=RELXW
      SOR2=(1.0-RELXW)
      DO 100 ISORP=1,ISOR
      RSD=0
C----------------------   POINT RELAXATION SCHEME  ---------------------
      DO 210 I=2,LI-1
      DO 212 J=2,LJ-1
      IF(ISKIP(J,I).GE.3) GO TO 212
      FNEW=SOR1*(RHS1(J,I)-AN(J,I)*W(J+1,I)-AS(J,I)*W(J-1,I)
     1 -AE(J,I)*W(J,I+1)-AW(J,I)*W(J,I-1))/AP(J,I)+SOR2*W(J,I)
      RSD=RSD+ABS(W(J,I)-FNEW)
      W(J,I)=FNEW
  212 CONTINUE
  210 CONTINUE
      IF((RSD/ENORM).LE.TOLRW) GO TO 101
  100 CONTINUE
      ISORP=ISORP-1
  101 CONTINUE
      ISOR=ISORP
      RESDW=0.0
      DO 400 I=1,LI
      DO J=1,LJ
      RESDW=RESDW+ABS(FPA(J,I)-W(J,I))
      ENDDO
  400 CONTINUE
      RETURN
      END