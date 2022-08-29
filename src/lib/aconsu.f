      SUBROUTINE ACONSU(SIGMA)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LIJ=LI+LJ,
     1     LPD=55*LE-14*LE-2*LIJ-4) 
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
      COMMON/DUMMY/ UTMP(LJ,LI),AP(LJ,LI),AEE(LJ,LI),AE(LJ,LI),
     1 AW(LJ,LI),AWW(LJ,LI),ANN(LJ,LI),AN(LJ,LI),AS(LJ,LI),ASS(LJ,LI),
     2 RHS1(LJ,LI),FP(LJ,LI),FPA(LJ+2,LI+2),EMU(LJ,LI),DUM(LPD)
      SYM=DFLOAT(ISYM)
      RINF=RREF/RSTR
      DO 10 I=1,LI
      DO J=1,LJ
      EMU(J,I)=REI*(VMU(J,I)+TMU(J,I)/SIGMA)
      FPA(J+1,I+1)=FP(J,I)
      ENDDO
   10 CONTINUE
      DO 12 J=1,LJ
      FPA(J+1,LI+2)=FBXP(2,J)
      FPA(J+1,2)=FBXM(2,J)
   12 CONTINUE
      DO 14 I=1,LI
      FPA(LJ+2,I+1)=FBYP(2,I)
      FPA(1,I+1)=FBYM(2,I)
   14 CONTINUE
C-QUICKEST--QUICKEST--QUICKEST--QUICKEST--QUICKEST--QUICKEST--QUICKEST--
C------------------------  ON WEST & EAST FACES  -----------------------
      DO 110 J=2,LJ-1
      DO I=3,LI
      UIP1=FPA(J+1,I+2)
      UIM2=FPA(J+1,I-1)
C
      DXE=XXC(I)
      DXW=XXC(I-1)
      DXWW=XXC(I-2)
      DXCP=XXS(I)
      DXCW=XXS(I-1)
C
      DYN=YYS(J+1)
      DYS=YYS(J)
      DYCP=YYC(J)
C
      IF(U(J,I-1).GE.0.0) THEN
         CURVI=((U(J,I)-U(J,I-1))/DXW-(U(J,I-1)-UIM2)/DXWW)/DXCW
         CURVSI=((U(J+1,I-1)-U(J,I-1))/DYN
     1         -(U(J,I-1)-U(J-1,I-1))/DYS)/DYCP
         ELSE
         CURVI=((UIP1-U(J,I))/DXE-(U(J,I)-U(J,I-1))/DXW)/DXCP
         CURVSI=((U(J+1,I)-U(J,I))/DYN-(U(J,I)-U(J-1,I))/DYS)/DYCP
         END IF
      UFACE=0.5*(U(J,I)+U(J,I-1))
     1     -DXW*DXW*CURVI/8.0+DYCP*DYCP*CURVSI/24.0
      AE(J,I-1)=UFACE
      AW(J,I)=UFACE
      AN(J,I)=0.5*(V(J+1,I-1)+V(J+1,I))
      AS(J,I)=0.5*(V(J,I-1)+V(J,I))
      ENDDO
  110 CONTINUE
      DO 112 I=3,LI-1
      DO J=2,LJ-1
      AEE(J,I)=0.0
      IF(AE(J,I).GE.0.0) AEE(J,I)=1.0
      AWW(J,I)=0.0
      IF(AW(J,I).GE.0.0) AWW(J,I)=1.0
      ANN(J,I)=0.0
      IF(AN(J,I).GE.0.0) ANN(J,I)=1.0
      ASS(J,I)=0.0
      IF(AS(J,I).GE.0.0) ASS(J,I)=1.0
      ENDDO
  112 CONTINUE
C----------------------------COEFFICIENTS-------------------------------
      DO 210 J=2,LJ-1
      DO I=3,LI-1
      JPP=J+2
      IF(JPP.GT.LJ) JPP=LJ
      DXEE=XXC(I+1)
      DXE=XXC(I)
      DXW=XXC(I-1)
      DXWW=XXC(I-2)
      DXCE=XXS(I+1)
      DXCP=XXS(I)
      DXCW=XXS(I-1)
C
      DYNN=YYS(JPP)
      DYN=YYS(J+1)
      DYS=YYS(J)
      DYSS=YYS(J-1)
      DYCN=YYC(J+1)
      DYCP=YYC(J)
      DYCS=YYC(J-1)
C
      RHOP=0.5*(RHO(J,I-1)+RHO(J,I))
      RHOPT=0.5*(RHONP(J,I-1)+RHONP(J,I))
      EMUP=0.5*(EMU(J,I-1)+EMU(J,I))
      VP=0.5*(V(J+1,I)+V(J,I))
      RHOUE=RHO(J,I)*AE(J,I)
      RHOUW=RHO(J,I-1)*AW(J,I)
      RHOVN=0.25*(RHO(J+1,I-1)+RHO(J+1,I)+RHO(J,I-1)+RHO(J,I))*AN(J,I)
      RHOVS=0.25*(RHO(J,I-1)+RHO(J,I)+RHO(J-1,I-1)+RHO(J-1,I))*AS(J,I)
      EMUN=0.25*(EMU(J+1,I-1)+EMU(J+1,I)+EMU(J,I-1)+EMU(J,I))
      EMUS=0.25*(EMU(J,I-1)+EMU(J,I)+EMU(J-1,I-1)+EMU(J-1,I))
      CNTE=AE(J,I)*DT/DXE
      CNTW=AW(J,I)*DT/DXW
      CNTN=AN(J,I)*DT/DYN
      CNTS=AS(J,I)*DT/DYS
      BY6E=1.0/6.0-EMU(J,I)*DT/DXE/DXE-CNTE*CNTE/6.0
      BY6W=1.0/6.0-EMU(J,I-1)*DT/DXW/DXW-CNTW*CNTW/6.0
      BY6N=1.0/6.0-EMUN*DT/DYN/DYN-CNTN*CNTN/6.0
      BY6S=1.0/6.0-EMUS*DT/DYS/DYS-CNTS*CNTS/6.0
C
      SW1=AEE(J,I)
      SW2=1.0-SW1
      A1=SW2/DXCE/DXEE
      A2=SW1/DXCP/DXE-(A1+SW2/DXCE/DXE)
      A4=SW1/DXCP/DXW
      A3=-(A1+A2+A4)
C
      SW1=AWW(J,I)
      SW2=1.0-SW1
      B1=SW2/DXCP/DXE
      B2=SW1/DXCW/DXW-(B1+SW2/DXCP/DXW)
      B4=SW1/DXCW/DXWW
      B3=-(B1+B2+B4)
C
      SW1=ANN(J,I)
      SW2=1.0-SW1
      C1=SW2/DYCN/DYNN
      C2=SW1/DYCP/DYN-(C1+SW2/DYCN/DYN)
      C4=SW1/DYCP/DYS
      C3=-(C1+C2+C4)
C
      SW1=ASS(J,I)
      SW2=1.0-SW1
      D1=SW2/DYCP/DYN
      D2=SW1/DYCS/DYS-(D1+SW2/DYCP/DYS)
      D4=SW1/DYCS/DYSS
      D3=-(D1+D2+D4)
      DTBRX=DT/DXCP
      DTBRY=DT/DYCP
      VARE=RHOUE*DXE*DXE*BY6E
      VARW=RHOUW*DXW*DXW*BY6W
      VARN=RHOVN*DYN*DYN*BY6N
      VARS=RHOVS*DYS*DYS*BY6S
      SYMT=SYM*DT/Y(J,I)
      AEE(J,I)=DTBRX*(-VARE*A1)
      AEP=DTBRX*(RHOUE/2.0-EMU(J,I)/DXE)
      AE(J,I)=AEP-DTBRX*(RHOUE*CNTE/2.0+VARE*A2-VARW*B1)
      AWP=DTBRX*(-RHOUW/2.0-EMU(J,I-1)/DXW)
      AW(J,I)=AWP-DTBRX*(VARE*A4+RHOUW*CNTW/2.0-VARW*B3)
      AWW(J,I)=DTBRX*VARW*B4
      ANN(J,I)=DTBRY*(-VARN*C1)
      ANP=DTBRY*(RHOVN/2.0-EMUN/DYN)-SYMT*0.5*EMUP/DYN
      AN(J,I)=ANP-DTBRY*(RHOVN*CNTN/2.0+VARN*C2-VARS*D1)
      ASP=DTBRY*(-RHOVS/2.0-EMUS/DYS)+SYMT*0.5*EMUP/DYS
      AS(J,I)=ASP-DTBRY*(VARN*C4+RHOVS*CNTS/2.0-VARS*D3)
      ASS(J,I)=DTBRY*VARS*D4
      AP(J,I)=RHOP-(AEP+AWP+ANP+ASP)
     1  +DTBRX*(RHOUE*CNTE/2.0+RHOUW*CNTW/2.0-VARE*A3+VARW*B2)
     2  +DTBRY*(RHOVN*CNTN/2.0+RHOVS*CNTS/2.0-VARN*C3+VARS*D2)
C-----------------------------SOURCE TERMS------------------------------
      SU=(EMU(J,I)*((U(J,I+1)-U(J,I))/DXE-2.0*(V(J+1,I)-V(J,I))/DYCP)
     1 -EMU(J,I-1)*((U(J,I)-U(J,I-1))/DXW
     2 -2.0*(V(J+1,I-1)-V(J,I-1))/DYCP))/3.0/DXCP
     3 +(EMUN*(V(J+1,I)-V(J+1,I-1))-EMUS*(V(J,I)-V(J,I-1)))/DYCP/DXCP
     4 +SYM*(EMUP*0.5*((V(J+1,I)-V(J+1,I-1))/DXE+(V(J,I)-V(J,I-1))/DXCP)
     5 -2.0*VP*(EMU(J,I)-EMU(J,I-1))/DXCP)/Y(J,I)/3.0
      IF(IGRAV.NE.0) SU=SU-BETA4*(0.5*(RHONP(J,I)+RHONP(J,I-1))-RINF)
      RHS1(J,I)=RHOP*U(J,I)+DT*SU
      ENDDO
  210 CONTINUE
C----------------------------- INSERT BODY -----------------------------
      IF(NBODY.GE.1) THEN
          DO 300 N=1,NBODY
          IBODYM=IBDM(N)
          IBODYP=IBDP(N)
          JBODYM=JBDM(N)
          JBODYP=JBDP(N)
          DO 310 I=IBODYM+1,IBODYP
          IF(ISKIP(JBODYP,I-1).GE.10) THEN
          JA=JBODYP+1
          RHS1(JA,I)=RHS1(JA,I)+ASS(JA,I)*U(JA,I)
          AS(JA,I)=0.0
          ASS(JA,I)=0.0
          END IF
          IF(ISKIP(JBODYM,I-1).GE.10) THEN
          JA=JBODYM-1
          RHS1(JA,I)=RHS1(JA,I)+ANN(JA,I)*U(JA,I)
          AN(JA,I)=0.0
          ANN(JA,I)=0.0
          END IF
  310     CONTINUE
          DO 312 J=JBODYM,JBODYP
          IF(ISKIP(J,IBODYP).GE.10) THEN
          IA=IBODYP+1
          RHS1(J,IA)=RHS1(J,IA)+AW(J,IA)*U(J,IA)+AWW(J,IA)*U(J,IA+1)
          AW(J,IA)=0.0
          AWW(J,IA)=0.0
          END IF
          IF(ISKIP(J,IBODYM).GE.10) THEN
          RHS1(J,IBODYM)=RHS1(J,IBODYM)+AE(J,IBODYM)*U(J,IBODYM)
     1        +AEE(J,IBODYM)*U(J,IBODYM-1)
          AE(J,IBODYM)=0.0
          AEE(J,IBODYM)=0.0
          END IF
  312     CONTINUE
  300     CONTINUE
          END IF
C-----------------------------------------------------------------------
      IF(IFLOW.EQ.2) THEN
C---------------------CORRECTIONS FOR TURBULENT FLOWS-------------------
      DO 410 I=3,LI-1
      IF(IBOT(I).NE.1) GO TO 410
      DYS=YYS(2)
      DELY=0.5*DYS
      DYCP=YYC(2)
      DTBRY=DT/DYCP
      WALLR=0.25*(RHO(2,I-1)+RHO(2,I)+RHO(1,I-1)+RHO(1,I))
      WALLK=0.25*(AK(2,I-1)+AK(2,I)+AK(1,I-1)+AK(1,I))
      WALLM=0.25*(VMU(2,I-1)+VMU(2,I)+VMU(1,I-1)+VMU(1,I))
      TERM1=WALLR*CMKEQ*DSQRT(WALLK)
      FINV=TERM1*AVON/DLOG(EE*DELY*TERM1/WALLM/REI)
      EMUS=0.25*(EMU(2,I-1)+EMU(2,I)+EMU(1,I-1)+EMU(1,I))
      AS(2,I)=AS(2,I)+DTBRY*EMUS/DYS
      AP(2,I)=AP(2,I)-DTBRY*EMUS/DYS+DTBRY*FINV
  410 CONTINUE
      DO 420 I=3,LI-1
      IF(ITOP(I).NE.1) GO TO 420
      DYN=YYS(LJ-1)
      DELY=0.5*DYN
      DYCP=YYC(LJ-1)
      DTBRY=DT/DYCP
      WALLR=0.25*(RHO(LJ,I-1)+RHO(LJ,I)+RHO(LJ-1,I-1)+RHO(LJ-1,I))
      WALLK=0.25*(AK(LJ,I-1)+AK(LJ,I)+AK(LJ-1,I-1)+AK(LJ-1,I))
      WALLM=0.25*(VMU(LJ,I-1)+VMU(LJ,I)+VMU(LJ-1,I-1)+VMU(LJ-1,I))
      TERM1=WALLR*CMKEQ*DSQRT(WALLK)
      FINV=TERM1*AVON/DLOG(EE*DELY*TERM1/WALLM/REI)
      EMUN=0.25*(EMU(LJ,I-1)+EMU(LJ,I)+EMU(LJ-1,I-1)+EMU(LJ-1,I))
      AN(LJ-1,I)=AN(LJ-1,I)+DTBRY*EMUN/DYN
      AP(LJ-1,I)=AP(LJ-1,I)-DTBRY*EMUN/DYN+DTBRY*FINV
  420 CONTINUE
      IF(NBODY.GE.1) THEN
        DO 500 N=1,NBODY
        IBODYM=IBDM(N)
        IBODYP=IBDP(N)
        JBODYM=JBDM(N)
        JBODYP=JBDP(N)
        JAM=JBODYM-1
        IF(JAM.LT.1) JAM=1
        JAP=JBODYP+1
        IF(JAP.GT.LJ) JAP=LJ
        DO 510 I=IBODYM+1,IBODYP
        IF(ISKIP(JBODYP,I-1).GE.10) THEN
        DYS=YYS(JAP)
        DELY=0.5*DYS
        DYCP=YYC(JAP)
        DTBRY=DT/DYCP
        WALLR=0.25*(RHO(JAP,I-1)+RHO(JAP,I)+RHO(JAP-1,I-1)+RHO(JAP-1,I))
        WALLK=0.25*(AK(JAP,I-1)+AK(JAP,I)+AK(JAP-1,I-1)+AK(JAP-1,I))
        WALLM=0.25*(VMU(JAP,I-1)+VMU(JAP,I)+VMU(JAP-1,I-1)+VMU(JAP-1,I))
        TERM1=WALLR*CMKEQ*DSQRT(WALLK)
        FINV=TERM1*AVON/DLOG(EE*DELY*TERM1/WALLM/REI)
        EMUS=0.25*(EMU(JAP,I-1)+EMU(JAP,I)+EMU(JAP-1,I-1)+EMU(JAP-1,I))
        AS(JAP,I)=AS(JAP,I)+DTBRY*EMUS/DYS
        AP(JAP,I)=AP(JAP,I)-DTBRY*EMUS/DYS+DTBRY*FINV
        END IF
        IF(ISKIP(JBODYM,I-1).GE.10) THEN
        DYN=YYS(JAM)
        DELY=0.5*DYN
        DYCP=YYC(JAM)
        DTBRY=DT/DYCP
        WALLR=0.25*(RHO(JAM+1,I-1)+RHO(JAM+1,I)+RHO(JAM,I-1)+RHO(JAM,I))
        WALLK=0.25*(AK(JAM+1,I-1)+AK(JAM+1,I)+AK(JAM,I-1)+AK(JAM,I))
        WALLM=0.25*(VMU(JAM+1,I-1)+VMU(JAM+1,I)+VMU(JAM,I-1)+VMU(JAM,I))
        TERM1=WALLR*CMKEQ*DSQRT(WALLK)
        FINV=TERM1*AVON/DLOG(EE*DELY*TERM1/WALLM/REI)
        EMUN=0.25*(EMU(JAM+1,I-1)+EMU(JAM+1,I)+EMU(JAM,I-1)+EMU(JAM,I))
        AN(JAM,I)=AN(JAM,I)+DTBRY*EMUN/DYN
        AP(JAM,I)=AP(JAM,I)-DTBRY*EMUN/DYN+DTBRY*FINV
        END IF
  510   CONTINUE
  500   CONTINUE
        END IF
      END IF
C------------------------ INSERT INJECTIONS ----------------------------
      IF(NFINJ.GE.1) THEN
          DO 600 N=1,NFINJ
          IFINJM=IFIM(N)
          IFINJP=IFIP(N)
          JFINJM=JFIM(N)
          JFINJP=JFIP(N)
          DO 610 I=IFINJM,IFINJP
          DO J=JFINJM,JFINJP
          AP(J,I)=0.5*(RHO(J,I-1)+RHO(J,I))
          ANN(J,I)=0.0
          AN(J,I)=0.0
          AS(J,I)=0.0
          ASS(J,I)=0.0
          AEE(J,I)=0.0
          AE(J,I)=0.0
          AW(J,I)=0.0
          AWW(J,I)=0.0
          RHS1(J,I)=AP(J,I)*U(J,I)
          ENDDO
  610     CONTINUE
  600     CONTINUE
          END IF
C-----------------------------------------------------------------------
      RETURN
      END SUBROUTINE ACONSU