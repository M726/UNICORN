      SUBROUTINE ACONSV(SIGMA)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LIJ=LI+LJ,
     1    LPD=55*LE-14*LE-2*LIJ-4)
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
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/FINJ/NFINJ,IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30)
      COMMON/KEPS/ AKREF,EPREF,C1KE,C2KE,CMKE,CMKEQ,EE,AVON,PCON
      COMMON/DUMMY/ UTMP(LJ,LI),AP(LJ,LI),AEE(LJ,LI),AE(LJ,LI),
     1 AW(LJ,LI),AWW(LJ,LI),ANN(LJ,LI),AN(LJ,LI),AS(LJ,LI),ASS(LJ,LI),
     2 RHS1(LJ,LI),FP(LJ,LI),FPA(LJ+2,LI+2),EMU(LJ,LI),DUM(LPD)
      SYM=DFLOAT(ISYM)
      DO 10 I=1,LI
      DO J=1,LJ
      EMU(J,I)=REI*(VMU(J,I)+TMU(J,I)/SIGMA)
      FPA(J+1,I+1)=FP(J,I)
      ENDDO
   10 CONTINUE
      DO 12 J=1,LJ
      FPA(J+1,LI+2)=FBXP(3,J)
      FPA(J+1,1)=FBXM(3,J)
   12 CONTINUE
      DO 14 I=1,LI
      FPA(LJ+2,I+1)=FBYP(3,I)
      FPA(2,I+1)=FBYM(3,I)
   14 CONTINUE
C-QUICKEST--QUICKEST--QUICKEST--QUICKEST--QUICKEST--QUICKEST--QUICKEST--
C----------------------  ON SOUTH & NORTH FACES  -----------------------
      DO 110 J=3,LJ
      DO I=2,LI-1
      VJP1=FPA(J+2,I+1)
      VJM2=FPA(J-1,I+1)
      DXE=XXS(I+1)
      DXW=XXS(I)
      DXCP=XXC(I)
C
      DYN=YYC(J)
      DYS=YYC(J-1)
      DYSS=YYC(J-2)
      DYCP=YYS(J)
      DYCS=YYS(J-1)
C
      IF(V(J-1,I).GE.0.0) THEN
         CURVJ=((V(J,I)-V(J-1,I))/DYS-(V(J-1,I)-VJM2)/DYSS)/DYCS
         CURVSJ=((V(J-1,I+1)-V(J-1,I))/DXE
     1         -(V(J-1,I)-V(J-1,I-1))/DXW)/DXCP
         ELSE
         CURVJ=((VJP1-V(J,I))/DYN-(V(J,I)-V(J-1,I))/DYS)/DYCP
         CURVSJ=((V(J,I+1)-V(J,I))/DXE-(V(J,I)-V(J,I-1))/DXW)/DXCP
         END IF
      VFACE=0.5*(V(J,I)+V(J-1,I))
     1     -DYS*DYS*CURVJ/8.0+DXCP*DXCP*CURVSJ/24.0
      AN(J-1,I)=VFACE
      AS(J,I)=VFACE
      AE(J,I)=0.5*(U(J-1,I+1)+U(J,I+1))
      AW(J,I)=0.5*(U(J-1,I)+U(J,I))
      ENDDO
  110 CONTINUE
      DO 112 J=3,LJ-1
      DO I=2,LI-1
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
      DO 210 J=3,LJ-1
      DO I=2,LI-1
      IPP=I+2
      IF(IPP.GT.LI) IPP=LI
      DXEE=XXS(IPP)
      DXE=XXS(I+1)
      DXW=XXS(I)
      DXWW=XXS(I-1)
      DXCE=XXC(I+1)
      DXCP=XXC(I)
      DXCW=XXC(I-1)
C
      DYNN=YYC(J+1)
      DYN=YYC(J)
      DYS=YYC(J-1)
      DYSS=YYC(J-2)
      DYCN=YYS(J+1)
      DYCP=YYS(J)
      DYCS=YYS(J-1)
C
      RHOP=0.5*(RHO(J-1,I)+RHO(J,I))
      RHOPT=0.5*(RHONP(J-1,I)+RHONP(J,I))
      EMUP=0.5*(EMU(J-1,I)+EMU(J,I))
      VP=V(J,I)
      WP=0.5*(W(J-1,I)+W(J,I))
      RHOVN=RHO(J,I)*AN(J,I)
      RHOVS=RHO(J-1,I)*AS(J,I)
      RHOUE=0.25*(RHO(J,I)+RHO(J,I+1)+RHO(J-1,I)+RHO(J-1,I+1))*AE(J,I)
      RHOUW=0.25*(RHO(J,I-1)+RHO(J,I)+RHO(J-1,I-1)+RHO(J-1,I))*AW(J,I)
      EMUE=0.25*(EMU(J,I)+EMU(J,I+1)+EMU(J-1,I)+EMU(J-1,I+1))
      EMUW=0.25*(EMU(J,I-1)+EMU(J,I)+EMU(J-1,I-1)+EMU(J-1,I))
      CNTE=AE(J,I)*DT/DXE
      CNTW=AW(J,I)*DT/DXW
      CNTN=AN(J,I)*DT/DYN
      CNTS=AS(J,I)*DT/DYS
      BY6E=1.0/6.0-EMUE*DT/DXE/DXE-CNTE*CNTE/6.0
      BY6W=1.0/6.0-EMUW*DT/DXW/DXW-CNTW*CNTW/6.0
      BY6N=1.0/6.0-EMU(J,I)*DT/DYN/DYN-CNTN*CNTN/6.0
      BY6S=1.0/6.0-EMU(J-1,I)*DT/DYS/DYS-CNTS*CNTS/6.0
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
C
      DTBRX=DT/DXCP
      DTBRY=DT/DYCP
      VARE=RHOUE*DXE*DXE*BY6E
      VARW=RHOUW*DXW*DXW*BY6W
      VARN=RHOVN*DYN*DYN*BY6N
      VARS=RHOVS*DYS*DYS*BY6S
      SYMT=SYM*DT/0.5/(Y(J-1,I)+Y(J,I))
      AEE(J,I)=DTBRX*(-VARE*A1)
      AEP=DTBRX*(RHOUE/2.0-EMUE/DXE)
      AE(J,I)=AEP-DTBRX*(RHOUE*CNTE/2.0+VARE*A2-VARW*B1)
      AWP=DTBRX*(-RHOUW/2.0-EMUW/DXW)
      AW(J,I)=AWP-DTBRX*(VARE*A4+RHOUW*CNTW/2.0-VARW*B3)
      AWW(J,I)=DTBRX*VARW*B4
      ANN(J,I)=DTBRY*(-VARN*C1)
      ANP=DTBRY*(RHOVN/2.0-EMU(J,I)/DYN)-SYMT*0.5*EMUP/DYN
      AN(J,I)=ANP-DTBRY*(RHOVN*CNTN/2.0+VARN*C2-VARS*D1)
      ASP=DTBRY*(-RHOVS/2.0-EMU(J-1,I)/DYS)+SYMT*0.5*EMUP/DYS
      AS(J,I)=ASP-DTBRY*(VARN*C4+RHOVS*CNTS/2.0-VARS*D3)
      ASS(J,I)=DTBRY*VARS*D4
      AP(J,I)=RHOP-(AEP+AWP+ANP+ASP)
     1  +DTBRX*(RHOUE*CNTE/2.0+RHOUW*CNTW/2.0-VARE*A3+VARW*B2)
     2  +DTBRY*(RHOVN*CNTN/2.0+RHOVS*CNTS/2.0-VARN*C3+VARS*D2)
C-----------------------------SOURCE TERMS------------------------------
      SV=(EMU(J,I)*((V(J+1,I)-V(J,I))/DYN-2.0*(U(J,I+1)-U(J,I))/DXCP)
     1 -EMU(J-1,I)*((V(J,I)-V(J-1,I))/DYS
     2 -2.0*(U(J-1,I+1)-U(J-1,I))/DXCP))/3.0/DYCP
     3 +(EMUE*(U(J,I+1)-U(J-1,I+1))-EMUW*(U(J,I)-U(J-1,I)))/DYCP/DXCP
     4 +SYM*(EMUP*0.5*((V(J+1,I)-V(J,I))/DYN+(V(J,I)-V(J-1,I))/DYS)
     5 -4.0*EMUP*VP/Y(J,I)-2.0*VP*(EMU(J,I)-EMU(J-1,I))/DYCP)/3.0/Y(J,I)
      RHS1(J,I)=RHOP*V(J,I)+DT*SV
      ENDDO
  210 CONTINUE
C----------------------------- INSERT BODY -----------------------------
      IF(NBODY.GE.1) THEN
          DO 300 N=1,NBODY
          IBODYM=IBDM(N)
          IBODYP=IBDP(N)
          JBODYM=JBDM(N)
          JBODYP=JBDP(N)
          DO 310 I=IBODYM,IBODYP
          IF(ISKIP(JBODYP,I).GE.10) THEN
          JA=JBODYP+1
          RHS1(JA,I)=RHS1(JA,I)+AS(JA,I)*V(JA,I)+ASS(JA,I)*V(JA+1,I)
          AS(JA,I)=0.0
          ASS(JA,I)=0.0
          END IF
          IF(ISKIP(JBODYM,I).GE.10) THEN
          RHS1(JBODYM,I)=RHS1(JBODYM,I)+AN(JBODYM,I)*V(JBODYM,I)
     1        +ANN(JBODYM,I)*V(JBODYM-1,I)
          AN(JBODYM,I)=0.0
          ANN(JBODYM,I)=0.0
          END IF
  310     CONTINUE
          DO 320 J=JBODYM+1,JBODYP
          IF(ISKIP(J-1,IBODYP).GE.10) THEN
          IA=IBODYP+1
          RHS1(J,IA)=RHS1(J,IA)+AWW(J,IA)*V(J,IA)
          AW(J,IA)=0.0
          AWW(J,IA)=0.0
          END IF
          IF(ISKIP(J-1,IBODYM).GE.10) THEN
          IA=IBODYM-1
          RHS1(J,IA)=RHS1(J,IA)+AEE(J,IA)*V(J,IA)
          AE(J,IA)=0.0
          AEE(J,IA)=0.0
          END IF
  320     CONTINUE
  300     CONTINUE
          END IF
C-----------------------------------------------------------------------
      IF(IFLOW.EQ.2) THEN
C------------------CORRECTIONS FOR THE TURBULENT FLOWS------------------
      DO 410 J=3,LJ-1
      IF(JLFT(J).NE.1) GO TO 410
      DXW=XXS(2)
      DELX=0.5*DXW
      DXCP=XXC(2)
      DTBRX=DT/DXCP
      WALLR=0.25*(RHO(J-1,2)+RHO(J,2)+RHO(J-1,1)+RHO(J,1))
      WALLK=0.25*(AK(J-1,2)+AK(J,2)+AK(J-1,1)+AK(J,1))
      WALLM=0.25*(VMU(J-1,2)+VMU(J,2)+VMU(J-1,1)+VMU(J,1))
      TERM1=WALLR*CMKEQ*DSQRT(WALLK)
      FINV=TERM1*AVON/DLOG(EE*DELX*TERM1/WALLM/REI)
      EMUW=0.25*(EMU(J-1,2)+EMU(J,2)+EMU(J-1,1)+EMU(J,1))
      AW(J,2)=AW(J,2)+DTBRX*EMUW/DXW
      AP(J,2)=AP(J,2)-DTBRX*EMUW/DXW+DTBRX*FINV
  410 CONTINUE
      DO 420 J=3,LJ-1
      IF(JRGT(J).NE.1) GO TO 420
      DXE=XXS(LI-1)
      DELX=0.5*DXE
      DXCP=XXC(LI-1)
      DTBRX=DT/DXCP
      WALLR=0.25*(RHO(J-1,LI)+RHO(J,LI)+RHO(J-1,LI-1)+RHO(J,LI-1))
      WALLK=0.25*(AK(J-1,LI)+AK(J,LI)+AK(J-1,LI-1)+AK(J,LI-1))
      WALLM=0.25*(VMU(J-1,LI)+VMU(J,LI)+VMU(J-1,LI-1)+VMU(J,LI-1))
      TERM1=WALLR*CMKEQ*DSQRT(WALLK)
      FINV=TERM1*AVON/DLOG(EE*DELX*TERM1/WALLM/REI)
      EMUE=0.25*(EMU(J-1,LI)+EMU(J,LI)+EMU(J-1,LI-1)+EMU(J,LI-1))
      AE(J,LI-1)=AE(J,LI-1)+DTBRX*EMUE/DXE
      AP(J,LI-1)=AP(J,LI-1)-DTBRX*EMUE/DXE+DTBRX*FINV
  420 CONTINUE
      IF(NBODY.GE.1) THEN
        DO 500 N=1,NBODY
        IBODYM=IBDM(N)
        IBODYP=IBDP(N)
        JBODYM=JBDM(N)
        JBODYP=JBDP(N)
        IAM=IBODYM-1
        IF(IAM.LT.1) IAM=1
        IAP=IBODYP+1
        IF(IAP.GT.LI) IAP=LI
        DO 510 J=JBODYM+1,JBODYP
        IF(ISKIP(J-1,IBODYP).GE.10) THEN
        DXW=XXS(IAP)
        DELX=0.5*DXW
        DXCP=XXC(IAP)
        DTBRX=DT/DXCP
        WALLR=0.25*(RHO(J-1,IAP)+RHO(J,IAP)+RHO(J-1,IAP-1)+RHO(J,IAP-1))
        WALLK=0.25*(AK(J-1,IAP)+AK(J,IAP)+AK(J-1,IAP-1)+AK(J,IAP-1))
        WALLM=0.25*(VMU(J-1,IAP)+VMU(J,IAP)+VMU(J-1,IAP-1)+VMU(J,IAP-1))
        TERM1=WALLR*CMKEQ*DSQRT(WALLK)
        FINV=TERM1*AVON/DLOG(EE*DELX*TERM1/WALLM/REI)
        EMUW=0.25*(EMU(J-1,IAP)+EMU(J,IAP)+EMU(J-1,IAP-1)+EMU(J,IAP-1))
        AW(J,IAP)=AW(J,IAP)+DTBRX*EMUW/DXW
        AP(J,IAP)=AP(J,IAP)-DTBRX*EMUW/DXW+DTBRX*FINV
        END IF
        IF(ISKIP(J-1,IBODYM).GE.10) THEN
        DXE=XXS(IAM)
        DELX=0.5*DXE
        DXCP=XXC(IAM)
        DTBRX=DT/DXCP
        WALLR=0.25*(RHO(J-1,IAM+1)+RHO(J,IAM+1)+RHO(J-1,IAM)+RHO(J,IAM))
        WALLK=0.25*(AK(J-1,IAM+1)+AK(J,IAM+1)+AK(J-1,IAM)+AK(J,IAM))
        WALLM=0.25*(VMU(J-1,IAM+1)+VMU(J,IAM+1)+VMU(J-1,IAM)+VMU(J,IAM))
        TERM1=WALLR*CMKEQ*DSQRT(WALLK)
        FINV=TERM1*AVON/DLOG(EE*DELX*TERM1/WALLM/REI)
        EMUE=0.25*(EMU(J-1,IAM+1)+EMU(J,IAM+1)+EMU(J-1,IAM)+EMU(J,IAM))
        AE(J,IAM)=AE(J,IAM)+DTBRX*EMUE/DXE
        AP(J,IAM)=AP(J,IAM)-DTBRX*EMUE/DXE+DTBRX*FINV
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
          AP(J,I)=0.5*(RHO(J-1,I)+RHO(J,I))
          ANN(J,I)=0.0
          AN(J,I)=0.0
          AS(J,I)=0.0
          ASS(J,I)=0.0
          AEE(J,I)=0.0
          AE(J,I)=0.0
          AW(J,I)=0.0
          AWW(J,I)=0.0
          RHS1(J,I)=AP(J,I)*V(J,I)
          ENDDO
  610     CONTINUE
  600     CONTINUE
          END IF
C-----------------------------------------------------------------------
      RETURN
      END