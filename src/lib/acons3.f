      SUBROUTINE ACONS3(SIGMA)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LIJ=LI+LJ,
     1    LPD=55*LE-16*LE-2*LIJ-4)
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
     2 RHS1(LJ,LI),FP(LJ,LI),FPA(LJ+2,LI+2),EMU(LJ,LI),
     3 GI(LJ,LI),GJ(LJ,LI),DUM(LPD)
      SYM=DFLOAT(ISYM)
      DO 10 J=1,LJ
      DO I=1,LI
      FPA(J+1,I+1)=W(J,I)
      EMU(J,I)=REI*(VMU(J,I)+TMU(J,I)/SIGMA)
      ENDDO
   10 CONTINUE
      DO 12 J=1,LJ
      FPA(J+1,LI+2)=FBXP(4,J)
      FPA(J+1,1)=FBXM(4,J)
   12 CONTINUE
      DO 14 I=1,LI
      FPA(LJ+2,I+1)=FBYP(4,I)
      FPA(1,I+1)=FBYM(4,I)
   14 CONTINUE
      DO 102 I=2,LI
      DO J=2,LJ
      GI(J,I)=0.5*(EMU(J,I)+EMU(J,I-1))
      GJ(J,I)=0.5*(EMU(J,I)+EMU(J-1,I))
      ENDDO
  102 CONTINUE
C-QUICKEST--QUICKEST--QUICKEST--QUICKEST--QUICKEST--QUICKEST--QUICKEST--
C------------------------  ON WEST & EAST FACES  -----------------------
      DO 110 J=2,LJ-1
      DO I=2,LI
      IF(I.LE.LI-1) THEN
        IP=I+1
        RIP1=RHO(J,I+1)
        GIP1=EMU(J,I+1)
        ELSE
        IP=LI
        RIP1=FBXP(1,J)
        GIP1=EMU(J,LI)
        END IF
      IF(I.GE.3) THEN
        RIM2=RHO(J,I-2)
        GIM2=EMU(J,I-2)
        ELSE
        RIM2=FBXM(1,J)
        GIM2=EMU(J,1)
        END IF
      DXE=XXS(IP)
      DXW=XXS(I)
      DXWW=XXS(I-1)
      DXCP=XXC(I)
      DXCW=XXC(I-1)
C
      DYN=YYS(J+1)
      DYS=YYS(J)
      DYCP=YYC(J)
C
      IF(U(J,I).GE.0.0) THEN
         CRVRI=((RHO(J,I)-RHO(J,I-1))/DXW
     1        -(RHO(J,I-1)-RIM2)/DXWW)/DXCW
         CRVGI=((EMU(J,I)-EMU(J,I-1))/DXW
     1        -(EMU(J,I-1)-GIM2)/DXWW)/DXCW
         CRVRSI=((RHO(J+1,I-1)-RHO(J,I-1))/DYN
     1         -(RHO(J,I-1)-RHO(J-1,I-1))/DYS)/DYCP
         CRVGSI=((EMU(J+1,I-1)-EMU(J,I-1))/DYN
     1         -(EMU(J,I-1)-EMU(J-1,I-1))/DYS)/DYCP
         ELSE
         CRVRI=((RIP1-RHO(J,I))/DXE
     1        -(RHO(J,I)-RHO(J,I-1))/DXW)/DXCP
         CRVGI=((GIP1-EMU(J,I))/DXE
     1        -(EMU(J,I)-EMU(J,I-1))/DXW)/DXCP
         CRVRSI=((RHO(J+1,I)-RHO(J,I))/DYN
     1         -(RHO(J,I)-RHO(J-1,I))/DYS)/DYCP
         CRVGSI=((EMU(J+1,I)-EMU(J,I))/DYN
     1         -(EMU(J,I)-EMU(J-1,I))/DYS)/DYCP
         END IF
      RHOFI=0.5*(RHO(J,I)+RHO(J,I-1))
     1     -DXW*DXW*CRVRI/8.0+DYCP*DYCP*CRVRSI/24.0
      GMAFI=0.5*(EMU(J,I)+EMU(J,I-1))
     1     -DXW*DXW*CRVGI/8.0+DYCP*DYCP*CRVGSI/24.0
      AE(J,I-1)=RHOFI
      AW(J,I)=RHOFI
      GI(J,I)=GMAFI
      ENDDO
  110 CONTINUE
C----------------------  ON SOUTH & NORTH FACES  -----------------------
      DO 112 J=2,LJ
      DO I=2,LI-1
      IF(J.LE.LJ-1) THEN
        JP=J+1
        RJP1=RHO(J+1,I)
        GJP1=EMU(J+1,I)
        ELSE
        JP=LJ
        RJP1=FBYP(1,I)
        GJP1=EMU(LJ,I)
        END IF
      IF(J.GE.3) THEN
        RJM2=RHO(J-2,I)
        GJM2=EMU(J-2,I)
        ELSE
        RJM2=FBYM(1,I)
        GJM2=EMU(1,I)
        END IF
      DXE=XXS(I+1)
      DXW=XXS(I)
      DXCP=XXC(I)
C
      DYN=YYS(JP)
      DYS=YYS(J)
      DYSS=YYS(J-1)
      DYCP=YYC(JP)
      DYCS=YYC(J-1)
C
      IF(V(J,I).GE.0.0) THEN
         CRVRJ=((RHO(J,I)-RHO(J-1,I))/DYS
     1        -(RHO(J-1,I)-RJM2)/DYSS)/DYCS
         CRVGJ=((EMU(J,I)-EMU(J-1,I))/DYS
     1        -(EMU(J-1,I)-GJM2)/DYSS)/DYCS
         CRVRSJ=((RHO(J-1,I+1)-RHO(J-1,I))/DXE
     1         -(RHO(J-1,I)-RHO(J-1,I-1))/DXW)/DXCP
         CRVGSJ=((EMU(J-1,I+1)-EMU(J-1,I))/DXE
     1         -(EMU(J-1,I)-EMU(J-1,I-1))/DXW)/DXCP
         ELSE
         CRVRJ=((RJP1-RHO(J,I))/DYN
     1        -(RHO(J,I)-RHO(J-1,I))/DYS)/DYCP
         CRVGJ=((GJP1-EMU(J,I))/DYN
     1        -(EMU(J,I)-EMU(J-1,I))/DYS)/DYCP
         CRVRSJ=((RHO(J,I+1)-RHO(J,I))/DXE
     1         -(RHO(J,I)-RHO(J,I-1))/DXW)/DXCP
         CRVGSJ=((EMU(J,I+1)-EMU(J,I))/DXE
     1         -(EMU(J,I)-EMU(J,I-1))/DXW)/DXCP
         END IF
      RHOFJ=0.5*(RHO(J,I)+RHO(J-1,I))
     1     -DYS*DYS*CRVRJ/8.0+DXCP*DXCP*CRVRSJ/24.0
      GMAFJ=0.5*(EMU(J,I)+EMU(J-1,I))
     1     -DYS*DYS*CRVGJ/8.0+DXCP*DXCP*CRVGSJ/24.0
      AN(J-1,I)=RHOFJ
      AS(J,I)=RHOFJ
      GJ(J,I)=GMAFJ
      ENDDO
  112 CONTINUE
      DO 114 I=2,LI-1
      DO J=2,LJ-1
      AEE(J,I)=0.0
      IF(U(J,I+1).GE.0.0) AEE(J,I)=1.0
      AWW(J,I)=0.0
      IF(U(J,I).GE.0.0) AWW(J,I)=1.0
      ANN(J,I)=0.0
      IF(V(J+1,I).GE.0.0) ANN(J,I)=1.0
      ASS(J,I)=0.0
      IF(V(J,I).GE.0.0) ASS(J,I)=1.0
      ENDDO
  114 CONTINUE
C----------------------------COEFFICIENTS-------------------------------
      DO 210 J=2,LJ-1
      DO I=2,LI-1
      IP=I+2
      IF(IP.GT.LI) IP=LI
      JP=J+2
      IF(JP.GT.LJ) JP=LJ
      DXEE=XXS(IP)
      DXE=XXS(I+1)
      DXW=XXS(I)
      DXWW=XXS(I-1)
      DXCE=XXC(I+1)
      DXCP=XXC(I)
      DXCW=XXC(I-1)
C
      DYNN=YYS(JP)
      DYN=YYS(J+1)
      DYS=YYS(J)
      DYSS=YYS(J-1)
      DYCN=YYC(J+1)
      DYCP=YYC(J)
      DYCS=YYC(J-1)
C
      RHOP=RHO(J,I)
      VP=0.5*(V(J,I)+V(J+1,I))
      RHOUE=AE(J,I)*U(J,I+1)
      RHOUW=AW(J,I)*U(J,I)
      RHOVN=AN(J,I)*V(J+1,I)
      RHOVS=AS(J,I)*V(J,I)
      CNTE=U(J,I+1)*DT/DXE
      CNTW=U(J,I)*DT/DXW
      CNTN=V(J+1,I)*DT/DYN
      CNTS=V(J,I)*DT/DYS
      BY6E=1.0/6.0-GI(J,I+1)*DT/DXE/DXE-CNTE*CNTE/6.0
      BY6W=1.0/6.0-GI(J,I)*DT/DXW/DXW-CNTW*CNTW/6.0
      BY6N=1.0/6.0-GJ(J+1,I)*DT/DYN/DYN-CNTN*CNTN/6.0
      BY6S=1.0/6.0-GJ(J,I)*DT/DYS/DYS-CNTS*CNTS/6.0
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
      SYMT=SYM*DT/Y(J,I)
      AEE(J,I)=DTBRX*(-VARE*A1)
      AEP=DTBRX*(RHOUE/2.0-GI(J,I+1)/DXE)
      AE(J,I)=AEP-DTBRX*(RHOUE*CNTE/2.0+VARE*A2-VARW*B1)
      AWP=DTBRX*(-RHOUW/2.0-GI(J,I)/DXW)
      AW(J,I)=AWP-DTBRX*(VARE*A4+RHOUW*CNTW/2.0-VARW*B3)
      AWW(J,I)=DTBRX*VARW*B4
      ANN(J,I)=DTBRY*(-VARN*C1)
      ANP=DTBRY*(RHOVN/2.0-GJ(J+1,I)/DYN)-SYMT*0.5*EMU(J,I)/DYN
      AN(J,I)=ANP-DTBRY*(RHOVN*CNTN/2.0+VARN*C2-VARS*D1)
      ASP=DTBRY*(-RHOVS/2.0-GJ(J,I)/DYS)+SYMT*0.5*EMU(J,I)/DYS
      AS(J,I)=ASP-DTBRY*(VARN*C4+RHOVS*CNTS/2.0-VARS*D3)
      ASS(J,I)=DTBRY*VARS*D4
      AP(J,I)=RHO(J,I)-(AEP+AWP+ANP+ASP)
     1  +DTBRX*(RHOUE*CNTE/2.0+RHOUW*CNTW/2.0-VARE*A3+VARW*B2)
     2  +DTBRY*(RHOVN*CNTN/2.0+RHOVS*CNTS/2.0-VARN*C3+VARS*D2)
C--------------------SOURCE TERMS IN SWIRL EQUATION---------------------
      SW=-(0.5*RHO(J,I)*(V(J,I)+V(J+1,I))
     1   +0.5*((EMU(J+1,I)-EMU(J,I))/DYN+(EMU(J,I)-EMU(J-1,I))/DYS)
     2   +EMU(J,I)/Y(J,I))*W(J,I)/Y(J,I)
      RHS1(J,I)=RHO(J,I)*W(J,I)+DT*SW
      ENDDO
  210 CONTINUE
C----------------------------- INSERT BODY -----------------------------
      IF(NBODY.GE.1) THEN
          DO 300 N=1,NBODY
          IBODYM=IBDM(N)
          IBODYP=IBDP(N)
          JBODYM=JBDM(N)
          JBODYP=JBDP(N)
          DO 312 I=IBODYM,IBODYP
          IF(ISKIP(JBODYP,I).GE.10) THEN
          JA=JBODYP+1
          RHS1(JA,I)=RHS1(JA,I)+ASS(JA,I)*FPA(JA+1,I+1)
          AS(JA,I)=0.0
          ASS(JA,I)=0.0
          END IF
          IF(ISKIP(JBODYM,I).GE.10) THEN
          JA=JBODYM-1
          RHS1(JA,I)=RHS1(JA,I)+ANN(JA,I)*FPA(JA+1,I+1)
          AN(JA,I)=0.0
          ANN(JA,I)=0.0
          END IF
  312     CONTINUE
          DO 314 J=JBODYM,JBODYP
          IF(ISKIP(J,IBODYP).GE.10) THEN
          IA=IBODYP+1
          RHS1(J,IA)=RHS1(J,IA)+AWW(J,IA)*FPA(J+1,IA+1)
          AW(J,IA)=0.0
          AWW(J,IA)=0.0
          END IF
          IF(ISKIP(J,IBODYM).GE.10) THEN
          IA=IBODYM-1
          RHS1(J,IA)=RHS1(J,IA)+AEE(J,IA)*FPA(J+1,IA+1)
          AE(J,IA)=0.0
          AEE(J,IA)=0.0
          END IF
  314     CONTINUE
  300     CONTINUE
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
          AP(J,I)=RHO(J,I)
          ANN(J,I)=0.0
          AN(J,I)=0.0
          AS(J,I)=0.0
          ASS(J,I)=0.0
          AEE(J,I)=0.0
          AE(J,I)=0.0
          AW(J,I)=0.0
          AWW(J,I)=0.0
          RHS1(J,I)=AP(J,I)*W(J,I)
          ENDDO
  610     CONTINUE
  600     CONTINUE
          END IF
C-----------------------------------------------------------------------
      RETURN
      END