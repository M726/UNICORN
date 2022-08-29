      SUBROUTINE KESOLV(ISOR,RELXKE,TOLRKE,SIGK,SIGE,RESDK,RESDE)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LIJ=LI+LJ,
     1    LPD=55*LE-15*LE)
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LE,LSP),U(LJ,LI),V(LJ,LI),
     1  W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDFSP(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/CB09/ IBOT(LI),ITOP(LI),JLFT(LJ),JRGT(LJ),
     1       FBOT(8+LSP,LI),FTOP(8+LSP,LI),FLFT(8+LSP,LJ),FRGT(8+LSP,LJ)
      COMMON/CB10/ FBXP(4,LJ),FBXM(4,LJ),FBYP(4,LI),FBYM(4,LI)
      COMMON/BODY/IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10),NBODY
      COMMON/KEPS/ AKREF,EPREF,C1KE,C2KE,CMKE,CMKEQ,EE,AVON,PCON
      COMMON/DUMMY/ AW(LJ,LI),AE(LJ,LI),AS(LJ,LI),AN(LJ,LI),AP(LJ,LI),
     1 BW(LJ,LI),BE(LJ,LI),BS(LJ,LI),BN(LJ,LI),BP(LJ,LI),EMU1(LJ,LI),
     2 EMU2(LJ,LI),VORT(LJ,LI),FPA(LJ,LI),FPB(LJ,LI),DUM(LPD)
      ENORM=DFLOAT(LE)
      SYM=DFLOAT(ISYM)
      DO 10 I=1,LI
      DO J=1,LJ
      EMU1(J,I)=REI*(VMU(J,I)+TMU(J,I)/SIGK)
      EMU2(J,I)=REI*(VMU(J,I)+TMU(J,I)/SIGE)
      FPA(J,I)=AK(J,I)
      FPB(J,I)=EPS(J,I)
      ENDDO
   10 CONTINUE
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
      SYM1=SYM*DT*0.5*EMU1(J,I)/Y(J,I)
      SYM2=SYM*DT*0.5*EMU2(J,I)/Y(J,I)
C
      ACONV=DABS(RHOUE)
      BCONV=ACONV
      ADIFU=(EMU1(J,I+1)+EMU1(J,I))/DXE
      BDIFU=(EMU2(J,I+1)+EMU2(J,I))/DXE
      IF(ADIFU.GE.ACONV) ACONV=ADIFU
      IF(BDIFU.GE.BCONV) BCONV=BDIFU
      AE(J,I)=DTBRX*0.5*(RHOUE-ACONV)
      BE(J,I)=DTBRX*0.5*(RHOUE-BCONV)
C
      ACONV=DABS(RHOUW)
      BCONV=ACONV
      ADIFU=(EMU1(J,I)+EMU1(J,I-1))/DXW
      BDIFU=(EMU2(J,I)+EMU2(J,I-1))/DXW
      IF(ADIFU.GE.ACONV) ACONV=ADIFU
      IF(BDIFU.GE.BCONV) BCONV=BDIFU
      AW(J,I)=DTBRX*0.5*(-RHOUW-ACONV)
      BW(J,I)=DTBRX*0.5*(-RHOUW-BCONV)
C
      ACONV=DABS(RHOVN)
      BCONV=ACONV
      ADIFU=(EMU1(J+1,I)+EMU1(J,I))/DYN
      BDIFU=(EMU2(J+1,I)+EMU2(J,I))/DYN
      IF(ADIFU.GE.ACONV) ACONV=ADIFU
      IF(BDIFU.GE.BCONV) BCONV=BDIFU
      AN(J,I)=DTBRY*0.5*(RHOVN-ACONV)-SYM1/DYN
      BN(J,I)=DTBRY*0.5*(RHOVN-BCONV)-SYM2/DYN
C
      ACONV=DABS(RHOVS)
      BCONV=ACONV
      ADIFU=(EMU1(J,I)+EMU1(J-1,I))/DYS
      BDIFU=(EMU2(J,I)+EMU2(J-1,I))/DYS
      IF(ADIFU.GE.ACONV) ACONV=ADIFU
      IF(BDIFU.GE.BCONV) BCONV=BDIFU
      AS(J,I)=DTBRY*0.5*(-RHOVS-ACONV)+SYM1/DYS
      BS(J,I)=DTBRY*0.5*(-RHOVS-BCONV)+SYM2/DYS
C
      AP(J,I)=RHO(J,I)-(AE(J,I)+AW(J,I)+AN(J,I)+AS(J,I))
      BP(J,I)=RHO(J,I)-(BE(J,I)+BW(J,I)+BN(J,I)+BS(J,I))
      UZ=(U(J,I+1)-U(J,I))/DXCP
      VZ=0.25*(V(J+1,I+1)+V(J,I+1)-V(J+1,I-1)-V(J,I-1))/DXCP
      UR=0.25*(U(J+1,I+1)+U(J+1,I)-U(J-1,I+1)-U(J-1,I))/DYCP
      VR=(V(J+1,I)-V(J,I))/DYCP
      VP=0.5*(V(J,I)+V(J+1,I))
      WZ=0.5*((W(J,I+1)-W(J,I))/DXE+(W(J,I)-W(J,I-1))/DXW)
      WR=0.5*((W(J+1,I)-W(J,I))/DYN+(W(J,I)-W(J-1,I))/DYS)
      VORT(J,I)=DT*REI*TMU(J,I)*(2.0*(UZ*UZ+VR*VR
     1  +SYM*VP*VP/Y(J,I)/Y(J,I))+(VZ+UR)**2+SYM*(WZ*WZ+WR*WR))
      ENDDO
   20 CONTINUE
      SOR1=RELXKE
      SOR2=(1.0-RELXKE)
      SOR3=RELXKE
      SOR4=(1.0-RELXKE)
      AKLMT=AKREF/AKSTR*1.0D-15
      EPLMT=EPREF/EPSTR*1.0D-15
      DO 100 ISORP=1,ISOR
      RSDK=0.0
      RSDE=0.0
      DO 102 I=2,LI-1
      IF(ISKIP(2,I).NE.1) GO TO 102
      DELY=0.5*YYC(2)
      AKRT=DSQRT(AK(2,I))
      TERM2=DLOG(EE*DELY*RHO(2,I)*CMKEQ*AKRT/VMU(2,I)/REI)
      EPS(2,I)=CMKE*AK(2,I)*AKRT*TERM2/DELY/AVON
  102 CONTINUE
      DO 104 I=2,LI-1
      IF(ISKIP(LJ-1,I).NE.1) GO TO 104
      DELY=0.5*YYC(LJ-1)
      AKRT=DSQRT(AK(LJ-1,I))
      TERM2=DLOG(EE*DELY*RHO(LJ-1,I)*CMKEQ*AKRT/VMU(LJ-1,I)/REI)
      EPS(LJ-1,I)=CMKE*AK(LJ-1,I)*AKRT*TERM2/DELY/AVON
  104 CONTINUE
      DO 106 J=2,LJ-1
      IF(ISKIP(J,2).NE.1) GO TO 106
      DELX=0.5*XXC(2)
      AKRT=DSQRT(AK(J,2))
      TERM2=DLOG(EE*DELX*RHO(J,2)*CMKEQ*AKRT/VMU(J,2)/REI)
      EPS(J,2)=CMKE*AK(J,2)*AKRT*TERM2/DELX/AVON
  106 CONTINUE
      DO 108 J=2,LJ-1
      IF(ISKIP(J,LI-1).NE.1) GO TO 108
      DELX=0.5*XXC(LI-1)
      AKRT=DSQRT(AK(J,LI-1))
      TERM2=DLOG(EE*DELX*RHO(J,LI-1)*CMKEQ*AKRT/VMU(J,LI-1)/REI)
      EPS(J,LI-1)=CMKE*AK(J,LI-1)*AKRT*TERM2/DELX/AVON
  108 CONTINUE
      IF(NBODY.GE.1) THEN
          DO 109 N=1,NBODY
          IBODYM=IBDM(N)
          IBODYP=IBDP(N)
          JBODYM=JBDM(N)
          JBODYP=JBDP(N)
          IAM=IBODYM-1
          IF(IAM.LT.1) IAM=1
          IAP=IBODYP+1
          IF(IAP.GT.LI) IAP=LI
          JAM=JBODYM-1
          IF(JAM.LT.1) JAM=1
          JAP=JBODYP+1
          IF(JAP.GT.LJ) JAP=LJ
          DO I=IBODYM,IBODYP
          IF(ISKIP(JBODYP,I).GE.11) THEN
          DELY=0.5*YYC(JAP)
          AKRT=DSQRT(AK(JAP,I))
          TERM2=DLOG(EE*DELY*RHO(JAP,I)*CMKEQ*AKRT/VMU(JAP,I)/REI)
          EPS(JAP,I)=CMKE*AK(JAP,I)*AKRT*TERM2/DELY/AVON
          AK(JBODYP,I)=AKREF/AKSTR
          END IF
          IF(ISKIP(JBODYM,I).GE.11) THEN
          DELY=0.5*YYC(JAM)
          AKRT=DSQRT(AK(JAM,I))
          TERM2=DLOG(EE*DELY*RHO(JAM,I)*CMKEQ*AKRT/VMU(JAM,I)/REI)
          EPS(JAM,I)=CMKE*AK(JAM,I)*AKRT*TERM2/DELY/AVON
          AK(JBODYM,I)=AKREF/AKSTR
          END IF
          ENDDO
          DO J=JBODYM,JBODYP
          IF(ISKIP(J,IBODYP).GE.11) THEN
          DELX=0.5*XXC(IAP)
          AKRT=DSQRT(AK(J,IAP))
          TERM2=DLOG(EE*DELX*RHO(J,IAP)*CMKEQ*AKRT/VMU(J,IAP)/REI)
          EPS(J,IAP)=CMKE*AK(J,IAP)*AKRT*TERM2/DELX/AVON
          AK(J,IBODYP)=AKREF/AKSTR
          END IF
          IF(ISKIP(J,IBODYM).GE.11) THEN
          DELX=0.5*XXC(IAM)
          AKRT=DSQRT(AK(J,IAM))
          TERM2=DLOG(EE*DELX*RHO(J,IAM)*CMKEQ*AKRT/VMU(J,IAM)/REI)
          EPS(J,IAM)=CMKE*AK(J,IAM)*AKRT*TERM2/DELX/AVON
          AK(J,IBODYM)=AKREF/AKSTR
          END IF
          ENDDO
  109     CONTINUE
          END IF
C----------------------   POINT RELAXATION SCHEME  ---------------------
      DO 210 I=2,LI-1
C---------------------     SOLVE K - EQUATION     ----------------------
      DO 212 J=2,LJ-1
      IF(ISKIP(J,I).GE.2) GO TO 212
      RHS1=RHO(J,I)*FPA(J,I)+VORT(J,I)
      RHS2=-DT*RHO(J,I)*EPS(J,I)/AK(J,I)
      FNEW=SOR1*(RHS1-AN(J,I)*AK(J+1,I)
     1 -AS(J,I)*AK(J-1,I)-AE(J,I)*AK(J,I+1)
     2 -AW(J,I)*AK(J,I-1))/(AP(J,I)-RHS2)+SOR2*AK(J,I)
      IF(FNEW.LT.AKLMT) FNEW=AKLMT
      RSDK=RSDK+DABS(AK(J,I)-FNEW)
      AK(J,I)=FNEW
  212 CONTINUE
C-----------------    SOLVE EPS - EQUATION        ----------------------
      DO 214 J=2,LJ-1
      IF(ISKIP(J,I).GE.1) GO TO 214
      RHS1=RHO(J,I)*FPB(J,I)+C1KE*VORT(J,I)*EPS(J,I)/AK(J,I)
      RHS2=-DT*RHO(J,I)*C2KE*EPS(J,I)/AK(J,I)
      FNEW=SOR3*(RHS1-BN(J,I)*EPS(J+1,I)
     1 -BS(J,I)*EPS(J-1,I)-BE(J,I)*EPS(J,I+1)
     2 -BW(J,I)*EPS(J,I-1))/(BP(J,I)-RHS2)+SOR4*EPS(J,I)
      IF(FNEW.LT.EPLMT) FNEW=EPLMT
      RSDE=RSDE+DABS(EPS(J,I)-FNEW)
      EPS(J,I)=FNEW
  214 CONTINUE
  210 CONTINUE
      IF((RSDK/ENORM).LE.TOLRKE) GO TO 101
  100 CONTINUE
      ISORP=ISORP-1
  101 CONTINUE
      ISOR=ISORP
      RESDK=0.0
      RESDE=0.0
      DO 400 I=1,LI
      DO J=1,LJ
      RESDK=RESDK+DABS(FPA(J,I)-AK(J,I))
      RESDE=RESDE+DABS(FPB(J,I)-EPS(J,I))
      ENDDO
  400 CONTINUE
      RESDK=RESDK/ENORM
      RESDE=RESDE/ENORM
      RETURN
      END                                                               