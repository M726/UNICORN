      SUBROUTINE BODYBC(INSP,PBACK)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)        
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),
     1  W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB05/ RHONP(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDFSP(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/CB08/ TPOL1,TPOL2,POLSP(14,LSP),CISP(12,LSP)
      COMMON/CB09/ IBOT(LI),ITOP(LI),JLFT(LJ),JRGT(LJ),
     1       FBOT(8+LSP,LI),FTOP(8+LSP,LI),FLFT(8+LSP,LJ),FRGT(8+LSP,LJ)
      COMMON/CB10/ FBXP(4,LJ),FBXM(4,LJ),FBYP(4,LI),FBYM(4,LI)
      COMMON/BODY/IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10),NBODY
      COMMON/FINJ/IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30),NFINJ
      COMMON/HFORM/ HFO(LSP)
      COMMON/KEPS/ AKREF,EPREF,C1KE,C2KE,CMKE,CMKEQ,EE,AVON,PCON
      COMMON/NOISE/ INOISE,IVXTYP,FNOISE,ANOISE(200),VXTYP(5),TTIME
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
C                         WALL CONDITIONS (J=LJ)
      HCONV=GASC/WMSTR/HSTR
      ININJ=06
C11111111111111111111111111111111111111111111111111111111111111111111111
C---- J= 1 ---- J= 1 ---- J= 1 ---- J= 1 ---- J= 1 ---- J= 1 ---- J= 1 -
C11111111111111111111111111111111111111111111111111111111111111111111111
      JJ=1
      JA=2
      JB=3
      DO 100 I=1,LI
      AL1=FBOT(1,I)
      AL2=1.0-AL1
      IF(IBOT(I).EQ.0) THEN
C----------------------  IMPOSE INFLOW CONDITIONS  ---------------------
      RHO(JJ,I)=FBOT(2,I)
      U(JJ,I)=FBOT(3,I)
      V(JA,I)=FBOT(4,I)
      W(JA,I)=FBOT(5,I)
      HT(JJ,I)=FBOT(6,I)
      AK(JJ,I)=FBOT(7,I)
      EPS(JJ,I)=FBOT(8,I)
      STMF(JJ,I)=0.0
      STND(JJ,I)=0.0
      EPS(JJ,I)=FBOT(8,I)
      DO 102 ISP=1,LSP
      FSP(JJ,I,ISP)=FBOT(8+ISP,I)
  102 CONTINUE
      FBYM(1,I)=3.0*RHO(JA,I)-8.0*FBOT(2,I)+6.0*RHO(JJ,I)
      FBYM(2,I)=3.0*U(JA,I)-8.0*FBOT(3,I)+6.0*U(JJ,I)
      FBYM(3,I)=V(JA,I)
      FBYM(4,I)=3.0*W(JA,I)-8.0*FBOT(5,I)+6.0*W(JJ,I)
C-------------------------  END INFLOW   -------------------------------
      END IF
      IF(IBOT(I).EQ.1) THEN
C----------------------  WALL CONDITIONS  ------------------------------
      ITCOND=INT(FBOT(2,I))
      TWALL=FBOT(3,I)
      QWALL=FBOT(4,I)
      U(JJ,I)=-U(JA,I)
      V(JA,I)=0.0
      W(JJ,I)=-W(JA,I)
      HT(JJ,I)=HT(JA,I)
C     AK(JJ,I)=0.0
C     EPS(JJ,I)=1.0D-10
      STMF(JJ,I)=STMF(JA,I)
      STND(JJ,I)=STND(JA,I)
      IF(ICHEM.GE.1) THEN
         DO 104 ISP=1,LSP
         FSP(JJ,I,ISP)=FSP(JA,I,ISP)
         IF(FSP(JJ,I,ISP).LT.0.0) FSP(JJ,I,ISP)=0.0
  104    CONTINUE
         END IF
      IF(ITCOND.GE.1) THEN
              TKD=TWALL*TSTR
              TKD2=TKD*TKD
              TKD3=TKD2*TKD
              TKD4=TKD3*TKD
              TKD5=TKD4*TKD
              IPOLY=1
              IF(TKD.GT.TPOL1) IPOLY=8
              HWALL=0.0
              CPWALL=0.0
              DO 106 ISP=1,LSP
              HWALL=HWALL+FSP(JJ,I,ISP)*((POLSP(IPOLY,ISP)*TKD
     1        +POLSP(IPOLY+1,ISP)*TKD2/2.0+POLSP(IPOLY+2,ISP)*TKD3/3.0
     2        +POLSP(IPOLY+3,ISP)*TKD4/4.0+POLSP(IPOLY+4,ISP)*TKD5/5.0
     3        +POLSP(IPOLY+5,ISP))*HCONV/WM(ISP)-HFO(ISP))
              CPWALL=CPWALL+FSP(JJ,I,ISP)*(POLSP(IPOLY,ISP)
     1        +POLSP(IPOLY+1,ISP)*TKD+POLSP(IPOLY+2,ISP)*TKD2
     2        +POLSP(IPOLY+3,ISP)*TKD3+POLSP(IPOLY+4,ISP)*TKD4)/WM(ISP)
  106         CONTINUE
              HT(JJ,I)=HWALL
              CPWALL=CPWALL/WMSTR*GASC*TSTR/HSTR
              END IF
      IF(ITCOND.EQ.2) THEN
         PRNL=VMU(JJ,I)*CPWALL/VTC(JJ,I)/BETA3
         DELHT=YYS(JA)*(PRNL/VMU(JJ,I)/REI)*QWALL 
         HT(JJ,I)=HT(JA,I)-0.5*(U(JA,I)**2+V(JA,I)**2)+DELHT
         END IF
      FBYM(1,I)=RHO(JB,I)
      FBYM(2,I)=-U(JB,I)
      FBYM(3,I)=-V(JB,I)
      FBYM(4,I)=-W(JB,I)
C-------------------------  END WALL  ----------------------------------
      END IF
      IF(IBOT(I).EQ.2) THEN
C--------------------  IMPOSE SYMMETRY CONDITIONS  ---------------------
      Y2S=(0.5*YYS(JA))**2
      Y3S=(YYS(JB)+0.5*YYS(JA))**2
      RHO(JJ,I)=RHO(JA,I)
      U(JJ,I)=U(JA,I)
      V(JA,I)=0.0
      W(JJ,I)=0.0
      HT(JJ,I)=HT(JA,I)
      AK(JJ,I)=AK(JA,I)
      EPS(JJ,I)=EPS(JA,I)
      STMF(JJ,I)=STMF(JA,I)
      STND(JJ,I)=STND(JA,I)
      FBYM(1,I)=RHO(JA,I)-8.0*Y2S*(RHO(JA,I)-RHO(JB,I))/(Y3S-Y2S)
      FBYM(2,I)=U(JA,I)-8.0*Y2S*(U(JA,I)-U(JB,I))/(Y3S-Y2S)
      FBYM(3,I)=-V(JB,I)
      FBYM(4,I)=-W(JA,I)
      IF(ICHEM.GE.1) THEN
         DO 108 ISP=1,LSP-1
         FSP(JJ,I,ISP)=FSP(JA,I,ISP)
         IF(FSP(JJ,I,ISP).LT.0.0) FSP(JJ,I,ISP)=0.0
  108    CONTINUE
         END IF
C-----------------------  END SYMMETRY   -------------------------------
      END IF
      IF(IBOT(I).EQ.3) THEN
C--------------------  IMPOSE FREE FLOW CONDITIONS  --------------------
      RHO(JJ,I)=AL1*RHO(JA,I)+AL2*RHO(JB,I)
      U(JJ,I)=AL1*U(JA,I)+AL2*U(JB,I)
      V(JA,I)=V(JB,I)
      W(JJ,I)=AL1*W(JA,I)+AL2*W(JB,I)
      HT(JJ,I)=AL1*HT(JA,I)+AL2*HT(JB,I)
      AK(JJ,I)=AL1*AK(JA,I)+AL2*AK(JB,I)
      EPS(JJ,I)=AL1*EPS(JA,I)+AL2*EPS(JB,I)
      STMF(JJ,I)=AL1*STMF(JA,I)+AL2*STMF(JB,I)
      STND(JJ,I)=AL1*STND(JA,I)+AL2*STND(JB,I)
      FBYM(1,I)=AL1*RHO(JJ,I)+AL2*RHO(JA,I)
      FBYM(2,I)=AL1*U(JJ,I)+AL2*U(JA,I)
      FBYM(3,I)=V(JA,I)
      FBYM(4,I)=AL1*W(JJ,I)+AL2*W(JA,I)
      IF(ICHEM.GE.1) THEN
         DO 110 ISP=1,LSP-1
         FSP(JJ,I,ISP)=AL1*FSP(JA,I,ISP)+AL2*FSP(JB,I,ISP)
         IF(FSP(JJ,I,ISP).LT.0.0) FSP(JJ,I,ISP)=0.0
  110    CONTINUE
         END IF
C-----------------------  END FREE FLOW   ------------------------------
      END IF
C--------------------- ENFORCE COMPATIBILITY  --------------------------
      FTOT=0.0
      DO 111 ISP=1,LSP-1
      FTOT=FTOT+FSP(JJ,I,ISP)
  111 CONTINUE
      FTOT=FTOT+STMF(JJ,I)
      FSP(JJ,I,LSP)=1.0-FTOT
      IF(FTOT.GT.1.0) THEN
                      DO 112 ISP=1,LSP-1
                      FSP(JJ,I,ISP)=FSP(JJ,I,ISP)/FTOT
  112                 CONTINUE
                      STMF(JJ,I)=STMF(JJ,I)/FTOT
                      FSP(JJ,I,LSP)=0.0
                      END IF
      IF(ITHRM.LE.0) GO TO 101
C----------------------  CALCULATE TEMPERATURE -------------------------
      ENTH=HT(JJ,I)-0.5*(U(JJ,I)**2+V(JJ,I)**2+W(JJ,I)**2)
      DO 113 ISP=1,LSP
      ENTH=ENTH+FSP(JJ,I,ISP)*HFO(ISP)
  113 CONTINUE
      ENTH=ENTH*HSTR
      IPL=1
      TPOLY=TPOL1
      TKD=TK(JJ,I)*TSTR
  114 CONTINUE
      IF(TKD.GT.TPOLY) TKD=TPOLY/2.0
      A1=0.0
      A2=0.0
      A3=0.0
      A4=0.0
      A5=0.0
      A6=0.0
      DO 115 ISP=1,LSP
      A1=A1+POLSP(IPL,ISP)*FSP(JJ,I,ISP)/WM(ISP)
      A2=A2+POLSP(IPL+1,ISP)*FSP(JJ,I,ISP)/WM(ISP)/2.0
      A3=A3+POLSP(IPL+2,ISP)*FSP(JJ,I,ISP)/WM(ISP)/3.0
      A4=A4+POLSP(IPL+3,ISP)*FSP(JJ,I,ISP)/WM(ISP)/4.0
      A5=A5+POLSP(IPL+4,ISP)*FSP(JJ,I,ISP)/WM(ISP)/5.0
      A6=A6+POLSP(IPL+5,ISP)*FSP(JJ,I,ISP)/WM(ISP)
  115 CONTINUE
      FUNT=((A1+(A2+(A3+(A4+A5*TPOLY)*TPOLY)*TPOLY)*TPOLY)*TPOLY
     1    +A6)*GASC/WMSTR-ENTH
      IF(FUNT.GE.1.0D-07) GO TO 116
      IPL=8
      IF(TKD.LE.TPOLY) TKD=TPOL2/2.0
      TPOLY=TPOL2
      GO TO 114
  116 CONTINUE
      DO 120 INEWT=1,100
      DTKD=((A1+(A2+(A3+(A4+A5*TKD)*TKD)*TKD)*TKD)*TKD+A6-ENTH*WMSTR
     1   /GASC)/(A1+(2.0*A2+(3.0*A3+(4.0*A4+5.0*A5*TKD)*TKD)*TKD)*TKD)
      IF(DABS(DTKD).LE.0.001) GO TO 122
      TKD=TKD-DTKD
  120 CONTINUE
  122 CONTINUE
      IF(TKD.GT.TPOL2) TKD=TPOL2
      TK(JJ,I)=TKD/TSTR
      RBAR=0.0
      DO 123 ISP=1,LSP
      RBAR=RBAR+FSP(JJ,I,ISP)/WM(ISP)
  123 CONTINUE
      P(JJ,I)=P(JA,I)
      RHO(JJ,I)=(PBACK)/(BETA1*TK(JJ,I)*RBAR)
  101 CONTINUE
  100 CONTINUE
C22222222222222222222222222222222222222222222222222222222222222222222222
C---- J=LJ ---- J=LJ ---- J=LJ ---- J=LJ ---- J=LJ ---- J=LJ ---- J=LJ -
C22222222222222222222222222222222222222222222222222222222222222222222222
      JJ=LJ
      JA=LJ-1
      JB=LJ-2
      DO 200 I=1,LI
      AL1=FTOP(1,I)
      AL2=1.0-AL1
      IF(ITOP(I).EQ.0) THEN
C----------------------  IMPOSE INFLOW CONDITIONS  ---------------------
      RHO(JJ,I)=FTOP(2,I)
      U(JJ,I)=FTOP(3,I)
      V(JJ,I)=FTOP(4,I)
      W(JJ,I)=FTOP(5,I)
      HT(JJ,I)=FTOP(6,I)
      AK(JJ,I)=FTOP(7,I)
      EPS(JJ,I)=FTOP(8,I)
      STMF(JJ,I)=0.0
      STND(JJ,I)=0.0
      DO 202 ISP=1,LSP
      FSP(JJ,I,ISP)=FTOP(8+ISP,I)
  202 CONTINUE
      FBYP(1,I)=3.0*RHO(JA,I)-8.0*FTOP(2,I)+6.0*RHO(JJ,I)
      FBYP(2,I)=3.0*U(JA,I)-8.0*FTOP(3,I)+6.0*U(JJ,I)
      FBYP(3,I)=V(JJ,I)
      FBYP(4,I)=3.0*W(JA,I)-8.0*FTOP(5,I)+6.0*W(JJ,I)
C-------------------------  END INFLOW   -------------------------------
      END IF
      IF(ITOP(I).EQ.1) THEN
C----------------------  WALL CONDITIONS  ------------------------------
      ITCOND=INT(FTOP(2,I))
      TWALL=FTOP(3,I)
      QWALL=FTOP(4,I)
      U(JJ,I)=-U(JA,I)
      V(JJ,I)=0.0
      W(JJ,I)=-W(JA,I)
C     AK(JJ,I)=0.0
C     EPS(JJ,I)=1.0D-10
      STMF(JJ,I)=STMF(JA,I)
      STND(JJ,I)=STND(JA,I)
      HT(JJ,I)=HT(JA,I)
      IF(ICHEM.GE.1) THEN
         DO 204 ISP=1,LSP
         FSP(JJ,I,ISP)=FSP(JA,I,ISP)
         IF(FSP(JJ,I,ISP).LT.0.0) FSP(JJ,I,ISP)=0.0
  204    CONTINUE
         END IF
      IF(ITCOND.GE.1) THEN
              TKD=TWALL*TSTR
              TKD2=TKD*TKD
              TKD3=TKD2*TKD
              TKD4=TKD3*TKD
              TKD5=TKD4*TKD
              IPOLY=1
              IF(TKD.GT.TPOL1) IPOLY=8
              HWALL=0.0
              CPWALL=0.0
              DO 206 ISP=1,LSP
              HWALL=HWALL+FSP(JJ,I,ISP)*((POLSP(IPOLY,ISP)*TKD
     1        +POLSP(IPOLY+1,ISP)*TKD2/2.0+POLSP(IPOLY+2,ISP)*TKD3/3.0
     2        +POLSP(IPOLY+3,ISP)*TKD4/4.0+POLSP(IPOLY+4,ISP)*TKD5/5.0
     3        +POLSP(IPOLY+5,ISP))*HCONV/WM(ISP)-HFO(ISP))
              CPWALL=CPWALL+FSP(JJ,I,ISP)*(POLSP(IPOLY,ISP)
     1        +POLSP(IPOLY+1,ISP)*TKD+POLSP(IPOLY+2,ISP)*TKD2
     2        +POLSP(IPOLY+3,ISP)*TKD3+POLSP(IPOLY+4,ISP)*TKD4)/WM(ISP)
  206         CONTINUE
              HT(JJ,I)=HWALL
              CPWALL=CPWALL/WMSTR*GASC*TSTR/HSTR
              END IF
      IF(ITCOND.EQ.2) THEN
         PRNL=VMU(JJ,I)*CPWALL/VTC(JJ,I)/BETA3
         DELHT=YYS(JA)*(PRNL/VMU(JJ,I)/REI)*QWALL 
         HT(JJ,I)=HT(JA,I)-0.5*(U(JA,I)**2+V(JA,I)**2)+DELHT
         END IF
      FBYP(1,I)=RHO(JB,I)
      FBYP(2,I)=-U(JB,I)
      FBYP(3,I)=-V(JA,I)
      FBYP(4,I)=-W(JB,I)
C-------------------------  END WALL  ----------------------------------
      END IF
      IF(ITOP(I).EQ.2) THEN
C--------------------  IMPOSE SYMMETRY CONDITIONS  ---------------------
      Y2S=(0.5*YYS(JA))**2
      Y3S=(YYS(JB)+0.5*YYS(JA))**2
      RHO(JJ,I)=RHO(JA,I)
      U(JJ,I)=U(JA,I)
      V(JJ,I)=0.0
      W(JJ,I)=0.0
      HT(JJ,I)=HT(JA,I)
      AK(JJ,I)=AK(JA,I)
      EPS(JJ,I)=EPS(JA,I)
      STMF(JJ,I)=STMF(JA,I)
      STND(JJ,I)=STND(JA,I)
      FBYP(1,I)=RHO(JA,I)-8.0*Y2S*(RHO(JA,I)-RHO(JB,I))/(Y3S-Y2S)
      FBYP(2,I)=U(JA,I)-8.0*Y2S*(U(JA,I)-U(JB,I))/(Y3S-Y2S)
      FBYP(3,I)=-V(JA,I)
      FBYP(4,I)=-W(JA,I)
      IF(ICHEM.GE.1) THEN
         DO 208 ISP=1,LSP-1
         FSP(JJ,I,ISP)=FSP(JA,I,ISP)
         IF(FSP(JJ,I,ISP).LT.0.0) FSP(JJ,I,ISP)=0.0
  208    CONTINUE
         END IF
C-----------------------  END SYMMETRY   -------------------------------
      END IF
      IF(ITOP(I).EQ.3) THEN
C--------------------  IMPOSE FREE FLOW CONDITIONS  --------------------
      RHO(JJ,I)=AL1*RHO(JA,I)+AL2*RHO(JB,I)
      U(JJ,I)=U(JA,I)
      V(JJ,I)=AL1*V(JA,I)+AL2*V(JB,I)
      W(JJ,I)=AL1*W(JA,I)+AL2*W(JB,I)
      HT(JJ,I)=AL1*HT(JA,I)+AL2*HT(JB,I)
      AK(JJ,I)=AL1*AK(JA,I)+AL2*AK(JB,I)
      EPS(JJ,I)=AL1*EPS(JA,I)+AL2*EPS(JB,I)
      STMF(JJ,I)=AL1*STMF(JA,I)+AL2*STMF(JB,I)
      STND(JJ,I)=AL1*STND(JA,I)+AL2*STND(JB,I)
      FBYP(1,I)=AL1*RHO(JJ,I)+AL2*RHO(JA,I)
      FBYP(2,I)=AL1*U(JJ,I)+AL2*U(JA,I)
      FBYP(3,I)=V(JJ,I)
      FBYP(4,I)=AL1*W(JJ,I)+AL2*W(JA,I)
      IF(ICHEM.GE.1) THEN
         DO 210 ISP=1,LSP-1
         FSP(JJ,I,ISP)=AL1*FSP(JA,I,ISP)+AL2*FSP(JB,I,ISP)
         IF(FSP(JJ,I,ISP).LT.0.0) FSP(JJ,I,ISP)=0.0
  210    CONTINUE
         END IF
C-----------------------  END FREE FLOW   ------------------------------
      END IF
C--------------------- ENFORCE COMPATIBILITY  --------------------------
      FTOT=0.0
      DO 211 ISP=1,LSP-1
      FTOT=FTOT+FSP(JJ,I,ISP)
  211 CONTINUE
      FTOT=FTOT+STMF(JJ,I)
      FSP(JJ,I,LSP)=1.0-FTOT
      IF(FTOT.GT.1.0) THEN
                      DO 212 ISP=1,LSP-1
                      FSP(JJ,I,ISP)=FSP(JJ,I,ISP)/FTOT
  212                 CONTINUE
                      STMF(JJ,I)=STMF(JJ,I)/FTOT
                      FSP(JJ,I,LSP)=0.0
                      END IF
      IF(ITHRM.LE.0) GO TO 201
C----------------------  CALCULATE TEMPERATURE -------------------------
      ENTH=HT(JJ,I)
     1    -0.5*(U(JJ,I)*U(JJ,I)+V(JJ,I)*V(JJ,I)+W(JJ,I)*W(JJ,I))
      DO 213 ISP=1,LSP
      ENTH=ENTH+FSP(JJ,I,ISP)*HFO(ISP)
  213 CONTINUE
      ENTH=ENTH*HSTR
      IPL=1
      TPOLY=TPOL1
      TKD=TK(JJ,I)*TSTR
  214 CONTINUE
      IF(TKD.GT.TPOLY) TKD=TPOLY/2.0
      A1=0.0
      A2=0.0
      A3=0.0
      A4=0.0
      A5=0.0
      A6=0.0
      DO 215 ISP=1,LSP
      A1=A1+POLSP(IPL,ISP)*FSP(JJ,I,ISP)/WM(ISP)
      A2=A2+POLSP(IPL+1,ISP)*FSP(JJ,I,ISP)/WM(ISP)/2.0
      A3=A3+POLSP(IPL+2,ISP)*FSP(JJ,I,ISP)/WM(ISP)/3.0
      A4=A4+POLSP(IPL+3,ISP)*FSP(JJ,I,ISP)/WM(ISP)/4.0
      A5=A5+POLSP(IPL+4,ISP)*FSP(JJ,I,ISP)/WM(ISP)/5.0
      A6=A6+POLSP(IPL+5,ISP)*FSP(JJ,I,ISP)/WM(ISP)
  215 CONTINUE
      FUNT=((A1+(A2+(A3+(A4+A5*TPOLY)*TPOLY)*TPOLY)*TPOLY)*TPOLY
     1    +A6)*GASC/WMSTR-ENTH
      IF(FUNT.GE.1.0D-07) GO TO 216
      IPL=8
      IF(TKD.LE.TPOLY) TKD=TPOL2/2.0
      TPOLY=TPOL2
      GO TO 214
  216 CONTINUE
      DO 220 INEWT=1,100
      DTKD=((A1+(A2+(A3+(A4+A5*TKD)*TKD)*TKD)*TKD)*TKD+A6-ENTH*WMSTR
     1   /GASC)/(A1+(2.0*A2+(3.0*A3+(4.0*A4+5.0*A5*TKD)*TKD)*TKD)*TKD)
      IF(DABS(DTKD).LE.0.001) GO TO 222
      TKD=TKD-DTKD
  220 CONTINUE
  222 CONTINUE
      IF(TKD.GT.TPOL2) TKD=TPOL2
      TK(JJ,I)=TKD/TSTR
      RBAR=0.0
      DO 223 ISP=1,LSP
      RBAR=RBAR+FSP(JJ,I,ISP)/WM(ISP)
  223 CONTINUE
      P(JJ,I)=P(JA,I)
      RHO(JJ,I)=(PBACK)/(BETA1*TK(JJ,I)*RBAR)
  201 CONTINUE
  200 CONTINUE
C33333333333333333333333333333333333333333333333333333333333333333333333
C---- I= 1 ---- I= 1 ---- I= 1 ---- I= 1 ---- I= 1 ---- I= 1 ---- I= 1 -
C33333333333333333333333333333333333333333333333333333333333333333333333
C--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE-
      UPERT=0.0
      IF(INOISE.LT.0) THEN
         IWIDTH=-INOISE-1
         IF(INOISE.LT.-100) IWIDTH=IWIDTH-100
         PI2=2.0*3.141592654
         FT=TTIME*DABS(FNOISE)
         ICYCLE=INT(FT)
         IF(FNOISE.LE.0.0) ICYCLE=0
         FTIME=(FT-DFLOAT(ICYCLE))
         IA=INT(FTIME*IWIDTH)+1
         IF(IA.GT.IWIDTH) IA=IWIDTH
         IB=IA+1
         UPERT=(ANOISE(IA)*DFLOAT(IB-1)-ANOISE(IB)*DFLOAT(IA-1)
     1        +FTIME*(ANOISE(IB)-ANOISE(IA))*DFLOAT(IWIDTH))/VSTR
         IF(FLFT(3,1).GT.0.0) UPERT=UPERT/FLFT(3,1)
         END IF
C--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE-
      II=1
      IA=2
      IB=3
      DO 300 J=1,LJ
      AL1=FLFT(1,J)
      AL2=1.0-AL1
      IF(JLFT(J).EQ.0) THEN
C----------------------  IMPOSE INFLOW CONDITIONS  ---------------------
      RHO(J,II)=FLFT(2,J)
      UU=0.0
C     IF(FLFT(7,J).EQ.FLFT(7,1)) UU=1.0
      IF(J.LE.21) UU=1.0
      U(J,IA)=FLFT(3,J)*(1.0+UU*UPERT)
      V(J,II)=FLFT(4,J)
      W(J,II)=FLFT(5,J)
      HT(J,II)=FLFT(6,J)
      AK(J,II)=FLFT(7,J)
      EPS(J,II)=FLFT(8,J)
      STMF(J,II)=0.0
      STND(J,II)=0.0
      FBXM(1,J)=3.0*RHO(J,IA)-8.0*FLFT(2,J)+6.0*RHO(J,II)
      FBXM(2,J)=U(J,IA)
      FBXM(3,J)=3.0*V(J,IA)-8.0*FLFT(4,J)+6.0*V(J,II)
      FBXM(4,J)=3.0*W(J,IA)-8.0*FLFT(5,J)+6.0*W(J,II)
      DO 302 ISP=1,LSP-1
      FSP(J,II,ISP)=0.0
  302 CONTINUE
      IF(UU.GE.1.0.AND.IVXTYP.EQ.1) THEN
          FSP(J,II,INSP)=VXTYP(1)
          FSP(J,II,2)=VXTYP(2)
          FSP(J,II,ININJ)=VXTYP(3)
          RHO(J,II)=VXTYP(4)
          HT(J,II)=VXTYP(5)
          FBXM(1,J)=RHO(J,II)
          ELSE
          DO 303 ISP=1,LSP
          FSP(J,II,ISP)=FLFT(8+ISP,J)
  303     CONTINUE
          FBXM(1,J)=RHO(J,II)
          END IF
C-------------------------  END INFLOW   -------------------------------
      END IF
      IF(JLFT(J).EQ.1) THEN
C----------------------  WALL CONDITIONS  ------------------------------
      ITCOND=INT(FLFT(2,J))
      TWALL=FLFT(3,J)
      QWALL=FLFT(4,J)
      U(J,IA)=0.0
      V(J,II)=-V(J,IA)
      W(J,II)=-W(J,IA)
      HT(J,II)=HT(J,IA)
C     AK(J,II)=0.0
C     EPS(J,II)=1.0D-10
      STMF(J,II)=0.0
      STND(J,II)=0.0
      IF(ICHEM.GE.1) THEN
         DO 304 ISP=1,LSP
         FSP(J,II,ISP)=FSP(J,IA,ISP)
         IF(FSP(J,II,ISP).LT.0.0) FSP(J,II,ISP)=0.0
  304    CONTINUE
         END IF
      IF(ITCOND.GE.1) THEN
              TKD=TWALL*TSTR
              TKD2=TKD*TKD
              TKD3=TKD2*TKD
              TKD4=TKD3*TKD
              TKD5=TKD4*TKD
              IPOLY=1
              IF(TKD.GT.TPOL1) IPOLY=8
              HWALL=0.0
              CPWALL=0.0
              DO 306 ISP=1,LSP
              HWALL=HWALL+FSP(J,II,ISP)*((POLSP(IPOLY,ISP)*TKD
     1        +POLSP(IPOLY+1,ISP)*TKD2/2.0+POLSP(IPOLY+2,ISP)*TKD3/3.0
     2        +POLSP(IPOLY+3,ISP)*TKD4/4.0+POLSP(IPOLY+4,ISP)*TKD5/5.0
     3        +POLSP(IPOLY+5,ISP))*HCONV/WM(ISP)-HFO(ISP))
              CPWALL=CPWALL+FSP(J,II,ISP)*(POLSP(IPOLY,ISP)
     1        +POLSP(IPOLY+1,ISP)*TKD+POLSP(IPOLY+2,ISP)*TKD2
     2        +POLSP(IPOLY+3,ISP)*TKD3+POLSP(IPOLY+4,ISP)*TKD4)/WM(ISP)
  306         CONTINUE
              HT(J,II)=HWALL
              CPWALL=CPWALL/WMSTR*GASC*TSTR/HSTR
              END IF
      IF(ITCOND.EQ.2) THEN
         PRNL=VMU(J,II)*CPWALL/VTC(J,II)/BETA3
         DELHT=XXS(IA)*(PRNL/VMU(J,II)/REI)*QWALL 
         HT(J,II)=HT(J,IA)-0.5*(U(J,IA)**2+V(J,IA)**2)+DELHT
         END IF
      FBXM(1,J)=RHO(J,IB)
      FBXM(2,J)=-U(J,IB)
      FBXM(3,J)=-V(J,IB)
      FBXM(4,J)=-W(J,IB)
C-------------------------  END WALL  ----------------------------------
      END IF
      IF(JLFT(J).EQ.2) THEN
C--------------------  IMPOSE SYMMETRY CONDITIONS  ---------------------
      X2S=(0.5*XXS(IA))**2
      X3S=(XXS(IB)+0.5*XXS(IA))**2
      RHO(J,II)=RHO(J,IA)
      U(J,IA)=0.0
      V(J,II)=V(J,IA)
      W(J,II)=W(J,IA)
      HT(J,II)=HT(J,IA)
      AK(J,II)=AK(J,IA)
      EPS(J,II)=EPS(J,IA)
      STMF(J,II)=STMF(J,IA)
      STND(J,II)=STND(J,IA)
      FBXM(1,J)=RHO(J,IA)-8.0*X2S*(RHO(J,IA)-RHO(J,IB))/(X3S-X2S)
      FBXM(2,J)=-U(J,IB)
      FBXM(3,J)=V(J,IA)-8.0*X2S*(V(J,IA)-V(J,IB))/(X3S-X2S)
      FBXM(4,J)=W(J,IA)-8.0*X2S*(W(J,IA)-W(J,IB))/(X3S-X2S)
      IF(ICHEM.GE.1) THEN
         DO 308 ISP=1,LSP-1
         FSP(J,II,ISP)=FSP(J,IA,ISP)
         IF(FSP(J,II,ISP).LT.0.0) FSP(J,II,ISP)=0.0
  308    CONTINUE
         END IF
C-----------------------  END SYMMETRY   -------------------------------
      END IF
      IF(JLFT(J).EQ.3) THEN
C--------------------  IMPOSE FREE FLOW CONDITIONS  --------------------
      RHO(J,II)=AL1*RHO(J,IA)+AL2*RHO(J,IB)
      U(J,IA)=U(J,IB)
      V(J,II)=AL1*V(J,IA)+AL2*V(J,IB)
      W(J,II)=AL1*W(J,IA)+AL2*W(J,IB)
      HT(J,II)=AL1*HT(J,IA)+AL2*HT(J,IB)
      AK(J,II)=AL1*AK(J,IA)+AL2*AK(J,IB)
      EPS(J,II)=AL1*EPS(J,IA)+AL2*EPS(J,IB)
      STMF(J,II)=AL1*STMF(J,IA)+AL2*STMF(J,IB)
      STND(J,II)=AL1*STND(J,IA)+AL2*STND(J,IB)
      FBXM(1,J)=AL1*RHO(J,II)+AL2*RHO(J,IA)
      FBXM(2,J)=AL1*U(J,II)+AL2*U(J,IA)
      FBXM(3,J)=V(J,IA)
      FBXM(4,J)=AL1*W(J,II)+AL2*W(J,IA)
      IF(ICHEM.GE.1) THEN
         DO 310 ISP=1,LSP-1
         FSP(J,II,ISP)=AL1*FSP(J,IA,ISP)+AL2*FSP(J,IB,ISP)
         IF(FSP(J,II,ISP).LT.0.0) FSP(J,II,ISP)=0.0
  310    CONTINUE
         END IF
C-----------------------  END FREE FLOW   ------------------------------
      END IF
C--------------------- ENFORCE COMPATIBILITY  --------------------------
      FTOT=0.0
      DO 311 ISP=1,LSP-1
      FTOT=FTOT+FSP(J,II,ISP)
  311 CONTINUE
      FTOT=FTOT+STMF(J,II)
      FSP(J,II,LSP)=1.0-FTOT
      IF(FTOT.GT.1.0) THEN
                      DO 312 ISP=1,LSP-1
                      FSP(J,II,ISP)=FSP(J,II,ISP)/FTOT
  312                 CONTINUE
                      STMF(J,II)=STMF(J,II)/FTOT
                      FSP(J,II,LSP)=0.0
                      END IF
      IF(ITHRM.LE.0) GO TO 301
C----------------------  CALCULATE TEMPERATURE -------------------------
      ENTH=HT(J,II)
     1    -0.5*(U(J,II)*U(J,II)+V(J,II)*V(J,II)+W(J,II)*W(J,II))
      DO 313 ISP=1,LSP
      ENTH=ENTH+FSP(J,II,ISP)*HFO(ISP)
  313 CONTINUE
      ENTH=ENTH*HSTR
      IPL=1
      TPOLY=TPOL1
      TKD=TK(J,II)*TSTR
  314 CONTINUE
      IF(TKD.GT.TPOLY) TKD=TPOLY/2.0
      A1=0.0
      A2=0.0
      A3=0.0
      A4=0.0
      A5=0.0
      A6=0.0
      DO 315 ISP=1,LSP
      A1=A1+POLSP(IPL,ISP)*FSP(J,II,ISP)/WM(ISP)
      A2=A2+POLSP(IPL+1,ISP)*FSP(J,II,ISP)/WM(ISP)/2.0
      A3=A3+POLSP(IPL+2,ISP)*FSP(J,II,ISP)/WM(ISP)/3.0
      A4=A4+POLSP(IPL+3,ISP)*FSP(J,II,ISP)/WM(ISP)/4.0
      A5=A5+POLSP(IPL+4,ISP)*FSP(J,II,ISP)/WM(ISP)/5.0
      A6=A6+POLSP(IPL+5,ISP)*FSP(J,II,ISP)/WM(ISP)
  315 CONTINUE
      FUNT=((A1+(A2+(A3+(A4+A5*TPOLY)*TPOLY)*TPOLY)*TPOLY)*TPOLY
     1    +A6)*GASC/WMSTR-ENTH
      IF(FUNT.GE.1.0D-07) GO TO 316
      IPL=8
      IF(TKD.LE.TPOLY) TKD=TPOL2/2.0
      TPOLY=TPOL2
      GO TO 314
  316 CONTINUE
      DO 320 INEWT=1,100
      DTKD=((A1+(A2+(A3+(A4+A5*TKD)*TKD)*TKD)*TKD)*TKD+A6-ENTH*WMSTR
     1   /GASC)/(A1+(2.0*A2+(3.0*A3+(4.0*A4+5.0*A5*TKD)*TKD)*TKD)*TKD)
      IF(DABS(DTKD).LE.0.001) GO TO 322
      TKD=TKD-DTKD
  320 CONTINUE
  322 CONTINUE
      IF(TKD.GT.TPOL2) TKD=TPOL2
      TK(J,II)=TKD/TSTR
      RBAR=0.0
      DO 323 ISP=1,LSP
      RBAR=RBAR+FSP(J,II,ISP)/WM(ISP)
  323 CONTINUE
      P(J,II)=P(J,IA)
      RHO(J,II)=(PBACK)/(BETA1*TK(J,II)*RBAR)
  301 CONTINUE
  300 CONTINUE
C44444444444444444444444444444444444444444444444444444444444444444444444
C---- I=LI ---- I=LI ---- I=LI ---- I=LI ---- I=LI ---- I=LI ---- I=LI -
C44444444444444444444444444444444444444444444444444444444444444444444444
C--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE-
      UPERT=0.0
      IF(INOISE.LT.-100) THEN
         IWIDTH=-INOISE-100-1
         PI2=2.0*3.141592654
         FT=TTIME*DABS(FNOISE)
         ICYCLE=INT(FT)
         IF(FNOISE.LE.0.0) ICYCLE=0
         FTIME=(FT-DFLOAT(ICYCLE))
         IA=INT(FTIME*IWIDTH)+1
         IF(IA.GT.IWIDTH) IA=IWIDTH
         IB=IA+1
         UPERT=(ANOISE(100+IA)*DFLOAT(IB-1)-ANOISE(100+IB)*DFLOAT(IA-1)
     1      +FTIME*(ANOISE(100+IB)-ANOISE(100+IA))*DFLOAT(IWIDTH))/VSTR
         IF(FRGT(3,1).GT.0.0) UPERT=UPERT/FRGT(3,1)
         END IF
C--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE--DRIVE-
      II=LI
      IA=LI-1
      IB=LI-2
      DO 400 J=1,LJ
      AL1=FRGT(1,J)
      AL2=1.0-AL1
      IF(JRGT(J).EQ.0) THEN
C----------------------  IMPOSE INFLOW CONDITIONS  ---------------------
      RHO(J,II)=FRGT(2,J)
      UU=0.0
C     IF(FLFT(7,J).EQ.FLFT(7,1)) UU=1.0
      IF(J.LE.21) UU=1.0
      U(J,II)=FRGT(3,J)*(1.0+UU*UPERT)
      V(J,II)=FRGT(4,J)
      W(J,II)=FRGT(5,J)
      HT(J,II)=FRGT(6,J)
      AK(J,II)=FRGT(7,J)
      EPS(J,II)=FRGT(8,J)
      STMF(J,II)=0.0
      STND(J,II)=0.0
      FBXP(1,J)=3.0*RHO(J,IA)-8.0*FRGT(2,J)+6.0*RHO(J,II)
      FBXP(2,J)=U(J,II)
      FBXP(3,J)=3.0*V(J,IA)-8.0*FRGT(4,J)+6.0*V(J,II)
      FBXP(4,J)=3.0*W(J,IA)-8.0*FRGT(5,J)+6.0*W(J,II)
      DO 402 ISP=1,LSP-1
      FSP(J,II,ISP)=0.0
  402 CONTINUE
      IF(UU.GE.1.0.AND.IVXTYP.EQ.1) THEN
          FSP(J,II,INSP)=VXTYP(1)
          FSP(J,II,2)=VXTYP(2)
          FSP(J,II,ININJ)=VXTYP(3)
          RHO(J,II)=VXTYP(4)
          HT(J,II)=VXTYP(5)
          FBXM(1,J)=RHO(J,II)
          ELSE
          DO 403 ISP=1,LSP
          FSP(J,II,ISP)=FRGT(8+ISP,J)
  403     CONTINUE
          END IF
C-------------------------  END INFLOW   -------------------------------
      END IF
      IF(JRGT(J).EQ.1) THEN
C----------------------  WALL CONDITIONS  ------------------------------
      ITCOND=INT(FRGT(2,J))
      TWALL=FRGT(3,J)
      QWALL=FRGT(4,J)
      U(J,II)=0.0
      V(J,II)=-V(J,IA)
      W(J,II)=-W(J,IA)
      HT(J,II)=HT(J,IA)
C     AK(J,II)=0.0
C     EPS(J,II)=1.0D-10
      STMF(J,II)=0.0
      STND(J,II)=0.0
      IF(ICHEM.GE.1) THEN
         DO 404 ISP=1,LSP
         FSP(J,II,ISP)=FSP(J,IA,ISP)
         IF(FSP(J,II,ISP).LT.0.0) FSP(J,II,ISP)=0.0
  404    CONTINUE
         END IF
      IF(ITCOND.GE.1) THEN
              TKD=TWALL*TSTR
              TKD2=TKD*TKD
              TKD3=TKD2*TKD
              TKD4=TKD3*TKD
              TKD5=TKD4*TKD
              IPOLY=1
              IF(TKD.GT.TPOL1) IPOLY=8
              HWALL=0.0
              CPWALL=0.0
              DO 406 ISP=1,LSP
              HWALL=HWALL+FSP(J,II,ISP)*((POLSP(IPOLY,ISP)*TKD
     1        +POLSP(IPOLY+1,ISP)*TKD2/2.0+POLSP(IPOLY+2,ISP)*TKD3/3.0
     2        +POLSP(IPOLY+3,ISP)*TKD4/4.0+POLSP(IPOLY+4,ISP)*TKD5/5.0
     3        +POLSP(IPOLY+5,ISP))*HCONV/WM(ISP)-HFO(ISP))
              CPWALL=CPWALL+FSP(J,II,ISP)*(POLSP(IPOLY,ISP)
     1        +POLSP(IPOLY+1,ISP)*TKD+POLSP(IPOLY+2,ISP)*TKD2
     2        +POLSP(IPOLY+3,ISP)*TKD3+POLSP(IPOLY+4,ISP)*TKD4)/WM(ISP)
  406         CONTINUE
              HT(J,II)=HWALL
              CPWALL=CPWALL/WMSTR*GASC*TSTR/HSTR
              END IF
      IF(ITCOND.EQ.2) THEN
         PRNL=VMU(J,II)*CPWALL/VTC(J,II)/BETA3
         DELHT=XXS(IA)*(PRNL/VMU(J,II)/REI)*QWALL 
         HT(J,II)=HT(J,IA)-0.5*(U(J,IA)**2+V(J,IA)**2)+DELHT
         END IF
      FBXP(1,J)=RHO(J,IB)
      FBXP(2,J)=-U(J,IA)
      FBXP(3,J)=-V(J,IB)
      FBXP(4,J)=-W(J,IB)
C-------------------------  END WALL  ----------------------------------
      END IF
      IF(JRGT(J).EQ.2) THEN
C--------------------  IMPOSE SYMMETRY CONDITIONS  ---------------------
      X2S=(0.5*XXS(IA))**2
      X3S=(XXS(IB)+0.5*XXS(IA))**2
      RHO(J,II)=RHO(J,IA)
      U(J,II)=0.0
      V(J,II)=V(J,IA)
      W(J,II)=W(J,IA)
      HT(J,II)=HT(J,IA)
      AK(J,II)=AK(J,IA)
      EPS(J,II)=EPS(J,IA)
      STMF(J,II)=STMF(J,IA)
      STND(J,II)=STND(J,IA)
      FBXP(1,J)=RHO(J,IA)-8.0*X2S*(RHO(J,IA)-RHO(J,IB))/(X3S-X2S)
      FBXP(2,J)=-U(J,IA)
      FBXP(3,J)=V(J,IA)-8.0*X2S*(V(J,IA)-V(J,IB))/(X3S-X2S)
      FBXP(4,J)=W(J,IA)-8.0*X2S*(W(J,IA)-W(J,IB))/(X3S-X2S)
      IF(ICHEM.GE.1) THEN
         DO 408 ISP=1,LSP-1
         FSP(J,II,ISP)=FSP(J,IA,ISP)
         IF(FSP(J,II,ISP).LT.0.0) FSP(J,II,ISP)=0.0
  408    CONTINUE
         END IF
C-----------------------  END SYMMETRY   -------------------------------
      END IF
      IF(JRGT(J).EQ.3) THEN
C--------------------  IMPOSE FREE FLOW CONDITIONS  --------------------
      RHO(J,II)=AL1*RHO(J,IA)+AL2*RHO(J,IB)
      U(J,II)=AL1*U(J,IA)+AL2*U(J,IB)
      V(J,II)=AL1*V(J,IA)+AL2*V(J,IB)
      W(J,II)=AL1*W(J,IA)+AL2*W(J,IB)
      HT(J,II)=AL1*HT(J,IA)+AL2*HT(J,IB)
      AK(J,II)=AL1*AK(J,IA)+AL2*AK(J,IB)
      EPS(J,II)=AL1*EPS(J,IA)+AL2*EPS(J,IB)
      STMF(J,II)=AL1*STMF(J,IA)+AL2*STMF(J,IB)
      STND(J,II)=AL1*STND(J,IA)+AL2*STND(J,IB)
      FBXP(1,J)=RHO(J,II)
      FBXP(2,J)=U(J,II)
      FBXP(3,J)=V(J,II)
      FBXP(4,J)=W(J,II)
      IF(ICHEM.GE.1) THEN
         DO 410 ISP=1,LSP-1
         FSP(J,II,ISP)=AL1*FSP(J,IA,ISP)+AL2*FSP(J,IB,ISP)
         IF(FSP(J,II,ISP).LT.0.0) FSP(J,II,ISP)=0.0
  410    CONTINUE
         END IF
C-----------------------  END FREE FLOW   ------------------------------
      END IF
C--------------------- ENFORCE COMPATIBILITY  --------------------------
      FTOT=0.0
      DO 411 ISP=1,LSP-1
      FTOT=FTOT+FSP(J,II,ISP)
  411 CONTINUE
      FTOT=FTOT+STMF(J,II)
      FSP(J,II,LSP)=1.0-FTOT
      IF(FTOT.GT.1.0) THEN
                      DO 412 ISP=1,LSP-1
                      FSP(J,II,ISP)=FSP(J,II,ISP)/FTOT
  412                 CONTINUE
                      STMF(J,II)=STMF(J,II)/FTOT
                      FSP(J,II,LSP)=0.0
                      END IF
      IF(ITHRM.LE.0) GO TO 401
C----------------------  CALCULATE TEMPERATURE -------------------------
      ENTH=HT(J,II)
     1    -0.5*(U(J,II)*U(J,II)+V(J,II)*V(J,II)+W(J,II)*W(J,II))
      DO 413 ISP=1,LSP
      ENTH=ENTH+FSP(J,II,ISP)*HFO(ISP)
  413 CONTINUE
      ENTH=ENTH*HSTR
      IPL=1
      TPOLY=TPOL1
      TKD=TK(J,II)*TSTR
  414 CONTINUE
      IF(TKD.GT.TPOLY) TKD=TPOLY/2.0
      A1=0.0
      A2=0.0
      A3=0.0
      A4=0.0
      A5=0.0
      A6=0.0
      DO 415 ISP=1,LSP
      A1=A1+POLSP(IPL,ISP)*FSP(J,II,ISP)/WM(ISP)
      A2=A2+POLSP(IPL+1,ISP)*FSP(J,II,ISP)/WM(ISP)/2.0
      A3=A3+POLSP(IPL+2,ISP)*FSP(J,II,ISP)/WM(ISP)/3.0
      A4=A4+POLSP(IPL+3,ISP)*FSP(J,II,ISP)/WM(ISP)/4.0
      A5=A5+POLSP(IPL+4,ISP)*FSP(J,II,ISP)/WM(ISP)/5.0
      A6=A6+POLSP(IPL+5,ISP)*FSP(J,II,ISP)/WM(ISP)
  415 CONTINUE
      FUNT=((A1+(A2+(A3+(A4+A5*TPOLY)*TPOLY)*TPOLY)*TPOLY)*TPOLY
     1    +A6)*GASC/WMSTR-ENTH
      IF(FUNT.GE.1.0D-07) GO TO 416
      IPL=8
      IF(TKD.LE.TPOLY) TKD=TPOL2/2.0
      TPOLY=TPOL2
      GO TO 414
  416 CONTINUE
      DO 420 INEWT=1,100
      DTKD=((A1+(A2+(A3+(A4+A5*TKD)*TKD)*TKD)*TKD)*TKD+A6-ENTH*WMSTR
     1   /GASC)/(A1+(2.0*A2+(3.0*A3+(4.0*A4+5.0*A5*TKD)*TKD)*TKD)*TKD)
      IF(DABS(DTKD).LE.0.001) GO TO 422
      TKD=TKD-DTKD
  420 CONTINUE
  422 CONTINUE
      IF(TKD.GT.TPOL2) TKD=TPOL2
      TK(J,II)=TKD/TSTR
      RBAR=0.0
      DO 423 ISP=1,LSP
      RBAR=RBAR+FSP(J,II,ISP)/WM(ISP)
  423 CONTINUE
      P(J,II)=P(J,IA)
      RHO(J,II)=(PBACK)/(BETA1*TK(J,II)*RBAR)
  401 CONTINUE
  400 CONTINUE
C--------------------       INSERT INNER BODY       --------------------
      IF(NBODY.GE.1) THEN
            DO 500 N=1,NBODY
            IBODYM=IBDM(N)
            IBODYP=IBDP(N)
            JBODYM=JBDM(N)
            JBODYP=JBDP(N)
            DO 502 I=IBODYM+1,IBODYP
            DO J=JBODYM,JBODYP
            U(J,I)=0.0
            ENDDO
  502       CONTINUE
            DO 504 I=IBODYM,IBODYP
            DO J=JBODYM+1,JBODYP
            V(J,I)=0.0
            ENDDO
  504       CONTINUE
C
C           DO 506 I=IBODYM+1,IBODYP
C           IF(ISKIP(JBODYP,I-1).GE.10) U(JBODYP-1,I)=-U(JBODYP+1,I)
C           IF(ISKIP(JBODYM,I-1).GE.10) U(JBODYM+1,I)=-U(JBODYM-1,I)
C 506       CONTINUE
C
C           DO 508 J=JBODYM+1,JBODYP
C           IF(ISKIP(J-1,IBODYP).GE.10) V(J,IBODYP-1)=-V(J,IBODYP+1)
C           IF(ISKIP(J-1,IBODYM).GE.10) V(J,IBODYM+1)=-V(J,IBODYM-1)
C 508       CONTINUE
C
            DO 510 I=IBODYM,IBODYP
            IF(ISKIP(JBODYP,I).GE.10) THEN
	       JJ=JBODYP+1
	       IF(JJ.GT.LJ) JJ=LJ
               W(JBODYP,I)=0.0
               STMF(JBODYP,I)=STMF(JJ,I)
               STND(JBODYP,I)=STND(JJ,I)
               IF(TBD(N).EQ.0.0) TK(JBODYP,I)=TK(JJ,I)
               IF(TBD(N).EQ.1.0) TK(JBODYP,I)=TK(JBODYP,I)
               IF(TBD(N).GT.1.0) TK(JBODYP,I)=TBD(N)/TSTR
               TKD=TK(JBODYP,I)*TSTR
               TKD2=TKD*TKD
               TKD3=TKD*TKD2
               TKD4=TKD*TKD3
               TKD5=TKD*TKD4
               IPOLY=1
               IF(TKD.GT.TPOL1) IPOLY=8
               RBAR=0.0
               ENTH=0.0
               DO 511 ISP=1,LSP
               FSP(JBODYP,I,ISP)=FSP(JJ,I,ISP)
               RBAR=RBAR+FSP(JBODYP,I,ISP)/WM(ISP)
               ENTH=ENTH+FSP(JBODYP,I,ISP)*((POLSP(IPOLY,ISP)*TKD
     1         +POLSP(IPOLY+1,ISP)*TKD2/2.0+POLSP(IPOLY+2,ISP)*TKD3/3.0
     2         +POLSP(IPOLY+3,ISP)*TKD4/4.0+POLSP(IPOLY+4,ISP)*TKD5/5.0
     3         +POLSP(IPOLY+5,ISP))*HCONV/WM(ISP)-HFO(ISP))
  511          CONTINUE
               RHO(JBODYP,I)=(PBACK)/(BETA1*TK(JBODYP,I)*RBAR)
               RHONP(JBODYP,I)=RHO(JBODYP,I)
               HT(JBODYP,I)=ENTH
               END IF
            IF(ISKIP(JBODYM,I).GE.10) THEN
	       JJ=JBODYM-1
	       IF(JJ.LT.1) JJ=1
               W(JBODYM,I)=0.0
               STMF(JBODYM,I)=STMF(JJ,I)
               STND(JBODYM,I)=STND(JJ,I)
	       IF(TBD(N).EQ.0.0) TK(JBODYM,I)=TK(JJ,I)
               IF(TBD(N).EQ.1.0) TK(JBODYM,I)=TK(JBODYM,I)
               IF(TBD(N).GT.1.0) TK(JBODYM,I)=TBD(N)/TSTR
               TKD=TK(JBODYM,I)*TSTR
               TKD2=TKD*TKD
               TKD3=TKD*TKD2
               TKD4=TKD*TKD3
               TKD5=TKD*TKD4
               IPOLY=1
               IF(TKD.GT.TPOL1) IPOLY=8
               RBAR=0.0
               ENTH=0.0
               DO 512 ISP=1,LSP
               FSP(JBODYM,I,ISP)=FSP(JJ,I,ISP)
               RBAR=RBAR+FSP(JBODYM,I,ISP)/WM(ISP)
               ENTH=ENTH+FSP(JBODYM,I,ISP)*((POLSP(IPOLY,ISP)*TKD
     1         +POLSP(IPOLY+1,ISP)*TKD2/2.0+POLSP(IPOLY+2,ISP)*TKD3/3.0
     2         +POLSP(IPOLY+3,ISP)*TKD4/4.0+POLSP(IPOLY+4,ISP)*TKD5/5.0
     3         +POLSP(IPOLY+5,ISP))*HCONV/WM(ISP)-HFO(ISP))
  512          CONTINUE
               RHO(JBODYM,I)=(PBACK)/(BETA1*TK(JBODYM,I)*RBAR)
               RHONP(JBODYM,I)=RHO(JBODYM,I)
               HT(JBODYM,I)=ENTH
               END IF
  510          CONTINUE
            DO 520 J=JBODYM,JBODYP
            IF(ISKIP(J,IBODYP).GE.10) THEN
	       II=IBODYP+1
	       IF(II.GT.LI) II=LI
               W(J,IBODYP)=0.0
               STMF(J,IBODYP)=STMF(J,II)
               STND(J,IBODYP)=STND(J,II)
	       IF(TBD(N).EQ.0.0) TK(J,IBODYP)=TK(J,II)
               IF(TBD(N).EQ.1.0) TK(J,IBODYP)=TK(J,IBODYP)
               IF(TBD(N).GT.1.0) TK(J,IBODYP)=TBD(N)/TSTR
               TKD=TK(J,IBODYP)*TSTR
               TKD2=TKD*TKD
               TKD3=TKD*TKD2
               TKD4=TKD*TKD3
               TKD5=TKD*TKD4
               IPOLY=1
               IF(TKD.GT.TPOL1) IPOLY=8
               RBAR=0.0
               ENTH=0.0
               DO 521 ISP=1,LSP
               FSP(J,IBODYP,ISP)=FSP(J,II,ISP)
               RBAR=RBAR+FSP(J,IBODYP,ISP)/WM(ISP)
               ENTH=ENTH+FSP(J,IBODYP,ISP)*((POLSP(IPOLY,ISP)*TKD
     1         +POLSP(IPOLY+1,ISP)*TKD2/2.0+POLSP(IPOLY+2,ISP)*TKD3/3.0
     2         +POLSP(IPOLY+3,ISP)*TKD4/4.0+POLSP(IPOLY+4,ISP)*TKD5/5.0
     3         +POLSP(IPOLY+5,ISP))*HCONV/WM(ISP)-HFO(ISP))
  521          CONTINUE
               RHO(J,IBODYP)=(PBACK)/(BETA1*TK(J,IBODYP)*RBAR)
               RHONP(J,IBODYP)=RHO(J,IBODYP)
               HT(J,IBODYP)=ENTH
               END IF
            IF(ISKIP(J,IBODYM).GE.10) THEN
	       II=IBODYM-1
	       IF(II.LT.1) II=1
               W(J,IBODYM)=0.0
               STMF(J,IBODYM)=STMF(J,II)
               STND(J,IBODYM)=STND(J,II)
               IF(TBD(N).EQ.0.0) TK(J,IBODYM)=TK(J,II)
               IF(TBD(N).EQ.1.0) TK(J,IBODYM)=TK(J,IBODYM)
               IF(TBD(N).GT.1.0) TK(J,IBODYM)=TBD(N)/TSTR
               TKD=TK(J,IBODYM)*TSTR
               TKD2=TKD*TKD
               TKD3=TKD*TKD2
               TKD4=TKD*TKD3
               TKD5=TKD*TKD4
               IPOLY=1
               IF(TKD.GT.TPOL1) IPOLY=8
               RBAR=0.0
               ENTH=0.0
               DO 522 ISP=1,LSP
               FSP(J,IBODYM,ISP)=FSP(J,II,ISP)
               RBAR=RBAR+FSP(J,IBODYM,ISP)/WM(ISP)
               ENTH=ENTH+FSP(J,IBODYM,ISP)*((POLSP(IPOLY,ISP)*TKD
     1         +POLSP(IPOLY+1,ISP)*TKD2/2.0+POLSP(IPOLY+2,ISP)*TKD3/3.0
     2         +POLSP(IPOLY+3,ISP)*TKD4/4.0+POLSP(IPOLY+4,ISP)*TKD5/5.0
     3         +POLSP(IPOLY+5,ISP))*HCONV/WM(ISP)-HFO(ISP))
  522          CONTINUE
               RHO(J,IBODYM)=(PBACK)/(BETA1*TK(J,IBODYM)*RBAR)
               RHONP(J,IBODYM)=RHO(J,IBODYM)
               HT(J,IBODYM)=ENTH
               END IF
  520          CONTINUE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
            DO 531 I=IBODYM+1,IBODYP-1
            DO 530 J=JBODYM+1,JBODYP-1
            RHONP(J,I)=1.0
            RHO(J,I)=1.0
C           TK(J,I)=TK(JBODYM,I)
C           HT(J,I)=HT(JBODYM,I)
            W(J,I)=0.0
            STMF(J,I)=0.0
            STND(J,I)=0.0
            VTC(J,I)=2000.0
	    IF(N.EQ.2) VTC(J,I)=2000.0
            TMU(J,I)=0.0
            DO 532 ISP=1,LSP-1
            FSP(J,I,ISP)=0.0
  532       CONTINUE
            FSP(J,I,LSP)=1.0
  530       CONTINUE
  531       CONTINUE
  500       CONTINUE
            END IF
C--------------------       INSERT INJECTIONS       --------------------
      IF(NFINJ.GE.1) THEN
            DO 600 N=1,NFINJ
            IFINJM=IFIM(N)
            IFINJP=IFIP(N)
            JFINJM=JFIM(N)
            JFINJP=JFIP(N)
            DO 602 I=IFINJM,IFINJP
            DO 608 J=JFINJM,JFINJP
            U(J,I)=FFI(N,1)
            V(J,I)=FFI(N,2)
            W(J,I)=FFI(N,3)
            TK(J,I)=FFI(N,4)
            AK(J,I)=FFI(N,5)
            EPS(J,I)=FFI(N,6)
            STMF(J,I)=0.0
            STND(J,I)=0.0
            DO 603 ISP=1,LSP-1
            FSP(J,I,ISP)=0.0
  603       CONTINUE
            FSP(J,I,INSP)=FFI(N,7)
            FSP(J,I,2)=FFI(N,8)
            FSP(J,I,ININJ)=FFI(N,9)
            FSP(J,I,LSP)=FFI(N,10)
            RHO(J,I)=FFI(N,11)
            P(J,I)=FFI(N,12)
            IF(FFI(N,1).GT.0.0) P(J,I)=P(J,IFINJP+1)
            IF(FFI(N,1).LT.0.0) P(J,I)=P(J,IFINJM-1)
            IF(FFI(N,2).GT.0.0) P(J,I)=P(JFINJP+1,I)
            IF(FFI(N,2).LT.0.0) P(J,I)=P(JFINJM-1,I)
            HT(J,I)=FFI(N,13)
  608       CONTINUE
  602       CONTINUE
C--------------------------for interior boundary conditions
C           IF(FFI(N,2).GT.0.0) THEN
C           DO 604 I=IFINJM,IFINJP
C           DO 604 J=JFINJM,JFINJP-1
C           RHO(J,I)=3.0*RHO(JFINJP+1,I)-2.0*RHO(JFINJP,I)
C           U(J,I)=3.0*U(JFINJP+1,I)-2.0*U(JFINJP,I)
C           V(J,I)=V(JFINJP,I)
C           W(J,I)=3.0*W(JFINJP+1,I)-2.0*W(JFINJP,I)
C 604       CONTINUE
C           END IF
C           IF(FFI(N,2).LT.0.0) THEN
C           DO 605 I=IFINJM,IFINJP
C           DO 605 J=JFINJM+1,JFINJP
C           RHO(J,I)=3.0*RHO(JFINJM-1,I)-2.0*RHO(JFINJM,I)
C           U(J,I)=3.0*U(JFINJM-1,I)-2.0*U(JFINJM,I)
C           V(J,I)=V(JFINJM,I)
C           W(J,I)=3.0*W(JFINJM-1,I)-2.0*W(JFINJM,I)
C 605       CONTINUE
C           END IF
C           IF(FFI(N,1).GT.0.0) THEN
C           DO 606 I=IFINJM,IFINJP-1
C           DO 606 J=JFINJM,JFINJP
C           RHO(J,I)=3.0*RHO(J,IFINJP+1)-2.0*RHO(J,IFINJP)
C           U(J,I)=U(J,IFINJP)
C           V(J,I)=3.0*V(J,IFINJP+1)-2.0*V(J,IFINJP)
C           W(J,I)=3.0*W(J,IFINJP+1)-2.0*W(J,IFINJP)
C 606       CONTINUE
C           END IF
C           IF(FFI(N,1).LT.0.0) THEN
C           DO 607 I=IFINJM+1,IFINJP
C           DO 607 J=JFINJM,JFINJP
C           RHO(J,I)=3.0*RHO(J,IFINJM-1)-2.0*RHO(J,IFINJM)
C           U(J,I)=U(J,IFINJM)
C           V(J,I)=3.0*V(J,IFINJM-1)-2.0*V(J,IFINJM)
C           W(J,I)=3.0*W(J,IFINJM-1)-2.0*W(J,IFINJM)
C 607       CONTINUE
C           END IF
  600       CONTINUE
            END IF
      RETURN                                                            
      END