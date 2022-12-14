      SUBROUTINE STDATA(TTIME)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)   
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LRX=544,
     1   IPMAX=200,JPMAX=600,LNPR=10,LNMX=10)
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),
     1  W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDF(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/CB08/ TPOL1,TPOL2,POLSP(14,LSP),CISP(12,LSP)
      COMMON/HFORM/ HFO(LSP)
      COMMON/HRLS/ QDOT(LJ,LI)
      COMMON/BODY/IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10),NBODY
      COMMON/DUMD/ QDSP(LJ,LI,LSP)
      COMMON/DUMR/ XRG(LJ,LI),YRG(LJ,LI),RHORG(LJ,LI),TKRG(LJ,LI),
     *             URG(LJ,LI),VRG(LJ,LI),QRG(LJ,LI),VMURG(LJ,LI),
     *             VTCRG(LJ,LI),FSPRG(LJ,LI,LSP),VDFRG(LJ,LI,LSP),
     *             QSPRG(LJ,LI,LSP),CP(LJ,LI),ENTH(LJ,LI,LSP),
     *             DUMRB(LRX*LE-218*LE)
C
C------CONVERT DATA FROM STAGGERED GRID TO REGULAR GRID
      XSHIFT=0.5*(X(1,1)+X(1,2))
      YSHIFT=0.5*(Y(1,1)+Y(2,1))
      DO 1320 I=1,LI-1
      DO 1326 J=1,LJ-1
      XRG(J,I)=(0.5*(X(J,I)+X(J,I+1))-XSHIFT)*ALSTR
      YRG(J,I)=(0.5*(Y(J,I)+Y(J+1,I))-YSHIFT)*ALSTR
      RHORG(J,I)=0.25*(RHO(J,I)+RHO(J+1,I)+RHO(J,I+1)+RHO(J+1,I+1))*RSTR
      TKRG(J,I)=0.25*(TK(J,I)+TK(J+1,I)+TK(J,I+1)+TK(J+1,I+1))*TSTR
      QRG(J,I)=0.25*(QDOT(J,I)+QDOT(J+1,I)+QDOT(J,I+1)+QDOT(J+1,I+1))
      VMURG(J,I)=0.25*(VMU(J,I)+VMU(J+1,I)
     1                +VMU(J,I+1)+VMU(J+1,I+1))*AMSTR
      VTCRG(J,I)=0.25*(VTC(J,I)+VTC(J+1,I)
     1                +VTC(J,I+1)+VTC(J+1,I+1))*ACSTR
      URG(J,I)=0.5*(U(J,I+1)+U(J+1,I+1))*VSTR
      VRG(J,I)=0.5*(V(J+1,I)+V(J+1,I+1))*VSTR
      FTOT=0.0
      DO 1321 ISP=1,LSP-1
      VDFRG(J,I,ISP)=0.25*(VDF(J,I,ISP)+VDF(J+1,I,ISP)
     1                  +VDF(J,I+1,ISP)+VDF(J+1,I+1,ISP))*DSTR
      QSPRG(J,I,ISP)=0.25*(QDSP(J,I,ISP)+QDSP(J+1,I,ISP)
     1                  +QDSP(J,I+1,ISP)+QDSP(J+1,I+1,ISP))
      FSPRG(J,I,ISP)=0.25*(FSP(J,I,ISP)+FSP(J+1,I,ISP)
     1                  +FSP(J,I+1,ISP)+FSP(J+1,I+1,ISP))
      FTOT=FTOT+FSPRG(J,I,ISP)
 1321 CONTINUE
      FSPRG(J,I,LSP)=1.0-FTOT
      VDFRG(J,I,LSP)=0.25*(VDF(J,I,LSP)+VDF(J+1,I,LSP)
     1                  +VDF(J,I+1,LSP)+VDF(J+1,I+1,LSP))*DSTR
      QSPRG(J,I,LSP)=0.25*(QDSP(J,I,LSP)+QDSP(J+1,I,LSP)
     1                  +QDSP(J,I+1,LSP)+QDSP(J+1,I+1,LSP))
 1326 CONTINUE
 1320 CONTINUE
      DO 1322 I=1,LI
      XRG(LJ,I)=XRG(LJ-1,I)
      YRG(LJ,I)=Y(LJ,I)*ALSTR
      RHORG(LJ,I)=RHORG(LJ-1,I)
      TKRG(LJ,I)=TKRG(LJ-1,I)
      QRG(LJ,I)=QRG(LJ-1,I)
      VMURG(LJ,I)=VMURG(LJ-1,I)
      VTCRG(LJ,I)=VTCRG(LJ-1,I)
      URG(LJ,I)=URG(LJ-1,I)
      VRG(LJ,I)=VRG(LJ-1,I)
      DO 1323 ISP=1,LSP
      FSPRG(LJ,I,ISP)=FSPRG(LJ-1,I,ISP)
      QSPRG(LJ,I,ISP)=QSPRG(LJ-1,I,ISP)
      VDFRG(LJ,I,ISP)=VDFRG(LJ-1,I,ISP)
 1323 CONTINUE
 1322 CONTINUE
      DO 1324 J=1,LJ
      XRG(J,LI)=X(J,LI)*ALSTR
      YRG(J,LI)=YRG(J,LI-1)
      RHORG(J,LI)=RHORG(J,LI-1)
      TKRG(J,LI)=TKRG(J,LI-1)
      QRG(J,LI)=QRG(J,LI-1)
      VMURG(J,LI)=VMURG(J,LI-1)
      VTCRG(J,LI)=VTCRG(J,LI-1)
      URG(J,LI)=URG(J,LI-1)
      VRG(J,LI)=VRG(J,LI-1)
      DO 1325 ISP=1,LSP
      FSPRG(J,LI,ISP)=FSPRG(J,LI-1,ISP)
      QSPRG(J,LI,ISP)=QSPRG(J,LI-1,ISP)
      VDFRG(J,LI,ISP)=VDFRG(J,LI-1,ISP)
 1325 CONTINUE
 1324 CONTINUE
C------------------CALCULATION OF TCP AND TENTH-------------------------
C--------GASC=8.3144, CP--J/g-mole/K, H--J/g-mole, WM---g/mole
      DO 10 I=1,LI
      DO 11 J=1,LJ
      TKD=TKRG(J,I)
      TKD2=TKD*TKD
      TKD3=TKD*TKD2
      TKD4=TKD*TKD3
      TKD5=TKD*TKD4
      IPOLY=1
      IF(TKD.GT.TPOL1) IPOLY=8
      TCP=0.0
      TENTH=0.0
      DO 12 ISP=1,LSP
C-----------------CP
      TCP=TCP+(POLSP(IPOLY,ISP)+POLSP(IPOLY+1,ISP)*TKD
     1 +POLSP(IPOLY+2,ISP)*TKD2+POLSP(IPOLY+3,ISP)*TKD3
     2 +POLSP(IPOLY+4,ISP)*TKD4)*FSPRG(J,I,ISP)*GASC/WM(ISP)/WMSTR
C------------------ENTHALPY
      ENTH(J,I,ISP)=((POLSP(IPOLY,ISP)*TKD
     1    +POLSP(IPOLY+1,ISP)*TKD2/2.0+POLSP(IPOLY+2,ISP)*TKD3/3.0
     2    +POLSP(IPOLY+3,ISP)*TKD4/4.0+POLSP(IPOLY+4,ISP)*TKD5/5.0
     3    +POLSP(IPOLY+5,ISP))*GASC/WM(ISP)/WMSTR)
C    4    -HFO(ISP)*HSTR)*FSPRG(J,I,ISP)
      TENTH=TENTH+ENTH(J,I,ISP)
   12 CONTINUE
      CP(J,I)=TCP
   11 CONTINUE
   10 CONTINUE
C----------------------------WRITING THE DATA---------------------------
      WRITE(13,102) ITR,TTIME
      WRITE(13,105)
  102 FORMAT(1X,I6,F15.10)
      DO 1340 I=1,LI
      DO 1341 J=1,LJ
C---I,J,XRG(J,I),YRG(J,I),
C---    RHO(J,I),U(J,I),V(J,I),T(J,I),QDOT(J,I),
C---    MU(J,I),TC(J,I),CP(J,I),
C---    FSP(J,I,ISP),DSP(J,I,ISP),QSP(J,I,ISP),ENTH(J,I,ISP)
      WRITE(13,104) I,J,XRG(J,I),YRG(J,I),RHORG(J,I),
     *          URG(J,I),VRG(J,I),
     *          TKRG(J,I),QRG(J,I),VMURG(J,I),VTCRG(J,I),CP(J,I),
     *         (FSPRG(J,I,ISP),VDFRG(J,I,ISP),QSPRG(J,I,ISP),
     *          ENTH(J,I,ISP),ISP=1,LSP)
 1341 CONTINUE
 1340 CONTINUE
  104 FORMAT(I4,',',I4,',',E15.7E3,217(',',E15.7E3))
  105 FORMAT('VARIABLES="i","j","z","r","RHO","v","u","T","QDOT",',
     * '"VMU","VTC","CP",',
     * '"Y-H2","Y-O2","Y-H","Y-O","Y-OH","Y-H2O","Y-HO2","Y-H2O2",',
     * '"Y-CO","Y-CO2","Y-HCO","Y-CH2O","Y-CH4","Y-CH3","Y-T-CH2",',
     * '"Y-S-CH2","Y-C2H4","Y-CH3O","Y-C2H5","Y-C2H6","Y-CH","Y-C2H2",',
     * '"Y-C2H3","Y-CH2CHO","Y-C2H4O","x-CH2CO","Y-HCCO","Y-C2H",',
     * '"Y-CH2OH","Y-CH3OH","Y-C2H5OH","Y-CH3CHO","Y-CH3CHOH",',
     * '"Y-CH2CH2OH","Y-CH3CO","Y-CH3CH2O","Y-C3H4","Y-C3H3","Y-C3H5",',
     * '"Y-C3H6","Y-C3H8","Y-I-C3H7","Y-N-C3H7","Y-C4H6","Y-CHO",',
     * '"Y-C5H8","Y-C7H16","Y-C4H8","Y-C5H10","Y-AR","Y-HE","Y-N2",',
     * '"VDF-H2","VDF-O2","VDF-H","VDF-O","VDF-OH","VDF-H2O",',
     * '"VDF-HO2","VDF-H2O2","VDF-CO","VDF-CO2","VDF-HCO","VDF-CH2O",',
     * '"VDF-CH4","VDF-CH3","VDF-T-CH2","VDF-S-CH2","VDF-C2H4",',
     * '"VDF-CH3O","VDF-C2H5","VDF-C2H6","VDF-CH","VDF-C2H2",',
     * '"VDF-C2H3","VDF-CH2CHO","VDF-C2H4O","VDF-CH2CO","VDF-HCCO",',
     * '"VDF-C2H","VDF-CH2OH","VDF-CH3OH","VDF-C2H5OH","VDF-CH3CHO",',
     * '"VDF-CH3CHOH","VDF-CH2CH2OH","VDF-CH3CO","VDF-CH3CH2O",',
     * '"VDF-C3H4","VDF-C3H3","VDF-C3H5","VDF-C3H6","VDF-C3H8",',
     * '"VDF-I-C3H7","VDF-N-C3H7","VDF-C4H6","VDF-CHO","VDF-C5H8",',
     * '"VDF-C7H16","VDF-C4H8","VDF-C5H10","VDF-AR","VDF-HE","VDF-N2"',
     * '"PRD-H2","PRD-O2","PRD-H","PRD-O","PRD-OH","PRD-H2O",',
     * '"PRD-HO2","PRD-H2O2","PRD-CO","PRD-CO2","PRD-HCO","PRD-CH2O",',
     * '"PRD-CH4","PRD-CH3","PRD-T-CH2","PRD-S-CH2","PRD-C2H4",',
     * '"PRD-CH3O","PRD-C2H5","PRD-C2H6","PRD-CH","PRD-C2H2",',
     * '"PRD-C2H3","PRD-CH2CHO","PRD-C2H4O","PRD-CH2CO","PRD-HCCO",',
     * '"PRD-C2H","PRD-CH2OH","PRD-CH3OH","PRD-C2H5OH","PRD-CH3CHO",',
     * '"PRD-CH3CHOH","PRD-CH2CH2OH","PRD-CH3CO","PRD-CH3CH2O",',
     * '"PRD-C3H4","PRD-C3H3","PRD-C3H5","PRD-C3H6","PRD-C3H8",',
     * '"PRD-I-C3H7","PRD-N-C3H7","PRD-C4H6","PRD-CHO","PRD-C5H8",',
     * '"PRD-C7H16","PRD-C4H8","PRD-C5H10","PRD-AR","PRD-HE","PRD-N2"',
     * '"H-H2","H-O2","H-H","H-O","H-OH","H-H2O","H-HO2","H-H2O2",',
     * '"H-CO","H-CO2","H-HCO","H-CH2O","H-CH4","H-CH3","H-T-CH2",',
     * '"H-S-CH2","H-C2H4","H-CH3O","H-C2H5","H-C2H6","H-CH","H-C2H2",',
     * '"H-C2H3","H-CH2CHO","H-C2H4O","x-CH2CO","H-HCCO","H-C2H",',
     * '"H-CH2OH","H-CH3OH","H-C2H5OH","H-CH3CHO","H-CH3CHOH",',
     * '"H-CH2CH2OH","H-CH3CO","H-CH3CH2O","H-C3H4","H-C3H3","H-C3H5",',
     * '"H-C3H6","H-C3H8","H-I-C3H7","H-N-C3H7","H-C4H6","H-CHO",',
     * '"H-C5H8","H-C7H16","H-C4H8","H-C5H10","H-AR","H-HE","H-N2"')      
      RETURN
      END