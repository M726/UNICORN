
      PROGRAM UNICORN
      IMPLICIT REAL *8 (A-H,O-Z)   
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LJ2=LJ*2,LSP=52,LRX=544,               
     1   LB=LJ-2,LL=(LI-2)*LB,LL1=LB*LL-(LB+1)*LB/2,LPD=43*LE,
     2   IPMAX=500,JPMAX=200,LNPR=200,LNMX=2000)
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LE),Y(LE),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LE),FSP(LE,LSP),U(LE),V(LE),W(LE),
     1             P(LE),HT(LE),TK(LE),AK(LE),EPS(LE)
      COMMON/CB05/ RHONP(LE)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1 HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LE),VTC(LE),VDFSP(LE,LSP),TMU(LE)
      COMMON/CB08/ TPOL1,TPOL2,POLSP(14,LSP),CISP(12,LSP)
      COMMON/CB09/ IBOT(LI),ITOP(LI),JLFT(LJ),JRGT(LJ),
     1       FBOT(8+LSP,LI),FTOP(8+LSP,LI),FLFT(8+LSP,LJ),FRGT(8+LSP,LJ)
      COMMON/CB10/ FBXP(4,LJ),FBXM(4,LJ),FBYP(4,LI),FBYM(4,LI)
      COMMON/BODY/IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10),NBODY
      COMMON/FINJ/IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30),NFINJ
      COMMON/LUMAT/ IDATLU(LB),PSL(LL1),PSD(LL)
      COMMON/REAC/ LREV,NRLIND,ITBEF(LRX),ALF(LRX),AKF(LRX),EXAF(LRX),
     1             ALLOW(LRX),AKLOW(LRX),EALOW(LRX),TROE(LRX,4)
      COMMON/HFORM/ HFO(LSP)
      COMMON/HRLS/ QDOT(LE)
      COMMON/KEPS/ AKREF,EPREF,C1KE,C2KE,CMKE,CMKEQ,EE,AVON,PCON
      COMMON/SOOT/ RSOOT,STMF(LE),STND(LE),STDF(LE)
      COMMON/DUMMY/ UTMP(LE),AP(LE),AEE(LE),AE(LE),AW(LE),AWW(LE),
     1  ANN(LE),AN(LE),AS(LE),ASS(LE),RHS1(LE),FP(LE),DUM(LPD)
      COMMON/DUM0/ DUM0(149*LE)
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
      COMMON/DUMD/ QDSP(LJ,LI,LSP)
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
      COMMON/DUMX/ DUMX(4*LE)
      COMMON/DUMY/ DUMY(LSP*LE+3*LE)
      COMMON/DUMR/ DUMR(LRX*LE)
      COMMON/DUMS/ DUMS(12*LE)
      COMMON/NOISE/ INOISE,IVXTYP,FNOISE,ANOISE(200),VXTYP(5),TTIME
      COMMON/TRACK/ IPICK(101),JPICK(101),XPR(LNPR),YPR(LNPR),
     1              PATHX(LNMX,LNPR),PATHY(LNMX,LNPR),
     2              PUVEL(LNMX,LNPR),PVVEL(LNMX,LNPR),
     3              GUVEL(LNMX,LNPR),GVVEL(LNMX,LNPR)
      COMMON/FLOW/ RASMN,RASMX,RAS2MN,RAS2MX,XS(LI),YS(LJ+LJ-1),KSYM,LIP,LJP,IOFF,JOFF,NCLR(LNPR)
      COMMON/APIXL/ IPAGE,NCOLOR,NBACK,NAA(IPMAX,JPMAX)
      COMMON/FRAME/ IPIXL,JPIXL,XMIN,YMIN,XMAX,YMAX,XPIXL,YPIXL
      COMMON/PEN/NPEN,NAPEN,NBPEN
      DIMENSION UAVE(LE),VAVE(LE),TAVE(LE),FSPAV(LE,LSP)
      DIMENSION ISIDE(10),ITYPE(10),BCVAL(8+LSP,10)
      DIMENSION IREGN(10),XREGN(10),JREGN(10),YREGN(10)
      DIMENSION IRSW(LRX),CI11(4,5),WSP(LSP),CDIA(LSP),TLJ(LSP),
     1 ETLJ1(LSP,LSP),ETLJ2(LSP,LSP),ETLJ3(LSP,LSP),ETLJ4(LSP,LSP),
     2 PHIK(LSP),PHDIK(LSP),EFDIK(LSP),WMCD(LSP,LSP),WMIK(LSP,LSP),
     3 DWMIK(LSP,LSP),XSP(LSP),AMUSP(LSP),ATCSP(LSP),DIFSP(LSP,LSP)
      DIMENSION ISPQ(LSP),ISPQR(LSP)
      DIMENSION WMQ(LSP),CDIAQ(LSP),ETLJ1Q(LSP,LSP),ETLJ2Q(LSP,LSP),
     1 ETLJ3Q(LSP,LSP),ETLJ4Q(LSP,LSP),WMCDQ(LSP,LSP),WMIKQ(LSP,LSP),
     2 DWMIKQ(LSP,LSP),VDFQ(LSP),POLQ(14,LSP),CIQ(12,LSP)
      DIMENSION XBDM(10),YBDM(10),XBDP(10),YBDP(10)
      DIMENSION XFIM(10),YFIM(10),XFIP(10),YFIP(10)
      DIMENSION IXEVOL(10),XEVOL(10),AL1E(10),FEVOL(LJ,4,10)
      DIMENSION IDRV(10),ADRV(10)
      CHARACTER *25 FILEIN,FILEOT,FILE6,FILETM,FILEPT,FILEDR,FILEDG,
     1              FILEMV,FILEAV
      CHARACTER *3 SCHMU,SCHMV,SCHMW,SCHMH,SCHMSP,SCHMKE,
     1             SCHM1,SCHM2
      CHARACTER *5 EJTYPE(5)
      CHARACTER *6 FUELTYPES(11)
      DATA SCHM1,SCHM2 /'ADI','PNT'/
      DATA EJTYPE /'FUEL ',' AIR ','  N2 ','SUPPR','LOCAL'/
C
C
      DATA WSP / 2.016D+00, 3.200D+01, 1.008D+00, 1.600D+01, 1.701D+01,
     *           1.802D+01, 3.301D+01, 3.401D+01, 2.801D+01, 4.401D+01,
     *           2.902D+01, 3.003D+01, 1.604D+01, 1.504D+01, 1.403D+01,
     *           1.403D+01, 2.805D+01, 3.103D+01, 2.906D+01, 3.007D+01,
     *           1.302D+01, 2.604D+01, 2.705D+01, 4.305D+01, 4.405D+01,
     *           4.204D+01, 4.103D+01, 2.503D+01, 3.103D+01, 3.204D+01,
     *           4.607D+01, 4.405D+01, 4.506D+01, 4.506D+01, 4.305D+01,
     *           4.506D+01, 4.007D+01, 3.906D+01, 4.107D+01, 4.208D+01,
     *           4.410D+01, 4.309D+01, 4.309D+01, 5.409D+01, 2.902D+01,
     *           6.812D+01, 1.002D+02, 5.611D+01, 7.014D+01, 3.995D+01,
     *           4.003D+00, 2.801D+01/
      DATA CDIA/ 2.920D+00, 3.458D+00, 2.050D+00, 2.750D+00, 2.750D+00,
     *           2.605D+00, 3.458D+00, 3.458D+00, 3.650D+00, 3.763D+00,
     *           3.590D+00, 3.590D+00, 3.746D+00, 3.800D+00, 3.800D+00,
     *           3.800D+00, 3.496D+00, 3.690D+00, 4.350D+00, 4.350D+00,
     *           2.750D+00, 3.721D+00, 3.721D+00, 3.970D+00, 3.970D+00,
     *           3.970D+00, 2.500D+00, 3.721D+00, 3.690D+00, 3.626D+00,
     *           4.410D+00, 3.970D+00, 4.530D+00, 4.530D+00, 3.970D+00,
     *           4.410D+00, 4.290D+00, 4.290D+00, 4.220D+00, 4.140D+00,
     *           4.810D+00, 4.810D+00, 4.810D+00, 4.720D+00, 3.590D+00,
     *           5.200D+00, 6.253D+00, 4.650D+00, 5.489D+00, 3.330D+00,
     *           2.576D+00, 3.621D+00/
      DATA TLJ / 3.800D+01, 1.074D+02, 1.450D+02, 8.000D+01, 8.000D+01,
     *           5.724D+02, 1.074D+02, 1.074D+02, 9.810D+01, 2.440D+02,
     *           4.980D+02, 4.980D+02, 1.414D+02, 1.440D+02, 1.440D+02,
     *           1.440D+02, 2.384D+02, 4.170D+02, 2.475D+02, 2.475D+02,
     *           8.000D+01, 2.653D+02, 2.653D+02, 4.360D+02, 4.360D+02,
     *           4.360D+02, 1.500D+02, 2.653D+02, 4.170D+02, 4.818D+02,
     *           4.706D+02, 4.360D+02, 3.626D+02, 3.626D+02, 4.360D+02,
     *           4.706D+02, 3.248D+02, 3.248D+02, 3.160D+02, 3.078D+02,
     *           3.034D+02, 3.034D+02, 3.034D+02, 3.570D+02, 4.980D+02,
     *           4.080D+02, 4.596D+02, 3.550D+02, 3.862D+02, 1.365D+02,
     *           1.020D+01, 9.753D+01/
      DATA CI11/ 4.00443D+00, 2.08578D+00, 1.04197D+00, 6.89404D-01,
     1          -5.20345D+00,-8.30180D-01,-2.98978D-02,-1.23293D-03,
     2           4.20192D+00, 2.54502D-01, 1.00205D-03, 2.81316D-06,
     3          -1.68184D+00,-3.79990D-02,-1.55597D-05,-3.09971D-09,
     4           2.64301D-01, 2.20751D-03, 8.74264D-08, 1.26015D-12/
      DATA FUELTYPES /'H2','CH4','CH3OH','C2H2','C2H4','C2H6','C3H8',
     1           'C3H6','CH2O','CO','C7H16'/    
      CALL RESET
CDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
C
C       Read input.uni file here
C
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
      MP=15
      OPEN(MP,FILE='cfg/input.uni',ACCESS='SEQUENTIAL',STATUS='OLD')
      READ(MP,1) ICAL
    1 FORMAT(I1)
      READ(MP,*) ISYM,IREAD,IGNIT
      READ(MP,*) ISTDY,INOISE,XNOISE,YNOISE,ANOISE(1),FNOISE,
     1           IVXTYP,VXTYP(1),VXTYP(2),VXTYP(3),VXTYP(4)
      NEJEC1=0
      IF(INOISE.EQ.99) THEN
          INOISE=0
          READ(MP,*) NEJEC1,FEJEC,XEJEC1,XEJEC2,YEJEC,
     1           AMEJEC,DTEJEC,MTEJEC
          END IF
      IF(INOISE.LT.0) THEN
          IWIDTH=-INOISE
          IF(INOISE.LT.-100) IWIDTH=IWIDTH-100
          NT1=INT((IWIDTH-1)/10)
          NT2=IWIDTH-NT1*10
          DO 1001 N=1,NT1
          NA=(N-1)*10+1
          READ(MP,*) (ANOISE(I),I=NA,NA+9)
 1001     CONTINUE
          NA=NT1*10+1
          IF(NT2.GE.1) READ(MP,*) (ANOISE(I),I=NA,NA+NT2)
          IF(INOISE.LT.-100) THEN
             DO 1002 N=1,NT1
             NA=100+(N-1)*10+1
             READ(MP,*) (ANOISE(I),I=NA,NA+9)
 1002        CONTINUE
             NA=100+NT1*10+1
             IF(NT2.GE.1) READ(MP,*) (ANOISE(I),I=NA,NA+NT2)
             END IF
          END IF
      READ(MP,*) RTIN,RTOT,ALENG                  
      READ(MP,*) VREF,TREF,PREF,RREF,AKREF,EPREF,FO2IN,INFUEL
      READ(MP,*) IFLOW,ISWRL,ITHRM,ICHEM,IPROP,IGRAV
      READ(MP,*) NSEG
      DO 1000 N=1,NSEG
      READ(MP,*) ISIDE(N),ITYPE(N),(BCVAL(I,N),I=1,LSP+7)
 1000 CONTINUE
      READ(MP,*) NBODY
      IF(NBODY.GT.0) THEN 
          DO 1003 N=1,NBODY
          READ(MP,*) XBDM(N),YBDM(N),XBDP(N),YBDP(N),TBD(N)
 1003     CONTINUE
          END IF
      READ(MP,*) NFINJ
      IF(NFINJ.GT.0) THEN 
          DO 1004 N=1,NFINJ
          READ(MP,*) XFIM(N),YFIM(N),XFIP(N),YFIP(N),(FFI(N,I),I=1,9)
 1004     CONTINUE
          END IF
      READ(MP,*) NIREGN,(IREGN(I),XREGN(I),I=1,NIREGN)
      READ(MP,*) NJREGN,(JREGN(J),YREGN(J),J=1,NJREGN)
      READ(MP,*) ITEND,ISECS,CFLNO,ISTORE,ISTB
      READ(MP,*) ITPRNT,IPRES
      READ(MP,*) SCHMU,SCHMV,SCHMW,SCHMH,SCHMSP,SCHMKE
      READ(MP,*) ISORU,ISORV,ISORW,ISORH,ISORSP,ISORKE
      READ(MP,*) RELXU,RELXV,RELXW,RELXH,RELXSP,RELXKE
      READ(MP,*) TOLRU,TOLRV,TOLRW,TOLRH,TOLRSP,TOLRKE
      READ(MP,*) EAG,AKG,(IRSW(I),I=1,LRX)
      READ(MP,*) IBEVOL,ISEVOL,NEVOL,(XEVOL(N),N=1,NEVOL)
      READ(MP,*) IBDRV,NDRV,(IDRV(I),I=1,10)
      READ(MP,*) IBDRG,ISDRG
      READ(MP,*) NOPR,IBINJ,ITINJ,IEINJ,PDIA,PDEN,PTHR,PVEL
          IF(NOPR.GT.0) THEN
          NT1=INT((NOPR-1)/4)
          NT2=NOPR-NT1*4
          DO 1005 N=1,NT1
          NA=(N-1)*4+1
          READ(MP,*) (XPR(I),YPR(I),I=NA,NA+3)
 1005     CONTINUE
          NA=NT1*4+1
          IF(NT2.GE.1) READ(MP,*) (XPR(I),YPR(I),I=NA,NA+NT2)
          END IF
          IF(NOPR.LT.0) THEN
          NOPR=-NOPR
          NNA=1
 1006     CONTINUE
          READ(MP,*) NN,XP1,YP1,XP2,YP2
          NNB=NN+NNA-1
          IF(NNB.GT.NOPR) NNB=NOPR
          XPR(NNA)=XP1
          YPR(NNA)=YP1
          XPR(NNB)=XP2
          YPR(NNB)=YP2
          XPRD=(XPR(NNB)-XPR(NNA))/FLOAT(NN-1)
          YPRD=(YPR(NNB)-YPR(NNA))/FLOAT(NN-1)
          DO 1007 N=NNA+1,NNB-1
          XPR(N)=XPR(NNA)+XPRD*FLOAT(N-NNA)
          YPR(N)=YPR(NNA)+YPRD*FLOAT(N-NNA)
 1007     CONTINUE
          NNA=NNB+1
          IF(NNB.LT.NOPR) GO TO 1006
          END IF
      READ(MP,*) IBANM,ISANM,KSYM,IPANM,X1ANM,X2ANM,Y1ANM,Y2ANM,
     1           NFSURF,KORNT
      READ(MP,*) NBAVE,NEAVE
      READ(MP,*) FILEIN
      READ(MP,*) FILEOT
      READ(MP,*) FILETM
      READ(MP,*) FILEDR
      READ(MP,*) FILEDG
      READ(MP,*) FILEPT
      READ(MP,*) FILEMV
      READ(MP,*) FILEAV
      READ(MP,*) FILE6
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
      IPASS=0
      IGKEEP=20
      IF(IFLOW.GE.10) THEN
                      IPASS=1
                      IFLOW=IFLOW/10
                      IF(ISTORE.GE.1) WRITE(*,802) IFLOW
  802 FORMAT(2X,'DO YOU WANT TO STORE THE DATA ?  IFLOW=',I3)
                      ISTORE=0
                      END IF
      IF(IPRES.EQ.0) IPRES=ITEND+1
      NTPRNT=DABS(DFLOAT(ITPRNT)/2.0)
      ITRANS=0
      IF(ITPRNT.LT.0) THEN
                      ITRANS=1
                      ITPRNT=-ITPRNT
                      END IF
      IF(ITPRNT.EQ.2) ITPRNT=1
      IF(ITPRNT.EQ.1.AND.ITEND.GE.1) ITPRNT=ITEND
      IF(ITPRNT.EQ.0) ITPRNT=ITEND+1
      NSTORE=2
      IF(ISTORE.LT.0) THEN
                      NSTORE=1
                      ISTORE=-ISTORE
                      END IF
      IF(ISTORE.EQ.1.AND.ITEND.GE.1) ISTORE=ITEND
      IF(ISTORE.EQ.0) ISTORE=ITEND+1
      IF(NEAVE.EQ.1) NEAVE=ITEND
C OPEN OPEN OPEN OPEN OPEN OPEN OPEN  OPEN OPEN OPEN OPEN OPEN OPEN OPEN
           OPEN(16,FILE=FILE6,ACCESS='SEQUENTIAL',STATUS='OLD')
      IF(NEVOL.GE.1)
     1     OPEN(14,FILE=FILETM,ACCESS='SEQUENTIAL',STATUS='OLD')
      IF(ISTORE.EQ.1.OR.ISTORE.LE.ITEND) 
     1     OPEN(12,FILE=FILEOT,ACCESS='SEQUENTIAL',STATUS='OLD')
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
      IF(ISTORE.EQ.1.OR.ISTORE.LE.ITEND) 
     1     OPEN(13,FILE='FDETAIL.DATA',ACCESS='SEQUENTIAL',
     2             STATUS='UNKNOWN')
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
      IF(NOPR.GT.0)
     1  OPEN(20,FILE=FILEPT,ACCESS='SEQUENTIAL',STATUS='OLD')
      IF(IREAD.GE.1) 
     1     OPEN(32,FILE=FILEIN,ACCESS='SEQUENTIAL',STATUS='OLD')
      IF(NDRV.GE.1) 
     1     OPEN(72,FILE=FILEDR,ACCESS='SEQUENTIAL',STATUS='OLD')
      IPAGE=1
      IF(ISANM.LT.0) THEN
                     ISANM=-ISANM
                     IPAGE=0
                     END IF
      IPIXL=IPMAX
      JPIXL=JPMAX
      IF((KORNT.EQ.2).OR.(KORNT.EQ.4)) THEN
                 JJ=FLOAT(JPIXL)/4.0
                 JPIXL=JJ*4
                 ELSE
                 JJ=FLOAT(IPIXL)/4.0
                 IPIXL=JJ*4
                 END IF
      IF(NBAVE.GE.1) 
     1      OPEN(13,FILE=FILEAV,ACCESS='SEQUENTIAL',STATUS='OLD')
C OVER OVER OVER OVER OVER OVER OVER  OVER OVER OVER OVER OVER OVER OVER
      INSP=1
      INSIMP=0
      IF(INFUEL.LT.0) THEN
                      INSIMP=1
                      INFUEL=-INFUEL
                      END IF
C-----------------------------------------------------------------------
C   FUEL: 1 - H2,  2--CH4,  3--CH3OH,  4--C2H2, 5--C2H4, 6--C2H6, 
C         7--C3H8, 8--C3H6, 9--CH2O,  10--CO,  11--C7H16
C-----------------------------------------------------------------------
      IF(INFUEL.LT.1.OR.INFUEL.GT.11) THEN
          ICODE=0
          WRITE(*,803) INFUEL
  803     FORMAT(/10X,'--- WHAT FUEL IS THIS ???? ---, INFUEL = ',I4)
          GO TO 102
          END IF
      IF(INFUEL.EQ. 1) INSP=01
      IF(INFUEL.EQ. 2) INSP=13
      IF(INFUEL.EQ. 3) INSP=30
      IF(INFUEL.EQ. 4) INSP=22
      IF(INFUEL.EQ. 5) INSP=17
      IF(INFUEL.EQ. 6) INSP=20
      IF(INFUEL.EQ. 7) INSP=41
      IF(INFUEL.EQ. 8) INSP=40
      IF(INFUEL.EQ. 9) INSP=12
      IF(INFUEL.EQ.10) INSP=09
      IF(INFUEL.EQ.11) INSP=47
      IF(ICHEM.EQ.2) ICHEM=1
      WRITE(*,804) LSP,LRX/LREV,LREV,FUELTYPES(INFUEL)
  804 FORMAT(
     1 /10X,'---FINITE-RATE MODEL FOR H2 AND HYROCARBON FUELS--',
     2 /10X,'        ----(SANDIEGO MECHANISM)----',
     3 /10X,'  ----(',I3,' SPECIES & ',I4,' X ',I1,' REACTIONS)----',
     4 /10X,'--- (FUEL USED IN THIS CALCULATION IS ',A6,') ---')
C
C--VALUES---VALUES---VALUES---VALUES---VALUES---VALUES---VALUES---VALUES
      FSECS=0.95*DFLOAT(ISECS)
      IF(FSECS.LE.0.0) FSECS=1.0D+10
C     CLOCK1=SECOND()
      CLOCK1=0.0
      PI=3.14159265358979312D+00
      GFORCE=9.80665
      GASC=8.3144
      RUATM=82.06
      RUCAL=1.9872
      GAMA=1.4
      C1KE=1.44
      C2KE=1.92
      CMKE=0.09
      CMKEQ=CMKE**(0.25)
      EE=9.0
      AVON=0.4 
      PCON=9.24
      SIGK=1.0
      SIGE=1.21
      SIGSP=0.9
      SIGH=0.9
      SIGW=1.0
      AVGNO=6.022169D+23
C New Const for when ICHEM == 0
      RSOOT=1.9000D+03
      SIG0=5.669D-08
      SIGB=1.38D-23
C--VALUES---VALUES---VALUES---VALUES---VALUES---VALUES---VALUES---VALUES
C--------NON-REVERSIBLE REACTIONS--------NON-REVERSIBLE REACTIONS-------
C        
C        ------NO. OF REVERSE REACTIONS SET TO ZERO --   0
C        
C---HIGH-RATE-COEFFICIENT REACTIONS----HIGH-RATE-COEFFICIENT REACTIONS--
C
C                          NONE
C
C-----------------------------------------------------------------------
      NRLIND=0
      NRTOT=0
      DO 10 I=1,LRX,LREV
      EXAF(I)=EXAF(I)/RUCAL
      IF(IRSW(I).EQ.0) THEN
                       ITBEF(I)=1
                       AKF(I)=0.0
                       EALOW(I)=0.0
                       AKLOW(I)=0.0
                       ALLOW(I)=0.0
                       TROE(I,1)=0.0
                       TROE(I,2)=0.0
                       TROE(I,3)=0.0
                       TROE(I,4)=0.0
                       ELSE
                       NRTOT=NRTOT+1
                       IF(AKF(I).LT.0) THEN
                                 NRLIND=NRLIND+1
                                 AKF(I)=-AKF(I)
                                 EALOW(I)=EALOW(I)/RUCAL
                                 ELSE
                                 ITBEF(I)=1
                                 EALOW(I)=0.0
                                 AKLOW(I)=0.0
                                 ALLOW(I)=0.0
                                 TROE(I,1)=0.0
                                 TROE(I,2)=0.0
                                 TROE(I,3)=0.0
                                 TROE(I,4)=0.0
                                 END IF
                       END IF
      IF(LREV.EQ.2) THEN
                    AKF(I+1)=1.0
                    IF(IRSW(I+1).EQ.0) THEN
                       AKF(I+1)=0.0
                       ELSE
                       NRTOT=NRTOT+1
                       END IF
                    ITBEF(I+1)=1
                    EXAF(I+1)=0.0
                    ALF(I+1)=0.0
                    EALOW(I+1)=0.0
                    AKLOW(I+1)=0.0
                    ALLOW(I+1)=0.0
                    TROE(I+1,1)=0.0
                    TROE(I+1,2)=0.0
                    TROE(I+1,3)=0.0
                    TROE(I+1,4)=0.0
                    END IF
   10 CONTINUE
      WRITE(*,807) NRTOT,NRLIND,LREV
  807 FORMAT(12X,'NO. OF REACTIONS USED IN THIS CALCULATION= ',I4/
     1       12X,'NO. OF REACTIONS IN FALL-OFF REGION = ',I3,' X ',I1/)
      EAGC=EAG/RUCAL
      TSTP=273.15+15.0
      SYM=DFLOAT(ISYM)
CFFFFFFFFFFFFFFFFFF    FREESTREAM PROPERTIES   FFFFFFFFFFFFFFFFFFFFFFFFF
      WMREF=1.0/(FO2IN/WSP(2)+(1.0-FO2IN)/WSP(LSP))/1000.0
      AMREF=5.69428D-07*(TREF**0.62245)
      ACREF=2.78763D-04*(TREF**0.796783)
      DREF=1.2D-04
      CPREF=GAMA*(GASC/WMREF)/(GAMA-1.0)
      PRREF=CPREF*AMREF/ACREF
      RREF=PREF/(GASC/WMREF)/TREF
      RSTP=PREF/(GASC/WMREF)/TSTP
      HREF=CPREF*TREF
CFFFFFFFFFFFFFFFFFFFFFF      END        FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
C**************  DEFINING NON-DIMENSIONALIZING FACTORS *****************
      ALSTR=RTOT-RTIN                                                         
      VSTR=VREF                                                         
      RSTR=RREF                                                        
      PSTR=RREF*VREF*VREF                                               
      PBACK=PREF/PSTR
      TSTR=TREF
      HSTR=VREF*VREF
      AKSTR=VREF*VREF
      EPSTR=VREF*VREF*VREF/ALSTR
      CPSTR=VREF*VREF/TREF
      ACSTR=ACREF
      TMSTR=ALSTR/VREF                                                  
      AMSTR=AMREF                                                       
      DSTR=DREF
      WMSTR=WMREF
      REI=AMREF/RREF/VREF/ALSTR
      BETA1=GASC*RSTR*TSTR/PSTR/WMSTR
      BETA2=TSTR*ACSTR/AMSTR/VSTR/VSTR
      BETA3=RSTR*DSTR/AMSTR
      BETA4=DFLOAT(IGRAV)*GFORCE*ALSTR/VSTR/VSTR
      IF(IGRAV.LE.-1) BETA4=GFORCE*ALSTR/VSTR/VSTR/DFLOAT(-IGRAV)
      HCONV=GASC/WMSTR/HSTR
C*************************      END     ********************************
C44444444444   NON-DIMENSIONALIZIG THE FLOW PARAMETERS     4444444444444
      CONVRC=418.68*1.9891D-04*DSQRT(TSTR/1000.0/WMSTR)/ACSTR
      CONVRM=0.1*2.6693D-05*DSQRT(1000.0*WMSTR*TSTR)/AMSTR
      CONVRD=0.02628*TSTR*DSQRT(TSTR)/PSTR/DSQRT(WMSTR*1000.0)/DSTR
C------------------------  Heats of Formation --------------------------
      DO 11 I=1,LSP
      WM(I)=WSP(I)/1000.0/WMSTR
      HFO(I)=HCONV*POLSP(6,I)/WM(I)
   11 CONTINUE
      HTLMT=HCONV*((POLSP(8,LSP)+(POLSP(9,LSP)/2.0+(POLSP(10,LSP)/3.0
     1 +(POLSP(11,LSP)/4.0+POLSP(12,LSP)/5.0*TPOL2)*TPOL2)*TPOL2)
     2 *TPOL2)*TPOL2+POLSP(13,LSP))/WM(LSP)
      HTLMT=5.0*HTLMT
      TKDMIN=TSTR-1.0
C--CONSTANTS---CONSTANTS---CONSTANTS---CONSTANTS---CONSTANTS---CONSTANTS
      DO I=1,LSP
		  DO J=1,LSP
              WMIK(I,J)=(WM(I)/WM(J))**0.25
			  IF(I.EQ.J) WMIK(I,J)=DSQRT(WM(I))
			  DWMIK(I,J)=2.828427125*DSQRT(1.0+WM(I)/WM(J))
			  ECDIA=0.5*(CDIA(I)+CDIA(J))
			  WMCD(I,J)=DSQRT((WM(I)+WM(J))/2.0/WM(I)/WM(J))/ECDIA/ECDIA
			  ETLJ1(I,J)=DSQRT(TLJ(I)*TLJ(J))
			  ETLJ2(I,J)=ETLJ1(I,J)*ETLJ1(I,J)
			  ETLJ3(I,J)=ETLJ2(I,J)*ETLJ1(I,J)
			  ETLJ4(I,J)=ETLJ3(I,J)*ETLJ1(I,J)
		  ENDDO
	  ENDDO
C--CONSTANTS---CONSTANTS---CONSTANTS---CONSTANTS---CONSTANTS---CONSTANTS
C-----------------------------------------------------------------------
      RENO=RREF*VREF*ALSTR/AMREF                                             
      PRNO=PRREF                                           
      WRITE(16,810) VREF,RREF,PREF,TREF,AMREF,RENO,PRNO       
  810 FORMAT(//1X,26('F'),'   Freestream Properties  ',26('F')//        
     1 2X,'Reference Velocity =',F9.4,'m/s;   Rho =',1PD11.4,
     3 ' Kg/m/m/m;'/2X,'Pressure=',1PD11.4,' N/m/m;' ,
     4 '  Temperature=',0PF7.2,' K;'/2X,'Mu =',1PD11.4,' Kg/m/s;'        
     5 '  Re No.=',D11.4,';  Pr No.=',0PF7.3///79('F'))    
      IF(ISYM.EQ.1) THEN
          IF(IGRAV.EQ.0) WRITE(16,812)  
          IF(IGRAV.GE.1) WRITE(16,813) IGRAV
          IF(IGRAV.LE.-1) WRITE(16,814) -1*IGRAV
  812     FORMAT(5X,'Cylindrical Tube-')            
  813     FORMAT(5X,'Cylindrical Vertical Pipe- Gravitational Force =',
     1         I2,'g s')            
  814     FORMAT(5X,'Cylindrical Vertical Pipe- Gravitational Force =',
     1         '1/',I2,'th g')            
          WRITE(16,815) RTIN*100.0,RTOT*100.0,ALENG*100.0
  815     FORMAT(//5X,'Inner Radius =',F6.2,'cm;  Outer Radius = ',            
     1    F6.2,'cm;'/5X,'Pipe Length=',F6.2,'cm.') 
          END IF
      IF(ISYM.EQ.0) THEN
          IF(IGRAV.EQ.0) WRITE(16,816)  
          IF(IGRAV.GE.1) WRITE(16,817) IGRAV
          IF(IGRAV.LE.-1) WRITE(16,818) -1*IGRAV
  816     FORMAT(5X,'Square Duct-')           
  817     FORMAT(5X,'Square Duct- Gravitational Force =',I2,'g s')
  818     FORMAT(5X,'Square Duct- Gravitational Force = 1/',I2,'th g')
          WRITE(16,819) RTIN*100.0,RTOT*100.0,ALENG*100.0
  819     FORMAT(//5X,'Height of the inner wall=',F6.2,'cm; Height of',            
     1    ' the Outer Wall=',F6.2,'cm;'/
     2    5X,'Channel Length=',F6.2,'cm.') 
          END IF
      WRITE(16,804) LSP,LRX/LREV,LREV,FUELTYPES(INFUEL)
      IF(IPROP.EQ.0) WRITE(16,822)
  822 FORMAT(/15X,'---CONSTANT TRANSPORT PROPERTIES---')
      IF(IPROP.EQ.1) WRITE(16,823)
  823 FORMAT(/15X,'---VARIABLE TRANSPORT PROPERTIES---')
      WRITE(16,824) LI,LJ
  824 FORMAT(//5X,'GRID SYSTEM:-  LI=',I3,';   LJ=',I3//79('G')//) 
      N=0                                                               
      ITR=0
      ICODE=0
      CALL PRESET(ICODE,RTIN,RTOT,ALENG,
     1            XBDM,YBDM,XBDP,YBDP,XFIM,YFIM,XFIP,YFIP,
     2            NIREGN,IREGN,XREGN,NJREGN,JREGN,YREGN)
      IF(ICODE.GE.1) GO TO 102
      CALL INFLOW(ICODE,INSP,
     1   NSEG,ISIDE,ITYPE,BCVAL,
     3   FO2IN)
      CALL SKIPDF
      IF(ICODE.GE.1) GO TO 102
      IF(IREAD.GE.1) THEN
         WRITE(*,826) IREAD,FILEIN
         WRITE(16,826) IREAD,FILEIN
  826    FORMAT(10X,'IREAD = ',I1,'; READING DATA FROM ',2X,A15)
         CALL READFF(32,IREAD,TTIME0)
         END IF
C--------     PREPARE THE LU MATRIX FOR PRESSURE SOLVER      -----------
      KBOUND=1
      CALL LUMATX(KBOUND)
C---------------------   FINDING THE TIME STEP  'DT'    ----------------
      DXMIN=1.0
      DYMIN=1.0
      DO 58 I=LJ+1,LE-LJ-1,LJ
      DX=X(I+LJ)-X(I)
      IF(DX.LE.DXMIN) DXMIN=DX
   58 CONTINUE
      DO 59 I=LJ+1,LJ+LJ-1
      DY=Y(I+1)-Y(I)
      IF(DY.LE.DYMIN) DYMIN=DY
   59 CONTINUE
      DTREF=DXMIN
      IF(DYMIN.LT.DXMIN) DTREF=DYMIN
      DT=CFLNO*DTREF
      DTS=DT*ALSTR/VSTR
      DTMS=DTS*1000.0
      DTPR=DT
      WRITE(*,721) LI,LJ,DXMIN*ALSTR*1000.0,DYMIN*ALSTR*1000.0,DTS
  721 FORMAT(
     1 /25X,'---GRID SYSTEM ',I3,'X',I3,'---',
     2 /18X,'---DXMIN = ',F6.3,' mm; DYMIN = ',F6.3,' mm---',
     3 /18X,'---Time Step = ',F15.12,' seconds---'/)
C----------------------------STORING THE EVOLUTION----------------------
      IF(NEVOL.GE.1) THEN
           WRITE(16,827) NEVOL
  827      FORMAT(1X,29('-'),'STORING THE EVOLUTION',29('-')//
     1        20X,'AT ',I2,' X-LOCATIONS'//
     2        4X,' N',3X,'  X(mm) ',3X,' I ',3X,'  AL1')
           DO 701 N=1,NEVOL
           IXLOC=0
           DO 702 I=1,LE,LJ
           IXLOC=I
           IF(X(I)*ALSTR.GE.XEVOL(N)) GO TO 703
  702      CONTINUE
  703      CONTINUE
           IF(IXLOC.LT.LJ) IXLOC=LJ+1
           AL1E(N)=(X(IXLOC)-XEVOL(N)/ALSTR)/(X(IXLOC)-X(IXLOC-LJ))
           IXEVOL(N)=IXLOC
           IXLOCI=(IXLOC-1)/LJ+1
           WRITE(16,828) N,XEVOL(N)*1000.0,IXLOCI,AL1E(N)
  828      FORMAT(4X,I2,3X,F8.2,3X,I3,3X,F8.5)
  701      CONTINUE
           WRITE(14,704) LJ,NEVOL
           WRITE(14,705) (XEVOL(N),N=1,NEVOL),(Y(J)*ALSTR,J=1,LJ)
           END IF
  704 FORMAT(2I6)
  705 FORMAT(1P9E14.7)
C-------------------------INTRODUCING A VORTEX--------------------------
      ININJ=06
      IF(INOISE.LT.0.AND.IVXTYP.EQ.1) THEN
           TKD=VXTYP(4)
           VXN2=1.0-VXTYP(1)-VXTYP(2)-VXTYP(3)
           VXTYP(4)=(PBACK)/(BETA1*(TKD/TSTR)*(VXTYP(1)/WM(INSP)
     1    +VXTYP(2)/WM(2)+VXTYP(3)/WM(ININJ)+VXN2/WM(LSP)))
           TKD2=TKD*TKD
           TKD3=TKD*TKD2
           TKD4=TKD*TKD3
           TKD5=TKD*TKD4
           IPOLY=1
           IF(TKD.GT.TPOL1) IPOLY=8
           CON1=GASC/WMSTR/HSTR
           HFU=(POLSP(IPOLY,INSP)*TKD+POLSP(IPOLY+1,INSP)*TKD2/2.0
     1    +POLSP(IPOLY+2,INSP)*TKD3/3.0+POLSP(IPOLY+3,INSP)*TKD4/4.0
     2    +POLSP(IPOLY+4,INSP)*TKD5/5.0
     3    +POLSP(IPOLY+5,INSP))*CON1/WM(INSP)
           HO2=(POLSP(IPOLY,2)*TKD+POLSP(IPOLY+1,2)*TKD2/2.0
     1    +POLSP(IPOLY+2,2)*TKD3/3.0+POLSP(IPOLY+3,2)*TKD4/4.0
     2    +POLSP(IPOLY+4,2)*TKD5/5.0+POLSP(IPOLY+5,2))*CON1/WM(2)
           HSU=(POLSP(IPOLY,ININJ)*TKD+POLSP(IPOLY+1,ININJ)*TKD2/2.0
     1    +POLSP(IPOLY+2,ININJ)*TKD3/3.0+POLSP(IPOLY+3,ININJ)*TKD4/4.0
     2    +POLSP(IPOLY+4,ININJ)*TKD5/5.0
     3    +POLSP(IPOLY+5,ININJ))*CON1/WM(ININJ)
           HN2=(POLSP(IPOLY,LSP)*TKD+POLSP(IPOLY+1,LSP)*TKD2/2.0
     1    +POLSP(IPOLY+2,LSP)*TKD3/3.0+POLSP(IPOLY+3,LSP)*TKD4/4.0
     2    +POLSP(IPOLY+4,LSP)*TKD5/5.+POLSP(IPOLY+5,LSP))*CON1/WM(LSP)
           VXTYP(5)=(HFU-HFO(INSP))*VXTYP(1)+(HO2-HFO(2))*VXTYP(2)
     1    +(HSU-HFO(ININJ))*VXTYP(3)
     2    +(HN2-HFO(LSP))*VXN2
          END IF
C-----------------------INTRODUCING RANDOM NOISE------------------------
      IF(INOISE.GT.0) THEN
           XNOISE=XNOISE/ALSTR
           YNOISE=YNOISE/ALSTR
           DO 720 J=1,LJ
           IF(Y(J).GT.YNOISE) GO TO 722
  720      CONTINUE
  722      CONTINUE
           JANOIS=J-INOISE
           JBNOIS=J+INOISE
           IF(JANOIS.LT.2) JANOIS=2
           IF(JBNOIS.GT.LJ-1) JBNOIS=LJ-1
           II=0
           DO 724 I=1,LE,LJ
           II=II+1
           IF(X(I).GT.XNOISE) GO TO 726
  724      CONTINUE
  726      CONTINUE
           IANOIS=II-INOISE
           IBNOIS=II+INOISE
           IF(IANOIS.LT.2) IANOIS=2
           IF(IBNOIS.GT.LI-1) IBNOIS=LI-1
           END IF
C-----------------PREPARING FOR THE PARTICLE TRACKING-------------------
      IF(NOPR.GT.0) THEN
           IF(ITINJ.LE.0) ITINJ=1
           IF(IEINJ.LE.1) IEINJ=ITEND
           NOINJ=0
           PDIA=PDIA*1.0D-06/ALSTR
           PDEN=PDEN/RSTR
           PTHR=PTHR/ACSTR
           DO 732 N=1,NOPR
           XPR(N)=XPR(N)/ALSTR
           YPR(N)=YPR(N)/ALSTR
           DO I=1,LNMX
           PATHX(I,N)=0.0
           PATHY(I,N)=0.0
           ENDDO
  732      CONTINUE
           DO 733 I=1,101
           AA=DFLOAT(I-1)*0.01
           DO 734 J=1,LI
           JJ=(J-1)*LJ+1
           XA=X(JJ)/X(LE)
           IF(XA.GT.AA) GO TO 735
  734      CONTINUE
  735      CONTINUE
           IPICK(I)=J-1
           DO 736 J=1,LJ
           YA=Y(J)
           IF(YA.GT.AA) GO TO 737
  736      CONTINUE
  737      CONTINUE
           JPICK(I)=J-1
  733      CONTINUE
           END IF
C-------------------- Particles along iso-temperature line -------------
      DO 8010 N=1,NOPR
      XX=XPR(N)+DFLOAT(N-1)*0.0004/ALSTR
      DO 8002 I=1,LI
      II=(I-1)*LJ+1
      IF(X(II).GE.XX) GO TO 8003
 8002 CONTINUE
 8003 CONTINUE
      XPR(N)=X(II)
      TTMAX=0.0
      DO 8001 J=1,LJ
      II=(I-1)*LJ+J
      IF(TK(II).GE.TTMAX) TTMAX=TK(II)
 8001 CONTINUE
      IF(TTMAX*TSTR.LT.1200.00) GO TO 8010
      TK1=1200.0/TSTR
          DO 8006 J=LJ-1,1,-1
          II=(I-1)*LJ+J
          IF(TK(II).GE.TK1) THEN
               YY=Y(II+1)-(Y(II+1)-Y(II))*(TK1-TK(II+1))
     1                   /(TK(II)-TK(II+1))
               GO TO 8005
               END IF
 8006     CONTINUE
          YY=YPR(N)
 8005 CONTINUE
      YPR(N)=YY
      WRITE(*,8009) N,XPR(N)*ALSTR,YPR(N)*ALSTR,TK1*TSTR,TK(II)*TSTR
 8009 FORMAT(1X,I3,2X,2F12.8,2X,2F8.2)
 8010 CONTINUE
C-----------------------INTRODUCING RADIAL EJECTION---------------------
      IF(NEJEC1.GT.0) THEN
           XEJEC1=XEJEC1/ALSTR
           XEJEC2=XEJEC2/ALSTR
           YEJEC=YEJEC/ALSTR
           DO 742 J=1,LJ
           IF(Y(J).LE.YEJEC) JEJEC=J
  742      CONTINUE
           IF(JEJEC.LT.2) JEJEC=2
           II=0
           DO 744 I=1,LE,LJ
           II=II+1
           IF(X(I).LE.XEJEC1) IEJEC1=II
           IF(X(I).LE.XEJEC2) IEJEC2=II
  744      CONTINUE
           IF(IEJEC1.LT.2) IEJEC1=2
           IF(IEJEC2.GT.LI-1) IEJEC2=LI-1
           IA1=(IEJEC1-1)*LJ+1
           IA2=(IEJEC2-1)*LJ+1
           DXEJ=0.5*(X(IA2+LJ)+X(IA2)-X(IA1)-X(IA1-LJ))
           VEJEC=AMEJEC/(6283100.0*Y(JEJEC)*DXEJ*ALSTR*ALSTR)
           NEJECD=INT(DTEJEC/DTMS)+1
           NEJEC2=NEJEC1+NEJECD
           IBEJEC=NEJEC1
           NCYCLE=0
           IF(FEJEC.GT.0.0) NCYCLE=INT(1.0/FEJEC/DTS)
           IF(MTEJEC.EQ.1) THEN
               EJFS1=FSP(1,INSP)
               EJFS2=FSP(1,2)
               EJFS3=0.0
               EJTK=TK(1)
               END IF
           IF(MTEJEC.EQ.2) THEN
               EJFS1=FSP(LJ,INSP)
               EJFS2=FSP(LJ,2)
               EJFS3=0.0
               EJTK=TK(LJ)
               END IF
           IF(MTEJEC.EQ.3) THEN
               EJFS1=0.0
               EJFS2=0.0
               EJFS3=0.0
               EJTK=1.0
               END IF
           IF(MTEJEC.EQ.4) THEN
               EJFS1=0.0
               EJFS2=0.0
               EJFS3=1.0
               EJTK=1.0
               END IF
C          IF(MTEJEC.EQ.5) THEN-----LOCAL FLUID
           EJN2=(1.0-EJFS1-EJFS2-EJFS3)
           EJRHO=(PBACK)/(BETA1*EJTK*(EJFS1/WM(INSP)
     1    +EJFS2/WM(2)+EJFS3/WM(ININJ)+EJN2/WM(LSP)))
           TKD=EJTK*TSTR
           TKD2=TKD*TKD
           TKD3=TKD*TKD2
           TKD4=TKD*TKD3
           TKD5=TKD*TKD4
           IPOLY=1
           IF(TKD.GT.TPOL1) IPOLY=8
           CON1=GASC/WMSTR/HSTR
           HFU=(POLSP(IPOLY,INSP)*TKD+POLSP(IPOLY+1,INSP)*TKD2/2.0
     1    +POLSP(IPOLY+2,INSP)*TKD3/3.0+POLSP(IPOLY+3,INSP)*TKD4/4.0
     2    +POLSP(IPOLY+4,INSP)*TKD5/5.0
     3    +POLSP(IPOLY+5,INSP))*CON1/WM(INSP)
           HO2=(POLSP(IPOLY,2)*TKD+POLSP(IPOLY+1,2)*TKD2/2.0
     1    +POLSP(IPOLY+2,2)*TKD3/3.0+POLSP(IPOLY+3,2)*TKD4/4.0
     2    +POLSP(IPOLY+4,2)*TKD5/5.0+POLSP(IPOLY+5,2))*CON1/WM(2)
           HSU=(POLSP(IPOLY,ININJ)*TKD+POLSP(IPOLY+1,ININJ)*TKD2/2.0
     1    +POLSP(IPOLY+2,ININJ)*TKD3/3.0+POLSP(IPOLY+3,ININJ)*TKD4/4.0
     2    +POLSP(IPOLY+4,ININJ)*TKD5/5.0
     3    +POLSP(IPOLY+5,ININJ))*CON1/WM(ININJ)
           HN2=(POLSP(IPOLY,LSP)*TKD+POLSP(IPOLY+1,LSP)*TKD2/2.0
     1    +POLSP(IPOLY+2,LSP)*TKD3/3.0+POLSP(IPOLY+3,LSP)*TKD4/4.0
     2    +POLSP(IPOLY+4,LSP)*TKD5/5.+POLSP(IPOLY+5,LSP))*CON1/WM(LSP)
           EJHT=(HFU-HFO(INSP))*EJFS1+(HO2-HFO(2))*EJFS2
     1    +(HSU-HFO(ININJ))*EJFS3
     1    +(HN2-HFO(LSP))*EJN2+0.5*VEJEC*VEJEC/VSTR/VSTR
           WRITE(*,745) NEJEC1,NEJEC2,VEJEC,EJTYPE(MTEJEC),NCYCLE,FEJEC
  745      FORMAT(2X,'RADIAL MASS EJECTION;  NA = ',I4,'  NB = ',I4/10X,
     1         '   VEJEC = ',F8.5,'m/s  ','   TYPE = ',A5/10X,
     2         '   AT A RATE',I4,' (',F6.1,' Hz.)')
           END IF
C---------------------   STORE THE DRIVING HISTORY     -----------------
      IF(NDRV.GE.1) THEN
                    WRITE(72,746) (IDRV(I),I=1,NDRV)
  746               FORMAT('* TIME (ms)',10(',J=',I3))
                    END IF
C-----------------------------------------------------------------------
      DO 52 I=1,LE
      RHONP(I)=RHO(I)
   52 CONTINUE
C-------------------       START ADVANCING THE SOLUTION      -----------
      TTIME=TTIME0
      CALL BODYBC(INSP,PBACK)
      IF(NTPRNT.GE.1) CALL PRNTFF(INSP,ITRANS)
	  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C-------------------START TIME STEP------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 100 ITR=1,ITEND
	  
C     CLOCK=SECOND()-CLOCK1
      CLOCK=0.0
      IF(CLOCK.GE.FSECS) THEN
                         ITEND=0
                         GO TO 101
                         END IF
C------------------  INITIALIZE THE NEW VARIABLES ----------------------
      DO 109 I=1,LE
      RHONP(I)=RHO(I)
  109 CONTINUE
C11111111111111111111  TRANSPORT PROPERTIES   11111111111111111111111111
      DO 110 I=1,LE
      TMASS=0.0
      DO 111 ISP=1,LSP-1
      IF(FSP(I,ISP).LT.0.0) FSP(I,ISP)=0.0
      TMASS=TMASS+FSP(I,ISP)
  111 CONTINUE
      FSP(I,LSP)=1.0-TMASS
  110 CONTINUE
      IF(IPROP.EQ.1) GO TO 112
      DO 113 I=1,LE
      VMU(I)=AMREF/AMSTR
      VTC(I)=ACREF/ACSTR
      DO 114 NP=1,LSP
      VDFSP(I,NP)=DREF/DSTR
  114 CONTINUE
  113 CONTINUE
      GO TO 129
  112 CONTINUE
C-----------------------------------------------------------------------
C-----------------LIMIT MASS FRACTION IN THE SUB-MIXTURE----------------
      FLMT=1.0D-03
C-----------------------------------------------------------------------
      DO 128 I=1,LI 
      DO J=1,LJ 
      II=(I-1)*LJ+J
      IA=0
      TMASS=0.0
      TMOLE=0.0
      DO 130 ISP=1,LSP-1
      ISPQ(ISP)=0
      IF(FSP(II,ISP).LT.FLMT) GO TO 130
      IA=IA+1
      ISPQ(ISP)=IA
      ISPQR(IA)=ISP
      WMQ(IA)=WM(ISP)
      CDIAQ(IA)=CDIA(ISP)
      DO 131 IB=1,IA
      ISPB=ISPQR(IB)
      WMIKQ(IA,IB)=WMIK(ISP,ISPB)
      WMCDQ(IA,IB)=WMCD(ISP,ISPB)
      DWMIKQ(IA,IB)=DWMIK(ISP,ISPB)
      ETLJ1Q(IA,IB)=ETLJ1(ISP,ISPB)
      ETLJ2Q(IA,IB)=ETLJ2(ISP,ISPB)
      ETLJ3Q(IA,IB)=ETLJ3(ISP,ISPB)
      ETLJ4Q(IA,IB)=ETLJ4(ISP,ISPB)
C
      WMIKQ(IB,IA)=WMIK(ISPB,ISP)
      WMCDQ(IB,IA)=WMCD(ISPB,ISP)
      DWMIKQ(IB,IA)=DWMIK(ISPB,ISP)
      ETLJ1Q(IB,IA)=ETLJ1(ISPB,ISP)
      ETLJ2Q(IB,IA)=ETLJ2(ISPB,ISP)
      ETLJ3Q(IB,IA)=ETLJ3(ISPB,ISP)
      ETLJ4Q(IB,IA)=ETLJ4(ISPB,ISP)
  131 CONTINUE
      DO 132 N=1,12
      POLQ(N,IA)=POLSP(N,ISP)
      CIQ(N,IA)=CISP(N,ISP)
  132 CONTINUE
      POLQ(13,IA)=POLSP(13,ISP)
      POLQ(14,IA)=POLSP(14,ISP)
      TMASS=TMASS+FSP(II,ISP)
      XSP(IA)=FSP(II,ISP)/WM(ISP)
      TMOLE=TMOLE+XSP(IA)
  130 CONTINUE
      LSPT=IA+1
      IF(LSPT.GE.LSP) LSPT=LSP
      ISPQR(LSPT)=LSP
      WMQ(LSPT)=WM(LSP)
      CDIAQ(LSPT)=CDIA(LSP)
      DO 133 IB=1,LSPT
      ISPB=ISPQR(IB)
      WMIKQ(LSPT,IB)=WMIK(LSP,ISPB)
      WMCDQ(LSPT,IB)=WMCD(LSP,ISPB)
      DWMIKQ(LSPT,IB)=DWMIK(LSP,ISPB)
      ETLJ1Q(LSPT,IB)=ETLJ1(LSP,ISPB)
      ETLJ2Q(LSPT,IB)=ETLJ2(LSP,ISPB)
      ETLJ3Q(LSPT,IB)=ETLJ3(LSP,ISPB)
      ETLJ4Q(LSPT,IB)=ETLJ4(LSP,ISPB)
C
      WMIKQ(IB,LSPT)=WMIK(ISPB,LSP)
      WMCDQ(IB,LSPT)=WMCD(ISPB,LSP)
      DWMIKQ(IB,LSPT)=DWMIK(ISPB,LSP)
      ETLJ1Q(IB,LSPT)=ETLJ1(ISPB,LSP)
      ETLJ2Q(IB,LSPT)=ETLJ2(ISPB,LSP)
      ETLJ3Q(IB,LSPT)=ETLJ3(ISPB,LSP)
      ETLJ4Q(IB,LSPT)=ETLJ4(ISPB,LSP)
  133 CONTINUE
      DO 134 N=1,12
      POLQ(N,LSPT)=POLSP(N,LSP)
      CIQ(N,LSPT)=CISP(N,LSP)
  134 CONTINUE
      POLQ(13,LSPT)=POLSP(13,LSP)
      POLQ(14,LSPT)=POLSP(14,LSP)
      XSP(LSPT)=(1.0-TMASS)/WMQ(LSPT)
      TMOLE=TMOLE+XSP(LSPT)
C
      TKRT=DSQRT(TK(II))
      TKD1=TK(II)*TSTR
      TKD2=TKD1*TKD1
      TKD3=TKD2*TKD1
      TKD4=TKD3*TKD1
      TKD5=TKD4*TKD1
      IPL1=1
      IPL2=1
      IF(TKD1.GE.TPOL1) THEN
                        IPL1=8
                        IPL2=7
                        END IF
      TKBP=TK(II)*TKRT/PBACK
C------Calculation of individual viscosity & thermal conductivity
      DO 122 ISP=1,LSPT
      EUKEN=0.115+0.354*(POLQ(IPL1,ISP)
     1   +POLQ(IPL1+1,ISP)*TKD1+POLQ(IPL1+2,ISP)*TKD2
     2   +POLQ(IPL1+3,ISP)*TKD3+POLQ(IPL1+4,ISP)*TKD4)
C     IF(ISP.EQ.25.OR.ISP.EQ.27) EUKEN=1.0
      SIOME=CDIAQ(ISP)*CDIAQ(ISP)*(CIQ(IPL2,ISP)
     1      +CIQ(IPL2+1,ISP)*TKD1+CIQ(IPL2+2,ISP)*TKD2
     2      +CIQ(IPL2+3,ISP)*TKD3+CIQ(IPL2+4,ISP)*TKD4
     3      +CIQ(IPL2+5,ISP)*TKD5)
      ATCSP(ISP)=CONVRC*EUKEN*TKRT/WMIKQ(ISP,ISP)/SIOME
      AMUSP(ISP)=CONVRM*TKRT*WMIKQ(ISP,ISP)/SIOME
C------Calculation of binary diffusion coefficients (m - n, m=n)
      DO 123 N=ISP+1,LSPT
      TKDS=TKD1/ETLJ1Q(ISP,N)
      K=1
      IF(TKDS.GE.2.0) K=2
      IF(TKDS.GT.5.0) K=3
      IF(TKDS.GT.70.0) K=4
      CI11T=CI11(K,1)+CI11(K,2)*TKDS+CI11(K,3)*TKD2/ETLJ2Q(ISP,N)
     1     +CI11(K,4)*TKD3/ETLJ3Q(ISP,N)+CI11(K,5)*TKD4/ETLJ4Q(ISP,N)
      DIFSP(ISP,N)=WMCDQ(ISP,N)*TKBP/CI11T
      DIFSP(N,ISP)=WMCDQ(N,ISP)*TKBP/CI11T
  123 CONTINUE
  122 CONTINUE
C-------------Calculating the mixture properties
      DO 124 ISP=1,LSPT
      PHIK(ISP)=0.0
      PHDIK(ISP)=0.0
      EFDIK(ISP)=0.0
      DO 125 JSP=1,LSPT
      IF(JSP.EQ.ISP) GO TO 125
      PHIK(ISP)=PHIK(ISP)+XSP(JSP)*((1.0
     1 +DSQRT(ATCSP(ISP)/ATCSP(JSP))*WMIKQ(ISP,JSP))**2)/DWMIKQ(ISP,JSP)
      PHDIK(ISP)=PHDIK(ISP)+XSP(JSP)*((1.0
     1 +DSQRT(AMUSP(ISP)/AMUSP(JSP))*WMIKQ(ISP,JSP))**2)/DWMIKQ(ISP,JSP)
      EFDIK(ISP)=EFDIK(ISP)+XSP(JSP)/DIFSP(ISP,JSP)
  125 CONTINUE
  124 CONTINUE
      VTC(II)=0.0
      VMU(II)=0.0
      DO 126 ISP=1,LSPT
      VTC(II)=VTC(II)+ATCSP(ISP)*XSP(ISP)/(XSP(ISP)+1.065*PHIK(ISP))
      VMU(II)=VMU(II)+AMUSP(ISP)*XSP(ISP)/(XSP(ISP)+PHDIK(ISP))
      AA=TMOLE-XSP(ISP)
      IF(AA.LE.1.0D-10) THEN
            VDFQ(ISP)=CONVRD*DIFSP(ISP,LSPT)
            ELSE
            VDFQ(ISP)=CONVRD*(TMOLE-XSP(ISP))/EFDIK(ISP)
            END IF
  126 CONTINUE
      DO 127 ISP=1,LSP
      IS=ISPQ(ISP)
      IF(IS.EQ.0) THEN
                  EFD=0.0
                  DO 135 IB=1,LSPT
                  N=ISPQR(IB)
                  TKDS=TKD1/ETLJ1(ISP,N)
                  K=1
                  IF(TKDS.GE.2.0) K=2
                  IF(TKDS.GT.5.0) K=3
                  IF(TKDS.GT.70.0) K=4
                  CI11T=CI11(K,1)+CI11(K,2)*TKDS
     1                 +CI11(K,3)*TKD2/ETLJ2(ISP,N)
     2                 +CI11(K,4)*TKD3/ETLJ3(ISP,N)
     3                 +CI11(K,5)*TKD4/ETLJ4(ISP,N)
                  EFD=EFD+XSP(IB)/(WMCD(ISP,N)*TKBP/CI11T)
  135             CONTINUE
                  VDFSP(II,ISP)=CONVRD*TMOLE/EFD
                  ELSE
                  VDFSP(II,ISP)=VDFQ(IS)
                  END IF
  127 CONTINUE
      STDF(II)=0.2*VDFSP(II,LSP)
      ENDDO
  128 CONTINUE
  129 CONTINUE
      IF(IFLOW.EQ.2) THEN
          DO 139 I=1,LE
          TMU(I)=CMKE*RHO(I)*AK(I)*AK(I)/EPS(I)/REI
  139     CONTINUE
          END IF
C11111111111111111111111111111111111111111111111111111111111111111111111
      RESID=0.0
      RESDM=0.0
      RESDP=0.0
      RESDK=0.0
      RESDE=0.0
      RESDH=0.0
      RESDC=0.0
      RESDW=0.0
      IRLXSP=0
      IRLXH=0
      IRLXU=0
      IRLXV=0
      IRLXW=0
      IRLXKE=0
      IF(ICHEM.EQ.0.AND.ITHRM.LE.0) GO TO 501
C-----------------             ADD IGNITION                 ------------
      IF(ITR.EQ.IGNIT) CALL FIRE(INSP,IGKEEP)
      IF(IGNIT.GE.1.AND.ITR.GE.IGNIT.AND.ITR.LE.(IGNIT+IGKEEP)) THEN
          DO 508 I=1,LE
          IF(ISKIP(I).EQ.-9) THEN
                             TK(I)=1700.0/TSTR
                             FSP(I, 3)=0.0001
                             FSP(I, 4)=0.0001
                             FSP(I, 5)=0.001
                             FSP(I, 6)=0.001
                             FSP(I, 9)=0.001
                             FSP(I,10)=0.001
                             END IF
  508     CONTINUE
          END IF
C-----------------   SOLVING SPECIES & ENTHALPY EQUATIONS   ------------
      IRLXSP=ISORSP
      IRLXH=ISORH
      RELX=RELXSP
      ISCHM1=1
      ISCHM2=1
      IF(SCHMSP.EQ.SCHM1) ISCHM1=2
      IF(SCHMH.EQ.SCHM1) ISCHM2=2
	  
      CALL SPGSOLV(IRLXSP,RELXSP,TOLRSP,RESDC,
     1  IRLXH,RELXH,TOLRH,RESDH,SIGSP,SIGH)
      IF(ITHRM.LE.0) GO TO 501
C--------------         CALCULATING THE DENSITY         ----------------
      DO 510 I=1,LE
      RBAR=0.0
      DO 511 ISP=1,LSP
      RBAR=RBAR+FSP(I,ISP)/WM(ISP)
  511 CONTINUE
      RHONP(I)=(PBACK)/(BETA1*TK(I)*RBAR)
  510 CONTINUE
      DO 532 I=LE-LJ+1,LE
      RHONP(I)=2.0*RHONP(I-LJ)-RHONP(I-LJ2)
  532 CONTINUE
      DO 534 I=1,LE,LJ
      RHONP(I)=RHONP(I+1)
  534 CONTINUE
      IF(ISTDY.EQ.1) THEN
         DO 535 I=1,LE
         RHO(I)=RHONP(I)
  535    CONTINUE
         END IF
  501 CONTINUE
      CALL BODYBC(INSP,PBACK)
      IF(IPASS.GE.1) GO TO 269
      IF(IFLOW.EQ.0) GO TO 269 
C-----------         SOLVING U-MOMENTUM EQUATION    --------------------
      RINF=RREF/RSTR
      IEQN=1
      DO 212 I=1,LE
      FP(I)=U(I)
  212 CONTINUE
      SIGMA=1.0
      CALL ACONSU(SIGMA)
      IRLXU=ISORU
      IF(SCHMU.EQ.SCHM1) CALL ADISOLV(IEQN,IRLXU,RELXU,TOLRU)
      IF(SCHMU.EQ.SCHM2) CALL PTRELX(IEQN,IRLXU,RELXU,TOLRU)
      DO 214 I=1,LE
      UTMP(I)=FP(I)
  214 CONTINUE
C-----------         SOLVING V-MOMENTUM EQUATION    --------------------
      IEQN=2
      DO 220 I=1,LE
      FP(I)=V(I)
  220 CONTINUE
      SIGMA=1.0
      CALL ACONSV(SIGMA)
      IRLXV=ISORV
      IF(SCHMV.EQ.SCHM1) CALL ADISOLV(IEQN,IRLXV,RELXV,TOLRV)
      IF(SCHMV.EQ.SCHM2) CALL PTRELX(IEQN,IRLXV,RELXV,TOLRV)
      DO 224 I=1,LE
      V(I)=FP(I)
      UALT=UTMP(I)
      UTMP(I)=U(I)
      U(I)=UALT
  224 CONTINUE
C------------------------------GENERATE NOISE---------------------------
      IF(INOISE.GT.0) THEN
         CALL PERTRB(IANOIS,IBNOIS,JANOIS,JBNOIS)
         END IF
C------------------------INPUT RADIAL EJECTION--------------------------
      IF(NEJEC1.GT.0) THEN
      IF(ITR.GE.NEJEC1.AND.ITR.LE.NEJEC2) THEN
         SLEJEC=0.0
         ZMOVE=DFLOAT(ITR-NEJEC1)*DTS*SLEJEC/ALSTR
         IA=(IEJEC1-1)*LJ+1
         DO 240 I=1,LE-IA,LJ
         Z1=X(IA+I-1)-X(IA)
         IF(Z1.GE.ZMOVE) GO TO 241
  240    CONTINUE
         I=1
  241    CONTINUE
         IEMOVE=I-1
         DO 242 I=IEJEC1,IEJEC2
         IN=IEMOVE+(I-1)*LJ+JEJEC
         DO 244 II=IN,IN+1
         V(II)=VEJEC/VSTR
         IF(MTEJEC.LE.4) THEN
            DO 245 ISP=1,LSP-1
            FSP(II,ISP)=0.0
  245       CONTINUE
            FSP(II,INSP)=EJFS1
            FSP(II,2)=EJFS2
            FSP(II,ININJ)=EJFS3
            FSP(II,LSP)=EJN2
            TK(II)=EJTK
            HT(II)=EJHT
            RHONP(II)=EJRHO
            RHO(II)=EJRHO
            END IF
  244    CONTINUE
  242    CONTINUE
         END IF
      ITEJEC=ITR-IBEJEC
      IF(ITEJEC.LT.1) ITEJEC=1
      IF(NCYCLE.GT.0.AND.MOD(ITEJEC,NCYCLE).EQ.0) THEN
         NEJEC1=ITR
         NEJEC2=NEJEC1+NEJECD
         END IF
         END IF
      CALL BODYBC(INSP,PBACK)
C-----------------     SOLVING PRESSURE POISSON EQUATION      ----------
      CALL PRDIRC(KBOUND,RESDM)
      CALL BODYBC(INSP,PBACK)
      DO 260 I=1,LE
      RESID=RESID+DABS(UTMP(I)-U(I))
  260 CONTINUE
      RESID=RESID/DFLOAT(LE)
C--------------         SOLVING SWIRL EQUATION    ----------------------
      IF(ISWRL.EQ.0) GO TO 239 
      IEQN=3
      DO 230 I=1,LE
      FP(I)=W(I)
  230 CONTINUE
      SIGMA=SIGW
      IRLXW=ISORW
C     CALL SWLSOLV(IRLXW,RELXW,TOLRW,SIGMA,RESDW)
      CALL ACONS3(SIGMA)
      IF(SCHMW.EQ.SCHM1) 
     1     CALL ADISOLV(IEQN,IRLXW,RELXW,TOLRW)
      IF(SCHMW.EQ.SCHM2) 
     1     CALL PTRELX(IEQN,IRLXW,RELXW,TOLRW)
      DO 232 I=1,LE
      RESDW=RESDW+ABS(W(I)-FP(I))
      W(I)=FP(I)
  232 CONTINUE
      RESDW=RESDW/DFLOAT(LE)
  239 CONTINUE
      IF(IFLOW.LE.1) GO TO 259 
C------------     SOLVING  K-EPS - ENERGY  EQUATIONS    ----------------     
      IRLXKE=ISORKE
      ISCHEM=1
      IF(SCHMKE.EQ.SCHM2) ISCHEM=2
      CALL KESOLV(IRLXKE,RELXKE,TOLRKE,SIGK,SIGE,RESDK,RESDE)
  259 CONTINUE
C-----------------------------------------------------------------------
  269 CONTINUE
      DO 267 I=1,LE
      RHO(I)=RHONP(I)
  267 CONTINUE
      TTIME=TTIME+DT*ALSTR/VSTR
      IA=ITR+IPRES*20-1
      IF(IA.GT.ITEND) IA=ITEND
      IF(MOD(ITR-1,IPRES*20).EQ.0) WRITE(*,754) ITR,IA,DTMS
  754 FORMAT(1X,13('-'),' ITERATIONS FROM ',I5,' TO ',I5,
     1       2X,'(DT = ',F8.6,' ms) ',13('-')/
     2       2X,'U-RESID            M-RESID   HT-RESID     CHEM-RES',
     3       6X,'K & EPS-RESID')
      IF(MOD(ITR,IPRES).EQ.0) WRITE(*,756) RESID,IRLXU,IRLXV,IRLXW,
     1        RESDM,RESDH,IRLXH,RESDC,IRLXSP,RESDK,RESDE,IRLXKE
      IF(ITR.EQ.1) WRITE(16,754) ITR,ITEND,DTS
      IF(MOD(ITR,100).EQ.0) WRITE(16,756) RESID,IRLXU,IRLXV,IRLXW,
     1        RESDM,RESDH,IRLXH,RESDC,IRLXSP,RESDK,RESDE,IRLXKE
  756 FORMAT(1X,1PE8.1,'(',I2,',',I2,',',I2,') ',1PE8.1,1X,
     1  1PE8.1,'(',I2,') ',1PE8.1,'(',I2,') ',
     1  1PE8.1,' &',1PE8.1,'(',I2,')')
      IF(MOD(ITR,ITPRNT).EQ.0) THEN
      WRITE(16,758) ITR,DTMS,TTIME,RESID,IRLXU,IRLXV,RESDW,IRLXW,
     1        RESDM,RESDH,IRLXH,RESDC,IRLXSP,RESDK,RESDE,IRLXKE
  758 FORMAT(1X,'N=',I4,', DT=',F8.6,'ms, TIME=',F9.1,'s'/
     1  'U-RES=',1PE9.2,'(',I2,',',I2,') ; ',
     2  'W-RES=',1PE9.2,'(',I2,') ','; M-RES=',1PE9.2/
     2  'H-RES=',1PE9.2,'(',I2,') ','; CHEM-RES=',1PE9.2,'(',I2,')'/
     3  'K & EPS-RES=',1PE9.2,' &',1PE9.2,'(',I2,')')
        CALL PRNTFF(INSP,ITRANS)
        END IF
      IF(ITR.GE.ISTB.AND.MOD(ITR-ISTB,ISTORE).EQ.0) THEN
        CALL PRNTDK(12,INFUEL)
        IF(NSTORE.EQ.1) REWIND(12)
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
        CALL STDATA(TTIME)
        IF(NSTORE.EQ.1) REWIND(13)
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
        END IF
C------------------------STORE THE EVOLUITION---------------------------
      IF(NEVOL.GE.1.AND.ITR.GE.IBEVOL) THEN
         IF(MOD(ITR,ISEVOL).EQ.0) THEN
           DO 707 N=1,NEVOL
           IXLOC=IXEVOL(N)
           DO 706 J=IXLOC,IXLOC+LJ-1
           FEVOL(J-IXLOC+1,1,N)=(TK(J-LJ)*AL1E(N)
     1        +TK(J)*(1.0-AL1E(N)))*TSTR
           FEVOL(J-IXLOC+1,2,N)=(U(J-LJ)*AL1E(N)
     1        +U(J)*(1.0-AL1E(N)))*VSTR
  706      CONTINUE
  707      CONTINUE
           WRITE(14,705) TTIME,((FEVOL(J,1,N),J=1,LJ),
     1                          (FEVOL(J,2,N),J=1,LJ),N=1,NEVOL)
           END IF
         END IF
C--------------------------PARTICLE TRACKING----------------------------
      IF(NOPR.GT.0.AND.ITR.GT.IBINJ) THEN
             NEWINJ=0
             IF(ITR.LE.IEINJ.AND.MOD(ITR-IBINJ,ITINJ).EQ.0) THEN
                NEWINJ=1
C-------------------- Particles along iso-temperature line -------------
                DO 9010 N=1,NOPR
                XX=XPR(N)
                DO 9002 I=1,LI
                II=(I-1)*LJ+1
                IF(X(II).GE.XX) GO TO 9003
 9002           CONTINUE
 9003           CONTINUE
                TK1=1200.0/TSTR
                DO 9006 J=LJ-1,1,-1
                II=(I-1)*LJ+J
                IF(TK(II).GE.TK1) THEN
                    YY=Y(II+1)-(Y(II+1)-Y(II))*(TK1-TK(II+1))
     1                   /(TK(II)-TK(II+1))
                    GO TO 9005
                    END IF
 9006           CONTINUE
                YY=0.0
 9005           CONTINUE
                YPR(N)=YY
 9010           CONTINUE
                END IF
C-----------------------------------------------------------------------
             SLPART=0.0
             PRMOVE=DFLOAT(ITR-IBINJ)*DTS*SLPART/ALSTR
             DTPR=0.02*DT
             DO 9100 IDP=1,50
             IF(IDP.GE.2) NEWINJ=0
             CALL PRTRACK(NOPR,NOINJ,NEWINJ,DTPR,PRMOVE,
     1                    PDIA,PDEN,PTHR,PVEL)
 9100 CONTINUE
             IF(ITR.GE.ISTB.AND.MOD(ITR-ISTB,ISTORE).EQ.0) THEN
              WRITE(20,738) LI,LJ,NOPR,NOINJ,ALSTR,ITR
              WRITE(20,739) (X(I),I=1,LE,LJ),(Y(J),J=1,LJ)
C             WRITE(20,739) ((PATHX(J,I),PATHY(J,I),I=1,NOPR),J=1,NOINJ)
              WRITE(20,739) ((PATHX(J,I),PATHY(J,I),
C    1                        GUVEL(J,I),GVVEL(J,I),PUVEL(J,I),
     2                        GUVEL(J,I),I=1,NOPR),J=1,NOINJ)
              IF(NSTORE.EQ.1) REWIND(20)
              END IF
             END IF
  738        FORMAT(4I6,F12.7,I6)
  739        FORMAT(1P10E12.5)
C---------------------------ANIMATION PLOTS-----------------------------
      IF((IBANM.GE.1).AND.(ITR.GE.IBANM).AND.
     1    (MOD(ITR-IBANM,ISANM).EQ.0)) 
     2    CALL PLOTS(NOPR,NOINJ,IPANM,
     3               X1ANM,Y1ANM,X2ANM,Y2ANM,NFSURF,KORNT)
C---------------------   STORE THE DRIVING HISTORY     -----------------
      IF((NDRV.GE.1).AND.(ITR.GE.IBDRV)) THEN
             DO 748 I=1,NDRV
             II=(IDRV(I)-1)*LJ+3
             ADRV(I)=U(II)*VSTR
  748        CONTINUE
             WRITE(72,749) TTIME*1000.0,(ADRV(I),I=1,NDRV)
  749        FORMAT(F10.6,10(',',F10.5))
             END IF
C----------------------------FLOW AVERAGE-------------------------------
      IF((ITR.GE.NBAVE).AND.(ITR.LE.NEAVE)) THEN
             IAVE=IAVE+1
             DO 763 I=1,LE
             UAVE(I)=UAVE(I)+U(I)
             VAVE(I)=VAVE(I)+V(I)
             TAVE(I)=TAVE(I)+TK(I)
             DO 764 ISP=1,LSP
             FSPAV(I,ISP)=FSPAV(I,ISP)+FSP(I,ISP)
  764        CONTINUE
  763        CONTINUE
             END IF
			 
C----------------------------END OF A TIME-STEP-------------------------
  100 CONTINUE
  101 CONTINUE
      ITR=ITR-1
      IF(NEVOL.GE.1) THEN
             TTIME=-1.0
             WRITE(14,705) TTIME,((FEVOL(J,1,N),J=1,LJ),
     1         (FEVOL(J,2,N),J=1,LJ),(FEVOL(J,3,N),J=1,LJ),N=1,NEVOL),
     2         (FEVOL(J,4,N),J=1,LJ)
             CLOSE(14)
             END IF
      IF(ITEND.EQ.0.AND.ITPRNT.GE.1) CALL PRNTFF(INSP,ITRANS)
      IF(ITEND.EQ.0.AND.ISTORE.GE.1) CALL PRNTDK(12,INFUEL)
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
      IF(ITEND.EQ.0.AND.ISTORE.GE.1) CALL STDATA(TTIME)
      IF(ISTORE.GE.1) CLOSE(13)
	  WRITE(*, 899)
  899 FORMAT('END')
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
C--------------------------WRITING AVERAGE DATA-------------------------
      IF(NBAVE.GE.1) THEN
             AAVE=FLOAT(IAVE)
             ITR=IAVE
             DO 766 I=1,LE
             U(I)=UAVE(I)/AAVE
             V(I)=VAVE(I)/AAVE
             TK(I)=TAVE(I)/AAVE
             DO ISP=1,LSP
             FSP(I,ISP)=FSPAV(I,ISP)/AAVE
             ENDDO
  766        CONTINUE
             CALL PRNTDK(13,INFUEL)
             CLOSE(13)
             END IF
C---CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE---
                      CLOSE(6)
      IF(ISTORE.GE.1) CLOSE(12)
      IF(NOINJ.GE.1)  CLOSE(20)
      IF(IREAD.GE.1)  CLOSE(32)
      IF(IBANM.GE.1)  CLOSE(82)
C---CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE-CLOSE---
  102 CONTINUE
      WRITE(16,852) ICODE
  852 FORMAT(5X,'ENDED WITH ERROR CODE',2X,I3)
      CLOSE(16)
      STOP
      END PROGRAM UNICORN



C***********************************************************************



      INCLUDE 'lib\acons3.f'
      INCLUDE 'lib\aconsu.f'
      INCLUDE 'lib\aconsv.f'
      INCLUDE 'lib\adisolv.f'
      INCLUDE 'lib\arrow.f'
      INCLUDE 'lib\bodybc.f'
      INCLUDE 'lib\ckrates.f'
      INCLUDE 'lib\dot.f'
      INCLUDE 'lib\dotpr.f'
      INCLUDE 'lib\drawby.f'
      INCLUDE 'lib\fire.f'
      INCLUDE 'lib\header.f'
      INCLUDE 'lib\inflow.f'
      INCLUDE 'lib\kesolv.f'
      INCLUDE 'lib\line.f'
      INCLUDE 'lib\lumatx.f'
      INCLUDE 'lib\num2bit4.f'
      INCLUDE 'lib\pertrb.f'
      INCLUDE 'lib\plots.f'
      INCLUDE 'lib\prdirc.f'
      INCLUDE 'lib\preset.f'
      INCLUDE 'lib\prntdk.f'
      INCLUDE 'lib\prntff.f'
      INCLUDE 'lib\prtrack.f'
      INCLUDE 'lib\ptrelx.f'
      INCLUDE 'lib\readff.f'
      INCLUDE 'lib\reset.f'
      INCLUDE 'lib\skipdf.f'
      INCLUDE 'lib\solveh.f'
      INCLUDE 'lib\solvesp.f'
      INCLUDE 'lib\spgsolv.f'
      INCLUDE 'lib\stdata.f'
      INCLUDE 'lib\swlsolv.f'

      
      INCLUDE 'lib\DATA.f'
