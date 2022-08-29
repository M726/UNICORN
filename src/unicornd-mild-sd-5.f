C**** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****@
C
C                     WEAK-COMPRESSIBLE FLOW SOLUTION
C            (UNSTEADY EQUATIONS FOR CHEMICALLY REACTING FLOWS)
C                                                      Viswanath R Katta
C
C***********************************************************************
C                                                                       
C                      Staggered Mesh: LI=711,LJ=131                               
C                                                                       
C***********************************************************************
C                                                                       
C       Combustion of Heptane with Sandiego Mechanism (2005/12/01)   
C          http://www-mae.ucsd.edu/~combustion/cermech/                                     
C
C   FUEL: 1 - H2,  2--CH4,  3--CH3OH,  4--C2H2, 5--C2H4, 6--C2H6, 
C         7--C3H8, 8--C3H6, 9--CH2O,  10--CO,  11--C7H16
C
C     1-H2          2-O2          3-H           4-O           5-OH        
C     6-H2O         7-HO2         8-H2O2        9-CO         10-CO2       
C    11-HCO        12-CH2O       13-CH4        14-CH3        15-T-CH2     
C    16-S-CH2      17-C2H4       18-CH3O       19-C2H5       20-C2H6      
C    21-CH         22-C2H2       23-C2H3       24-CH2CHO     25-C2H4O     
C    26-CH2CO      27-HCCO       28-C2H        29-CH2OH      30-CH3OH     
C    31-C2H5OH     32-CH3CHO     33-CH3CHOH    34-CH2CH2OH   35-CH3CO     
C    36-CH3CH2O    37-C3H4       38-C3H3       39-C3H5       40-C3H6      
C    41-C3H8       42-I-C3H7     43-N-C3H7     44-C4H6       45-CHO       
C    46-C5H8       47-C7H16      48-C4H8       49-C5H10      50-AR        
C    51-HE         52-N2        
C
C                        FINTE-RATE CHEMISTRY
C                                                                       
C***********************************************************************
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
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/FINJ/NFINJ,IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30)
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
      COMMON/FLOW/ KSYM,LIP,LJP,IOFF,JOFF,NCLR(LNPR),
     1  XS(LI),YS(LJ+LJ-1),RASMN,RASMX,RAS2MN,RAS2MX
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
      OPEN(MP,FILE='input.uni',ACCESS='SEQUENTIAL',STATUS='OLD')
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
   12 CONTINUE
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
     1   RTIN,RTOT,ALENG,NSEG,ISIDE,ITYPE,BCVAL,
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
         CALL PERTRB(IANOIS,IBNOIS,JANOIS,JBNOIS,IRAND)
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
  249 CONTINUE
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
      CALL KESOLV(IRLXKE,RELXKE,TOLRKE,ISCHEM,SIGK,SIGE,RESDK,RESDE)
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
      END
C***********************************************************************
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
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/FINJ/NFINJ,IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30)
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
      END
C***********************************************************************
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
C***********************************************************************
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
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/FINJ/NFINJ,IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30)
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
C***********************************************************************
      SUBROUTINE PTRELX(IEQN,ISOR,RELX,TOLR)                   
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LIJ=LI+LJ,
     1   LPD=55*LE-13*LE-2*LIJ-4) 
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/DUMMY/ UTMP(LJ,LI),AP(LJ,LI),AEE(LJ,LI),AE(LJ,LI),
     1 AW(LJ,LI),AWW(LJ,LI),ANN(LJ,LI),AN(LJ,LI),AS(LJ,LI),ASS(LJ,LI),
     2 RHS1(LJ,LI),FP(LJ,LI),FPA(LJ+2,LI+2),DUM(LPD)
      ENORM=DFLOAT(LE)
C                      RELAX AT POINT BY POINT
      SOR1=DABS(RELX)
      SOR2=(1.0-RELX)
      JBIG=2
      IBIG=2
      IF(IEQN.EQ.1) IBIG=3
      IF(IEQN.EQ.2) JBIG=3
      DO 100 ISORP=1,ISOR
      RSD=0.0
      DO 20 J=JBIG,LJ-1
      DO I=IBIG,LI-1
      JJ=J+1
      II=I+1
      FPAO=FPA(JJ,II)
      FPA(JJ,II)=SOR2*FPAO+(RHS1(J,I)
     1  -ANN(J,I)*FPA(JJ+2,II)-ASS(J,I)*FPA(JJ-2,II)
     2  -AN(J,I)*FPA(JJ+1,II)-AS(J,I)*FPA(JJ-1,II)
     3  -AEE(J,I)*FPA(JJ,II+2)-AWW(J,I)*FPA(JJ,II-2)
     4  -AE(J,I)*FPA(JJ,II+1)-AW(J,I)*FPA(JJ,II-1))*SOR1/AP(J,I)
      RSD=RSD+DABS(FPA(JJ,II)-FPAO)
      ENDDO
   20 CONTINUE
      IF((RSD/ENORM).LE.TOLR) GO TO 101
  100 CONTINUE
      ISOR=ISORP-1
  101 CONTINUE
      ISOR=ISORP
      DO 200 J=JBIG,LJ-1 
      DO I=IBIG,LI-1 
      FP(J,I)=FPA(J+1,I+1)
      ENDDO
  200 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE ADISOLV(IEQN,ISOR,RELX,TOLR)                   
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52,LIJ=LI+LJ,
     1     LPD=55*LE-15*LE-4*LIJ-4)
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/DUMMY/ UTMP(LJ,LI),AP(LJ,LI),AEE(LJ,LI),AE(LJ,LI),
     1 AW(LJ,LI),AWW(LJ,LI),ANN(LJ,LI),AN(LJ,LI),AS(LJ,LI),ASS(LJ,LI),
     2 RHS1(LJ,LI),RHS(LJ,LI),FPA(LJ+2,LI+2),AU(LJ,LI),BU(LJ,LI),
     3 AWM(LJ),ASM(LI),APM(LIJ),DUM(LPD)
      ENORM=DFLOAT(LE)
C                      RELAX BY ADI SCHEME
      SOR1=RELX
      SOR2=(1.0-RELX)
      JBIG=2
      IBIG=2
      IF(IEQN.EQ.1) IBIG=3
      IF(IEQN.EQ.2) JBIG=3
      DO 100 ISORP=1,ISOR
      RSD=0.0
C                      SWEEP IN X - DIRECTION
      DO 20 I=IBIG,LI-1
      II=I+1
      DO J=JBIG,LJ-1
      JJ=J+1
      RHS(J,I)=SOR1*(RHS1(J,I)-(ANN(J,I)*FPA(JJ+2,II)
     1 +AN(J,I)*FPA(JJ+1,II)+AS(J,I)*FPA(JJ-1,II)
     2 +ASS(J,I)*FPA(JJ-2,II)))+SOR2*(AP(J,I)*FPA(JJ,II)
     3 +AEE(J,I)*FPA(JJ,II+2)+AE(J,I)*FPA(JJ,II+1)
     4 +AW(J,I)*FPA(JJ,II-1)+AWW(J,I)*FPA(JJ,II-2))
      ENDDO
   20 CONTINUE
      DO 22 J=JBIG,LJ-1
      JJ=J+1
      A1=AP(J,IBIG)
      AU(J,IBIG)=AE(J,IBIG)/A1
      BU(J,IBIG)=AEE(J,IBIG)/A1
      RHS(J,IBIG)=(RHS(J,IBIG)-FPA(JJ,IBIG)*AW(J,IBIG)
     1         -FPA(JJ,IBIG-1)*AWW(J,IBIG))/A1
      RHS(J,IBIG+1)=(RHS(J,IBIG+1)-FPA(JJ,IBIG)*AWW(J,IBIG+1))
      RHS(J,LI-2)=RHS(J,LI-2)-AEE(J,LI-2)*FPA(JJ,LI+1)
      AWM(J)=AW(J,IBIG+1)
      APM(J)=AP(J,IBIG+1)
   22 CONTINUE
      DO 24 I=IBIG+1,LI-2
      II=I+1
      DO J=JBIG,LJ-1
      JJ=J+1
      A1=APM(J)-AU(J,I-1)*AWM(J)
      AU(J,I)=(AE(J,I)-BU(J,I-1)*AWM(J))/A1
      BU(J,I)=AEE(J,I)/A1
      RHS(J,I)=(RHS(J,I)-RHS(J,I-1)*AWM(J))/A1
      AWM(J)=AW(J,I+1)-AU(J,I-1)*AWW(J,I+1)
      APM(J)=AP(J,I+1)-BU(J,I-1)*AWW(J,I+1)
      RHS(J,I+1)=RHS(J,I+1)-RHS(J,I-1)*AWW(J,I+1)
      ENDDO
   24 CONTINUE
      DO 26 J=JBIG,LJ-1
      JJ=J+1
      A1=APM(J)-AU(J,LI-2)*AWM(J)
      RHS(J,LI-1)=(RHS(J,LI-1)-AEE(J,LI-1)*FPA(JJ,LI+2)
     1         -AE(J,LI-1)*FPA(JJ,LI+1)-RHS(J,LI-2)*AWM(J))/A1
      RHS(J,LI-2)=RHS(J,LI-2)-RHS(J,LI-1)*AU(J,LI-2)
      FPA(JJ,LI)=RHS(J,LI-1)
      FPA(JJ,LI-1)=RHS(J,LI-2)
   26 CONTINUE
      DO 28 I=LI-3,IBIG,-1
      DO J=JBIG,LJ-1
      RHS(J,I)=RHS(J,I)-RHS(J,I+1)*AU(J,I)-RHS(J,I+2)*BU(J,I)
      FPA(J+1,I+1)=RHS(J,I)
      ENDDO
   28 CONTINUE
C                      SWEEP IN Y - DIRECTION
      DO 30 J=JBIG,LJ-1
      JJ=J+1
      DO I=IBIG,LI-1
      II=I+1
      RHS(J,I)=SOR1*(RHS1(J,I)-(AEE(J,I)*FPA(JJ,II+2)
     1 +AE(J,I)*FPA(JJ,II+1)+AW(J,I)*FPA(JJ,II-1)
     2 +AWW(J,I)*FPA(JJ,II-2)))+SOR2*(AP(J,I)*FPA(JJ,II)
     3 +ANN(J,I)*FPA(JJ+2,II)+AN(J,I)*FPA(JJ+1,II)
     4 +AS(J,I)*FPA(JJ-1,II)+ASS(J,I)*FPA(JJ-2,II))
      ENDDO
   30 CONTINUE
      DO 32 I=IBIG,LI-1
      II=I+1
      A1=AP(JBIG,I)
      AU(JBIG,I)=AN(JBIG,I)/A1
      BU(JBIG,I)=ANN(JBIG,I)/A1
      RHS(JBIG,I)=(RHS(JBIG,I)-FPA(JBIG,II)*AS(JBIG,I)
     1         -FPA(JBIG-1,II)*ASS(JBIG,I))/A1
      RHS(JBIG+1,I)=(RHS(JBIG+1,I)-FPA(JBIG,II)*ASS(JBIG+1,I))
      RHS(LJ-2,I)=RHS(LJ-2,I)-ANN(LJ-2,I)*FPA(LJ+1,II)
      ASM(I)=AS(JBIG+1,I)
      APM(I)=AP(JBIG+1,I)
   32 CONTINUE
      DO 34 J=JBIG+1,LJ-2
      JJ=J+1
      DO I=IBIG,LI-1
      II=I+1
      A1=APM(I)-AU(J-1,I)*ASM(I)
      AU(J,I)=(AN(J,I)-BU(J-1,I)*ASM(I))/A1
      BU(J,I)=ANN(J,I)/A1
      RHS(J,I)=(RHS(J,I)-RHS(J-1,I)*ASM(I))/A1
      ASM(I)=AS(J+1,I)-AU(J-1,I)*ASS(J+1,I)
      APM(I)=AP(J+1,I)-BU(J-1,I)*ASS(J+1,I)
      RHS(J+1,I)=RHS(J+1,I)-RHS(J-1,I)*ASS(J+1,I)
      ENDDO
   34 CONTINUE
      DO 36 I=IBIG,LI-1
      II=I+1
      A1=APM(I)-AU(LJ-2,I)*ASM(I)
      RHS(LJ-1,I)=(RHS(LJ-1,I)-ANN(LJ-1,I)*FPA(LJ+2,II)
     1         -AN(LJ-1,I)*FPA(LJ+1,II)-RHS(LJ-2,I)*ASM(I))/A1
      RHS(LJ-2,I)=RHS(LJ-2,I)-RHS(LJ-1,I)*AU(LJ-2,I)
   36 CONTINUE
      DO 38 J=LJ-3,JBIG,-1
      DO I=IBIG,LI-1
      RHS(J,I)=RHS(J,I)-RHS(J+1,I)*AU(J,I)-RHS(J+2,I)*BU(J,I)
      ENDDO
   38 CONTINUE
      DO 40 J=JBIG,LJ-1 
      DO I=IBIG,LI-1
      RSD=RSD+DABS(FPA(J+1,I+1)-RHS(J,I))
      FPA(J+1,I+1)=RHS(J,I)
      ENDDO
   40 CONTINUE
      IF((RSD/ENORM).LE.TOLR) GO TO 101
  100 CONTINUE
      ISOR=ISORP-1
  101 CONTINUE
      ISOR=ISORP
      DO 200 J=1,LJ
      DO I=1,LI
      RHS(J,I)=FPA(J+1,I+1)
      ENDDO
  200 CONTINUE
      RETURN
      END
C***********************************************************************
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
      COMMON/LUMAT/ IBAND(LB),PSL(LL1),PSD(LL)
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
C***********************************************************************
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
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/FINJ/NFINJ,IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30)
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
  200 CONTINUE
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
C***********************************************************************
      SUBROUTINE SPGSOLV(ISOR1,RELXSP,TOLRSP,RESDC,
     1  ISOR2,RELXH,TOLRH,RESDH,SIGH,SIGSP)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LJ2=LJ*2,LSP=52,LRX=544)
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FA0(LJ,LI),FA1(LJ,LI),FA2(LJ,LI),
     *  FA3(LJ,LI),FA4(LJ,LI),FA5(LJ,LI),FA6(LJ,LI),FA7(LJ,LI),
     *  FA8(LJ,LI),FA9(LJ,LI),FAA(LJ,LI),FAB(LJ,LI),FAC(LJ,LI),
     *  FAD(LJ,LI),FAE(LJ,LI),FAF(LJ,LI),FAG(LJ,LI),FAH(LJ,LI),
     *  FAI(LJ,LI),FAJ(LJ,LI),FB0(LJ,LI),FB1(LJ,LI),FB2(LJ,LI),
     *  FB3(LJ,LI),FB4(LJ,LI),FB5(LJ,LI),FB6(LJ,LI),FB7(LJ,LI),
     *  FB8(LJ,LI),FB9(LJ,LI),FBA(LJ,LI),FBB(LJ,LI),FBC(LJ,LI),
     *  FBD(LJ,LI),FBE(LJ,LI),FBF(LJ,LI),FBG(LJ,LI),FBH(LJ,LI),
     *  FBI(LJ,LI),FBJ(LJ,LI),FC0(LJ,LI),FC1(LJ,LI),FC2(LJ,LI),
     *  FC3(LJ,LI),FC4(LJ,LI),FC5(LJ,LI),FC6(LJ,LI),FC7(LJ,LI),
     *  FC8(LJ,LI),FC9(LJ,LI),FCA(LJ,LI),FCB(LJ,LI),
     *  U(LJ,LI),V(LJ,LI),W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),
     *  AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDF(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/CB08/ TPOL1,TPOL2,POLSP(14,LSP),CISP(12,LSP)
      COMMON/CB09/ IBOT(LI),ITOP(LI),JLFT(LJ),JRGT(LJ),
     1       FBOT(8+LSP,LI),FTOP(8+LSP,LI),FLFT(8+LSP,LJ),FRGT(8+LSP,LJ)
      COMMON/CB10/ FBXP(4,LJ),FBXM(4,LJ),FBYP(4,LI),FBYM(4,LI)
      COMMON/REAC/ LREV,NRLIND,ITBEF(LRX),ALF(LRX),AKF(LRX),EXAF(LRX),
     1             ALLOW(LRX),AKLOW(LRX),EALOW(LRX),TROE(LRX,4)
      COMMON/HFORM/ HFO(LSP)
      COMMON/KEPS/ AKREF,EPREF,C1KE,C2KE,CMKE,CMKEQ,EE,AVON,PCON
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
      COMMON/DUMMY/
     *  EA0(LJ,LI),WA0(LJ,LI),PA0(LJ,LI),QA0(LJ,LI),EA1(LJ,LI),
     *  WA1(LJ,LI),PA1(LJ,LI),QA1(LJ,LI),EA2(LJ,LI),WA2(LJ,LI),
     *  PA2(LJ,LI),QA2(LJ,LI),EA3(LJ,LI),WA3(LJ,LI),PA3(LJ,LI),
     *  QA3(LJ,LI),EA4(LJ,LI),WA4(LJ,LI),PA4(LJ,LI),QA4(LJ,LI),
     *  EA5(LJ,LI),WA5(LJ,LI),PA5(LJ,LI),QA5(LJ,LI),EA6(LJ,LI),
     *  WA6(LJ,LI),PA6(LJ,LI),QA6(LJ,LI),EA7(LJ,LI),WA7(LJ,LI),
     *  PA7(LJ,LI),QA7(LJ,LI),EA8(LJ,LI),WA8(LJ,LI),PA8(LJ,LI),
     *  QA8(LJ,LI),EA9(LJ,LI),WA9(LJ,LI),PA9(LJ,LI),QA9(LJ,LI),
     *  EAA(LJ,LI),WAA(LJ,LI),PAA(LJ,LI),QAA(LJ,LI),EAB(LJ,LI),
     *  WAB(LJ,LI),PAB(LJ,LI),QAB(LJ,LI),EAC(LJ,LI),WAC(LJ,LI),
     *  PAC(LJ,LI),QAC(LJ,LI),EAD(LJ,LI),WAD(LJ,LI),PAD(LJ,LI)
      COMMON/DUM0/
     *  QAD(LJ,LI),EAE(LJ,LI),WAE(LJ,LI),PAE(LJ,LI),QAE(LJ,LI),
     *  EAF(LJ,LI),WAF(LJ,LI),PAF(LJ,LI),QAF(LJ,LI),EAG(LJ,LI),
     *  WAG(LJ,LI),PAG(LJ,LI),QAG(LJ,LI),EAH(LJ,LI),WAH(LJ,LI),
     *  PAH(LJ,LI),QAH(LJ,LI),EAI(LJ,LI),WAI(LJ,LI),PAI(LJ,LI),
     *  QAI(LJ,LI),EAJ(LJ,LI),WAJ(LJ,LI),PAJ(LJ,LI),QAJ(LJ,LI),
     *  EB0(LJ,LI),WB0(LJ,LI),PB0(LJ,LI),QB0(LJ,LI),EB1(LJ,LI),
     *  WB1(LJ,LI),PB1(LJ,LI),QB1(LJ,LI),EB2(LJ,LI),WB2(LJ,LI),
     *  PB2(LJ,LI),QB2(LJ,LI),EB3(LJ,LI),WB3(LJ,LI),PB3(LJ,LI),
     *  QB3(LJ,LI),EB4(LJ,LI),WB4(LJ,LI),PB4(LJ,LI),QB4(LJ,LI),
     *  EB5(LJ,LI),WB5(LJ,LI),PB5(LJ,LI),QB5(LJ,LI),EB6(LJ,LI),
     *  WB6(LJ,LI),PB6(LJ,LI),QB6(LJ,LI),EB7(LJ,LI),WB7(LJ,LI),
     *  PB7(LJ,LI),QB7(LJ,LI),EB8(LJ,LI),WB8(LJ,LI),PB8(LJ,LI),
     *  QB8(LJ,LI),EB9(LJ,LI),WB9(LJ,LI),PB9(LJ,LI),QB9(LJ,LI),
     *  EBA(LJ,LI),WBA(LJ,LI),PBA(LJ,LI),QBA(LJ,LI),EBB(LJ,LI),
     *  WBB(LJ,LI),PBB(LJ,LI),QBB(LJ,LI),EBC(LJ,LI),WBC(LJ,LI),
     *  PBC(LJ,LI),QBC(LJ,LI),EBD(LJ,LI),WBD(LJ,LI),PBD(LJ,LI),
     *  QBD(LJ,LI),EBE(LJ,LI),WBE(LJ,LI),PBE(LJ,LI),QBE(LJ,LI),
     *  EBF(LJ,LI),WBF(LJ,LI),PBF(LJ,LI),QBF(LJ,LI),EBG(LJ,LI),
     *  WBG(LJ,LI),PBG(LJ,LI),QBG(LJ,LI),EBH(LJ,LI),WBH(LJ,LI),
     *  PBH(LJ,LI),QBH(LJ,LI),EBI(LJ,LI),WBI(LJ,LI),PBI(LJ,LI),
     *  QBI(LJ,LI),EBJ(LJ,LI),WBJ(LJ,LI),PBJ(LJ,LI),QBJ(LJ,LI),
     *  EC0(LJ,LI),WC0(LJ,LI),PC0(LJ,LI),QC0(LJ,LI),EC1(LJ,LI),
     *  WC1(LJ,LI),PC1(LJ,LI),QC1(LJ,LI),EC2(LJ,LI),WC2(LJ,LI),
     *  PC2(LJ,LI),QC2(LJ,LI),EC3(LJ,LI),WC3(LJ,LI),PC3(LJ,LI),
     *  QC3(LJ,LI),EC4(LJ,LI),WC4(LJ,LI),PC4(LJ,LI),QC4(LJ,LI),
     *  EC5(LJ,LI),WC5(LJ,LI),PC5(LJ,LI),QC5(LJ,LI),EC6(LJ,LI),
     *  WC6(LJ,LI),PC6(LJ,LI),QC6(LJ,LI),EC7(LJ,LI),WC7(LJ,LI),
     *  PC7(LJ,LI),QC7(LJ,LI),EC8(LJ,LI),WC8(LJ,LI),PC8(LJ,LI),
     *  QC8(LJ,LI),EC9(LJ,LI),WC9(LJ,LI),PC9(LJ,LI),QC9(LJ,LI),
     *  ECA(LJ,LI),WCA(LJ,LI),PCA(LJ,LI),QCA(LJ,LI)
      COMMON/DUMX/ ZE(LJ,LI),ZW(LJ,LI),ZN(LJ,LI),ZS(LJ,LI)
      COMMON/DUMY/            SA0(LJ,LI),SA1(LJ,LI),SA2(LJ,LI),
     *  SA3(LJ,LI),SA4(LJ,LI),SA5(LJ,LI),SA6(LJ,LI),SA7(LJ,LI),
     *  SA8(LJ,LI),SA9(LJ,LI),SAA(LJ,LI),SAB(LJ,LI),SAC(LJ,LI),
     *  SAD(LJ,LI),SAE(LJ,LI),SAF(LJ,LI),SAG(LJ,LI),SAH(LJ,LI),
     *  SAI(LJ,LI),SAJ(LJ,LI),SB0(LJ,LI),SB1(LJ,LI),SB2(LJ,LI),
     *  SB3(LJ,LI),SB4(LJ,LI),SB5(LJ,LI),SB6(LJ,LI),SB7(LJ,LI),
     *  SB8(LJ,LI),SB9(LJ,LI),SBA(LJ,LI),SBB(LJ,LI),SBC(LJ,LI),
     *  SBD(LJ,LI),SBE(LJ,LI),SBF(LJ,LI),SBG(LJ,LI),SBH(LJ,LI),
     *  SBI(LJ,LI),SBJ(LJ,LI),SC0(LJ,LI),SC1(LJ,LI),SC2(LJ,LI),
     *  SC3(LJ,LI),SC4(LJ,LI),SC5(LJ,LI),SC6(LJ,LI),SC7(LJ,LI),
     *  SC8(LJ,LI),SC9(LJ,LI),SCA(LJ,LI),SCB(LJ,LI),
     *  FPZ(LJ,LI),EMU(LJ,LI),EKT(LJ,LI)
      COMMON/DUMS/ STE(LJ,LI),STW(LJ,LI),STN(LJ,LI),STS(LJ,LI),
     1  SST1(LJ,LI),SST2(LJ,LI),SOUR(LJ,LI,6)
      DIMENSION TCPSP(LSP)
      SYM=DFLOAT(ISYM)
      HCONVR=GASC/WMSTR/HSTR
      HNORM=DFLOAT(LE)
      SPNORM=DFLOAT(LE)*DFLOAT(LSP-1)
      DO 10 I=1,LI
      DO J=1,LJ
      TKD=TK(J,I)*TSTR
      TKD2=TKD*TKD
      TKD3=TKD*TKD2
      TKD4=TKD*TKD3
      TKD5=TKD*TKD4
      IPOLY=1
      IF(TKD.GT.TPOL1) IPOLY=8
      TCP=0.0
      DO 11 ISP=1,LSP
      TCPSP(ISP)=(POLSP(IPOLY,ISP)+POLSP(IPOLY+1,ISP)*TKD
     1 +POLSP(IPOLY+2,ISP)*TKD2+POLSP(IPOLY+3,ISP)*TKD3
     2 +POLSP(IPOLY+4,ISP)*TKD4)/WM(ISP)
   11 CONTINUE
      TCPSP(  1)=TCPSP(  1)*FA0(J,I)
      TCPSP(  2)=TCPSP(  2)*FA1(J,I)
      TCPSP(  3)=TCPSP(  3)*FA2(J,I)
      TCPSP(  4)=TCPSP(  4)*FA3(J,I)
      TCPSP(  5)=TCPSP(  5)*FA4(J,I)
      TCPSP(  6)=TCPSP(  6)*FA5(J,I)
      TCPSP(  7)=TCPSP(  7)*FA6(J,I)
      TCPSP(  8)=TCPSP(  8)*FA7(J,I)
      TCPSP(  9)=TCPSP(  9)*FA8(J,I)
      TCPSP( 10)=TCPSP( 10)*FA9(J,I)
      TCPSP( 11)=TCPSP( 11)*FAA(J,I)
      TCPSP( 12)=TCPSP( 12)*FAB(J,I)
      TCPSP( 13)=TCPSP( 13)*FAC(J,I)
      TCPSP( 14)=TCPSP( 14)*FAD(J,I)
      TCPSP( 15)=TCPSP( 15)*FAE(J,I)
      TCPSP( 16)=TCPSP( 16)*FAF(J,I)
      TCPSP( 17)=TCPSP( 17)*FAG(J,I)
      TCPSP( 18)=TCPSP( 18)*FAH(J,I)
      TCPSP( 19)=TCPSP( 19)*FAI(J,I)
      TCPSP( 20)=TCPSP( 20)*FAJ(J,I)
      TCPSP( 21)=TCPSP( 21)*FB0(J,I)
      TCPSP( 22)=TCPSP( 22)*FB1(J,I)
      TCPSP( 23)=TCPSP( 23)*FB2(J,I)
      TCPSP( 24)=TCPSP( 24)*FB3(J,I)
      TCPSP( 25)=TCPSP( 25)*FB4(J,I)
      TCPSP( 26)=TCPSP( 26)*FB5(J,I)
      TCPSP( 27)=TCPSP( 27)*FB6(J,I)
      TCPSP( 28)=TCPSP( 28)*FB7(J,I)
      TCPSP( 29)=TCPSP( 29)*FB8(J,I)
      TCPSP( 30)=TCPSP( 30)*FB9(J,I)
      TCPSP( 31)=TCPSP( 31)*FBA(J,I)
      TCPSP( 32)=TCPSP( 32)*FBB(J,I)
      TCPSP( 33)=TCPSP( 33)*FBC(J,I)
      TCPSP( 34)=TCPSP( 34)*FBD(J,I)
      TCPSP( 35)=TCPSP( 35)*FBE(J,I)
      TCPSP( 36)=TCPSP( 36)*FBF(J,I)
      TCPSP( 37)=TCPSP( 37)*FBG(J,I)
      TCPSP( 38)=TCPSP( 38)*FBH(J,I)
      TCPSP( 39)=TCPSP( 39)*FBI(J,I)
      TCPSP( 40)=TCPSP( 40)*FBJ(J,I)
      TCPSP( 41)=TCPSP( 41)*FC0(J,I)
      TCPSP( 42)=TCPSP( 42)*FC1(J,I)
      TCPSP( 43)=TCPSP( 43)*FC2(J,I)
      TCPSP( 44)=TCPSP( 44)*FC3(J,I)
      TCPSP( 45)=TCPSP( 45)*FC4(J,I)
      TCPSP( 46)=TCPSP( 46)*FC5(J,I)
      TCPSP( 47)=TCPSP( 47)*FC6(J,I)
      TCPSP( 48)=TCPSP( 48)*FC7(J,I)
      TCPSP( 49)=TCPSP( 49)*FC8(J,I)
      TCPSP( 50)=TCPSP( 50)*FC9(J,I)
      TCPSP( 51)=TCPSP( 51)*FCA(J,I)
      TCPSP( 52)=TCPSP( 52)*FCB(J,I)
      TCP=0.0
      DO 12 ISP=1,LSP
      TCP=TCP+TCPSP(ISP)
   12 CONTINUE
      TCP=TCP/WMSTR*GASC*TSTR/HSTR
      EKT(J,I)=REI*(BETA2*VTC(J,I)/TCP+TMU(J,I)/SIGH)
      EMU(J,I)=REI*BETA3*RHO(J,I)
      FPZ(J,I)=HT(J,I)
      SST1(J,I)=STMF(J,I)
      SST2(J,I)=STND(J,I)
      SA0(J,I)=FA0(J,I)
      SA1(J,I)=FA1(J,I)
      SA2(J,I)=FA2(J,I)
      SA3(J,I)=FA3(J,I)
      SA4(J,I)=FA4(J,I)
      SA5(J,I)=FA5(J,I)
      SA6(J,I)=FA6(J,I)
      SA7(J,I)=FA7(J,I)
      SA8(J,I)=FA8(J,I)
      SA9(J,I)=FA9(J,I)
      SAA(J,I)=FAA(J,I)
      SAB(J,I)=FAB(J,I)
      SAC(J,I)=FAC(J,I)
      SAD(J,I)=FAD(J,I)
      SAE(J,I)=FAE(J,I)
      SAF(J,I)=FAF(J,I)
      SAG(J,I)=FAG(J,I)
      SAH(J,I)=FAH(J,I)
      SAI(J,I)=FAI(J,I)
      SAJ(J,I)=FAJ(J,I)
      SB0(J,I)=FB0(J,I)
      SB1(J,I)=FB1(J,I)
      SB2(J,I)=FB2(J,I)
      SB3(J,I)=FB3(J,I)
      SB4(J,I)=FB4(J,I)
      SB5(J,I)=FB5(J,I)
      SB6(J,I)=FB6(J,I)
      SB7(J,I)=FB7(J,I)
      SB8(J,I)=FB8(J,I)
      SB9(J,I)=FB9(J,I)
      SBA(J,I)=FBA(J,I)
      SBB(J,I)=FBB(J,I)
      SBC(J,I)=FBC(J,I)
      SBD(J,I)=FBD(J,I)
      SBE(J,I)=FBE(J,I)
      SBF(J,I)=FBF(J,I)
      SBG(J,I)=FBG(J,I)
      SBH(J,I)=FBH(J,I)
      SBI(J,I)=FBI(J,I)
      SBJ(J,I)=FBJ(J,I)
      SC0(J,I)=FC0(J,I)
      SC1(J,I)=FC1(J,I)
      SC2(J,I)=FC2(J,I)
      SC3(J,I)=FC3(J,I)
      SC4(J,I)=FC4(J,I)
      SC5(J,I)=FC5(J,I)
      SC6(J,I)=FC6(J,I)
      SC7(J,I)=FC7(J,I)
      SC8(J,I)=FC8(J,I)
      SC9(J,I)=FC9(J,I)
      SCA(J,I)=FCA(J,I)
      SCB(J,I)=FCB(J,I)
      ENDDO
   10 CONTINUE
      DO 14 I=2,LI
      DO J=2,LJ
      M=I-1
      N=J-1
      TURS=REI*(TMU(J,I)+TMU(N,I))/SIGSP
      TURW=REI*(TMU(J,I)+TMU(J,M))/SIGSP
C-----------------------------------------------------------------------
      ZS(J,I)=(EKT(J,I)+EKT(N,I))/YYS(J)
      STS(J,I)=(EMU(J,I)*STDF(J,I)+EMU(N,I)*STDF(N,I)+TURS)/YYS(J)
      QA0(J,I)=(EMU(J,I)*VDF(J,I,  1)+EMU(N,I)*VDF(N,I,  1)+TURS)/YYS(J)
      QA1(J,I)=(EMU(J,I)*VDF(J,I,  2)+EMU(N,I)*VDF(N,I,  2)+TURS)/YYS(J)
      QA2(J,I)=(EMU(J,I)*VDF(J,I,  3)+EMU(N,I)*VDF(N,I,  3)+TURS)/YYS(J)
      QA3(J,I)=(EMU(J,I)*VDF(J,I,  4)+EMU(N,I)*VDF(N,I,  4)+TURS)/YYS(J)
      QA4(J,I)=(EMU(J,I)*VDF(J,I,  5)+EMU(N,I)*VDF(N,I,  5)+TURS)/YYS(J)
      QA5(J,I)=(EMU(J,I)*VDF(J,I,  6)+EMU(N,I)*VDF(N,I,  6)+TURS)/YYS(J)
      QA6(J,I)=(EMU(J,I)*VDF(J,I,  7)+EMU(N,I)*VDF(N,I,  7)+TURS)/YYS(J)
      QA7(J,I)=(EMU(J,I)*VDF(J,I,  8)+EMU(N,I)*VDF(N,I,  8)+TURS)/YYS(J)
      QA8(J,I)=(EMU(J,I)*VDF(J,I,  9)+EMU(N,I)*VDF(N,I,  9)+TURS)/YYS(J)
      QA9(J,I)=(EMU(J,I)*VDF(J,I, 10)+EMU(N,I)*VDF(N,I, 10)+TURS)/YYS(J)
      QAA(J,I)=(EMU(J,I)*VDF(J,I, 11)+EMU(N,I)*VDF(N,I, 11)+TURS)/YYS(J)
      QAB(J,I)=(EMU(J,I)*VDF(J,I, 12)+EMU(N,I)*VDF(N,I, 12)+TURS)/YYS(J)
      QAC(J,I)=(EMU(J,I)*VDF(J,I, 13)+EMU(N,I)*VDF(N,I, 13)+TURS)/YYS(J)
      QAD(J,I)=(EMU(J,I)*VDF(J,I, 14)+EMU(N,I)*VDF(N,I, 14)+TURS)/YYS(J)
      QAE(J,I)=(EMU(J,I)*VDF(J,I, 15)+EMU(N,I)*VDF(N,I, 15)+TURS)/YYS(J)
      QAF(J,I)=(EMU(J,I)*VDF(J,I, 16)+EMU(N,I)*VDF(N,I, 16)+TURS)/YYS(J)
      QAG(J,I)=(EMU(J,I)*VDF(J,I, 17)+EMU(N,I)*VDF(N,I, 17)+TURS)/YYS(J)
      QAH(J,I)=(EMU(J,I)*VDF(J,I, 18)+EMU(N,I)*VDF(N,I, 18)+TURS)/YYS(J)
      QAI(J,I)=(EMU(J,I)*VDF(J,I, 19)+EMU(N,I)*VDF(N,I, 19)+TURS)/YYS(J)
      QAJ(J,I)=(EMU(J,I)*VDF(J,I, 20)+EMU(N,I)*VDF(N,I, 20)+TURS)/YYS(J)
      QB0(J,I)=(EMU(J,I)*VDF(J,I, 21)+EMU(N,I)*VDF(N,I, 21)+TURS)/YYS(J)
      QB1(J,I)=(EMU(J,I)*VDF(J,I, 22)+EMU(N,I)*VDF(N,I, 22)+TURS)/YYS(J)
      QB2(J,I)=(EMU(J,I)*VDF(J,I, 23)+EMU(N,I)*VDF(N,I, 23)+TURS)/YYS(J)
      QB3(J,I)=(EMU(J,I)*VDF(J,I, 24)+EMU(N,I)*VDF(N,I, 24)+TURS)/YYS(J)
      QB4(J,I)=(EMU(J,I)*VDF(J,I, 25)+EMU(N,I)*VDF(N,I, 25)+TURS)/YYS(J)
      QB5(J,I)=(EMU(J,I)*VDF(J,I, 26)+EMU(N,I)*VDF(N,I, 26)+TURS)/YYS(J)
      QB6(J,I)=(EMU(J,I)*VDF(J,I, 27)+EMU(N,I)*VDF(N,I, 27)+TURS)/YYS(J)
      QB7(J,I)=(EMU(J,I)*VDF(J,I, 28)+EMU(N,I)*VDF(N,I, 28)+TURS)/YYS(J)
      QB8(J,I)=(EMU(J,I)*VDF(J,I, 29)+EMU(N,I)*VDF(N,I, 29)+TURS)/YYS(J)
      QB9(J,I)=(EMU(J,I)*VDF(J,I, 30)+EMU(N,I)*VDF(N,I, 30)+TURS)/YYS(J)
      QBA(J,I)=(EMU(J,I)*VDF(J,I, 31)+EMU(N,I)*VDF(N,I, 31)+TURS)/YYS(J)
      QBB(J,I)=(EMU(J,I)*VDF(J,I, 32)+EMU(N,I)*VDF(N,I, 32)+TURS)/YYS(J)
      QBC(J,I)=(EMU(J,I)*VDF(J,I, 33)+EMU(N,I)*VDF(N,I, 33)+TURS)/YYS(J)
      QBD(J,I)=(EMU(J,I)*VDF(J,I, 34)+EMU(N,I)*VDF(N,I, 34)+TURS)/YYS(J)
      QBE(J,I)=(EMU(J,I)*VDF(J,I, 35)+EMU(N,I)*VDF(N,I, 35)+TURS)/YYS(J)
      QBF(J,I)=(EMU(J,I)*VDF(J,I, 36)+EMU(N,I)*VDF(N,I, 36)+TURS)/YYS(J)
      QBG(J,I)=(EMU(J,I)*VDF(J,I, 37)+EMU(N,I)*VDF(N,I, 37)+TURS)/YYS(J)
      QBH(J,I)=(EMU(J,I)*VDF(J,I, 38)+EMU(N,I)*VDF(N,I, 38)+TURS)/YYS(J)
      QBI(J,I)=(EMU(J,I)*VDF(J,I, 39)+EMU(N,I)*VDF(N,I, 39)+TURS)/YYS(J)
      QBJ(J,I)=(EMU(J,I)*VDF(J,I, 40)+EMU(N,I)*VDF(N,I, 40)+TURS)/YYS(J)
      QC0(J,I)=(EMU(J,I)*VDF(J,I, 41)+EMU(N,I)*VDF(N,I, 41)+TURS)/YYS(J)
      QC1(J,I)=(EMU(J,I)*VDF(J,I, 42)+EMU(N,I)*VDF(N,I, 42)+TURS)/YYS(J)
      QC2(J,I)=(EMU(J,I)*VDF(J,I, 43)+EMU(N,I)*VDF(N,I, 43)+TURS)/YYS(J)
      QC3(J,I)=(EMU(J,I)*VDF(J,I, 44)+EMU(N,I)*VDF(N,I, 44)+TURS)/YYS(J)
      QC4(J,I)=(EMU(J,I)*VDF(J,I, 45)+EMU(N,I)*VDF(N,I, 45)+TURS)/YYS(J)
      QC5(J,I)=(EMU(J,I)*VDF(J,I, 46)+EMU(N,I)*VDF(N,I, 46)+TURS)/YYS(J)
      QC6(J,I)=(EMU(J,I)*VDF(J,I, 47)+EMU(N,I)*VDF(N,I, 47)+TURS)/YYS(J)
      QC7(J,I)=(EMU(J,I)*VDF(J,I, 48)+EMU(N,I)*VDF(N,I, 48)+TURS)/YYS(J)
      QC8(J,I)=(EMU(J,I)*VDF(J,I, 49)+EMU(N,I)*VDF(N,I, 49)+TURS)/YYS(J)
      QC9(J,I)=(EMU(J,I)*VDF(J,I, 50)+EMU(N,I)*VDF(N,I, 50)+TURS)/YYS(J)
      QCA(J,I)=(EMU(J,I)*VDF(J,I, 51)+EMU(N,I)*VDF(N,I, 51)+TURS)/YYS(J)
C-----------------------------------------------------------------------
      ZW(J,I)=(EKT(J,I)+EKT(J,M))/XXS(I)
      STW(J,I)=(EMU(J,I)*STDF(J,I)+EMU(J,M)*STDF(J,M)+TURW)/XXS(I)
      WA0(J,I)=(EMU(J,I)*VDF(J,I,  1)+EMU(J,M)*VDF(J,M,  1)+TURW)/XXS(I)
      WA1(J,I)=(EMU(J,I)*VDF(J,I,  2)+EMU(J,M)*VDF(J,M,  2)+TURW)/XXS(I)
      WA2(J,I)=(EMU(J,I)*VDF(J,I,  3)+EMU(J,M)*VDF(J,M,  3)+TURW)/XXS(I)
      WA3(J,I)=(EMU(J,I)*VDF(J,I,  4)+EMU(J,M)*VDF(J,M,  4)+TURW)/XXS(I)
      WA4(J,I)=(EMU(J,I)*VDF(J,I,  5)+EMU(J,M)*VDF(J,M,  5)+TURW)/XXS(I)
      WA5(J,I)=(EMU(J,I)*VDF(J,I,  6)+EMU(J,M)*VDF(J,M,  6)+TURW)/XXS(I)
      WA6(J,I)=(EMU(J,I)*VDF(J,I,  7)+EMU(J,M)*VDF(J,M,  7)+TURW)/XXS(I)
      WA7(J,I)=(EMU(J,I)*VDF(J,I,  8)+EMU(J,M)*VDF(J,M,  8)+TURW)/XXS(I)
      WA8(J,I)=(EMU(J,I)*VDF(J,I,  9)+EMU(J,M)*VDF(J,M,  9)+TURW)/XXS(I)
      WA9(J,I)=(EMU(J,I)*VDF(J,I, 10)+EMU(J,M)*VDF(J,M, 10)+TURW)/XXS(I)
      WAA(J,I)=(EMU(J,I)*VDF(J,I, 11)+EMU(J,M)*VDF(J,M, 11)+TURW)/XXS(I)
      WAB(J,I)=(EMU(J,I)*VDF(J,I, 12)+EMU(J,M)*VDF(J,M, 12)+TURW)/XXS(I)
      WAC(J,I)=(EMU(J,I)*VDF(J,I, 13)+EMU(J,M)*VDF(J,M, 13)+TURW)/XXS(I)
      WAD(J,I)=(EMU(J,I)*VDF(J,I, 14)+EMU(J,M)*VDF(J,M, 14)+TURW)/XXS(I)
      WAE(J,I)=(EMU(J,I)*VDF(J,I, 15)+EMU(J,M)*VDF(J,M, 15)+TURW)/XXS(I)
      WAF(J,I)=(EMU(J,I)*VDF(J,I, 16)+EMU(J,M)*VDF(J,M, 16)+TURW)/XXS(I)
      WAG(J,I)=(EMU(J,I)*VDF(J,I, 17)+EMU(J,M)*VDF(J,M, 17)+TURW)/XXS(I)
      WAH(J,I)=(EMU(J,I)*VDF(J,I, 18)+EMU(J,M)*VDF(J,M, 18)+TURW)/XXS(I)
      WAI(J,I)=(EMU(J,I)*VDF(J,I, 19)+EMU(J,M)*VDF(J,M, 19)+TURW)/XXS(I)
      WAJ(J,I)=(EMU(J,I)*VDF(J,I, 20)+EMU(J,M)*VDF(J,M, 20)+TURW)/XXS(I)
      WB0(J,I)=(EMU(J,I)*VDF(J,I, 21)+EMU(J,M)*VDF(J,M, 21)+TURW)/XXS(I)
      WB1(J,I)=(EMU(J,I)*VDF(J,I, 22)+EMU(J,M)*VDF(J,M, 22)+TURW)/XXS(I)
      WB2(J,I)=(EMU(J,I)*VDF(J,I, 23)+EMU(J,M)*VDF(J,M, 23)+TURW)/XXS(I)
      WB3(J,I)=(EMU(J,I)*VDF(J,I, 24)+EMU(J,M)*VDF(J,M, 24)+TURW)/XXS(I)
      WB4(J,I)=(EMU(J,I)*VDF(J,I, 25)+EMU(J,M)*VDF(J,M, 25)+TURW)/XXS(I)
      WB5(J,I)=(EMU(J,I)*VDF(J,I, 26)+EMU(J,M)*VDF(J,M, 26)+TURW)/XXS(I)
      WB6(J,I)=(EMU(J,I)*VDF(J,I, 27)+EMU(J,M)*VDF(J,M, 27)+TURW)/XXS(I)
      WB7(J,I)=(EMU(J,I)*VDF(J,I, 28)+EMU(J,M)*VDF(J,M, 28)+TURW)/XXS(I)
      WB8(J,I)=(EMU(J,I)*VDF(J,I, 29)+EMU(J,M)*VDF(J,M, 29)+TURW)/XXS(I)
      WB9(J,I)=(EMU(J,I)*VDF(J,I, 30)+EMU(J,M)*VDF(J,M, 30)+TURW)/XXS(I)
      WBA(J,I)=(EMU(J,I)*VDF(J,I, 31)+EMU(J,M)*VDF(J,M, 31)+TURW)/XXS(I)
      WBB(J,I)=(EMU(J,I)*VDF(J,I, 32)+EMU(J,M)*VDF(J,M, 32)+TURW)/XXS(I)
      WBC(J,I)=(EMU(J,I)*VDF(J,I, 33)+EMU(J,M)*VDF(J,M, 33)+TURW)/XXS(I)
      WBD(J,I)=(EMU(J,I)*VDF(J,I, 34)+EMU(J,M)*VDF(J,M, 34)+TURW)/XXS(I)
      WBE(J,I)=(EMU(J,I)*VDF(J,I, 35)+EMU(J,M)*VDF(J,M, 35)+TURW)/XXS(I)
      WBF(J,I)=(EMU(J,I)*VDF(J,I, 36)+EMU(J,M)*VDF(J,M, 36)+TURW)/XXS(I)
      WBG(J,I)=(EMU(J,I)*VDF(J,I, 37)+EMU(J,M)*VDF(J,M, 37)+TURW)/XXS(I)
      WBH(J,I)=(EMU(J,I)*VDF(J,I, 38)+EMU(J,M)*VDF(J,M, 38)+TURW)/XXS(I)
      WBI(J,I)=(EMU(J,I)*VDF(J,I, 39)+EMU(J,M)*VDF(J,M, 39)+TURW)/XXS(I)
      WBJ(J,I)=(EMU(J,I)*VDF(J,I, 40)+EMU(J,M)*VDF(J,M, 40)+TURW)/XXS(I)
      WC0(J,I)=(EMU(J,I)*VDF(J,I, 41)+EMU(J,M)*VDF(J,M, 41)+TURW)/XXS(I)
      WC1(J,I)=(EMU(J,I)*VDF(J,I, 42)+EMU(J,M)*VDF(J,M, 42)+TURW)/XXS(I)
      WC2(J,I)=(EMU(J,I)*VDF(J,I, 43)+EMU(J,M)*VDF(J,M, 43)+TURW)/XXS(I)
      WC3(J,I)=(EMU(J,I)*VDF(J,I, 44)+EMU(J,M)*VDF(J,M, 44)+TURW)/XXS(I)
      WC4(J,I)=(EMU(J,I)*VDF(J,I, 45)+EMU(J,M)*VDF(J,M, 45)+TURW)/XXS(I)
      WC5(J,I)=(EMU(J,I)*VDF(J,I, 46)+EMU(J,M)*VDF(J,M, 46)+TURW)/XXS(I)
      WC6(J,I)=(EMU(J,I)*VDF(J,I, 47)+EMU(J,M)*VDF(J,M, 47)+TURW)/XXS(I)
      WC7(J,I)=(EMU(J,I)*VDF(J,I, 48)+EMU(J,M)*VDF(J,M, 48)+TURW)/XXS(I)
      WC8(J,I)=(EMU(J,I)*VDF(J,I, 49)+EMU(J,M)*VDF(J,M, 49)+TURW)/XXS(I)
      WC9(J,I)=(EMU(J,I)*VDF(J,I, 50)+EMU(J,M)*VDF(J,M, 50)+TURW)/XXS(I)
      WCA(J,I)=(EMU(J,I)*VDF(J,I, 51)+EMU(J,M)*VDF(J,M, 51)+TURW)/XXS(I)
      ENDDO
   14 CONTINUE
C-----------------------------COEFFICIENTS------------------------------
      DO 20 I=2,LI-1
      DO J=2,LJ-1
      IF(ISKIP(J,I).GE.10) GO TO 21
      DXE=XXS(I+1)
      DXW=XXS(I)
      DXCP=XXC(I)
C
      DYN=YYS(J+1)
      DYS=YYS(J)
      DYCP=YYC(J)
C
      DTBRX=DT/DXCP
      DTBRY=DT/DYCP
      RHOUE=0.5*(RHO(J,I+1)+RHO(J,I))*U(J,I+1)
      RHOUW=0.5*(RHO(J,I)+RHO(J,I-1))*U(J,I)
      RHOVN=0.5*(RHO(J+1,I)+RHO(J,I))*V(J+1,I)
      RHOVS=0.5*(RHO(J,I)+RHO(J-1,I))*V(J,I)
C-----------------------------------------------------------------------
      SYMZ=SYM*DT*0.5*EKT(J,I)/Y(J,I)
      SYM1=SYM*DT*0.5*EMU(J,I)/Y(J,I)
      SYM2=SYM*DT*0.5*REI*TMU(J,I)/SIGSP/Y(J,I)
      SYMST=SYM1*STDF(J,I)+SYM2
      SYMA0=SYM1*VDF(J,I,  1)+SYM2
      SYMA1=SYM1*VDF(J,I,  2)+SYM2
      SYMA2=SYM1*VDF(J,I,  3)+SYM2
      SYMA3=SYM1*VDF(J,I,  4)+SYM2
      SYMA4=SYM1*VDF(J,I,  5)+SYM2
      SYMA5=SYM1*VDF(J,I,  6)+SYM2
      SYMA6=SYM1*VDF(J,I,  7)+SYM2
      SYMA7=SYM1*VDF(J,I,  8)+SYM2
      SYMA8=SYM1*VDF(J,I,  9)+SYM2
      SYMA9=SYM1*VDF(J,I, 10)+SYM2
      SYMAA=SYM1*VDF(J,I, 11)+SYM2
      SYMAB=SYM1*VDF(J,I, 12)+SYM2
      SYMAC=SYM1*VDF(J,I, 13)+SYM2
      SYMAD=SYM1*VDF(J,I, 14)+SYM2
      SYMAE=SYM1*VDF(J,I, 15)+SYM2
      SYMAF=SYM1*VDF(J,I, 16)+SYM2
      SYMAG=SYM1*VDF(J,I, 17)+SYM2
      SYMAH=SYM1*VDF(J,I, 18)+SYM2
      SYMAI=SYM1*VDF(J,I, 19)+SYM2
      SYMAJ=SYM1*VDF(J,I, 20)+SYM2
      SYMB0=SYM1*VDF(J,I, 21)+SYM2
      SYMB1=SYM1*VDF(J,I, 22)+SYM2
      SYMB2=SYM1*VDF(J,I, 23)+SYM2
      SYMB3=SYM1*VDF(J,I, 24)+SYM2
      SYMB4=SYM1*VDF(J,I, 25)+SYM2
      SYMB5=SYM1*VDF(J,I, 26)+SYM2
      SYMB6=SYM1*VDF(J,I, 27)+SYM2
      SYMB7=SYM1*VDF(J,I, 28)+SYM2
      SYMB8=SYM1*VDF(J,I, 29)+SYM2
      SYMB9=SYM1*VDF(J,I, 30)+SYM2
      SYMBA=SYM1*VDF(J,I, 31)+SYM2
      SYMBB=SYM1*VDF(J,I, 32)+SYM2
      SYMBC=SYM1*VDF(J,I, 33)+SYM2
      SYMBD=SYM1*VDF(J,I, 34)+SYM2
      SYMBE=SYM1*VDF(J,I, 35)+SYM2
      SYMBF=SYM1*VDF(J,I, 36)+SYM2
      SYMBG=SYM1*VDF(J,I, 37)+SYM2
      SYMBH=SYM1*VDF(J,I, 38)+SYM2
      SYMBI=SYM1*VDF(J,I, 39)+SYM2
      SYMBJ=SYM1*VDF(J,I, 40)+SYM2
      SYMC0=SYM1*VDF(J,I, 41)+SYM2
      SYMC1=SYM1*VDF(J,I, 42)+SYM2
      SYMC2=SYM1*VDF(J,I, 43)+SYM2
      SYMC3=SYM1*VDF(J,I, 44)+SYM2
      SYMC4=SYM1*VDF(J,I, 45)+SYM2
      SYMC5=SYM1*VDF(J,I, 46)+SYM2
      SYMC6=SYM1*VDF(J,I, 47)+SYM2
      SYMC7=SYM1*VDF(J,I, 48)+SYM2
      SYMC8=SYM1*VDF(J,I, 49)+SYM2
      SYMC9=SYM1*VDF(J,I, 50)+SYM2
      SYMCA=SYM1*VDF(J,I, 51)+SYM2
C-----------------------------------------------------------------------
      ZCONV=DABS(RHOUE)
      STCONV=ZCONV
      A0CONV=ZCONV
      A1CONV=ZCONV
      A2CONV=ZCONV
      A3CONV=ZCONV
      A4CONV=ZCONV
      A5CONV=ZCONV
      A6CONV=ZCONV
      A7CONV=ZCONV
      A8CONV=ZCONV
      A9CONV=ZCONV
      AACONV=ZCONV
      ABCONV=ZCONV
      ACCONV=ZCONV
      ADCONV=ZCONV
      AECONV=ZCONV
      AFCONV=ZCONV
      AGCONV=ZCONV
      AHCONV=ZCONV
      AICONV=ZCONV
      AJCONV=ZCONV
      B0CONV=ZCONV
      B1CONV=ZCONV
      B2CONV=ZCONV
      B3CONV=ZCONV
      B4CONV=ZCONV
      B5CONV=ZCONV
      B6CONV=ZCONV
      B7CONV=ZCONV
      B8CONV=ZCONV
      B9CONV=ZCONV
      BACONV=ZCONV
      BBCONV=ZCONV
      BCCONV=ZCONV
      BDCONV=ZCONV
      BECONV=ZCONV
      BFCONV=ZCONV
      BGCONV=ZCONV
      BHCONV=ZCONV
      BICONV=ZCONV
      BJCONV=ZCONV
      C0CONV=ZCONV
      C1CONV=ZCONV
      C2CONV=ZCONV
      C3CONV=ZCONV
      C4CONV=ZCONV
      C5CONV=ZCONV
      C6CONV=ZCONV
      C7CONV=ZCONV
      C8CONV=ZCONV
      C9CONV=ZCONV
      CACONV=ZCONV
C-----------------------------------------------------------------------
      ZDIFU=ZW(J,I+1)
      STDIFU=STW(J,I+1)
      A0DIFU=WA0(J,I+1)
      A1DIFU=WA1(J,I+1)
      A2DIFU=WA2(J,I+1)
      A3DIFU=WA3(J,I+1)
      A4DIFU=WA4(J,I+1)
      A5DIFU=WA5(J,I+1)
      A6DIFU=WA6(J,I+1)
      A7DIFU=WA7(J,I+1)
      A8DIFU=WA8(J,I+1)
      A9DIFU=WA9(J,I+1)
      AADIFU=WAA(J,I+1)
      ABDIFU=WAB(J,I+1)
      ACDIFU=WAC(J,I+1)
      ADDIFU=WAD(J,I+1)
      AEDIFU=WAE(J,I+1)
      AFDIFU=WAF(J,I+1)
      AGDIFU=WAG(J,I+1)
      AHDIFU=WAH(J,I+1)
      AIDIFU=WAI(J,I+1)
      AJDIFU=WAJ(J,I+1)
      B0DIFU=WB0(J,I+1)
      B1DIFU=WB1(J,I+1)
      B2DIFU=WB2(J,I+1)
      B3DIFU=WB3(J,I+1)
      B4DIFU=WB4(J,I+1)
      B5DIFU=WB5(J,I+1)
      B6DIFU=WB6(J,I+1)
      B7DIFU=WB7(J,I+1)
      B8DIFU=WB8(J,I+1)
      B9DIFU=WB9(J,I+1)
      BADIFU=WBA(J,I+1)
      BBDIFU=WBB(J,I+1)
      BCDIFU=WBC(J,I+1)
      BDDIFU=WBD(J,I+1)
      BEDIFU=WBE(J,I+1)
      BFDIFU=WBF(J,I+1)
      BGDIFU=WBG(J,I+1)
      BHDIFU=WBH(J,I+1)
      BIDIFU=WBI(J,I+1)
      BJDIFU=WBJ(J,I+1)
      C0DIFU=WC0(J,I+1)
      C1DIFU=WC1(J,I+1)
      C2DIFU=WC2(J,I+1)
      C3DIFU=WC3(J,I+1)
      C4DIFU=WC4(J,I+1)
      C5DIFU=WC5(J,I+1)
      C6DIFU=WC6(J,I+1)
      C7DIFU=WC7(J,I+1)
      C8DIFU=WC8(J,I+1)
      C9DIFU=WC9(J,I+1)
      CADIFU=WCA(J,I+1)
C-----------------------------------------------------------------------
      IF(ZDIFU.GE.ZCONV) ZCONV=ZDIFU
      IF(STDIFU.GE.STCONV) STCONV=STDIFU
      IF(A0DIFU.GE.A0CONV) A0CONV=A0DIFU
      IF(A1DIFU.GE.A1CONV) A1CONV=A1DIFU
      IF(A2DIFU.GE.A2CONV) A2CONV=A2DIFU
      IF(A3DIFU.GE.A3CONV) A3CONV=A3DIFU
      IF(A4DIFU.GE.A4CONV) A4CONV=A4DIFU
      IF(A5DIFU.GE.A5CONV) A5CONV=A5DIFU
      IF(A6DIFU.GE.A6CONV) A6CONV=A6DIFU
      IF(A7DIFU.GE.A7CONV) A7CONV=A7DIFU
      IF(A8DIFU.GE.A8CONV) A8CONV=A8DIFU
      IF(A9DIFU.GE.A9CONV) A9CONV=A9DIFU
      IF(AADIFU.GE.AACONV) AACONV=AADIFU
      IF(ABDIFU.GE.ABCONV) ABCONV=ABDIFU
      IF(ACDIFU.GE.ACCONV) ACCONV=ACDIFU
      IF(ADDIFU.GE.ADCONV) ADCONV=ADDIFU
      IF(AEDIFU.GE.AECONV) AECONV=AEDIFU
      IF(AFDIFU.GE.AFCONV) AFCONV=AFDIFU
      IF(AGDIFU.GE.AGCONV) AGCONV=AGDIFU
      IF(AHDIFU.GE.AHCONV) AHCONV=AHDIFU
      IF(AIDIFU.GE.AICONV) AICONV=AIDIFU
      IF(AJDIFU.GE.AJCONV) AJCONV=AJDIFU
      IF(B0DIFU.GE.B0CONV) B0CONV=B0DIFU
      IF(B1DIFU.GE.B1CONV) B1CONV=B1DIFU
      IF(B2DIFU.GE.B2CONV) B2CONV=B2DIFU
      IF(B3DIFU.GE.B3CONV) B3CONV=B3DIFU
      IF(B4DIFU.GE.B4CONV) B4CONV=B4DIFU
      IF(B5DIFU.GE.B5CONV) B5CONV=B5DIFU
      IF(B6DIFU.GE.B6CONV) B6CONV=B6DIFU
      IF(B7DIFU.GE.B7CONV) B7CONV=B7DIFU
      IF(B8DIFU.GE.B8CONV) B8CONV=B8DIFU
      IF(B9DIFU.GE.B9CONV) B9CONV=B9DIFU
      IF(BADIFU.GE.BACONV) BACONV=BADIFU
      IF(BBDIFU.GE.BBCONV) BBCONV=BBDIFU
      IF(BCDIFU.GE.BCCONV) BCCONV=BCDIFU
      IF(BDDIFU.GE.BDCONV) BDCONV=BDDIFU
      IF(BEDIFU.GE.BECONV) BECONV=BEDIFU
      IF(BFDIFU.GE.BFCONV) BFCONV=BFDIFU
      IF(BGDIFU.GE.BGCONV) BGCONV=BGDIFU
      IF(BHDIFU.GE.BHCONV) BHCONV=BHDIFU
      IF(BIDIFU.GE.BICONV) BICONV=BIDIFU
      IF(BJDIFU.GE.BJCONV) BJCONV=BJDIFU
      IF(C0DIFU.GE.C0CONV) C0CONV=C0DIFU
      IF(C1DIFU.GE.C1CONV) C1CONV=C1DIFU
      IF(C2DIFU.GE.C2CONV) C2CONV=C2DIFU
      IF(C3DIFU.GE.C3CONV) C3CONV=C3DIFU
      IF(C4DIFU.GE.C4CONV) C4CONV=C4DIFU
      IF(C5DIFU.GE.C5CONV) C5CONV=C5DIFU
      IF(C6DIFU.GE.C6CONV) C6CONV=C6DIFU
      IF(C7DIFU.GE.C7CONV) C7CONV=C7DIFU
      IF(C8DIFU.GE.C8CONV) C8CONV=C8DIFU
      IF(C9DIFU.GE.C9CONV) C9CONV=C9DIFU
      IF(CADIFU.GE.CACONV) CACONV=CADIFU
C-----------------------------------------------------------------------
      ZE(J,I)=DTBRX*0.5*(RHOUE-ZCONV)
      STE(J,I)=DTBRX*0.5*(RHOUE-STCONV)
      EA0(J,I)=DTBRX*0.5*(RHOUE-A0CONV)
      EA1(J,I)=DTBRX*0.5*(RHOUE-A1CONV)
      EA2(J,I)=DTBRX*0.5*(RHOUE-A2CONV)
      EA3(J,I)=DTBRX*0.5*(RHOUE-A3CONV)
      EA4(J,I)=DTBRX*0.5*(RHOUE-A4CONV)
      EA5(J,I)=DTBRX*0.5*(RHOUE-A5CONV)
      EA6(J,I)=DTBRX*0.5*(RHOUE-A6CONV)
      EA7(J,I)=DTBRX*0.5*(RHOUE-A7CONV)
      EA8(J,I)=DTBRX*0.5*(RHOUE-A8CONV)
      EA9(J,I)=DTBRX*0.5*(RHOUE-A9CONV)
      EAA(J,I)=DTBRX*0.5*(RHOUE-AACONV)
      EAB(J,I)=DTBRX*0.5*(RHOUE-ABCONV)
      EAC(J,I)=DTBRX*0.5*(RHOUE-ACCONV)
      EAD(J,I)=DTBRX*0.5*(RHOUE-ADCONV)
      EAE(J,I)=DTBRX*0.5*(RHOUE-AECONV)
      EAF(J,I)=DTBRX*0.5*(RHOUE-AFCONV)
      EAG(J,I)=DTBRX*0.5*(RHOUE-AGCONV)
      EAH(J,I)=DTBRX*0.5*(RHOUE-AHCONV)
      EAI(J,I)=DTBRX*0.5*(RHOUE-AICONV)
      EAJ(J,I)=DTBRX*0.5*(RHOUE-AJCONV)
      EB0(J,I)=DTBRX*0.5*(RHOUE-B0CONV)
      EB1(J,I)=DTBRX*0.5*(RHOUE-B1CONV)
      EB2(J,I)=DTBRX*0.5*(RHOUE-B2CONV)
      EB3(J,I)=DTBRX*0.5*(RHOUE-B3CONV)
      EB4(J,I)=DTBRX*0.5*(RHOUE-B4CONV)
      EB5(J,I)=DTBRX*0.5*(RHOUE-B5CONV)
      EB6(J,I)=DTBRX*0.5*(RHOUE-B6CONV)
      EB7(J,I)=DTBRX*0.5*(RHOUE-B7CONV)
      EB8(J,I)=DTBRX*0.5*(RHOUE-B8CONV)
      EB9(J,I)=DTBRX*0.5*(RHOUE-B9CONV)
      EBA(J,I)=DTBRX*0.5*(RHOUE-BACONV)
      EBB(J,I)=DTBRX*0.5*(RHOUE-BBCONV)
      EBC(J,I)=DTBRX*0.5*(RHOUE-BCCONV)
      EBD(J,I)=DTBRX*0.5*(RHOUE-BDCONV)
      EBE(J,I)=DTBRX*0.5*(RHOUE-BECONV)
      EBF(J,I)=DTBRX*0.5*(RHOUE-BFCONV)
      EBG(J,I)=DTBRX*0.5*(RHOUE-BGCONV)
      EBH(J,I)=DTBRX*0.5*(RHOUE-BHCONV)
      EBI(J,I)=DTBRX*0.5*(RHOUE-BICONV)
      EBJ(J,I)=DTBRX*0.5*(RHOUE-BJCONV)
      EC0(J,I)=DTBRX*0.5*(RHOUE-C0CONV)
      EC1(J,I)=DTBRX*0.5*(RHOUE-C1CONV)
      EC2(J,I)=DTBRX*0.5*(RHOUE-C2CONV)
      EC3(J,I)=DTBRX*0.5*(RHOUE-C3CONV)
      EC4(J,I)=DTBRX*0.5*(RHOUE-C4CONV)
      EC5(J,I)=DTBRX*0.5*(RHOUE-C5CONV)
      EC6(J,I)=DTBRX*0.5*(RHOUE-C6CONV)
      EC7(J,I)=DTBRX*0.5*(RHOUE-C7CONV)
      EC8(J,I)=DTBRX*0.5*(RHOUE-C8CONV)
      EC9(J,I)=DTBRX*0.5*(RHOUE-C9CONV)
      ECA(J,I)=DTBRX*0.5*(RHOUE-CACONV)
C-----------------------------------------------------------------------
      ZCONV=DABS(RHOUW)
      STCONV=ZCONV
      A0CONV=ZCONV
      A1CONV=ZCONV
      A2CONV=ZCONV
      A3CONV=ZCONV
      A4CONV=ZCONV
      A5CONV=ZCONV
      A6CONV=ZCONV
      A7CONV=ZCONV
      A8CONV=ZCONV
      A9CONV=ZCONV
      AACONV=ZCONV
      ABCONV=ZCONV
      ACCONV=ZCONV
      ADCONV=ZCONV
      AECONV=ZCONV
      AFCONV=ZCONV
      AGCONV=ZCONV
      AHCONV=ZCONV
      AICONV=ZCONV
      AJCONV=ZCONV
      B0CONV=ZCONV
      B1CONV=ZCONV
      B2CONV=ZCONV
      B3CONV=ZCONV
      B4CONV=ZCONV
      B5CONV=ZCONV
      B6CONV=ZCONV
      B7CONV=ZCONV
      B8CONV=ZCONV
      B9CONV=ZCONV
      BACONV=ZCONV
      BBCONV=ZCONV
      BCCONV=ZCONV
      BDCONV=ZCONV
      BECONV=ZCONV
      BFCONV=ZCONV
      BGCONV=ZCONV
      BHCONV=ZCONV
      BICONV=ZCONV
      BJCONV=ZCONV
      C0CONV=ZCONV
      C1CONV=ZCONV
      C2CONV=ZCONV
      C3CONV=ZCONV
      C4CONV=ZCONV
      C5CONV=ZCONV
      C6CONV=ZCONV
      C7CONV=ZCONV
      C8CONV=ZCONV
      C9CONV=ZCONV
      CACONV=ZCONV
C-----------------------------------------------------------------------
      ZDIFU=ZW(J,I)
      STDIFU=STW(J,I)
      A0DIFU=WA0(J,I)
      A1DIFU=WA1(J,I)
      A2DIFU=WA2(J,I)
      A3DIFU=WA3(J,I)
      A4DIFU=WA4(J,I)
      A5DIFU=WA5(J,I)
      A6DIFU=WA6(J,I)
      A7DIFU=WA7(J,I)
      A8DIFU=WA8(J,I)
      A9DIFU=WA9(J,I)
      AADIFU=WAA(J,I)
      ABDIFU=WAB(J,I)
      ACDIFU=WAC(J,I)
      ADDIFU=WAD(J,I)
      AEDIFU=WAE(J,I)
      AFDIFU=WAF(J,I)
      AGDIFU=WAG(J,I)
      AHDIFU=WAH(J,I)
      AIDIFU=WAI(J,I)
      AJDIFU=WAJ(J,I)
      B0DIFU=WB0(J,I)
      B1DIFU=WB1(J,I)
      B2DIFU=WB2(J,I)
      B3DIFU=WB3(J,I)
      B4DIFU=WB4(J,I)
      B5DIFU=WB5(J,I)
      B6DIFU=WB6(J,I)
      B7DIFU=WB7(J,I)
      B8DIFU=WB8(J,I)
      B9DIFU=WB9(J,I)
      BADIFU=WBA(J,I)
      BBDIFU=WBB(J,I)
      BCDIFU=WBC(J,I)
      BDDIFU=WBD(J,I)
      BEDIFU=WBE(J,I)
      BFDIFU=WBF(J,I)
      BGDIFU=WBG(J,I)
      BHDIFU=WBH(J,I)
      BIDIFU=WBI(J,I)
      BJDIFU=WBJ(J,I)
      C0DIFU=WC0(J,I)
      C1DIFU=WC1(J,I)
      C2DIFU=WC2(J,I)
      C3DIFU=WC3(J,I)
      C4DIFU=WC4(J,I)
      C5DIFU=WC5(J,I)
      C6DIFU=WC6(J,I)
      C7DIFU=WC7(J,I)
      C8DIFU=WC8(J,I)
      C9DIFU=WC9(J,I)
      CADIFU=WCA(J,I)
C-----------------------------------------------------------------------
      IF(ZDIFU.GE.ZCONV) ZCONV=ZDIFU
      IF(STDIFU.GE.ZCONV) STCONV=STDIFU
      IF(A0DIFU.GE.A0CONV) A0CONV=A0DIFU
      IF(A1DIFU.GE.A1CONV) A1CONV=A1DIFU
      IF(A2DIFU.GE.A2CONV) A2CONV=A2DIFU
      IF(A3DIFU.GE.A3CONV) A3CONV=A3DIFU
      IF(A4DIFU.GE.A4CONV) A4CONV=A4DIFU
      IF(A5DIFU.GE.A5CONV) A5CONV=A5DIFU
      IF(A6DIFU.GE.A6CONV) A6CONV=A6DIFU
      IF(A7DIFU.GE.A7CONV) A7CONV=A7DIFU
      IF(A8DIFU.GE.A8CONV) A8CONV=A8DIFU
      IF(A9DIFU.GE.A9CONV) A9CONV=A9DIFU
      IF(AADIFU.GE.AACONV) AACONV=AADIFU
      IF(ABDIFU.GE.ABCONV) ABCONV=ABDIFU
      IF(ACDIFU.GE.ACCONV) ACCONV=ACDIFU
      IF(ADDIFU.GE.ADCONV) ADCONV=ADDIFU
      IF(AEDIFU.GE.AECONV) AECONV=AEDIFU
      IF(AFDIFU.GE.AFCONV) AFCONV=AFDIFU
      IF(AGDIFU.GE.AGCONV) AGCONV=AGDIFU
      IF(AHDIFU.GE.AHCONV) AHCONV=AHDIFU
      IF(AIDIFU.GE.AICONV) AICONV=AIDIFU
      IF(AJDIFU.GE.AJCONV) AJCONV=AJDIFU
      IF(B0DIFU.GE.B0CONV) B0CONV=B0DIFU
      IF(B1DIFU.GE.B1CONV) B1CONV=B1DIFU
      IF(B2DIFU.GE.B2CONV) B2CONV=B2DIFU
      IF(B3DIFU.GE.B3CONV) B3CONV=B3DIFU
      IF(B4DIFU.GE.B4CONV) B4CONV=B4DIFU
      IF(B5DIFU.GE.B5CONV) B5CONV=B5DIFU
      IF(B6DIFU.GE.B6CONV) B6CONV=B6DIFU
      IF(B7DIFU.GE.B7CONV) B7CONV=B7DIFU
      IF(B8DIFU.GE.B8CONV) B8CONV=B8DIFU
      IF(B9DIFU.GE.B9CONV) B9CONV=B9DIFU
      IF(BADIFU.GE.BACONV) BACONV=BADIFU
      IF(BBDIFU.GE.BBCONV) BBCONV=BBDIFU
      IF(BCDIFU.GE.BCCONV) BCCONV=BCDIFU
      IF(BDDIFU.GE.BDCONV) BDCONV=BDDIFU
      IF(BEDIFU.GE.BECONV) BECONV=BEDIFU
      IF(BFDIFU.GE.BFCONV) BFCONV=BFDIFU
      IF(BGDIFU.GE.BGCONV) BGCONV=BGDIFU
      IF(BHDIFU.GE.BHCONV) BHCONV=BHDIFU
      IF(BIDIFU.GE.BICONV) BICONV=BIDIFU
      IF(BJDIFU.GE.BJCONV) BJCONV=BJDIFU
      IF(C0DIFU.GE.C0CONV) C0CONV=C0DIFU
      IF(C1DIFU.GE.C1CONV) C1CONV=C1DIFU
      IF(C2DIFU.GE.C2CONV) C2CONV=C2DIFU
      IF(C3DIFU.GE.C3CONV) C3CONV=C3DIFU
      IF(C4DIFU.GE.C4CONV) C4CONV=C4DIFU
      IF(C5DIFU.GE.C5CONV) C5CONV=C5DIFU
      IF(C6DIFU.GE.C6CONV) C6CONV=C6DIFU
      IF(C7DIFU.GE.C7CONV) C7CONV=C7DIFU
      IF(C8DIFU.GE.C8CONV) C8CONV=C8DIFU
      IF(C9DIFU.GE.C9CONV) C9CONV=C9DIFU
      IF(CADIFU.GE.CACONV) CACONV=CADIFU
C-----------------------------------------------------------------------
      ZW(J,I)=DTBRX*0.5*(-RHOUW-ZCONV)
      STW(J,I)=DTBRX*0.5*(-RHOUW-STCONV)
      WA0(J,I)=DTBRX*0.5*(-RHOUW-A0CONV)
      WA1(J,I)=DTBRX*0.5*(-RHOUW-A1CONV)
      WA2(J,I)=DTBRX*0.5*(-RHOUW-A2CONV)
      WA3(J,I)=DTBRX*0.5*(-RHOUW-A3CONV)
      WA4(J,I)=DTBRX*0.5*(-RHOUW-A4CONV)
      WA5(J,I)=DTBRX*0.5*(-RHOUW-A5CONV)
      WA6(J,I)=DTBRX*0.5*(-RHOUW-A6CONV)
      WA7(J,I)=DTBRX*0.5*(-RHOUW-A7CONV)
      WA8(J,I)=DTBRX*0.5*(-RHOUW-A8CONV)
      WA9(J,I)=DTBRX*0.5*(-RHOUW-A9CONV)
      WAA(J,I)=DTBRX*0.5*(-RHOUW-AACONV)
      WAB(J,I)=DTBRX*0.5*(-RHOUW-ABCONV)
      WAC(J,I)=DTBRX*0.5*(-RHOUW-ACCONV)
      WAD(J,I)=DTBRX*0.5*(-RHOUW-ADCONV)
      WAE(J,I)=DTBRX*0.5*(-RHOUW-AECONV)
      WAF(J,I)=DTBRX*0.5*(-RHOUW-AFCONV)
      WAG(J,I)=DTBRX*0.5*(-RHOUW-AGCONV)
      WAH(J,I)=DTBRX*0.5*(-RHOUW-AHCONV)
      WAI(J,I)=DTBRX*0.5*(-RHOUW-AICONV)
      WAJ(J,I)=DTBRX*0.5*(-RHOUW-AJCONV)
      WB0(J,I)=DTBRX*0.5*(-RHOUW-B0CONV)
      WB1(J,I)=DTBRX*0.5*(-RHOUW-B1CONV)
      WB2(J,I)=DTBRX*0.5*(-RHOUW-B2CONV)
      WB3(J,I)=DTBRX*0.5*(-RHOUW-B3CONV)
      WB4(J,I)=DTBRX*0.5*(-RHOUW-B4CONV)
      WB5(J,I)=DTBRX*0.5*(-RHOUW-B5CONV)
      WB6(J,I)=DTBRX*0.5*(-RHOUW-B6CONV)
      WB7(J,I)=DTBRX*0.5*(-RHOUW-B7CONV)
      WB8(J,I)=DTBRX*0.5*(-RHOUW-B8CONV)
      WB9(J,I)=DTBRX*0.5*(-RHOUW-B9CONV)
      WBA(J,I)=DTBRX*0.5*(-RHOUW-BACONV)
      WBB(J,I)=DTBRX*0.5*(-RHOUW-BBCONV)
      WBC(J,I)=DTBRX*0.5*(-RHOUW-BCCONV)
      WBD(J,I)=DTBRX*0.5*(-RHOUW-BDCONV)
      WBE(J,I)=DTBRX*0.5*(-RHOUW-BECONV)
      WBF(J,I)=DTBRX*0.5*(-RHOUW-BFCONV)
      WBG(J,I)=DTBRX*0.5*(-RHOUW-BGCONV)
      WBH(J,I)=DTBRX*0.5*(-RHOUW-BHCONV)
      WBI(J,I)=DTBRX*0.5*(-RHOUW-BICONV)
      WBJ(J,I)=DTBRX*0.5*(-RHOUW-BJCONV)
      WC0(J,I)=DTBRX*0.5*(-RHOUW-C0CONV)
      WC1(J,I)=DTBRX*0.5*(-RHOUW-C1CONV)
      WC2(J,I)=DTBRX*0.5*(-RHOUW-C2CONV)
      WC3(J,I)=DTBRX*0.5*(-RHOUW-C3CONV)
      WC4(J,I)=DTBRX*0.5*(-RHOUW-C4CONV)
      WC5(J,I)=DTBRX*0.5*(-RHOUW-C5CONV)
      WC6(J,I)=DTBRX*0.5*(-RHOUW-C6CONV)
      WC7(J,I)=DTBRX*0.5*(-RHOUW-C7CONV)
      WC8(J,I)=DTBRX*0.5*(-RHOUW-C8CONV)
      WC9(J,I)=DTBRX*0.5*(-RHOUW-C9CONV)
      WCA(J,I)=DTBRX*0.5*(-RHOUW-CACONV)
C-----------------------------------------------------------------------
      ZCONV=DABS(RHOVN)
      STCONV=ZCONV
      A0CONV=ZCONV
      A1CONV=ZCONV
      A2CONV=ZCONV
      A3CONV=ZCONV
      A4CONV=ZCONV
      A5CONV=ZCONV
      A6CONV=ZCONV
      A7CONV=ZCONV
      A8CONV=ZCONV
      A9CONV=ZCONV
      AACONV=ZCONV
      ABCONV=ZCONV
      ACCONV=ZCONV
      ADCONV=ZCONV
      AECONV=ZCONV
      AFCONV=ZCONV
      AGCONV=ZCONV
      AHCONV=ZCONV
      AICONV=ZCONV
      AJCONV=ZCONV
      B0CONV=ZCONV
      B1CONV=ZCONV
      B2CONV=ZCONV
      B3CONV=ZCONV
      B4CONV=ZCONV
      B5CONV=ZCONV
      B6CONV=ZCONV
      B7CONV=ZCONV
      B8CONV=ZCONV
      B9CONV=ZCONV
      BACONV=ZCONV
      BBCONV=ZCONV
      BCCONV=ZCONV
      BDCONV=ZCONV
      BECONV=ZCONV
      BFCONV=ZCONV
      BGCONV=ZCONV
      BHCONV=ZCONV
      BICONV=ZCONV
      BJCONV=ZCONV
      C0CONV=ZCONV
      C1CONV=ZCONV
      C2CONV=ZCONV
      C3CONV=ZCONV
      C4CONV=ZCONV
      C5CONV=ZCONV
      C6CONV=ZCONV
      C7CONV=ZCONV
      C8CONV=ZCONV
      C9CONV=ZCONV
      CACONV=ZCONV
C-----------------------------------------------------------------------
      ZDIFU=ZS(J+1,I)
      STDIFU=STS(J+1,I)
      A0DIFU=QA0(J+1,I)
      A1DIFU=QA1(J+1,I)
      A2DIFU=QA2(J+1,I)
      A3DIFU=QA3(J+1,I)
      A4DIFU=QA4(J+1,I)
      A5DIFU=QA5(J+1,I)
      A6DIFU=QA6(J+1,I)
      A7DIFU=QA7(J+1,I)
      A8DIFU=QA8(J+1,I)
      A9DIFU=QA9(J+1,I)
      AADIFU=QAA(J+1,I)
      ABDIFU=QAB(J+1,I)
      ACDIFU=QAC(J+1,I)
      ADDIFU=QAD(J+1,I)
      AEDIFU=QAE(J+1,I)
      AFDIFU=QAF(J+1,I)
      AGDIFU=QAG(J+1,I)
      AHDIFU=QAH(J+1,I)
      AIDIFU=QAI(J+1,I)
      AJDIFU=QAJ(J+1,I)
      B0DIFU=QB0(J+1,I)
      B1DIFU=QB1(J+1,I)
      B2DIFU=QB2(J+1,I)
      B3DIFU=QB3(J+1,I)
      B4DIFU=QB4(J+1,I)
      B5DIFU=QB5(J+1,I)
      B6DIFU=QB6(J+1,I)
      B7DIFU=QB7(J+1,I)
      B8DIFU=QB8(J+1,I)
      B9DIFU=QB9(J+1,I)
      BADIFU=QBA(J+1,I)
      BBDIFU=QBB(J+1,I)
      BCDIFU=QBC(J+1,I)
      BDDIFU=QBD(J+1,I)
      BEDIFU=QBE(J+1,I)
      BFDIFU=QBF(J+1,I)
      BGDIFU=QBG(J+1,I)
      BHDIFU=QBH(J+1,I)
      BIDIFU=QBI(J+1,I)
      BJDIFU=QBJ(J+1,I)
      C0DIFU=QC0(J+1,I)
      C1DIFU=QC1(J+1,I)
      C2DIFU=QC2(J+1,I)
      C3DIFU=QC3(J+1,I)
      C4DIFU=QC4(J+1,I)
      C5DIFU=QC5(J+1,I)
      C6DIFU=QC6(J+1,I)
      C7DIFU=QC7(J+1,I)
      C8DIFU=QC8(J+1,I)
      C9DIFU=QC9(J+1,I)
      CADIFU=QCA(J+1,I)
C-----------------------------------------------------------------------
      IF(ZDIFU.GE.ZCONV) ZCONV=ZDIFU
      IF(STDIFU.GE.STCONV) STCONV=STDIFU
      IF(A0DIFU.GE.A0CONV) A0CONV=A0DIFU
      IF(A1DIFU.GE.A1CONV) A1CONV=A1DIFU
      IF(A2DIFU.GE.A2CONV) A2CONV=A2DIFU
      IF(A3DIFU.GE.A3CONV) A3CONV=A3DIFU
      IF(A4DIFU.GE.A4CONV) A4CONV=A4DIFU
      IF(A5DIFU.GE.A5CONV) A5CONV=A5DIFU
      IF(A6DIFU.GE.A6CONV) A6CONV=A6DIFU
      IF(A7DIFU.GE.A7CONV) A7CONV=A7DIFU
      IF(A8DIFU.GE.A8CONV) A8CONV=A8DIFU
      IF(A9DIFU.GE.A9CONV) A9CONV=A9DIFU
      IF(AADIFU.GE.AACONV) AACONV=AADIFU
      IF(ABDIFU.GE.ABCONV) ABCONV=ABDIFU
      IF(ACDIFU.GE.ACCONV) ACCONV=ACDIFU
      IF(ADDIFU.GE.ADCONV) ADCONV=ADDIFU
      IF(AEDIFU.GE.AECONV) AECONV=AEDIFU
      IF(AFDIFU.GE.AFCONV) AFCONV=AFDIFU
      IF(AGDIFU.GE.AGCONV) AGCONV=AGDIFU
      IF(AHDIFU.GE.AHCONV) AHCONV=AHDIFU
      IF(AIDIFU.GE.AICONV) AICONV=AIDIFU
      IF(AJDIFU.GE.AJCONV) AJCONV=AJDIFU
      IF(B0DIFU.GE.B0CONV) B0CONV=B0DIFU
      IF(B1DIFU.GE.B1CONV) B1CONV=B1DIFU
      IF(B2DIFU.GE.B2CONV) B2CONV=B2DIFU
      IF(B3DIFU.GE.B3CONV) B3CONV=B3DIFU
      IF(B4DIFU.GE.B4CONV) B4CONV=B4DIFU
      IF(B5DIFU.GE.B5CONV) B5CONV=B5DIFU
      IF(B6DIFU.GE.B6CONV) B6CONV=B6DIFU
      IF(B7DIFU.GE.B7CONV) B7CONV=B7DIFU
      IF(B8DIFU.GE.B8CONV) B8CONV=B8DIFU
      IF(B9DIFU.GE.B9CONV) B9CONV=B9DIFU
      IF(BADIFU.GE.BACONV) BACONV=BADIFU
      IF(BBDIFU.GE.BBCONV) BBCONV=BBDIFU
      IF(BCDIFU.GE.BCCONV) BCCONV=BCDIFU
      IF(BDDIFU.GE.BDCONV) BDCONV=BDDIFU
      IF(BEDIFU.GE.BECONV) BECONV=BEDIFU
      IF(BFDIFU.GE.BFCONV) BFCONV=BFDIFU
      IF(BGDIFU.GE.BGCONV) BGCONV=BGDIFU
      IF(BHDIFU.GE.BHCONV) BHCONV=BHDIFU
      IF(BIDIFU.GE.BICONV) BICONV=BIDIFU
      IF(BJDIFU.GE.BJCONV) BJCONV=BJDIFU
      IF(C0DIFU.GE.C0CONV) C0CONV=C0DIFU
      IF(C1DIFU.GE.C1CONV) C1CONV=C1DIFU
      IF(C2DIFU.GE.C2CONV) C2CONV=C2DIFU
      IF(C3DIFU.GE.C3CONV) C3CONV=C3DIFU
      IF(C4DIFU.GE.C4CONV) C4CONV=C4DIFU
      IF(C5DIFU.GE.C5CONV) C5CONV=C5DIFU
      IF(C6DIFU.GE.C6CONV) C6CONV=C6DIFU
      IF(C7DIFU.GE.C7CONV) C7CONV=C7DIFU
      IF(C8DIFU.GE.C8CONV) C8CONV=C8DIFU
      IF(C9DIFU.GE.C9CONV) C9CONV=C9DIFU
      IF(CADIFU.GE.CACONV) CACONV=CADIFU
C-----------------------------------------------------------------------
      ZN(J,I)=DTBRY*0.5*(RHOVN-ZCONV)-SYMZ/DYN
      STN(J,I)=DTBRY*0.5*(RHOVN-STCONV)-SYMST/DYN
      PA0(J,I)=DTBRY*0.5*(RHOVN-A0CONV)-SYMA0/DYN
      PA1(J,I)=DTBRY*0.5*(RHOVN-A1CONV)-SYMA1/DYN
      PA2(J,I)=DTBRY*0.5*(RHOVN-A2CONV)-SYMA2/DYN
      PA3(J,I)=DTBRY*0.5*(RHOVN-A3CONV)-SYMA3/DYN
      PA4(J,I)=DTBRY*0.5*(RHOVN-A4CONV)-SYMA4/DYN
      PA5(J,I)=DTBRY*0.5*(RHOVN-A5CONV)-SYMA5/DYN
      PA6(J,I)=DTBRY*0.5*(RHOVN-A6CONV)-SYMA6/DYN
      PA7(J,I)=DTBRY*0.5*(RHOVN-A7CONV)-SYMA7/DYN
      PA8(J,I)=DTBRY*0.5*(RHOVN-A8CONV)-SYMA8/DYN
      PA9(J,I)=DTBRY*0.5*(RHOVN-A9CONV)-SYMA9/DYN
      PAA(J,I)=DTBRY*0.5*(RHOVN-AACONV)-SYMAA/DYN
      PAB(J,I)=DTBRY*0.5*(RHOVN-ABCONV)-SYMAB/DYN
      PAC(J,I)=DTBRY*0.5*(RHOVN-ACCONV)-SYMAC/DYN
      PAD(J,I)=DTBRY*0.5*(RHOVN-ADCONV)-SYMAD/DYN
      PAE(J,I)=DTBRY*0.5*(RHOVN-AECONV)-SYMAE/DYN
      PAF(J,I)=DTBRY*0.5*(RHOVN-AFCONV)-SYMAF/DYN
      PAG(J,I)=DTBRY*0.5*(RHOVN-AGCONV)-SYMAG/DYN
      PAH(J,I)=DTBRY*0.5*(RHOVN-AHCONV)-SYMAH/DYN
      PAI(J,I)=DTBRY*0.5*(RHOVN-AICONV)-SYMAI/DYN
      PAJ(J,I)=DTBRY*0.5*(RHOVN-AJCONV)-SYMAJ/DYN
      PB0(J,I)=DTBRY*0.5*(RHOVN-B0CONV)-SYMB0/DYN
      PB1(J,I)=DTBRY*0.5*(RHOVN-B1CONV)-SYMB1/DYN
      PB2(J,I)=DTBRY*0.5*(RHOVN-B2CONV)-SYMB2/DYN
      PB3(J,I)=DTBRY*0.5*(RHOVN-B3CONV)-SYMB3/DYN
      PB4(J,I)=DTBRY*0.5*(RHOVN-B4CONV)-SYMB4/DYN
      PB5(J,I)=DTBRY*0.5*(RHOVN-B5CONV)-SYMB5/DYN
      PB6(J,I)=DTBRY*0.5*(RHOVN-B6CONV)-SYMB6/DYN
      PB7(J,I)=DTBRY*0.5*(RHOVN-B7CONV)-SYMB7/DYN
      PB8(J,I)=DTBRY*0.5*(RHOVN-B8CONV)-SYMB8/DYN
      PB9(J,I)=DTBRY*0.5*(RHOVN-B9CONV)-SYMB9/DYN
      PBA(J,I)=DTBRY*0.5*(RHOVN-BACONV)-SYMBA/DYN
      PBB(J,I)=DTBRY*0.5*(RHOVN-BBCONV)-SYMBB/DYN
      PBC(J,I)=DTBRY*0.5*(RHOVN-BCCONV)-SYMBC/DYN
      PBD(J,I)=DTBRY*0.5*(RHOVN-BDCONV)-SYMBD/DYN
      PBE(J,I)=DTBRY*0.5*(RHOVN-BECONV)-SYMBE/DYN
      PBF(J,I)=DTBRY*0.5*(RHOVN-BFCONV)-SYMBF/DYN
      PBG(J,I)=DTBRY*0.5*(RHOVN-BGCONV)-SYMBG/DYN
      PBH(J,I)=DTBRY*0.5*(RHOVN-BHCONV)-SYMBH/DYN
      PBI(J,I)=DTBRY*0.5*(RHOVN-BICONV)-SYMBI/DYN
      PBJ(J,I)=DTBRY*0.5*(RHOVN-BJCONV)-SYMBJ/DYN
      PC0(J,I)=DTBRY*0.5*(RHOVN-C0CONV)-SYMC0/DYN
      PC1(J,I)=DTBRY*0.5*(RHOVN-C1CONV)-SYMC1/DYN
      PC2(J,I)=DTBRY*0.5*(RHOVN-C2CONV)-SYMC2/DYN
      PC3(J,I)=DTBRY*0.5*(RHOVN-C3CONV)-SYMC3/DYN
      PC4(J,I)=DTBRY*0.5*(RHOVN-C4CONV)-SYMC4/DYN
      PC5(J,I)=DTBRY*0.5*(RHOVN-C5CONV)-SYMC5/DYN
      PC6(J,I)=DTBRY*0.5*(RHOVN-C6CONV)-SYMC6/DYN
      PC7(J,I)=DTBRY*0.5*(RHOVN-C7CONV)-SYMC7/DYN
      PC8(J,I)=DTBRY*0.5*(RHOVN-C8CONV)-SYMC8/DYN
      PC9(J,I)=DTBRY*0.5*(RHOVN-C9CONV)-SYMC9/DYN
      PCA(J,I)=DTBRY*0.5*(RHOVN-CACONV)-SYMCA/DYN
C-----------------------------------------------------------------------
      ZCONV=DABS(RHOVS)
      STCONV=ZCONV
      A0CONV=ZCONV
      A1CONV=ZCONV
      A2CONV=ZCONV
      A3CONV=ZCONV
      A4CONV=ZCONV
      A5CONV=ZCONV
      A6CONV=ZCONV
      A7CONV=ZCONV
      A8CONV=ZCONV
      A9CONV=ZCONV
      AACONV=ZCONV
      ABCONV=ZCONV
      ACCONV=ZCONV
      ADCONV=ZCONV
      AECONV=ZCONV
      AFCONV=ZCONV
      AGCONV=ZCONV
      AHCONV=ZCONV
      AICONV=ZCONV
      AJCONV=ZCONV
      B0CONV=ZCONV
      B1CONV=ZCONV
      B2CONV=ZCONV
      B3CONV=ZCONV
      B4CONV=ZCONV
      B5CONV=ZCONV
      B6CONV=ZCONV
      B7CONV=ZCONV
      B8CONV=ZCONV
      B9CONV=ZCONV
      BACONV=ZCONV
      BBCONV=ZCONV
      BCCONV=ZCONV
      BDCONV=ZCONV
      BECONV=ZCONV
      BFCONV=ZCONV
      BGCONV=ZCONV
      BHCONV=ZCONV
      BICONV=ZCONV
      BJCONV=ZCONV
      C0CONV=ZCONV
      C1CONV=ZCONV
      C2CONV=ZCONV
      C3CONV=ZCONV
      C4CONV=ZCONV
      C5CONV=ZCONV
      C6CONV=ZCONV
      C7CONV=ZCONV
      C8CONV=ZCONV
      C9CONV=ZCONV
      CACONV=ZCONV
C-----------------------------------------------------------------------
      ZDIFU=ZS(J,I)
      STDIFU=STS(J,I)
      A0DIFU=QA0(J,I)
      A1DIFU=QA1(J,I)
      A2DIFU=QA2(J,I)
      A3DIFU=QA3(J,I)
      A4DIFU=QA4(J,I)
      A5DIFU=QA5(J,I)
      A6DIFU=QA6(J,I)
      A7DIFU=QA7(J,I)
      A8DIFU=QA8(J,I)
      A9DIFU=QA9(J,I)
      AADIFU=QAA(J,I)
      ABDIFU=QAB(J,I)
      ACDIFU=QAC(J,I)
      ADDIFU=QAD(J,I)
      AEDIFU=QAE(J,I)
      AFDIFU=QAF(J,I)
      AGDIFU=QAG(J,I)
      AHDIFU=QAH(J,I)
      AIDIFU=QAI(J,I)
      AJDIFU=QAJ(J,I)
      B0DIFU=QB0(J,I)
      B1DIFU=QB1(J,I)
      B2DIFU=QB2(J,I)
      B3DIFU=QB3(J,I)
      B4DIFU=QB4(J,I)
      B5DIFU=QB5(J,I)
      B6DIFU=QB6(J,I)
      B7DIFU=QB7(J,I)
      B8DIFU=QB8(J,I)
      B9DIFU=QB9(J,I)
      BADIFU=QBA(J,I)
      BBDIFU=QBB(J,I)
      BCDIFU=QBC(J,I)
      BDDIFU=QBD(J,I)
      BEDIFU=QBE(J,I)
      BFDIFU=QBF(J,I)
      BGDIFU=QBG(J,I)
      BHDIFU=QBH(J,I)
      BIDIFU=QBI(J,I)
      BJDIFU=QBJ(J,I)
      C0DIFU=QC0(J,I)
      C1DIFU=QC1(J,I)
      C2DIFU=QC2(J,I)
      C3DIFU=QC3(J,I)
      C4DIFU=QC4(J,I)
      C5DIFU=QC5(J,I)
      C6DIFU=QC6(J,I)
      C7DIFU=QC7(J,I)
      C8DIFU=QC8(J,I)
      C9DIFU=QC9(J,I)
      CADIFU=QCA(J,I)
C-----------------------------------------------------------------------
      IF(ZDIFU.GE.ZCONV) ZCONV=ZDIFU
      IF(STDIFU.GE.STCONV) STCONV=STDIFU
      IF(A0DIFU.GE.A0CONV) A0CONV=A0DIFU
      IF(A1DIFU.GE.A1CONV) A1CONV=A1DIFU
      IF(A2DIFU.GE.A2CONV) A2CONV=A2DIFU
      IF(A3DIFU.GE.A3CONV) A3CONV=A3DIFU
      IF(A4DIFU.GE.A4CONV) A4CONV=A4DIFU
      IF(A5DIFU.GE.A5CONV) A5CONV=A5DIFU
      IF(A6DIFU.GE.A6CONV) A6CONV=A6DIFU
      IF(A7DIFU.GE.A7CONV) A7CONV=A7DIFU
      IF(A8DIFU.GE.A8CONV) A8CONV=A8DIFU
      IF(A9DIFU.GE.A9CONV) A9CONV=A9DIFU
      IF(AADIFU.GE.AACONV) AACONV=AADIFU
      IF(ABDIFU.GE.ABCONV) ABCONV=ABDIFU
      IF(ACDIFU.GE.ACCONV) ACCONV=ACDIFU
      IF(ADDIFU.GE.ADCONV) ADCONV=ADDIFU
      IF(AEDIFU.GE.AECONV) AECONV=AEDIFU
      IF(AFDIFU.GE.AFCONV) AFCONV=AFDIFU
      IF(AGDIFU.GE.AGCONV) AGCONV=AGDIFU
      IF(AHDIFU.GE.AHCONV) AHCONV=AHDIFU
      IF(AIDIFU.GE.AICONV) AICONV=AIDIFU
      IF(AJDIFU.GE.AJCONV) AJCONV=AJDIFU
      IF(B0DIFU.GE.B0CONV) B0CONV=B0DIFU
      IF(B1DIFU.GE.B1CONV) B1CONV=B1DIFU
      IF(B2DIFU.GE.B2CONV) B2CONV=B2DIFU
      IF(B3DIFU.GE.B3CONV) B3CONV=B3DIFU
      IF(B4DIFU.GE.B4CONV) B4CONV=B4DIFU
      IF(B5DIFU.GE.B5CONV) B5CONV=B5DIFU
      IF(B6DIFU.GE.B6CONV) B6CONV=B6DIFU
      IF(B7DIFU.GE.B7CONV) B7CONV=B7DIFU
      IF(B8DIFU.GE.B8CONV) B8CONV=B8DIFU
      IF(B9DIFU.GE.B9CONV) B9CONV=B9DIFU
      IF(BADIFU.GE.BACONV) BACONV=BADIFU
      IF(BBDIFU.GE.BBCONV) BBCONV=BBDIFU
      IF(BCDIFU.GE.BCCONV) BCCONV=BCDIFU
      IF(BDDIFU.GE.BDCONV) BDCONV=BDDIFU
      IF(BEDIFU.GE.BECONV) BECONV=BEDIFU
      IF(BFDIFU.GE.BFCONV) BFCONV=BFDIFU
      IF(BGDIFU.GE.BGCONV) BGCONV=BGDIFU
      IF(BHDIFU.GE.BHCONV) BHCONV=BHDIFU
      IF(BIDIFU.GE.BICONV) BICONV=BIDIFU
      IF(BJDIFU.GE.BJCONV) BJCONV=BJDIFU
      IF(C0DIFU.GE.C0CONV) C0CONV=C0DIFU
      IF(C1DIFU.GE.C1CONV) C1CONV=C1DIFU
      IF(C2DIFU.GE.C2CONV) C2CONV=C2DIFU
      IF(C3DIFU.GE.C3CONV) C3CONV=C3DIFU
      IF(C4DIFU.GE.C4CONV) C4CONV=C4DIFU
      IF(C5DIFU.GE.C5CONV) C5CONV=C5DIFU
      IF(C6DIFU.GE.C6CONV) C6CONV=C6DIFU
      IF(C7DIFU.GE.C7CONV) C7CONV=C7DIFU
      IF(C8DIFU.GE.C8CONV) C8CONV=C8DIFU
      IF(C9DIFU.GE.C9CONV) C9CONV=C9DIFU
      IF(CADIFU.GE.CACONV) CACONV=CADIFU
C-----------------------------------------------------------------------
      ZS(J,I)=DTBRY*0.5*(-RHOVS-ZCONV)+SYMZ/DYS
      STS(J,I)=DTBRY*0.5*(-RHOVS-STCONV)+SYMST/DYS
      QA0(J,I)=DTBRY*0.5*(-RHOVS-A0CONV)+SYMA0/DYS
      QA1(J,I)=DTBRY*0.5*(-RHOVS-A1CONV)+SYMA1/DYS
      QA2(J,I)=DTBRY*0.5*(-RHOVS-A2CONV)+SYMA2/DYS
      QA3(J,I)=DTBRY*0.5*(-RHOVS-A3CONV)+SYMA3/DYS
      QA4(J,I)=DTBRY*0.5*(-RHOVS-A4CONV)+SYMA4/DYS
      QA5(J,I)=DTBRY*0.5*(-RHOVS-A5CONV)+SYMA5/DYS
      QA6(J,I)=DTBRY*0.5*(-RHOVS-A6CONV)+SYMA6/DYS
      QA7(J,I)=DTBRY*0.5*(-RHOVS-A7CONV)+SYMA7/DYS
      QA8(J,I)=DTBRY*0.5*(-RHOVS-A8CONV)+SYMA8/DYS
      QA9(J,I)=DTBRY*0.5*(-RHOVS-A9CONV)+SYMA9/DYS
      QAA(J,I)=DTBRY*0.5*(-RHOVS-AACONV)+SYMAA/DYS
      QAB(J,I)=DTBRY*0.5*(-RHOVS-ABCONV)+SYMAB/DYS
      QAC(J,I)=DTBRY*0.5*(-RHOVS-ACCONV)+SYMAC/DYS
      QAD(J,I)=DTBRY*0.5*(-RHOVS-ADCONV)+SYMAD/DYS
      QAE(J,I)=DTBRY*0.5*(-RHOVS-AECONV)+SYMAE/DYS
      QAF(J,I)=DTBRY*0.5*(-RHOVS-AFCONV)+SYMAF/DYS
      QAG(J,I)=DTBRY*0.5*(-RHOVS-AGCONV)+SYMAG/DYS
      QAH(J,I)=DTBRY*0.5*(-RHOVS-AHCONV)+SYMAH/DYS
      QAI(J,I)=DTBRY*0.5*(-RHOVS-AICONV)+SYMAI/DYS
      QAJ(J,I)=DTBRY*0.5*(-RHOVS-AJCONV)+SYMAJ/DYS
      QB0(J,I)=DTBRY*0.5*(-RHOVS-B0CONV)+SYMB0/DYS
      QB1(J,I)=DTBRY*0.5*(-RHOVS-B1CONV)+SYMB1/DYS
      QB2(J,I)=DTBRY*0.5*(-RHOVS-B2CONV)+SYMB2/DYS
      QB3(J,I)=DTBRY*0.5*(-RHOVS-B3CONV)+SYMB3/DYS
      QB4(J,I)=DTBRY*0.5*(-RHOVS-B4CONV)+SYMB4/DYS
      QB5(J,I)=DTBRY*0.5*(-RHOVS-B5CONV)+SYMB5/DYS
      QB6(J,I)=DTBRY*0.5*(-RHOVS-B6CONV)+SYMB6/DYS
      QB7(J,I)=DTBRY*0.5*(-RHOVS-B7CONV)+SYMB7/DYS
      QB8(J,I)=DTBRY*0.5*(-RHOVS-B8CONV)+SYMB8/DYS
      QB9(J,I)=DTBRY*0.5*(-RHOVS-B9CONV)+SYMB9/DYS
      QBA(J,I)=DTBRY*0.5*(-RHOVS-BACONV)+SYMBA/DYS
      QBB(J,I)=DTBRY*0.5*(-RHOVS-BBCONV)+SYMBB/DYS
      QBC(J,I)=DTBRY*0.5*(-RHOVS-BCCONV)+SYMBC/DYS
      QBD(J,I)=DTBRY*0.5*(-RHOVS-BDCONV)+SYMBD/DYS
      QBE(J,I)=DTBRY*0.5*(-RHOVS-BECONV)+SYMBE/DYS
      QBF(J,I)=DTBRY*0.5*(-RHOVS-BFCONV)+SYMBF/DYS
      QBG(J,I)=DTBRY*0.5*(-RHOVS-BGCONV)+SYMBG/DYS
      QBH(J,I)=DTBRY*0.5*(-RHOVS-BHCONV)+SYMBH/DYS
      QBI(J,I)=DTBRY*0.5*(-RHOVS-BICONV)+SYMBI/DYS
      QBJ(J,I)=DTBRY*0.5*(-RHOVS-BJCONV)+SYMBJ/DYS
      QC0(J,I)=DTBRY*0.5*(-RHOVS-C0CONV)+SYMC0/DYS
      QC1(J,I)=DTBRY*0.5*(-RHOVS-C1CONV)+SYMC1/DYS
      QC2(J,I)=DTBRY*0.5*(-RHOVS-C2CONV)+SYMC2/DYS
      QC3(J,I)=DTBRY*0.5*(-RHOVS-C3CONV)+SYMC3/DYS
      QC4(J,I)=DTBRY*0.5*(-RHOVS-C4CONV)+SYMC4/DYS
      QC5(J,I)=DTBRY*0.5*(-RHOVS-C5CONV)+SYMC5/DYS
      QC6(J,I)=DTBRY*0.5*(-RHOVS-C6CONV)+SYMC6/DYS
      QC7(J,I)=DTBRY*0.5*(-RHOVS-C7CONV)+SYMC7/DYS
      QC8(J,I)=DTBRY*0.5*(-RHOVS-C8CONV)+SYMC8/DYS
      QC9(J,I)=DTBRY*0.5*(-RHOVS-C9CONV)+SYMC9/DYS
      QCA(J,I)=DTBRY*0.5*(-RHOVS-CACONV)+SYMCA/DYS
C
   21 CONTINUE
      ENDDO
   20 CONTINUE
      IF(ICHEM.LE.0) GO TO 399
C--------------------------FIND REACTION RATES--------------------------
      CALL CKRATES
C-----------------SOLVE THE SPECIES-ALGEBRAIC EQUATIONS-----------------
      CALL SOLVESP(ISOR1,RELXSP,TOLRSP,SPNORM)
  399 CONTINUE
C------------------------SOLVING ENTHALPY EQUATION-----------------------
      IF(ITHRM.GT.0) CALL SOLVEH(ISOR2,RELXH,TOLRH,SIGH,SIGSP,HNORM)
      RESDH=0.0
      RESDC=0.0
      DO 600 I=2,LI-1
      DO J=2,LJ-1
      IF(ISKIP(J,I).GE.3) GO TO 601
      RESDH=RESDH+DABS(FPZ(J,I)-HT(J,I))
      RESDC=RESDC+DABS(SST1(J,I)-STMF(J,I))
C     RESDC=RESDC+DABS(SST2(J,I)-STND(J,I))
      RESDC=RESDC+DABS(SA0(J,I)-FA0(J,I))
      RESDC=RESDC+DABS(SA1(J,I)-FA1(J,I))
      RESDC=RESDC+DABS(SA2(J,I)-FA2(J,I))
      RESDC=RESDC+DABS(SA3(J,I)-FA3(J,I))
      RESDC=RESDC+DABS(SA4(J,I)-FA4(J,I))
      RESDC=RESDC+DABS(SA5(J,I)-FA5(J,I))
      RESDC=RESDC+DABS(SA6(J,I)-FA6(J,I))
      RESDC=RESDC+DABS(SA7(J,I)-FA7(J,I))
      RESDC=RESDC+DABS(SA8(J,I)-FA8(J,I))
      RESDC=RESDC+DABS(SA9(J,I)-FA9(J,I))
      RESDC=RESDC+DABS(SAA(J,I)-FAA(J,I))
      RESDC=RESDC+DABS(SAB(J,I)-FAB(J,I))
      RESDC=RESDC+DABS(SAC(J,I)-FAC(J,I))
      RESDC=RESDC+DABS(SAD(J,I)-FAD(J,I))
      RESDC=RESDC+DABS(SAE(J,I)-FAE(J,I))
      RESDC=RESDC+DABS(SAF(J,I)-FAF(J,I))
      RESDC=RESDC+DABS(SAG(J,I)-FAG(J,I))
      RESDC=RESDC+DABS(SAH(J,I)-FAH(J,I))
      RESDC=RESDC+DABS(SAI(J,I)-FAI(J,I))
      RESDC=RESDC+DABS(SAJ(J,I)-FAJ(J,I))
      RESDC=RESDC+DABS(SB0(J,I)-FB0(J,I))
      RESDC=RESDC+DABS(SB1(J,I)-FB1(J,I))
      RESDC=RESDC+DABS(SB2(J,I)-FB2(J,I))
      RESDC=RESDC+DABS(SB3(J,I)-FB3(J,I))
      RESDC=RESDC+DABS(SB4(J,I)-FB4(J,I))
      RESDC=RESDC+DABS(SB5(J,I)-FB5(J,I))
      RESDC=RESDC+DABS(SB6(J,I)-FB6(J,I))
      RESDC=RESDC+DABS(SB7(J,I)-FB7(J,I))
      RESDC=RESDC+DABS(SB8(J,I)-FB8(J,I))
      RESDC=RESDC+DABS(SB9(J,I)-FB9(J,I))
      RESDC=RESDC+DABS(SBA(J,I)-FBA(J,I))
      RESDC=RESDC+DABS(SBB(J,I)-FBB(J,I))
      RESDC=RESDC+DABS(SBC(J,I)-FBC(J,I))
      RESDC=RESDC+DABS(SBD(J,I)-FBD(J,I))
      RESDC=RESDC+DABS(SBE(J,I)-FBE(J,I))
      RESDC=RESDC+DABS(SBF(J,I)-FBF(J,I))
      RESDC=RESDC+DABS(SBG(J,I)-FBG(J,I))
      RESDC=RESDC+DABS(SBH(J,I)-FBH(J,I))
      RESDC=RESDC+DABS(SBI(J,I)-FBI(J,I))
      RESDC=RESDC+DABS(SBJ(J,I)-FBJ(J,I))
      RESDC=RESDC+DABS(SC0(J,I)-FC0(J,I))
      RESDC=RESDC+DABS(SC1(J,I)-FC1(J,I))
      RESDC=RESDC+DABS(SC2(J,I)-FC2(J,I))
      RESDC=RESDC+DABS(SC3(J,I)-FC3(J,I))
      RESDC=RESDC+DABS(SC4(J,I)-FC4(J,I))
      RESDC=RESDC+DABS(SC5(J,I)-FC5(J,I))
      RESDC=RESDC+DABS(SC6(J,I)-FC6(J,I))
      RESDC=RESDC+DABS(SC7(J,I)-FC7(J,I))
      RESDC=RESDC+DABS(SC8(J,I)-FC8(J,I))
      RESDC=RESDC+DABS(SC9(J,I)-FC9(J,I))
      RESDC=RESDC+DABS(SCA(J,I)-FCA(J,I))
      HT(J,I)=FPZ(J,I)
      STMF(J,I)=SST1(J,I)
      STND(J,I)=SST2(J,I)
      FA0(J,I)=SA0(J,I)
      FA1(J,I)=SA1(J,I)
      FA2(J,I)=SA2(J,I)
      FA3(J,I)=SA3(J,I)
      FA4(J,I)=SA4(J,I)
      FA5(J,I)=SA5(J,I)
      FA6(J,I)=SA6(J,I)
      FA7(J,I)=SA7(J,I)
      FA8(J,I)=SA8(J,I)
      FA9(J,I)=SA9(J,I)
      FAA(J,I)=SAA(J,I)
      FAB(J,I)=SAB(J,I)
      FAC(J,I)=SAC(J,I)
      FAD(J,I)=SAD(J,I)
      FAE(J,I)=SAE(J,I)
      FAF(J,I)=SAF(J,I)
      FAG(J,I)=SAG(J,I)
      FAH(J,I)=SAH(J,I)
      FAI(J,I)=SAI(J,I)
      FAJ(J,I)=SAJ(J,I)
      FB0(J,I)=SB0(J,I)
      FB1(J,I)=SB1(J,I)
      FB2(J,I)=SB2(J,I)
      FB3(J,I)=SB3(J,I)
      FB4(J,I)=SB4(J,I)
      FB5(J,I)=SB5(J,I)
      FB6(J,I)=SB6(J,I)
      FB7(J,I)=SB7(J,I)
      FB8(J,I)=SB8(J,I)
      FB9(J,I)=SB9(J,I)
      FBA(J,I)=SBA(J,I)
      FBB(J,I)=SBB(J,I)
      FBC(J,I)=SBC(J,I)
      FBD(J,I)=SBD(J,I)
      FBE(J,I)=SBE(J,I)
      FBF(J,I)=SBF(J,I)
      FBG(J,I)=SBG(J,I)
      FBH(J,I)=SBH(J,I)
      FBI(J,I)=SBI(J,I)
      FBJ(J,I)=SBJ(J,I)
      FC0(J,I)=SC0(J,I)
      FC1(J,I)=SC1(J,I)
      FC2(J,I)=SC2(J,I)
      FC3(J,I)=SC3(J,I)
      FC4(J,I)=SC4(J,I)
      FC5(J,I)=SC5(J,I)
      FC6(J,I)=SC6(J,I)
      FC7(J,I)=SC7(J,I)
      FC8(J,I)=SC8(J,I)
      FC9(J,I)=SC9(J,I)
      FCA(J,I)=SCA(J,I)
  601 CONTINUE
      ENDDO
  600 CONTINUE
      RESDH=RESDH/HNORM
      RESDC=RESDC/SPNORM
      RETURN
      END
C***********************************************************************
      SUBROUTINE CKRATES
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LJ2=LJ*2,LSP=52,LRX=544)
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FA0(LJ,LI),FA1(LJ,LI),FA2(LJ,LI),
     *  FA3(LJ,LI),FA4(LJ,LI),FA5(LJ,LI),FA6(LJ,LI),FA7(LJ,LI),
     *  FA8(LJ,LI),FA9(LJ,LI),FAA(LJ,LI),FAB(LJ,LI),FAC(LJ,LI),
     *  FAD(LJ,LI),FAE(LJ,LI),FAF(LJ,LI),FAG(LJ,LI),FAH(LJ,LI),
     *  FAI(LJ,LI),FAJ(LJ,LI),FB0(LJ,LI),FB1(LJ,LI),FB2(LJ,LI),
     *  FB3(LJ,LI),FB4(LJ,LI),FB5(LJ,LI),FB6(LJ,LI),FB7(LJ,LI),
     *  FB8(LJ,LI),FB9(LJ,LI),FBA(LJ,LI),FBB(LJ,LI),FBC(LJ,LI),
     *  FBD(LJ,LI),FBE(LJ,LI),FBF(LJ,LI),FBG(LJ,LI),FBH(LJ,LI),
     *  FBI(LJ,LI),FBJ(LJ,LI),FC0(LJ,LI),FC1(LJ,LI),FC2(LJ,LI),
     *  FC3(LJ,LI),FC4(LJ,LI),FC5(LJ,LI),FC6(LJ,LI),FC7(LJ,LI),
     *  FC8(LJ,LI),FC9(LJ,LI),FCA(LJ,LI),FCB(LJ,LI),
     *  U(LJ,LI),V(LJ,LI),W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),
     *  AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB08/ TPOL1,TPOL2,POLSP(14,LSP),CISP(12,LSP)
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
      COMMON/REAC/ LREV,NRLIND,ITBEF(LRX),ALF(LRX),AKF(LRX),EXAF(LRX),
     1             ALLOW(LRX),AKLOW(LRX),EALOW(LRX),TROE(LRX,4)
      COMMON/DUMR/ RC(LE,LRX)
      COMMON/DUMS/ STE(LJ,LI),STW(LJ,LI),STN(LJ,LI),STS(LJ,LI),
     1  SST1(LJ,LI),SST2(LJ,LI),SOUR(LJ,LI,6)
      DIMENSION SPRF(LRX),GSP(LSP)
C---------------------------Soot Model----------------------------------
C
C         STMF - Soot Mass Fraction (kg. of soot/Kg)
C         STND - Soot Number Density (No. of Particles/Kg)
C         
C ** For not including soot set SOUR(I,J,n) = 0.0 at the end of this sub       
C         
C---  Based on H. Guo, F. liu and J. Smallwood
C Soot and NO formation in counterflow ethylene/oxygen/nitrogen
C diffusion flames; Combustion Theory Modeling Vol. 8, pp. 475-489, 2004
C 
C   1) C2H2 --> 2 C(S) + H2  (a = 1.0D+04; L = 0; E = 20632)
C   2) C2H2 + n C(S) --> (n+2) C(S) + H2 (a = 3.468D+03; L = 0; E=12077)
C   3) 0.5 O2 + C(S) --> CO (a = 1.0D+04; L = 0.5; E = 19626)
C   4) OH + C(S) --> CO +H (a = 1.06D+02; L = -0.5; E = 0)
C   5) O + C(S) --> CO (a = 5.54D+01; L = -0.5; E = 0)
C                     ki = a T^L exp(-E/T)
C      r1 = k1 [C2H2]
C      r2 = k2 f(As) [C2H2]        f(As) = As
C      r3 = k3 As [O2]          
C      r4 = k4 Poh As Xoh        Poh = 0.13, Xoh - mole fraction
C      r5 = k5 Po As Xo          Po = 0.50, Xo - mole fraction
C-------------------------Constants in Soot Model-----------------------
      PI=3.1415927D+00
      RSOOT=1.9000D+03
      SIG0=5.669D-08
      SIGB=1.38D-23
C           (Avegadro Number 6.0232D+23 particles/g-mole)
      AVNO=6.0232D+23
      CMIN=100.0
      CAGL=9.0
C
      AKST1=1.000D+04
      ALST1=0.0
      EAST1=2.0632D+04
C
      AKST2=3.468D+03
      ALST2=0.0
      EAST2=1.2077D+04
C
      AKST3=1.000D+04
      ALST3=0.5
      EAST3=1.9626D+04
C
      AKST4=1.060D+02
      ALST4=-0.5
      EAST4=0.0
C
      AKST5=5.54D+00
      ALST5=-0.5
      EAST5=0.0
C-----------------------------------------------------------------------
      C1ST=ALSTR/VSTR
      C2ND=ALSTR/VSTR*(RSTR/WMSTR*1.0D-06)
      C3RD=ALSTR/VSTR*(RSTR*RSTR/WMSTR/WMSTR*1.0D-12)
C---------------------CALCULATING THE REACTION RATES--------------------
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
C                                                                     RR
C              Heptane-Air Kinetics (Sandiego Mechanism)              RR
C                                                                     RR
CRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
      DO 202 I=2,LI-1
      DO J=2,LJ-1
      IF(ISKIP(J,I).GE.3) GO TO 203
      TKD=TK(J,I)*TSTR
      ELOGT=DLOG(TKD)
C--------------------- FORWARD SPECIFIC REACTION RATES  ----------------
      DO 204 IR=1,LRX,LREV
      ACTV=ALF(IR)*ELOGT-EXAF(IR)/TKD
      SPRF(IR)=AKF(IR)*DEXP(ACTV)
  204 CONTINUE
C--------------------- REVERSE SPECIFIC REACTION RATES  ----------------
      IF(LREV.EQ.1) GO TO 206
      TKD1=TKD/2.0
      TKD2=TKD*TKD/6.0
      TKD3=TKD1*TKD2
      TKD4=1.2*TKD1*TKD3
      IPOLY=1
      IF(TKD.GT.TPOL1) IPOLY=8
      DO 205 IS=1,LSP
      GSP0=POLSP(IPOLY,IS)*(ELOGT-1.0)+POLSP(IPOLY+1,IS)*TKD1
     1 +POLSP(IPOLY+2,IS)*TKD2+POLSP(IPOLY+3,IS)*TKD3
     2 +POLSP(IPOLY+4,IS)*TKD4-POLSP(IPOLY+5,IS)/TKD+POLSP(IPOLY+6,IS)
      GSP(IS)=DEXP(GSP0)
  205 CONTINUE
      RUT=82.06*TKD
C----------BACKWARD RATES--------
      SPRF(   2)=AKF(   2)*SPRF(   1)
     1           /(GSP(  5)*GSP(  4)/GSP(  3)/GSP(  2))
      SPRF(   4)=AKF(   4)*SPRF(   3)
     1           /(GSP(  5)*GSP(  3)/GSP(  1)/GSP(  4))
      SPRF(   6)=AKF(   6)*SPRF(   5)
     1           /(GSP(  6)*GSP(  3)/GSP(  1)/GSP(  5))
      SPRF(   8)=AKF(   8)*SPRF(   7)
     1           /(GSP(  5)*GSP(  5)/GSP(  6)/GSP(  4))
      SPRF(  10)=AKF(  10)*SPRF(   9)
     1           /RUT/(GSP(  1)/GSP(  3)/GSP(  3))
      SPRF(  12)=AKF(  12)*SPRF(  11)
     1           /RUT/(GSP(  6)/GSP(  3)/GSP(  5))
      SPRF(  14)=AKF(  14)*SPRF(  13)
     1           /RUT/(GSP(  2)/GSP(  4)/GSP(  4))
      SPRF(  16)=AKF(  16)*SPRF(  15)
     1           /RUT/(GSP(  5)/GSP(  3)/GSP(  4))
      SPRF(  18)=AKF(  18)*SPRF(  17)
     1           /RUT/(GSP(  7)/GSP(  4)/GSP(  5))
      SPRF(  20)=AKF(  20)*SPRF(  19)
     1           /RUT/(GSP(  7)/GSP(  3)/GSP(  2))
      SPRF(  22)=AKF(  22)*SPRF(  21)
     1           /(GSP(  5)*GSP(  5)/GSP(  7)/GSP(  3))
      SPRF(  24)=AKF(  24)*SPRF(  23)
     1           /(GSP(  1)*GSP(  2)/GSP(  7)/GSP(  3))
      SPRF(  26)=AKF(  26)*SPRF(  25)
     1           /(GSP(  6)*GSP(  4)/GSP(  7)/GSP(  3))
      SPRF(  28)=AKF(  28)*SPRF(  27)
     1           /(GSP(  5)*GSP(  2)/GSP(  7)/GSP(  4))
      SPRF(  30)=AKF(  30)*SPRF(  29)
     1           /(GSP(  6)*GSP(  2)/GSP(  7)/GSP(  5))
      SPRF(  32)=AKF(  32)*SPRF(  31)
     1           /RUT/(GSP(  8)/GSP(  5)/GSP(  5))
      SPRF(  34)=AKF(  34)*SPRF(  33)
     1           /(GSP(  8)*GSP(  2)/GSP(  7)/GSP(  7))
      SPRF(  36)=AKF(  36)*SPRF(  35)
     1           /(GSP(  7)*GSP(  1)/GSP(  8)/GSP(  3))
      SPRF(  38)=AKF(  38)*SPRF(  37)
     1           /(GSP(  6)*GSP(  5)/GSP(  8)/GSP(  3))
      SPRF(  40)=AKF(  40)*SPRF(  39)
     1           /(GSP(  6)*GSP(  7)/GSP(  8)/GSP(  5))
      SPRF(  42)=AKF(  42)*SPRF(  41)
     1           /(GSP(  7)*GSP(  5)/GSP(  8)/GSP(  4))
      SPRF(  44)=AKF(  44)*SPRF(  43)
     1           /(GSP( 10)*GSP(  3)/GSP(  9)/GSP(  5))
      SPRF(  46)=AKF(  46)*SPRF(  45)
     1           /(GSP( 10)*GSP(  5)/GSP(  9)/GSP(  7))
      SPRF(  48)=AKF(  48)*SPRF(  47)
     1           /(GSP( 10)*GSP(  4)/GSP(  9)/GSP(  2))
      SPRF(  50)=AKF(  50)*SPRF(  49)
     1           *RUT/(GSP(  9)*GSP(  3)/GSP( 11))
      SPRF(  52)=AKF(  52)*SPRF(  51)
     1           /(GSP(  9)*GSP(  1)/GSP( 11)/GSP(  3))
      SPRF(  54)=AKF(  54)*SPRF(  53)
     1           /(GSP(  9)*GSP(  5)/GSP( 11)/GSP(  4))
      SPRF(  56)=AKF(  56)*SPRF(  55)
     1           /(GSP( 10)*GSP(  3)/GSP( 11)/GSP(  4))
      SPRF(  58)=AKF(  58)*SPRF(  57)
     1           /(GSP(  9)*GSP(  6)/GSP( 11)/GSP(  5))
      SPRF(  60)=AKF(  60)*SPRF(  59)
     1           /(GSP(  9)*GSP(  7)/GSP( 11)/GSP(  2))
      SPRF(  62)=AKF(  62)*SPRF(  61)
     1           /(GSP(  9)*GSP( 13)/GSP( 11)/GSP( 14))
      SPRF(  64)=AKF(  64)*SPRF(  63)
     1           /RUT/(GSP( 12)/GSP(  3)/GSP( 11))
      SPRF(  66)=AKF(  66)*SPRF(  65)
     1           /(GSP( 11)*GSP(  1)/GSP( 12)/GSP(  3))
      SPRF(  68)=AKF(  68)*SPRF(  67)
     1           /(GSP( 11)*GSP(  5)/GSP( 12)/GSP(  4))
      SPRF(  70)=AKF(  70)*SPRF(  69)
     1           /(GSP( 11)*GSP(  6)/GSP( 12)/GSP(  5))
      SPRF(  72)=AKF(  72)*SPRF(  71)
     1           /(GSP( 11)*GSP(  7)/GSP( 12)/GSP(  2))
      SPRF(  74)=AKF(  74)*SPRF(  73)
     1           /(GSP( 11)*GSP(  8)/GSP( 12)/GSP(  7))
      SPRF(  76)=AKF(  76)*SPRF(  75)
     1           /(GSP(  1)*GSP( 14)/GSP( 13)/GSP(  3))
      SPRF(  78)=AKF(  78)*SPRF(  77)
     1           /(GSP(  6)*GSP( 14)/GSP( 13)/GSP(  5))
      SPRF(  80)=AKF(  80)*SPRF(  79)
     1           /(GSP( 14)*GSP(  5)/GSP( 13)/GSP(  4))
      SPRF(  82)=AKF(  82)*SPRF(  81)
     1           /(GSP( 14)*GSP(  7)/GSP( 13)/GSP(  2))
      SPRF(  84)=AKF(  84)*SPRF(  83)
     1           /(GSP( 14)*GSP(  8)/GSP( 13)/GSP(  7))
      SPRF(  86)=AKF(  86)*SPRF(  85)
     1           /(GSP( 15)*GSP(  1)/GSP( 14)/GSP(  3))
      SPRF(  88)=AKF(  88)*SPRF(  87)
     1           /(GSP( 16)*GSP(  1)/GSP( 14)/GSP(  3))
      SPRF(  90)=AKF(  90)*SPRF(  89)
     1           /(GSP( 16)*GSP(  6)/GSP( 14)/GSP(  5))
      SPRF(  92)=AKF(  92)*SPRF(  91)
     1           /(GSP( 12)*GSP(  3)/GSP( 14)/GSP(  4))
      SPRF(  94)=AKF(  94)*SPRF(  93)
     1           /(GSP( 17)*GSP(  3)/GSP( 14)/GSP( 15))
      SPRF(  96)=AKF(  96)*SPRF(  95)
     1           /(GSP( 18)*GSP(  5)/GSP( 14)/GSP(  7))
      SPRF(  98)=AKF(  98)*SPRF(  97)
     1           /(GSP( 12)*GSP(  5)/GSP( 14)/GSP(  2))
      SPRF( 100)=AKF( 100)*SPRF(  99)
     1           /(GSP( 18)*GSP(  4)/GSP( 14)/GSP(  2))
      SPRF( 102)=AKF( 102)*SPRF( 101)
     1           /(GSP( 17)*GSP(  1)/GSP( 14)/GSP( 14))
      SPRF( 104)=AKF( 104)*SPRF( 103)
     1           /(GSP( 19)*GSP(  3)/GSP( 14)/GSP( 14))
      SPRF( 106)=AKF( 106)*SPRF( 105)
     1           /RUT/(GSP( 13)/GSP(  3)/GSP( 14))
      SPRF( 108)=AKF( 108)*SPRF( 107)
     1           /RUT/(GSP( 20)/GSP( 14)/GSP( 14))
      SPRF( 110)=AKF( 110)*SPRF( 109)
     1           /(GSP( 12)*GSP(  3)/GSP( 16)/GSP(  5))
      SPRF( 112)=AKF( 112)*SPRF( 111)
     1           *RUT/(GSP(  9)*GSP(  5)*GSP(  3)/GSP( 16)/GSP(  2))
      SPRF( 114)=AKF( 114)*SPRF( 113)
     1           /(GSP(  9)*GSP( 12)/GSP( 16)/GSP( 10))
      SPRF( 118)=AKF( 118)*SPRF( 117)
     1           /(GSP( 21)*GSP(  1)/GSP( 15)/GSP(  3))
      SPRF( 120)=AKF( 120)*SPRF( 119)
     1           /(GSP( 12)*GSP(  3)/GSP( 15)/GSP(  5))
      SPRF( 122)=AKF( 122)*SPRF( 121)
     1           /(GSP( 21)*GSP(  6)/GSP( 15)/GSP(  5))
      SPRF( 124)=AKF( 124)*SPRF( 123)
     1           *RUT/(GSP(  9)*GSP(  3)*GSP(  3)/GSP( 15)/GSP(  4))
      SPRF( 126)=AKF( 126)*SPRF( 125)
     1           /(GSP(  9)*GSP(  1)/GSP( 15)/GSP(  4))
      SPRF( 128)=AKF( 128)*SPRF( 127)
     1           /(GSP( 10)*GSP(  1)/GSP( 15)/GSP(  2))
      SPRF( 130)=AKF( 130)*SPRF( 129)
     1           *RUT/(GSP(  9)*GSP(  5)*GSP(  3)/GSP( 15)/GSP(  2))
      SPRF( 132)=AKF( 132)*SPRF( 131)
     1           *RUT/(GSP( 22)*GSP(  3)*GSP(  3)/GSP( 15)/GSP( 15))
      SPRF( 134)=AKF( 134)*SPRF( 133)
     1           /(GSP(  9)*GSP(  3)/GSP( 21)/GSP(  4))
      SPRF( 136)=AKF( 136)*SPRF( 135)
     1           /(GSP( 11)*GSP(  4)/GSP( 21)/GSP(  2))
      SPRF( 138)=AKF( 138)*SPRF( 137)
     1           /(GSP( 12)*GSP(  3)/GSP( 21)/GSP(  6))
      SPRF( 140)=AKF( 140)*SPRF( 139)
     1           /(GSP( 11)*GSP(  9)/GSP( 21)/GSP( 10))
      SPRF( 142)=AKF( 142)*SPRF( 141)
     1           /(GSP( 12)*GSP(  1)/GSP( 18)/GSP(  3))
      SPRF( 144)=AKF( 144)*SPRF( 143)
     1           /(GSP( 16)*GSP(  6)/GSP( 18)/GSP(  3))
      SPRF( 146)=AKF( 146)*SPRF( 145)
     1           /(GSP( 12)*GSP(  6)/GSP( 18)/GSP(  5))
      SPRF( 148)=AKF( 148)*SPRF( 147)
     1           /(GSP(  5)*GSP( 12)/GSP( 18)/GSP(  4))
      SPRF( 150)=AKF( 150)*SPRF( 149)
     1           /(GSP( 12)*GSP(  7)/GSP( 18)/GSP(  2))
      SPRF( 152)=AKF( 152)*SPRF( 151)
     1           *RUT/(GSP( 12)*GSP(  3)/GSP( 18))
      SPRF( 154)=AKF( 154)*SPRF( 153)
     1           /(GSP( 19)*GSP(  1)/GSP( 20)/GSP(  3))
      SPRF( 156)=AKF( 156)*SPRF( 155)
     1           /(GSP( 19)*GSP(  5)/GSP( 20)/GSP(  4))
      SPRF( 158)=AKF( 158)*SPRF( 157)
     1           /(GSP( 19)*GSP(  6)/GSP( 20)/GSP(  5))
      SPRF( 160)=AKF( 160)*SPRF( 159)
     1           /(GSP( 19)*GSP( 13)/GSP( 20)/GSP( 14))
      SPRF( 162)=AKF( 162)*SPRF( 161)
     1           *RUT/(GSP( 19)*GSP(  3)/GSP( 20))
      SPRF( 164)=AKF( 164)*SPRF( 163)
     1           /(GSP( 19)*GSP(  8)/GSP( 20)/GSP(  7))
      SPRF( 166)=AKF( 166)*SPRF( 165)
     1           /(GSP( 17)*GSP(  1)/GSP( 19)/GSP(  3))
      SPRF( 168)=AKF( 168)*SPRF( 167)
     1           /(GSP( 17)*GSP(  5)/GSP( 19)/GSP(  4))
      SPRF( 170)=AKF( 170)*SPRF( 169)
     1           /(GSP( 14)*GSP( 12)/GSP( 19)/GSP(  4))
      SPRF( 172)=AKF( 172)*SPRF( 171)
     1           /(GSP( 17)*GSP(  7)/GSP( 19)/GSP(  2))
      SPRF( 174)=AKF( 174)*SPRF( 173)
     1           *RUT/(GSP( 17)*GSP(  3)/GSP( 19))
      SPRF( 176)=AKF( 176)*SPRF( 175)
     1           /(GSP( 23)*GSP(  1)/GSP( 17)/GSP(  3))
      SPRF( 178)=AKF( 178)*SPRF( 177)
     1           /(GSP( 23)*GSP(  6)/GSP( 17)/GSP(  5))
      SPRF( 180)=AKF( 180)*SPRF( 179)
     1           /(GSP( 14)*GSP( 11)/GSP( 17)/GSP(  4))
      SPRF( 182)=AKF( 182)*SPRF( 181)
     1           /(GSP( 24)*GSP(  3)/GSP( 17)/GSP(  4))
      SPRF( 184)=AKF( 184)*SPRF( 183)
     1           /(GSP( 23)*GSP( 19)/GSP( 17)/GSP( 17))
      SPRF( 186)=AKF( 186)*SPRF( 185)
     1           /(GSP( 23)*GSP(  7)/GSP( 17)/GSP(  2))
      SPRF( 188)=AKF( 188)*SPRF( 187)
     1           /(GSP( 25)*GSP(  5)/GSP( 17)/GSP(  7))
      SPRF( 190)=AKF( 190)*SPRF( 189)
     1           *RUT/(GSP( 14)*GSP(  9)*GSP(  8)/GSP( 25)/GSP(  7))
      SPRF( 192)=AKF( 192)*SPRF( 191)
     1           *RUT/(GSP( 23)*GSP(  3)/GSP( 17))
      SPRF( 194)=AKF( 194)*SPRF( 193)
     1           *RUT/(GSP( 22)*GSP(  1)/GSP( 17))
      SPRF( 196)=AKF( 196)*SPRF( 195)
     1           /(GSP( 22)*GSP(  1)/GSP( 23)/GSP(  3))
      SPRF( 198)=AKF( 198)*SPRF( 197)
     1           *RUT/(GSP( 22)*GSP(  3)/GSP( 23))
      SPRF( 200)=AKF( 200)*SPRF( 199)
     1           /(GSP( 12)*GSP( 11)/GSP( 23)/GSP(  2))
      SPRF( 202)=AKF( 202)*SPRF( 201)
     1           /(GSP( 24)*GSP(  4)/GSP( 23)/GSP(  2))
      SPRF( 204)=AKF( 204)*SPRF( 203)
     1           /(GSP( 22)*GSP(  7)/GSP( 23)/GSP(  2))
      SPRF( 206)=AKF( 206)*SPRF( 205)
     1           /(GSP( 27)*GSP(  3)/GSP( 22)/GSP(  4))
      SPRF( 208)=AKF( 208)*SPRF( 207)
     1           /(GSP( 15)*GSP(  9)/GSP( 22)/GSP(  4))
      SPRF( 210)=AKF( 210)*SPRF( 209)
     1           /(GSP( 12)*GSP(  9)/GSP( 22)/GSP(  2))
      SPRF( 212)=AKF( 212)*SPRF( 211)
     1           /(GSP( 26)*GSP(  3)/GSP( 22)/GSP(  5))
      SPRF( 214)=AKF( 214)*SPRF( 213)
     1           /(GSP( 28)*GSP(  6)/GSP( 22)/GSP(  5))
      SPRF( 216)=AKF( 216)*SPRF( 215)
     1           /(GSP( 14)*GSP(  9)/GSP( 26)/GSP(  3))
      SPRF( 218)=AKF( 218)*SPRF( 217)
     1           /(GSP( 15)*GSP( 10)/GSP( 26)/GSP(  4))
      SPRF( 220)=AKF( 220)*SPRF( 219)
     1           /(GSP( 27)*GSP(  5)/GSP( 26)/GSP(  4))
      SPRF( 222)=AKF( 222)*SPRF( 221)
     1           /(GSP( 19)*GSP(  9)/GSP( 26)/GSP( 14))
      SPRF( 224)=AKF( 224)*SPRF( 223)
     1           /(GSP( 16)*GSP(  9)/GSP( 27)/GSP(  3))
      SPRF( 226)=AKF( 226)*SPRF( 225)
     1           *RUT/(GSP( 11)*GSP(  9)*GSP(  3)/GSP( 27)/GSP(  5))
      SPRF( 228)=AKF( 228)*SPRF( 227)
     1           *RUT/(GSP(  9)*GSP(  9)*GSP(  3)/GSP( 27)/GSP(  4))
      SPRF( 230)=AKF( 230)*SPRF( 229)
     1           *RUT/(GSP(  9)*GSP(  9)*GSP(  5)/GSP( 27)/GSP(  2))
      SPRF( 232)=AKF( 232)*SPRF( 231)
     1           *RUT/(GSP( 10)*GSP(  9)*GSP(  3)/GSP( 27)/GSP(  2))
      SPRF( 234)=AKF( 234)*SPRF( 233)
     1           /(GSP( 27)*GSP(  3)/GSP( 28)/GSP(  5))
      SPRF( 236)=AKF( 236)*SPRF( 235)
     1           /(GSP(  9)*GSP( 21)/GSP( 28)/GSP(  4))
      SPRF( 238)=AKF( 238)*SPRF( 237)
     1           /(GSP( 27)*GSP(  4)/GSP( 28)/GSP(  2))
      SPRF( 240)=AKF( 240)*SPRF( 239)
     1           /(GSP( 21)*GSP( 10)/GSP( 28)/GSP(  2))
      SPRF( 242)=AKF( 242)*SPRF( 241)
     1           /(GSP( 11)*GSP(  9)/GSP( 28)/GSP(  2))
      SPRF( 244)=AKF( 244)*SPRF( 243)
     1           /(GSP( 12)*GSP(  1)/GSP( 29)/GSP(  3))
      SPRF( 246)=AKF( 246)*SPRF( 245)
     1           /(GSP( 14)*GSP(  5)/GSP( 29)/GSP(  3))
      SPRF( 248)=AKF( 248)*SPRF( 247)
     1           /(GSP( 12)*GSP(  6)/GSP( 29)/GSP(  5))
      SPRF( 250)=AKF( 250)*SPRF( 249)
     1           /(GSP( 12)*GSP(  7)/GSP( 29)/GSP(  2))
      SPRF( 252)=AKF( 252)*SPRF( 251)
     1           *RUT/(GSP( 12)*GSP(  3)/GSP( 29))
      SPRF( 256)=AKF( 256)*SPRF( 255)
     1           /(GSP( 29)*GSP(  9)/GSP( 26)/GSP(  5))
      SPRF( 258)=AKF( 258)*SPRF( 257)
     1           /(GSP( 29)*GSP(  6)/GSP( 30)/GSP(  5))
      SPRF( 260)=AKF( 260)*SPRF( 259)
     1           /(GSP( 18)*GSP(  6)/GSP( 30)/GSP(  5))
      SPRF( 262)=AKF( 262)*SPRF( 261)
     1           /(GSP( 29)*GSP(  1)/GSP( 30)/GSP(  3))
      SPRF( 264)=AKF( 264)*SPRF( 263)
     1           /(GSP( 18)*GSP(  1)/GSP( 30)/GSP(  3))
      SPRF( 266)=AKF( 266)*SPRF( 265)
     1           /(GSP( 29)*GSP(  5)/GSP( 30)/GSP(  4))
      SPRF( 268)=AKF( 268)*SPRF( 267)
     1           /(GSP( 29)*GSP(  8)/GSP( 30)/GSP(  7))
      SPRF( 270)=AKF( 270)*SPRF( 269)
     1           /(GSP( 29)*GSP(  7)/GSP( 30)/GSP(  2))
      SPRF( 272)=AKF( 272)*SPRF( 271)
     1           *RUT/(GSP( 14)*GSP(  5)/GSP( 30))
      SPRF( 274)=AKF( 274)*SPRF( 273)
     1           *RUT/(GSP( 26)*GSP(  3)/GSP( 24))
      SPRF( 276)=AKF( 276)*SPRF( 275)
     1           /(GSP( 14)*GSP( 11)/GSP( 24)/GSP(  3))
      SPRF( 278)=AKF( 278)*SPRF( 277)
     1           /(GSP( 26)*GSP(  1)/GSP( 24)/GSP(  3))
      SPRF( 280)=AKF( 280)*SPRF( 279)
     1           /(GSP( 12)*GSP( 11)/GSP( 24)/GSP(  4))
      SPRF( 282)=AKF( 282)*SPRF( 281)
     1           /(GSP( 26)*GSP(  6)/GSP( 24)/GSP(  5))
      SPRF( 284)=AKF( 284)*SPRF( 283)
     1           *RUT/(GSP( 12)*GSP(  9)*GSP(  5)/GSP( 24)/GSP(  2))
      SPRF( 286)=AKF( 286)*SPRF( 285)
     1           *RUT/(GSP( 19)*GSP(  9)*GSP(  3)/GSP( 24)/GSP( 14))
      SPRF( 288)=AKF( 288)*SPRF( 287)
     1           *RUT/(GSP( 12)*GSP( 11)*GSP(  5)/GSP( 24)/GSP(  7))
      SPRF( 290)=AKF( 290)*SPRF( 289)
     1           /(GSP( 32)*GSP(  2)/GSP( 24)/GSP(  7))
      SPRF( 292)=AKF( 292)*SPRF( 291)
     1           *RUT/(GSP( 14)*GSP(  9)/GSP( 24))
      SPRF( 294)=AKF( 294)*SPRF( 293)
     1           *RUT/(GSP( 14)*GSP( 11)/GSP( 32))
      SPRF( 296)=AKF( 296)*SPRF( 295)
     1           *RUT/(GSP( 14)*GSP(  9)/GSP( 35))
      SPRF( 298)=AKF( 298)*SPRF( 297)
     1           /(GSP( 35)*GSP(  6)/GSP( 32)/GSP(  5))
      SPRF( 300)=AKF( 300)*SPRF( 299)
     1           /(GSP( 24)*GSP(  6)/GSP( 32)/GSP(  5))
      SPRF( 302)=AKF( 302)*SPRF( 301)
     1           /(GSP( 35)*GSP(  5)/GSP( 32)/GSP(  4))
      SPRF( 304)=AKF( 304)*SPRF( 303)
     1           /(GSP( 24)*GSP(  5)/GSP( 32)/GSP(  4))
      SPRF( 306)=AKF( 306)*SPRF( 305)
     1           /(GSP( 35)*GSP(  1)/GSP( 32)/GSP(  3))
      SPRF( 308)=AKF( 308)*SPRF( 307)
     1           /(GSP( 24)*GSP(  1)/GSP( 32)/GSP(  3))
      SPRF( 310)=AKF( 310)*SPRF( 309)
     1           /(GSP( 35)*GSP( 13)/GSP( 32)/GSP( 14))
      SPRF( 312)=AKF( 312)*SPRF( 311)
     1           /(GSP( 24)*GSP( 13)/GSP( 32)/GSP( 14))
      SPRF( 314)=AKF( 314)*SPRF( 313)
     1           /(GSP( 35)*GSP(  8)/GSP( 32)/GSP(  7))
      SPRF( 316)=AKF( 316)*SPRF( 315)
     1           /(GSP( 24)*GSP(  8)/GSP( 32)/GSP(  7))
      SPRF( 318)=AKF( 318)*SPRF( 317)
     1           /(GSP( 35)*GSP(  7)/GSP( 32)/GSP(  2))
      SPRF( 320)=AKF( 320)*SPRF( 319)
     1           *RUT/(GSP( 14)*GSP( 29)/GSP( 31))
      SPRF( 322)=AKF( 322)*SPRF( 321)
     1           *RUT/(GSP( 17)*GSP(  6)/GSP( 31))
      SPRF( 324)=AKF( 324)*SPRF( 323)
     1           /(GSP( 34)*GSP(  6)/GSP( 31)/GSP(  5))
      SPRF( 326)=AKF( 326)*SPRF( 325)
     1           /(GSP( 33)*GSP(  6)/GSP( 31)/GSP(  5))
      SPRF( 328)=AKF( 328)*SPRF( 327)
     1           /(GSP( 36)*GSP(  6)/GSP( 31)/GSP(  5))
      SPRF( 330)=AKF( 330)*SPRF( 329)
     1           /(GSP( 34)*GSP(  1)/GSP( 31)/GSP(  3))
      SPRF( 332)=AKF( 332)*SPRF( 331)
     1           /(GSP( 33)*GSP(  1)/GSP( 31)/GSP(  3))
      SPRF( 334)=AKF( 334)*SPRF( 333)
     1           /(GSP( 36)*GSP(  1)/GSP( 31)/GSP(  3))
      SPRF( 336)=AKF( 336)*SPRF( 335)
     1           /(GSP( 34)*GSP(  5)/GSP( 31)/GSP(  4))
      SPRF( 338)=AKF( 338)*SPRF( 337)
     1           /(GSP( 33)*GSP(  5)/GSP( 31)/GSP(  4))
      SPRF( 340)=AKF( 340)*SPRF( 339)
     1           /(GSP( 36)*GSP(  5)/GSP( 31)/GSP(  4))
      SPRF( 342)=AKF( 342)*SPRF( 341)
     1           /(GSP( 34)*GSP( 13)/GSP( 31)/GSP( 14))
      SPRF( 344)=AKF( 344)*SPRF( 343)
     1           /(GSP( 33)*GSP( 13)/GSP( 31)/GSP( 14))
      SPRF( 346)=AKF( 346)*SPRF( 345)
     1           /(GSP( 36)*GSP( 13)/GSP( 31)/GSP( 14))
      SPRF( 348)=AKF( 348)*SPRF( 347)
     1           /(GSP( 33)*GSP(  8)/GSP( 31)/GSP(  7))
      SPRF( 350)=AKF( 350)*SPRF( 349)
     1           /(GSP( 34)*GSP(  8)/GSP( 31)/GSP(  7))
      SPRF( 352)=AKF( 352)*SPRF( 351)
     1           /(GSP( 36)*GSP(  8)/GSP( 31)/GSP(  7))
      SPRF( 354)=AKF( 354)*SPRF( 353)
     1           /RUT/(GSP( 34)/GSP( 17)/GSP(  5))
      SPRF( 356)=AKF( 356)*SPRF( 355)
     1           /(GSP( 36)*GSP(  5)/GSP( 19)/GSP(  7))
      SPRF( 358)=AKF( 358)*SPRF( 357)
     1           *RUT/(GSP( 32)*GSP(  3)/GSP( 36))
      SPRF( 360)=AKF( 360)*SPRF( 359)
     1           *RUT/(GSP( 14)*GSP( 12)/GSP( 36))
      SPRF( 362)=AKF( 362)*SPRF( 361)
     1           /(GSP( 32)*GSP(  7)/GSP( 36)/GSP(  2))
      SPRF( 364)=AKF( 364)*SPRF( 363)
     1           /(GSP( 19)*GSP( 10)/GSP( 36)/GSP(  9))
      SPRF( 366)=AKF( 366)*SPRF( 365)
     1           /(GSP( 14)*GSP( 29)/GSP( 36)/GSP(  3))
      SPRF( 368)=AKF( 368)*SPRF( 367)
     1           /(GSP( 17)*GSP(  6)/GSP( 36)/GSP(  3))
      SPRF( 370)=AKF( 370)*SPRF( 369)
     1           /(GSP( 32)*GSP(  6)/GSP( 36)/GSP(  5))
      SPRF( 372)=AKF( 372)*SPRF( 371)
     1           /(GSP( 32)*GSP(  7)/GSP( 33)/GSP(  2))
      SPRF( 374)=AKF( 374)*SPRF( 373)
     1           /(GSP( 32)*GSP(  5)/GSP( 33)/GSP(  4))
      SPRF( 376)=AKF( 376)*SPRF( 375)
     1           /(GSP( 17)*GSP(  6)/GSP( 33)/GSP(  3))
      SPRF( 378)=AKF( 378)*SPRF( 377)
     1           /(GSP( 14)*GSP( 29)/GSP( 33)/GSP(  3))
      SPRF( 380)=AKF( 380)*SPRF( 379)
     1           *RUT/(GSP( 32)*GSP(  5)*GSP(  5)/GSP( 33)/GSP(  7))
      SPRF( 382)=AKF( 382)*SPRF( 381)
     1           /(GSP( 32)*GSP(  6)/GSP( 33)/GSP(  5))
      SPRF( 384)=AKF( 384)*SPRF( 383)
     1           *RUT/(GSP( 32)*GSP(  3)/GSP( 33))
      SPRF( 386)=AKF( 386)*SPRF( 385)
     1           /(GSP( 17)*GSP(  9)/GSP( 37)/GSP(  4))
      SPRF( 388)=AKF( 388)*SPRF( 387)
     1           /(GSP( 37)*GSP(  3)/GSP( 14)/GSP( 22))
      SPRF( 390)=AKF( 390)*SPRF( 389)
     1           /(GSP( 27)*GSP( 14)/GSP( 37)/GSP(  4))
      SPRF( 392)=AKF( 392)*SPRF( 391)
     1           /RUT/(GSP( 37)/GSP( 38)/GSP(  3))
      SPRF( 394)=AKF( 394)*SPRF( 393)
     1           /(GSP( 37)*GSP(  2)/GSP( 38)/GSP(  7))
      SPRF( 396)=AKF( 396)*SPRF( 395)
     1           /(GSP( 38)*GSP(  6)/GSP( 37)/GSP(  5))
      SPRF( 398)=AKF( 398)*SPRF( 397)
     1           /(GSP( 26)*GSP( 11)/GSP( 38)/GSP(  2))
      SPRF( 400)=AKF( 400)*SPRF( 399)
     1           /RUT/(GSP( 39)/GSP( 37)/GSP(  3))
      SPRF( 402)=AKF( 402)*SPRF( 401)
     1           /(GSP( 37)*GSP(  1)/GSP( 39)/GSP(  3))
      SPRF( 404)=AKF( 404)*SPRF( 403)
     1           /(GSP( 37)*GSP(  7)/GSP( 39)/GSP(  2))
      SPRF( 406)=AKF( 406)*SPRF( 405)
     1           /(GSP( 37)*GSP( 13)/GSP( 39)/GSP( 14))
      SPRF( 408)=AKF( 408)*SPRF( 407)
     1           /RUT/(GSP( 39)/GSP( 22)/GSP( 14))
      SPRF( 410)=AKF( 410)*SPRF( 409)
     1           /(GSP( 37)*GSP(  6)/GSP( 39)/GSP(  5))
      SPRF( 412)=AKF( 412)*SPRF( 411)
     1           /(GSP( 37)*GSP(  9)/GSP( 38)/GSP( 11))
      SPRF( 414)=AKF( 414)*SPRF( 413)
     1           *RUT/(GSP(  5)*GSP(  9)*GSP( 23)/GSP( 38)/GSP(  7))
      SPRF( 416)=AKF( 416)*SPRF( 415)
     1           *RUT/(GSP( 14)*GSP( 11)*GSP(  9)/GSP( 37)/GSP(  2))
      SPRF( 418)=AKF( 418)*SPRF( 417)
     1           /(GSP( 19)*GSP( 11)/GSP( 40)/GSP(  4))
      SPRF( 420)=AKF( 420)*SPRF( 419)
     1           /(GSP( 39)*GSP(  6)/GSP( 40)/GSP(  5))
      SPRF( 422)=AKF( 422)*SPRF( 421)
     1           *RUT/(GSP( 26)*GSP( 14)*GSP(  3)/GSP( 40)/GSP(  4))
      SPRF( 424)=AKF( 424)*SPRF( 423)
     1           /(GSP( 39)*GSP(  1)/GSP( 40)/GSP(  3))
      SPRF( 426)=AKF( 426)*SPRF( 425)
     1           /RUT/(GSP( 40)/GSP( 39)/GSP(  3))
      SPRF( 428)=AKF( 428)*SPRF( 427)
     1           /(GSP( 40)*GSP(  2)/GSP( 39)/GSP(  7))
      SPRF( 430)=AKF( 430)*SPRF( 429)
     1           *RUT/(GSP(  5)*GSP( 23)*GSP( 12)/GSP( 39)/GSP(  7))
      SPRF( 432)=AKF( 432)*SPRF( 431)
     1           /RUT/(GSP( 40)/GSP( 23)/GSP( 14))
      SPRF( 434)=AKF( 434)*SPRF( 433)
     1           /(GSP( 17)*GSP( 14)/GSP( 40)/GSP(  3))
      SPRF( 436)=AKF( 436)*SPRF( 435)
     1           /(GSP( 39)*GSP(  3)/GSP( 14)/GSP( 23))
      SPRF( 438)=AKF( 438)*SPRF( 437)
     1           *RUT/(GSP( 14)*GSP( 19)/GSP( 41))
      SPRF( 440)=AKF( 440)*SPRF( 439)
     1           /(GSP( 42)*GSP(  7)/GSP( 41)/GSP(  2))
      SPRF( 442)=AKF( 442)*SPRF( 441)
     1           /(GSP( 43)*GSP(  7)/GSP( 41)/GSP(  2))
      SPRF( 444)=AKF( 444)*SPRF( 443)
     1           /(GSP( 42)*GSP(  1)/GSP( 41)/GSP(  3))
      SPRF( 446)=AKF( 446)*SPRF( 445)
     1           /(GSP( 43)*GSP(  1)/GSP( 41)/GSP(  3))
      SPRF( 448)=AKF( 448)*SPRF( 447)
     1           /(GSP( 42)*GSP(  5)/GSP( 41)/GSP(  4))
      SPRF( 450)=AKF( 450)*SPRF( 449)
     1           /(GSP( 43)*GSP(  5)/GSP( 41)/GSP(  4))
      SPRF( 452)=AKF( 452)*SPRF( 451)
     1           /(GSP( 43)*GSP(  6)/GSP( 41)/GSP(  5))
      SPRF( 454)=AKF( 454)*SPRF( 453)
     1           /(GSP( 42)*GSP(  6)/GSP( 41)/GSP(  5))
      SPRF( 456)=AKF( 456)*SPRF( 455)
     1           /(GSP( 42)*GSP(  8)/GSP( 41)/GSP(  7))
      SPRF( 458)=AKF( 458)*SPRF( 457)
     1           /(GSP( 43)*GSP(  8)/GSP( 41)/GSP(  7))
      SPRF( 460)=AKF( 460)*SPRF( 459)
     1           /(GSP( 43)*GSP( 41)/GSP( 42)/GSP( 41))
      SPRF( 462)=AKF( 462)*SPRF( 461)
     1           /RUT/(GSP( 42)/GSP( 40)/GSP(  3))
      SPRF( 464)=AKF( 464)*SPRF( 463)
     1           /(GSP( 40)*GSP(  7)/GSP( 42)/GSP(  2))
      SPRF( 466)=AKF( 466)*SPRF( 465)
     1           *RUT/(GSP( 14)*GSP( 17)/GSP( 43))
      SPRF( 468)=AKF( 468)*SPRF( 467)
     1           /RUT/(GSP( 43)/GSP(  3)/GSP( 40))
      SPRF( 470)=AKF( 470)*SPRF( 469)
     1           /(GSP( 40)*GSP(  7)/GSP( 43)/GSP(  2))
      SPRF( 474)=AKF( 474)*SPRF( 473)
     1           *RUT/(GSP( 23)*GSP( 23)/GSP( 44))
      SPRF( 476)=AKF( 476)*SPRF( 475)
     1           /RUT/(GSP( 44)/GSP( 23)/GSP( 23))
      SPRF( 478)=AKF( 478)*SPRF( 477)
     1           /(GSP( 23)*GSP( 17)/GSP( 44)/GSP(  3))
      SPRF( 480)=AKF( 480)*SPRF( 479)
     1           *RUT/(GSP(  1)*GSP( 22)*GSP( 23)/GSP( 44)/GSP(  3))
      SPRF( 482)=AKF( 482)*SPRF( 481)
     1           *RUT/(GSP( 45)*GSP(  3)*GSP( 39)/GSP( 44)/GSP(  5))
      SPRF( 484)=AKF( 484)*SPRF( 483)
     1           *RUT/(GSP( 13)*GSP( 22)*GSP( 23)/GSP( 44)/GSP( 14))
      SPRF( 486)=AKF( 486)*SPRF( 485)
     1           /RUT/(GSP( 44)/GSP( 38)/GSP( 14))
      SPRF( 488)=AKF( 488)*SPRF( 487)
     1           *RUT/(GSP( 40)*GSP( 22)/GSP( 46))
      SPRF( 490)=AKF( 490)*SPRF( 489)
     1           *RUT/(GSP( 37)*GSP( 17)/GSP( 46))
      SPRF( 492)=AKF( 492)*SPRF( 491)
     1           *RUT/(GSP( 39)*GSP( 23)/GSP( 46))
      SPRF( 494)=AKF( 494)*SPRF( 493)
     1           *RUT/(GSP( 22)*GSP( 39)*GSP(  7)/GSP( 46)/GSP(  2))
      SPRF( 496)=AKF( 496)*SPRF( 495)
     1           *RUT/(GSP( 23)*GSP( 37)*GSP(  7)/GSP( 46)/GSP(  2))
      SPRF( 498)=AKF( 498)*SPRF( 497)
     1           *RUT/(GSP( 22)*GSP( 39)*GSP(  8)/GSP( 46)/GSP(  7))
      SPRF( 500)=AKF( 500)*SPRF( 499)
     1           *RUT/(GSP( 23)*GSP( 37)*GSP(  8)/GSP( 46)/GSP(  7))
      SPRF( 510)=AKF( 510)*SPRF( 509)
     1           *RUT/(GSP(  1)*GSP( 49)*GSP( 19)/GSP( 47)/GSP(  3))
      SPRF( 512)=AKF( 512)*SPRF( 511)
     1           *RUT/(GSP(  6)*GSP( 49)*GSP( 19)/GSP( 47)/GSP(  5))
      SPRF( 520)=AKF( 520)*SPRF( 519)
     1           *RUT/(GSP( 39)*GSP( 19)/GSP( 49))
      SPRF( 522)=AKF( 522)*SPRF( 521)
     1           *RUT/(GSP( 40)*GSP( 17)/GSP( 49))
      SPRF( 524)=AKF( 524)*SPRF( 523)
     1           *RUT/(GSP(  6)*GSP( 40)*GSP( 23)/GSP( 49)/GSP(  5))
      SPRF( 526)=AKF( 526)*SPRF( 525)
     1           *RUT/(GSP(  1)*GSP( 40)*GSP( 23)/GSP( 49)/GSP(  3))
      SPRF( 528)=AKF( 528)*SPRF( 527)
     1           *RUT/(GSP( 17)*GSP( 17)*GSP( 14)/GSP( 49)/GSP(  3))
      SPRF( 530)=AKF( 530)*SPRF( 529)
     1           /(GSP( 40)*GSP( 19)/GSP( 49)/GSP(  3))
      SPRF( 532)=AKF( 532)*SPRF( 531)
     1           *RUT/(GSP(  1)*GSP( 17)*GSP( 39)/GSP( 49)/GSP(  3))
      SPRF( 534)=AKF( 534)*SPRF( 533)
     1           *RUT/(GSP(  1)*GSP( 44)*GSP( 14)/GSP( 49)/GSP(  3))
      SPRF( 536)=AKF( 536)*SPRF( 535)
     1           *RUT/(GSP( 39)*GSP( 14)/GSP( 48))
      SPRF( 538)=AKF( 538)*SPRF( 537)
     1           /(GSP( 17)*GSP( 19)/GSP( 48)/GSP(  3))
      SPRF( 540)=AKF( 540)*SPRF( 539)
     1           /(GSP( 40)*GSP( 14)/GSP( 48)/GSP(  3))
      SPRF( 542)=AKF( 542)*SPRF( 541)
     1           *RUT/(GSP(  1)*GSP( 23)*GSP( 17)/GSP( 48)/GSP(  3))
      SPRF( 544)=AKF( 544)*SPRF( 543)
     1           *RUT/(GSP(  6)*GSP( 44)*GSP(  3)/GSP( 48)/GSP(  5))
  206 CONTINUE
C--------------------------END BACKWARD RATES---------------------------
C------------------------FORWARD REACTION RATES-------------------------
      RHO1=RHO(J,I)
      RHO2=RHO(J,I)*RHO(J,I)
      RHO3=RHO(J,I)*RHO(J,I)*RHO(J,I)
      TBODY=FA0(J,I)/WM(1)
      TBODY=TBODY+FA1(J,I)/WM(  2)
      TBODY=TBODY+FA2(J,I)/WM(  3)
      TBODY=TBODY+FA3(J,I)/WM(  4)
      TBODY=TBODY+FA4(J,I)/WM(  5)
      TBODY=TBODY+FA5(J,I)/WM(  6)
      TBODY=TBODY+FA6(J,I)/WM(  7)
      TBODY=TBODY+FA7(J,I)/WM(  8)
      TBODY=TBODY+FA8(J,I)/WM(  9)
      TBODY=TBODY+FA9(J,I)/WM( 10)
      TBODY=TBODY+FAA(J,I)/WM( 11)
      TBODY=TBODY+FAB(J,I)/WM( 12)
      TBODY=TBODY+FAC(J,I)/WM( 13)
      TBODY=TBODY+FAD(J,I)/WM( 14)
      TBODY=TBODY+FAE(J,I)/WM( 15)
      TBODY=TBODY+FAF(J,I)/WM( 16)
      TBODY=TBODY+FAG(J,I)/WM( 17)
      TBODY=TBODY+FAH(J,I)/WM( 18)
      TBODY=TBODY+FAI(J,I)/WM( 19)
      TBODY=TBODY+FAJ(J,I)/WM( 20)
      TBODY=TBODY+FB0(J,I)/WM( 21)
      TBODY=TBODY+FB1(J,I)/WM( 22)
      TBODY=TBODY+FB2(J,I)/WM( 23)
      TBODY=TBODY+FB3(J,I)/WM( 24)
      TBODY=TBODY+FB4(J,I)/WM( 25)
      TBODY=TBODY+FB5(J,I)/WM( 26)
      TBODY=TBODY+FB6(J,I)/WM( 27)
      TBODY=TBODY+FB7(J,I)/WM( 28)
      TBODY=TBODY+FB8(J,I)/WM( 29)
      TBODY=TBODY+FB9(J,I)/WM( 30)
      TBODY=TBODY+FBA(J,I)/WM( 31)
      TBODY=TBODY+FBB(J,I)/WM( 32)
      TBODY=TBODY+FBC(J,I)/WM( 33)
      TBODY=TBODY+FBD(J,I)/WM( 34)
      TBODY=TBODY+FBE(J,I)/WM( 35)
      TBODY=TBODY+FBF(J,I)/WM( 36)
      TBODY=TBODY+FBG(J,I)/WM( 37)
      TBODY=TBODY+FBH(J,I)/WM( 38)
      TBODY=TBODY+FBI(J,I)/WM( 39)
      TBODY=TBODY+FBJ(J,I)/WM( 40)
      TBODY=TBODY+FC0(J,I)/WM( 41)
      TBODY=TBODY+FC1(J,I)/WM( 42)
      TBODY=TBODY+FC2(J,I)/WM( 43)
      TBODY=TBODY+FC3(J,I)/WM( 44)
      TBODY=TBODY+FC4(J,I)/WM( 45)
      TBODY=TBODY+FC5(J,I)/WM( 46)
      TBODY=TBODY+FC6(J,I)/WM( 47)
      TBODY=TBODY+FC7(J,I)/WM( 48)
      TBODY=TBODY+FC8(J,I)/WM( 49)
      TBODY=TBODY+FC9(J,I)/WM( 50)
      TBODY=TBODY+FCA(J,I)/WM( 51)
      TBODY=TBODY+FCB(J,I)/WM( 52)
      TB001=TBODY+1.50*FA0(J,I)/WM(01)+11.0*FA5(J,I)/WM(06)
     1           +0.90*FA8(J,I)/WM(09)+2.80*FA9(J,I)/WM(10)
     2           -0.50*FC9(J,I)/WM(50)-0.50*FCA(J,I)/WM(51)
      TB002=TBODY+1.50*FA0(J,I)/WM(01)+11.0*FA5(J,I)/WM(06)
     1           +0.90*FA8(J,I)/WM(09)+2.80*FA9(J,I)/WM(10)
     2           -0.62*FC9(J,I)/WM(50)-0.62*FCA(J,I)/WM(51)
      TB003=TBODY+1.50*FA0(J,I)/WM(01)+11.0*FA5(J,I)/WM(06)
     1           +0.90*FA8(J,I)/WM(09)+2.80*FA9(J,I)/WM(10)
     2           -0.80*FC9(J,I)/WM(50)-0.80*FCA(J,I)/WM(51)
      TB004=TBODY+1.50*FA0(J,I)/WM(01)+11.0*FA5(J,I)/WM(06)
     1           +0.90*FA8(J,I)/WM(09)+2.80*FA9(J,I)/WM(10)
     2           -0.25*FC9(J,I)/WM(50)-0.25*FCA(J,I)/WM(51)
      TB005=TBODY+1.50*FA0(J,I)/WM(01)+15.0*FA5(J,I)/WM(06)
     1           +0.20*FA8(J,I)/WM(09)+1.40*FA9(J,I)/WM(10)
     2                                +0.50*FAJ(J,I)/WM(20)
     3           -0.30*FC9(J,I)/WM(50)-0.30*FCA(J,I)/WM(51)
      TB006=TBODY+1.00*FA0(J,I)/WM(01)+5.00*FA5(J,I)/WM(06)
     1           +0.50*FA8(J,I)/WM(09)+1.00*FA9(J,I)/WM(10)
     2           +1.00*FAC(J,I)/WM(13)+2.00*FAJ(J,I)/WM(20)
     3           -0.60*FC9(J,I)/WM(50)-0.60*FCA(J,I)/WM(51)
      TB007=TBODY+0.90*FA0(J,I)/WM(01)+11.0*FA5(J,I)/WM(06)
     1           +1.50*FA8(J,I)/WM(09)+1.50*FA9(J,I)/WM(10)
      TB008=TBODY+1.00*FA0(J,I)/WM(01)+5.00*FA5(J,I)/WM(06)
     1           +0.50*FA8(J,I)/WM(09)+1.00*FA9(J,I)/WM(10)
     2           +1.00*FAC(J,I)/WM(13)+2.00*FAJ(J,I)/WM(20)
     3           -0.30*FC9(J,I)/WM(50)
      TB009=TBODY+1.00*FA0(J,I)/WM(01)+5.00*FA5(J,I)/WM(06)
     1           +0.50*FA8(J,I)/WM(09)+1.00*FA9(J,I)/WM(10)
     2           +1.00*FAC(J,I)/WM(13)
     3           -0.30*FC9(J,I)/WM(50)
      TB010=TBODY+1.40*FA0(J,I)/WM(01)+14.4*FA5(J,I)/WM(06)
     1           +0.80*FA8(J,I)/WM(09)+2.60*FA9(J,I)/WM(10)
C-----------------------------------------------------------------------
C-------Specific Reaction Rates for the Reactions in Fall-Off Form------
C            TROE(IR,1) = 0.0 ==> Lindemann Form (RFACT=1.0)
C            TROE(IR,1) = +a  ==> Troe Form
C            TROE(IR,1) = -a  ==> SRI Form
C            TROE(IR,2) = 0.0  ==> Linear Form
C            AKLOW(IR) = 0.0  ==> TSENG Form (K0=Kinf*polynomial)
C-----------------------------------------------------------------------
      IF(NRLIND.EQ.0) GO TO 208 
      DMOLEA=RHO(J,I)*RSTR/WMSTR*1.0D-06
      DO 207 IR=1,LRX,LREV
      DMOLE=TBODY*DMOLEA
      IF(IR.EQ.019) DMOLE=TB005*DMOLEA
      IF(IR.EQ.031) DMOLE=TB006*DMOLEA
      IF(IR.EQ.063) DMOLE=TB008*DMOLEA
      IF(IR.EQ.105) DMOLE=TB009*DMOLEA
      IF(IR.EQ.107) DMOLE=TB008*DMOLEA
      IF(IR.EQ.161) DMOLE=TB008*DMOLEA
      IF(IR.EQ.173) DMOLE=TB009*DMOLEA
      IF(IR.EQ.197) DMOLE=TB009*DMOLEA
      IF(IR.EQ.271) DMOLE=TB009*DMOLEA
      IF(IR.EQ.295) DMOLE=TB009*DMOLEA
      IF(IR.EQ.319) DMOLE=TB009*DMOLEA
      IF(IR.EQ.321) DMOLE=TB009*DMOLEA
      IF(IR.EQ.425) DMOLE=TB008*DMOLEA
      IF(IR.EQ.431) DMOLE=TB008*DMOLEA
      IF(IR.EQ.461) DMOLE=TB008*DMOLEA
      IF(IR.EQ.467) DMOLE=TB008*DMOLEA
      IF(AKLOW(IR).EQ.0.0) GO TO 207
      IF(AKLOW(IR).LT.0.0) THEN
           AK0M=(10**(TROE(IR,1)+(TROE(IR,2)+(TROE(IR,3)
     1          +TROE(IR,4)*TKD)*TKD)*TKD))*DMOLE
           AFP=(AK0M/(1.0+AK0M))
           SPRF(IR)=AFP*SPRF(IR)
           IF(LREV.EQ.2) SPRF(IR+1)=AFP*SPRF(IR+1)
           GO TO 207
           END IF
      ACTV=ALLOW(IR)*ELOGT-EALOW(IR)/TKD
      AK0M=AKLOW(IR)*DEXP(ACTV)*DMOLE/(SPRF(IR)+1.0D-20)
      IF(TROE(IR,1).GT.0.0) THEN
           IF(TROE(IR,2).EQ.0.0) THEN
                FCLOG=DLOG10(TROE(IR,1))
                ELSE
                FCLOG=DLOG10((1.0-TROE(IR,1))*DEXP(-TKD/TROE(IR,2))
     1          +TROE(IR,1)*DEXP(-TKD/TROE(IR,3))+DEXP(-TROE(IR,4)/TKD))
                END IF
           CPLOG=DLOG10(AK0M)-(0.4+0.67*FCLOG)
           RFACT=10**(FCLOG/(1.0
     1          +(CPLOG/(0.75-1.27*FCLOG-0.14*CPLOG))**2))
           ELSE
           IF(TROE(IR,1).EQ.0.0) THEN
                RFACT=1.0
                ELSE
                XFACT=1.0/(1.0+DLOG10(AK0M)**2)
                RFACT=(-TROE(IR,1)*DEXP(-TROE(IR,2)/TKD)
     1               +DEXP(-TKD/TROE(IR,3)))**XFACT
                END IF
           END IF
      AFP=(AK0M/(1.0+AK0M))*RFACT
      SPRF(IR)=AFP*SPRF(IR)
      IF(LREV.EQ.2) SPRF(IR+1)=AFP*SPRF(IR+1)
  207 CONTINUE
  208 CONTINUE
C-----------------------------------------------------------------
      IA=(I-1)*LJ+J
      RC(IA,   1)=C2ND*RHO2*SPRF(   1)/WM(  3)/WM(  2)
      RC(IA,   2)=C2ND*RHO2*SPRF(   2)/WM(  5)/WM(  4)
      RC(IA,   3)=C2ND*RHO2*SPRF(   3)/WM(  1)/WM(  4)
      RC(IA,   4)=C2ND*RHO2*SPRF(   4)/WM(  5)/WM(  3)
      RC(IA,   5)=C2ND*RHO2*SPRF(   5)/WM(  1)/WM(  5)
      RC(IA,   6)=C2ND*RHO2*SPRF(   6)/WM(  6)/WM(  3)
      RC(IA,   7)=C2ND*RHO2*SPRF(   7)/WM(  6)/WM(  4)
      RC(IA,   8)=C2ND*RHO2*SPRF(   8)/WM(  5)/WM(  5)
      RC(IA,   9)=C3RD*RHO3*TB001*SPRF(   9)/WM(  3)/WM(  3)
      RC(IA,  10)=C2ND*RHO2*TB001*SPRF(  10)/WM(  1)
      RC(IA,  11)=C3RD*RHO3*TB002*SPRF(  11)/WM(  3)/WM(  5)
      RC(IA,  12)=C2ND*RHO2*TB002*SPRF(  12)/WM(  6)
      RC(IA,  13)=C3RD*RHO3*TB003*SPRF(  13)/WM(  4)/WM(  4)
      RC(IA,  14)=C2ND*RHO2*TB003*SPRF(  14)/WM(  2)
      RC(IA,  15)=C3RD*RHO3*TB004*SPRF(  15)/WM(  3)/WM(  4)
      RC(IA,  16)=C2ND*RHO2*TB004*SPRF(  16)/WM(  5)
      RC(IA,  17)=C3RD*RHO3*TB004*SPRF(  17)/WM(  4)/WM(  5)
      RC(IA,  18)=C2ND*RHO2*TB004*SPRF(  18)/WM(  7)
      RC(IA,  19)=C2ND*RHO2*SPRF(  19)/WM(  3)/WM(  2)
      RC(IA,  20)=C1ST*RHO1*SPRF(  20)/WM(  7)
      RC(IA,  21)=C2ND*RHO2*SPRF(  21)/WM(  7)/WM(  3)
      RC(IA,  22)=C2ND*RHO2*SPRF(  22)/WM(  5)/WM(  5)
      RC(IA,  23)=C2ND*RHO2*SPRF(  23)/WM(  7)/WM(  3)
      RC(IA,  24)=C2ND*RHO2*SPRF(  24)/WM(  1)/WM(  2)
      RC(IA,  25)=C2ND*RHO2*SPRF(  25)/WM(  7)/WM(  3)
      RC(IA,  26)=C2ND*RHO2*SPRF(  26)/WM(  6)/WM(  4)
      RC(IA,  27)=C2ND*RHO2*SPRF(  27)/WM(  7)/WM(  4)
      RC(IA,  28)=C2ND*RHO2*SPRF(  28)/WM(  5)/WM(  2)
      RC(IA,  29)=C2ND*RHO2*SPRF(  29)/WM(  7)/WM(  5)
      RC(IA,  30)=C2ND*RHO2*SPRF(  30)/WM(  6)/WM(  2)
      RC(IA,  31)=C2ND*RHO2*SPRF(  31)/WM(  5)/WM(  5)
      RC(IA,  32)=C1ST*RHO1*SPRF(  32)/WM(  8)
      RC(IA,  33)=C2ND*RHO2*SPRF(  33)/WM(  7)/WM(  7)
      RC(IA,  34)=C2ND*RHO2*SPRF(  34)/WM(  8)/WM(  2)
      RC(IA,  35)=C2ND*RHO2*SPRF(  35)/WM(  8)/WM(  3)
      RC(IA,  36)=C2ND*RHO2*SPRF(  36)/WM(  7)/WM(  1)
      RC(IA,  37)=C2ND*RHO2*SPRF(  37)/WM(  8)/WM(  3)
      RC(IA,  38)=C2ND*RHO2*SPRF(  38)/WM(  6)/WM(  5)
      RC(IA,  39)=C2ND*RHO2*SPRF(  39)/WM(  8)/WM(  5)
      RC(IA,  40)=C2ND*RHO2*SPRF(  40)/WM(  6)/WM(  7)
      RC(IA,  41)=C2ND*RHO2*SPRF(  41)/WM(  8)/WM(  4)
      RC(IA,  42)=C2ND*RHO2*SPRF(  42)/WM(  7)/WM(  5)
      RC(IA,  43)=C2ND*RHO2*SPRF(  43)/WM(  9)/WM(  5)
      RC(IA,  44)=C2ND*RHO2*SPRF(  44)/WM( 10)/WM(  3)
      RC(IA,  45)=C2ND*RHO2*SPRF(  45)/WM(  9)/WM(  7)
      RC(IA,  46)=C2ND*RHO2*SPRF(  46)/WM( 10)/WM(  5)
      RC(IA,  47)=C2ND*RHO2*SPRF(  47)/WM(  9)/WM(  2)
      RC(IA,  48)=C2ND*RHO2*SPRF(  48)/WM( 10)/WM(  4)
      RC(IA,  49)=C2ND*RHO2*TB007*SPRF(  49)/WM( 11)
      RC(IA,  50)=C3RD*RHO3*TB007*SPRF(  50)/WM(  9)/WM(  3)
      RC(IA,  51)=C2ND*RHO2*SPRF(  51)/WM( 11)/WM(  3)
      RC(IA,  52)=C2ND*RHO2*SPRF(  52)/WM(  9)/WM(  1)
      RC(IA,  53)=C2ND*RHO2*SPRF(  53)/WM( 11)/WM(  4)
      RC(IA,  54)=C2ND*RHO2*SPRF(  54)/WM(  9)/WM(  5)
      RC(IA,  55)=C2ND*RHO2*SPRF(  55)/WM( 11)/WM(  4)
      RC(IA,  56)=C2ND*RHO2*SPRF(  56)/WM( 10)/WM(  3)
      RC(IA,  57)=C2ND*RHO2*SPRF(  57)/WM( 11)/WM(  5)
      RC(IA,  58)=C2ND*RHO2*SPRF(  58)/WM(  9)/WM(  6)
      RC(IA,  59)=C2ND*RHO2*SPRF(  59)/WM( 11)/WM(  2)
      RC(IA,  60)=C2ND*RHO2*SPRF(  60)/WM(  9)/WM(  7)
      RC(IA,  61)=C2ND*RHO2*SPRF(  61)/WM( 11)/WM( 14)
      RC(IA,  62)=C2ND*RHO2*SPRF(  62)/WM(  9)/WM( 13)
      RC(IA,  63)=C2ND*RHO2*SPRF(  63)/WM(  3)/WM( 11)
      RC(IA,  64)=C1ST*RHO1*SPRF(  64)/WM( 12)
      RC(IA,  65)=C2ND*RHO2*SPRF(  65)/WM( 12)/WM(  3)
      RC(IA,  66)=C2ND*RHO2*SPRF(  66)/WM( 11)/WM(  1)
      RC(IA,  67)=C2ND*RHO2*SPRF(  67)/WM( 12)/WM(  4)
      RC(IA,  68)=C2ND*RHO2*SPRF(  68)/WM( 11)/WM(  5)
      RC(IA,  69)=C2ND*RHO2*SPRF(  69)/WM( 12)/WM(  5)
      RC(IA,  70)=C2ND*RHO2*SPRF(  70)/WM( 11)/WM(  6)
      RC(IA,  71)=C2ND*RHO2*SPRF(  71)/WM( 12)/WM(  2)
      RC(IA,  72)=C2ND*RHO2*SPRF(  72)/WM( 11)/WM(  7)
      RC(IA,  73)=C2ND*RHO2*SPRF(  73)/WM( 12)/WM(  7)
      RC(IA,  74)=C2ND*RHO2*SPRF(  74)/WM( 11)/WM(  8)
      RC(IA,  75)=C2ND*RHO2*SPRF(  75)/WM( 13)/WM(  3)
      RC(IA,  76)=C2ND*RHO2*SPRF(  76)/WM(  1)/WM( 14)
      RC(IA,  77)=C2ND*RHO2*SPRF(  77)/WM( 13)/WM(  5)
      RC(IA,  78)=C2ND*RHO2*SPRF(  78)/WM(  6)/WM( 14)
      RC(IA,  79)=C2ND*RHO2*SPRF(  79)/WM( 13)/WM(  4)
      RC(IA,  80)=C2ND*RHO2*SPRF(  80)/WM( 14)/WM(  5)
      RC(IA,  81)=C2ND*RHO2*SPRF(  81)/WM( 13)/WM(  2)
      RC(IA,  82)=C2ND*RHO2*SPRF(  82)/WM( 14)/WM(  7)
      RC(IA,  83)=C2ND*RHO2*SPRF(  83)/WM( 13)/WM(  7)
      RC(IA,  84)=C2ND*RHO2*SPRF(  84)/WM( 14)/WM(  8)
      RC(IA,  85)=C2ND*RHO2*SPRF(  85)/WM( 14)/WM(  3)
      RC(IA,  86)=C2ND*RHO2*SPRF(  86)/WM( 15)/WM(  1)
      RC(IA,  87)=C2ND*RHO2*SPRF(  87)/WM( 14)/WM(  3)
      RC(IA,  88)=C2ND*RHO2*SPRF(  88)/WM( 16)/WM(  1)
      RC(IA,  89)=C2ND*RHO2*SPRF(  89)/WM( 14)/WM(  5)
      RC(IA,  90)=C2ND*RHO2*SPRF(  90)/WM( 16)/WM(  6)
      RC(IA,  91)=C2ND*RHO2*SPRF(  91)/WM( 14)/WM(  4)
      RC(IA,  92)=C2ND*RHO2*SPRF(  92)/WM( 12)/WM(  3)
      RC(IA,  93)=C2ND*RHO2*SPRF(  93)/WM( 14)/WM( 15)
      RC(IA,  94)=C2ND*RHO2*SPRF(  94)/WM( 17)/WM(  3)
      RC(IA,  95)=C2ND*RHO2*SPRF(  95)/WM( 14)/WM(  7)
      RC(IA,  96)=C2ND*RHO2*SPRF(  96)/WM( 18)/WM(  5)
      RC(IA,  97)=C2ND*RHO2*SPRF(  97)/WM( 14)/WM(  2)
      RC(IA,  98)=C2ND*RHO2*SPRF(  98)/WM( 12)/WM(  5)
      RC(IA,  99)=C2ND*RHO2*SPRF(  99)/WM( 14)/WM(  2)
      RC(IA, 100)=C2ND*RHO2*SPRF( 100)/WM( 18)/WM(  4)
      RC(IA, 101)=C2ND*RHO2*SPRF( 101)/WM( 14)/WM( 14)
      RC(IA, 102)=C2ND*RHO2*SPRF( 102)/WM( 17)/WM(  1)
      RC(IA, 103)=C2ND*RHO2*SPRF( 103)/WM( 14)/WM( 14)
      RC(IA, 104)=C2ND*RHO2*SPRF( 104)/WM( 19)/WM(  3)
      RC(IA, 105)=C2ND*RHO2*SPRF( 105)/WM(  3)/WM( 14)
      RC(IA, 106)=C1ST*RHO1*SPRF( 106)/WM( 13)
      RC(IA, 107)=C2ND*RHO2*SPRF( 107)/WM( 14)/WM( 14)
      RC(IA, 108)=C1ST*RHO1*SPRF( 108)/WM( 20)
      RC(IA, 109)=C2ND*RHO2*SPRF( 109)/WM( 16)/WM(  5)
      RC(IA, 110)=C2ND*RHO2*SPRF( 110)/WM( 12)/WM(  3)
      RC(IA, 111)=C2ND*RHO2*SPRF( 111)/WM( 16)/WM(  2)
      RC(IA, 112)=C3RD*RHO3*SPRF( 112)/WM(  9)/WM(  5)/WM(  3)
      RC(IA, 113)=C2ND*RHO2*SPRF( 113)/WM( 16)/WM( 10)
      RC(IA, 114)=C2ND*RHO2*SPRF( 114)/WM(  9)/WM( 12)
      RC(IA, 115)=C2ND*RHO2*TB010*SPRF( 115)/WM( 16)
      RC(IA, 116)=C2ND*RHO2*TB010*SPRF( 116)/WM( 15)
      RC(IA, 117)=C2ND*RHO2*SPRF( 117)/WM( 15)/WM(  3)
      RC(IA, 118)=C2ND*RHO2*SPRF( 118)/WM( 21)/WM(  1)
      RC(IA, 119)=C2ND*RHO2*SPRF( 119)/WM( 15)/WM(  5)
      RC(IA, 120)=C2ND*RHO2*SPRF( 120)/WM( 12)/WM(  3)
      RC(IA, 121)=C2ND*RHO2*SPRF( 121)/WM( 15)/WM(  5)
      RC(IA, 122)=C2ND*RHO2*SPRF( 122)/WM( 21)/WM(  6)
      RC(IA, 123)=C2ND*RHO2*SPRF( 123)/WM( 15)/WM(  4)
      RC(IA, 124)=C3RD*RHO3*SPRF( 124)/WM(  9)/WM(  3)/WM(  3)
      RC(IA, 125)=C2ND*RHO2*SPRF( 125)/WM( 15)/WM(  4)
      RC(IA, 126)=C2ND*RHO2*SPRF( 126)/WM(  9)/WM(  1)
      RC(IA, 127)=C2ND*RHO2*SPRF( 127)/WM( 15)/WM(  2)
      RC(IA, 128)=C2ND*RHO2*SPRF( 128)/WM( 10)/WM(  1)
      RC(IA, 129)=C2ND*RHO2*SPRF( 129)/WM( 15)/WM(  2)
      RC(IA, 130)=C3RD*RHO3*SPRF( 130)/WM(  9)/WM(  5)/WM(  3)
      RC(IA, 131)=C2ND*RHO2*SPRF( 131)/WM( 15)/WM( 15)
      RC(IA, 132)=C3RD*RHO3*SPRF( 132)/WM( 22)/WM(  3)/WM(  3)
      RC(IA, 133)=C2ND*RHO2*SPRF( 133)/WM( 21)/WM(  4)
      RC(IA, 134)=C2ND*RHO2*SPRF( 134)/WM(  9)/WM(  3)
      RC(IA, 135)=C2ND*RHO2*SPRF( 135)/WM( 21)/WM(  2)
      RC(IA, 136)=C2ND*RHO2*SPRF( 136)/WM( 11)/WM(  4)
      RC(IA, 137)=C2ND*RHO2*SPRF( 137)/WM( 21)/WM(  6)
      RC(IA, 138)=C2ND*RHO2*SPRF( 138)/WM( 12)/WM(  3)
      RC(IA, 139)=C2ND*RHO2*SPRF( 139)/WM( 21)/WM( 10)
      RC(IA, 140)=C2ND*RHO2*SPRF( 140)/WM( 11)/WM(  9)
      RC(IA, 141)=C2ND*RHO2*SPRF( 141)/WM( 18)/WM(  3)
      RC(IA, 142)=C2ND*RHO2*SPRF( 142)/WM( 12)/WM(  1)
      RC(IA, 143)=C2ND*RHO2*SPRF( 143)/WM( 18)/WM(  3)
      RC(IA, 144)=C2ND*RHO2*SPRF( 144)/WM( 16)/WM(  6)
      RC(IA, 145)=C2ND*RHO2*SPRF( 145)/WM( 18)/WM(  5)
      RC(IA, 146)=C2ND*RHO2*SPRF( 146)/WM( 12)/WM(  6)
      RC(IA, 147)=C2ND*RHO2*SPRF( 147)/WM( 18)/WM(  4)
      RC(IA, 148)=C2ND*RHO2*SPRF( 148)/WM(  5)/WM( 12)
      RC(IA, 149)=C2ND*RHO2*SPRF( 149)/WM( 18)/WM(  2)
      RC(IA, 150)=C2ND*RHO2*SPRF( 150)/WM( 12)/WM(  7)
      RC(IA, 151)=C2ND*RHO2*TB009*SPRF( 151)/WM( 18)
      RC(IA, 152)=C3RD*RHO3*TB009*SPRF( 152)/WM( 12)/WM(  3)
      RC(IA, 153)=C2ND*RHO2*SPRF( 153)/WM( 20)/WM(  3)
      RC(IA, 154)=C2ND*RHO2*SPRF( 154)/WM( 19)/WM(  1)
      RC(IA, 155)=C2ND*RHO2*SPRF( 155)/WM( 20)/WM(  4)
      RC(IA, 156)=C2ND*RHO2*SPRF( 156)/WM( 19)/WM(  5)
      RC(IA, 157)=C2ND*RHO2*SPRF( 157)/WM( 20)/WM(  5)
      RC(IA, 158)=C2ND*RHO2*SPRF( 158)/WM( 19)/WM(  6)
      RC(IA, 159)=C2ND*RHO2*SPRF( 159)/WM( 20)/WM( 14)
      RC(IA, 160)=C2ND*RHO2*SPRF( 160)/WM( 19)/WM( 13)
      RC(IA, 161)=C1ST*RHO1*SPRF( 161)/WM( 20)
      RC(IA, 162)=C2ND*RHO2*SPRF( 162)/WM( 19)/WM(  3)
      RC(IA, 163)=C2ND*RHO2*SPRF( 163)/WM( 20)/WM(  7)
      RC(IA, 164)=C2ND*RHO2*SPRF( 164)/WM( 19)/WM(  8)
      RC(IA, 165)=C2ND*RHO2*SPRF( 165)/WM( 19)/WM(  3)
      RC(IA, 166)=C2ND*RHO2*SPRF( 166)/WM( 17)/WM(  1)
      RC(IA, 167)=C2ND*RHO2*SPRF( 167)/WM( 19)/WM(  4)
      RC(IA, 168)=C2ND*RHO2*SPRF( 168)/WM( 17)/WM(  5)
      RC(IA, 169)=C2ND*RHO2*SPRF( 169)/WM( 19)/WM(  4)
      RC(IA, 170)=C2ND*RHO2*SPRF( 170)/WM( 14)/WM( 12)
      RC(IA, 171)=C2ND*RHO2*SPRF( 171)/WM( 19)/WM(  2)
      RC(IA, 172)=C2ND*RHO2*SPRF( 172)/WM( 17)/WM(  7)
      RC(IA, 173)=C1ST*RHO1*SPRF( 173)/WM( 19)
      RC(IA, 174)=C2ND*RHO2*SPRF( 174)/WM( 17)/WM(  3)
      RC(IA, 175)=C2ND*RHO2*SPRF( 175)/WM( 17)/WM(  3)
      RC(IA, 176)=C2ND*RHO2*SPRF( 176)/WM( 23)/WM(  1)
      RC(IA, 177)=C2ND*RHO2*SPRF( 177)/WM( 17)/WM(  5)
      RC(IA, 178)=C2ND*RHO2*SPRF( 178)/WM( 23)/WM(  6)
      RC(IA, 179)=C2ND*RHO2*SPRF( 179)/WM( 17)/WM(  4)
      RC(IA, 180)=C2ND*RHO2*SPRF( 180)/WM( 14)/WM( 11)
      RC(IA, 181)=C2ND*RHO2*SPRF( 181)/WM( 17)/WM(  4)
      RC(IA, 182)=C2ND*RHO2*SPRF( 182)/WM( 24)/WM(  3)
      RC(IA, 183)=C2ND*RHO2*SPRF( 183)/WM( 17)/WM( 17)
      RC(IA, 184)=C2ND*RHO2*SPRF( 184)/WM( 23)/WM( 19)
      RC(IA, 185)=C2ND*RHO2*SPRF( 185)/WM( 17)/WM(  2)
      RC(IA, 186)=C2ND*RHO2*SPRF( 186)/WM( 23)/WM(  7)
      RC(IA, 187)=C2ND*RHO2*SPRF( 187)/WM( 17)/WM(  7)
      RC(IA, 188)=C2ND*RHO2*SPRF( 188)/WM( 25)/WM(  5)
      RC(IA, 189)=C2ND*RHO2*SPRF( 189)/WM( 25)/WM(  7)
      RC(IA, 190)=C3RD*RHO3*SPRF( 190)/WM( 14)/WM(  9)/WM(  8)
      RC(IA, 191)=C2ND*RHO2*TB009*SPRF( 191)/WM( 17)
      RC(IA, 192)=C3RD*RHO3*TB009*SPRF( 192)/WM( 23)/WM(  3)
      RC(IA, 193)=C2ND*RHO2*TB009*SPRF( 193)/WM( 17)
      RC(IA, 194)=C3RD*RHO3*TB009*SPRF( 194)/WM( 22)/WM(  1)
      RC(IA, 195)=C2ND*RHO2*SPRF( 195)/WM( 23)/WM(  3)
      RC(IA, 196)=C2ND*RHO2*SPRF( 196)/WM( 22)/WM(  1)
      RC(IA, 197)=C1ST*RHO1*SPRF( 197)/WM( 23)
      RC(IA, 198)=C2ND*RHO2*SPRF( 198)/WM( 22)/WM(  3)
      RC(IA, 199)=C2ND*RHO2*SPRF( 199)/WM( 23)/WM(  2)
      RC(IA, 200)=C2ND*RHO2*SPRF( 200)/WM( 12)/WM( 11)
      RC(IA, 201)=C2ND*RHO2*SPRF( 201)/WM( 23)/WM(  2)
      RC(IA, 202)=C2ND*RHO2*SPRF( 202)/WM( 24)/WM(  4)
      RC(IA, 203)=C2ND*RHO2*SPRF( 203)/WM( 23)/WM(  2)
      RC(IA, 204)=C2ND*RHO2*SPRF( 204)/WM( 22)/WM(  7)
      RC(IA, 205)=C2ND*RHO2*SPRF( 205)/WM( 22)/WM(  4)
      RC(IA, 206)=C2ND*RHO2*SPRF( 206)/WM( 27)/WM(  3)
      RC(IA, 207)=C2ND*RHO2*SPRF( 207)/WM( 22)/WM(  4)
      RC(IA, 208)=C2ND*RHO2*SPRF( 208)/WM( 15)/WM(  9)
      RC(IA, 209)=C2ND*RHO2*SPRF( 209)/WM( 22)/WM(  2)
      RC(IA, 210)=C2ND*RHO2*SPRF( 210)/WM( 12)/WM(  9)
      RC(IA, 211)=C2ND*RHO2*SPRF( 211)/WM( 22)/WM(  5)
      RC(IA, 212)=C2ND*RHO2*SPRF( 212)/WM( 26)/WM(  3)
      RC(IA, 213)=C2ND*RHO2*SPRF( 213)/WM( 22)/WM(  5)
      RC(IA, 214)=C2ND*RHO2*SPRF( 214)/WM( 28)/WM(  6)
      RC(IA, 215)=C2ND*RHO2*SPRF( 215)/WM( 26)/WM(  3)
      RC(IA, 216)=C2ND*RHO2*SPRF( 216)/WM( 14)/WM(  9)
      RC(IA, 217)=C2ND*RHO2*SPRF( 217)/WM( 26)/WM(  4)
      RC(IA, 218)=C2ND*RHO2*SPRF( 218)/WM( 15)/WM( 10)
      RC(IA, 219)=C2ND*RHO2*SPRF( 219)/WM( 26)/WM(  4)
      RC(IA, 220)=C2ND*RHO2*SPRF( 220)/WM( 27)/WM(  5)
      RC(IA, 221)=C2ND*RHO2*SPRF( 221)/WM( 26)/WM( 14)
      RC(IA, 222)=C2ND*RHO2*SPRF( 222)/WM( 19)/WM(  9)
      RC(IA, 223)=C2ND*RHO2*SPRF( 223)/WM( 27)/WM(  3)
      RC(IA, 224)=C2ND*RHO2*SPRF( 224)/WM( 16)/WM(  9)
      RC(IA, 225)=C2ND*RHO2*SPRF( 225)/WM( 27)/WM(  5)
      RC(IA, 226)=C3RD*RHO3*SPRF( 226)/WM( 11)/WM(  9)/WM(  3)
      RC(IA, 227)=C2ND*RHO2*SPRF( 227)/WM( 27)/WM(  4)
      RC(IA, 228)=C3RD*RHO3*SPRF( 228)/WM(  9)/WM(  9)/WM(  3)
      RC(IA, 229)=C2ND*RHO2*SPRF( 229)/WM( 27)/WM(  2)
      RC(IA, 230)=C3RD*RHO3*SPRF( 230)/WM(  9)/WM(  9)/WM(  5)
      RC(IA, 231)=C2ND*RHO2*SPRF( 231)/WM( 27)/WM(  2)
      RC(IA, 232)=C3RD*RHO3*SPRF( 232)/WM( 10)/WM(  9)/WM(  3)
      RC(IA, 233)=C2ND*RHO2*SPRF( 233)/WM( 28)/WM(  5)
      RC(IA, 234)=C2ND*RHO2*SPRF( 234)/WM( 27)/WM(  3)
      RC(IA, 235)=C2ND*RHO2*SPRF( 235)/WM( 28)/WM(  4)
      RC(IA, 236)=C2ND*RHO2*SPRF( 236)/WM(  9)/WM( 21)
      RC(IA, 237)=C2ND*RHO2*SPRF( 237)/WM( 28)/WM(  2)
      RC(IA, 238)=C2ND*RHO2*SPRF( 238)/WM( 27)/WM(  4)
      RC(IA, 239)=C2ND*RHO2*SPRF( 239)/WM( 28)/WM(  2)
      RC(IA, 240)=C2ND*RHO2*SPRF( 240)/WM( 21)/WM( 10)
      RC(IA, 241)=C2ND*RHO2*SPRF( 241)/WM( 28)/WM(  2)
      RC(IA, 242)=C2ND*RHO2*SPRF( 242)/WM( 11)/WM(  9)
      RC(IA, 243)=C2ND*RHO2*SPRF( 243)/WM( 29)/WM(  3)
      RC(IA, 244)=C2ND*RHO2*SPRF( 244)/WM( 12)/WM(  1)
      RC(IA, 245)=C2ND*RHO2*SPRF( 245)/WM( 29)/WM(  3)
      RC(IA, 246)=C2ND*RHO2*SPRF( 246)/WM( 14)/WM(  5)
      RC(IA, 247)=C2ND*RHO2*SPRF( 247)/WM( 29)/WM(  5)
      RC(IA, 248)=C2ND*RHO2*SPRF( 248)/WM( 12)/WM(  6)
      RC(IA, 249)=C2ND*RHO2*SPRF( 249)/WM( 29)/WM(  2)
      RC(IA, 250)=C2ND*RHO2*SPRF( 250)/WM( 12)/WM(  7)
      RC(IA, 251)=C2ND*RHO2*TB009*SPRF( 251)/WM( 29)
      RC(IA, 252)=C3RD*RHO3*TB009*SPRF( 252)/WM( 12)/WM(  3)
      RC(IA, 253)=C2ND*RHO2*TB009*SPRF( 253)/WM( 18)
      RC(IA, 254)=C2ND*RHO2*TB009*SPRF( 254)/WM( 29)
      RC(IA, 255)=C2ND*RHO2*SPRF( 255)/WM( 26)/WM(  5)
      RC(IA, 256)=C2ND*RHO2*SPRF( 256)/WM( 29)/WM(  9)
      RC(IA, 257)=C2ND*RHO2*SPRF( 257)/WM( 30)/WM(  5)
      RC(IA, 258)=C2ND*RHO2*SPRF( 258)/WM( 29)/WM(  6)
      RC(IA, 259)=C2ND*RHO2*SPRF( 259)/WM( 30)/WM(  5)
      RC(IA, 260)=C2ND*RHO2*SPRF( 260)/WM( 18)/WM(  6)
      RC(IA, 261)=C2ND*RHO2*SPRF( 261)/WM( 30)/WM(  3)
      RC(IA, 262)=C2ND*RHO2*SPRF( 262)/WM( 29)/WM(  1)
      RC(IA, 263)=C2ND*RHO2*SPRF( 263)/WM( 30)/WM(  3)
      RC(IA, 264)=C2ND*RHO2*SPRF( 264)/WM( 18)/WM(  1)
      RC(IA, 265)=C2ND*RHO2*SPRF( 265)/WM( 30)/WM(  4)
      RC(IA, 266)=C2ND*RHO2*SPRF( 266)/WM( 29)/WM(  5)
      RC(IA, 267)=C2ND*RHO2*SPRF( 267)/WM( 30)/WM(  7)
      RC(IA, 268)=C2ND*RHO2*SPRF( 268)/WM( 29)/WM(  8)
      RC(IA, 269)=C2ND*RHO2*SPRF( 269)/WM( 30)/WM(  2)
      RC(IA, 270)=C2ND*RHO2*SPRF( 270)/WM( 29)/WM(  7)
      RC(IA, 271)=C1ST*RHO1*SPRF( 271)/WM( 30)
      RC(IA, 272)=C2ND*RHO2*SPRF( 272)/WM( 14)/WM(  5)
      RC(IA, 273)=C1ST*RHO1*SPRF( 273)/WM( 24)
      RC(IA, 274)=C2ND*RHO2*SPRF( 274)/WM( 26)/WM(  3)
      RC(IA, 275)=C2ND*RHO2*SPRF( 275)/WM( 24)/WM(  3)
      RC(IA, 276)=C2ND*RHO2*SPRF( 276)/WM( 14)/WM( 11)
      RC(IA, 277)=C2ND*RHO2*SPRF( 277)/WM( 24)/WM(  3)
      RC(IA, 278)=C2ND*RHO2*SPRF( 278)/WM( 26)/WM(  1)
      RC(IA, 279)=C2ND*RHO2*SPRF( 279)/WM( 24)/WM(  4)
      RC(IA, 280)=C2ND*RHO2*SPRF( 280)/WM( 12)/WM( 11)
      RC(IA, 281)=C2ND*RHO2*SPRF( 281)/WM( 24)/WM(  5)
      RC(IA, 282)=C2ND*RHO2*SPRF( 282)/WM( 26)/WM(  6)
      RC(IA, 283)=C2ND*RHO2*SPRF( 283)/WM( 24)/WM(  2)
      RC(IA, 284)=C3RD*RHO3*SPRF( 284)/WM( 12)/WM(  9)/WM(  5)
      RC(IA, 285)=C2ND*RHO2*SPRF( 285)/WM( 24)/WM( 14)
      RC(IA, 286)=C3RD*RHO3*SPRF( 286)/WM( 19)/WM(  9)/WM(  3)
      RC(IA, 287)=C2ND*RHO2*SPRF( 287)/WM( 24)/WM(  7)
      RC(IA, 288)=C3RD*RHO3*SPRF( 288)/WM( 12)/WM( 11)/WM(  5)
      RC(IA, 289)=C2ND*RHO2*SPRF( 289)/WM( 24)/WM(  7)
      RC(IA, 290)=C2ND*RHO2*SPRF( 290)/WM( 32)/WM(  2)
      RC(IA, 291)=C1ST*RHO1*SPRF( 291)/WM( 24)
      RC(IA, 292)=C2ND*RHO2*SPRF( 292)/WM( 14)/WM(  9)
      RC(IA, 293)=C1ST*RHO1*SPRF( 293)/WM( 32)
      RC(IA, 294)=C2ND*RHO2*SPRF( 294)/WM( 14)/WM( 11)
      RC(IA, 295)=C1ST*RHO1*SPRF( 295)/WM( 35)
      RC(IA, 296)=C2ND*RHO2*SPRF( 296)/WM( 14)/WM(  9)
      RC(IA, 297)=C2ND*RHO2*SPRF( 297)/WM( 32)/WM(  5)
      RC(IA, 298)=C2ND*RHO2*SPRF( 298)/WM( 35)/WM(  6)
      RC(IA, 299)=C2ND*RHO2*SPRF( 299)/WM( 32)/WM(  5)
      RC(IA, 300)=C2ND*RHO2*SPRF( 300)/WM( 24)/WM(  6)
      RC(IA, 301)=C2ND*RHO2*SPRF( 301)/WM( 32)/WM(  4)
      RC(IA, 302)=C2ND*RHO2*SPRF( 302)/WM( 35)/WM(  5)
      RC(IA, 303)=C2ND*RHO2*SPRF( 303)/WM( 32)/WM(  4)
      RC(IA, 304)=C2ND*RHO2*SPRF( 304)/WM( 24)/WM(  5)
      RC(IA, 305)=C2ND*RHO2*SPRF( 305)/WM( 32)/WM(  3)
      RC(IA, 306)=C2ND*RHO2*SPRF( 306)/WM( 35)/WM(  1)
      RC(IA, 307)=C2ND*RHO2*SPRF( 307)/WM( 32)/WM(  3)
      RC(IA, 308)=C2ND*RHO2*SPRF( 308)/WM( 24)/WM(  1)
      RC(IA, 309)=C2ND*RHO2*SPRF( 309)/WM( 32)/WM( 14)
      RC(IA, 310)=C2ND*RHO2*SPRF( 310)/WM( 35)/WM( 13)
      RC(IA, 311)=C2ND*RHO2*SPRF( 311)/WM( 32)/WM( 14)
      RC(IA, 312)=C2ND*RHO2*SPRF( 312)/WM( 24)/WM( 13)
      RC(IA, 313)=C2ND*RHO2*SPRF( 313)/WM( 32)/WM(  7)
      RC(IA, 314)=C2ND*RHO2*SPRF( 314)/WM( 35)/WM(  8)
      RC(IA, 315)=C2ND*RHO2*SPRF( 315)/WM( 32)/WM(  7)
      RC(IA, 316)=C2ND*RHO2*SPRF( 316)/WM( 24)/WM(  8)
      RC(IA, 317)=C2ND*RHO2*SPRF( 317)/WM( 32)/WM(  2)
      RC(IA, 318)=C2ND*RHO2*SPRF( 318)/WM( 35)/WM(  7)
      RC(IA, 319)=C1ST*RHO1*SPRF( 319)/WM( 31)
      RC(IA, 320)=C2ND*RHO2*SPRF( 320)/WM( 14)/WM( 29)
      RC(IA, 321)=C1ST*RHO1*SPRF( 321)/WM( 31)
      RC(IA, 322)=C2ND*RHO2*SPRF( 322)/WM( 17)/WM(  6)
      RC(IA, 323)=C2ND*RHO2*SPRF( 323)/WM( 31)/WM(  5)
      RC(IA, 324)=C2ND*RHO2*SPRF( 324)/WM( 34)/WM(  6)
      RC(IA, 325)=C2ND*RHO2*SPRF( 325)/WM( 31)/WM(  5)
      RC(IA, 326)=C2ND*RHO2*SPRF( 326)/WM( 33)/WM(  6)
      RC(IA, 327)=C2ND*RHO2*SPRF( 327)/WM( 31)/WM(  5)
      RC(IA, 328)=C2ND*RHO2*SPRF( 328)/WM( 36)/WM(  6)
      RC(IA, 329)=C2ND*RHO2*SPRF( 329)/WM( 31)/WM(  3)
      RC(IA, 330)=C2ND*RHO2*SPRF( 330)/WM( 34)/WM(  1)
      RC(IA, 331)=C2ND*RHO2*SPRF( 331)/WM( 31)/WM(  3)
      RC(IA, 332)=C2ND*RHO2*SPRF( 332)/WM( 33)/WM(  1)
      RC(IA, 333)=C2ND*RHO2*SPRF( 333)/WM( 31)/WM(  3)
      RC(IA, 334)=C2ND*RHO2*SPRF( 334)/WM( 36)/WM(  1)
      RC(IA, 335)=C2ND*RHO2*SPRF( 335)/WM( 31)/WM(  4)
      RC(IA, 336)=C2ND*RHO2*SPRF( 336)/WM( 34)/WM(  5)
      RC(IA, 337)=C2ND*RHO2*SPRF( 337)/WM( 31)/WM(  4)
      RC(IA, 338)=C2ND*RHO2*SPRF( 338)/WM( 33)/WM(  5)
      RC(IA, 339)=C2ND*RHO2*SPRF( 339)/WM( 31)/WM(  4)
      RC(IA, 340)=C2ND*RHO2*SPRF( 340)/WM( 36)/WM(  5)
      RC(IA, 341)=C2ND*RHO2*SPRF( 341)/WM( 31)/WM( 14)
      RC(IA, 342)=C2ND*RHO2*SPRF( 342)/WM( 34)/WM( 13)
      RC(IA, 343)=C2ND*RHO2*SPRF( 343)/WM( 31)/WM( 14)
      RC(IA, 344)=C2ND*RHO2*SPRF( 344)/WM( 33)/WM( 13)
      RC(IA, 345)=C2ND*RHO2*SPRF( 345)/WM( 31)/WM( 14)
      RC(IA, 346)=C2ND*RHO2*SPRF( 346)/WM( 36)/WM( 13)
      RC(IA, 347)=C2ND*RHO2*SPRF( 347)/WM( 31)/WM(  7)
      RC(IA, 348)=C2ND*RHO2*SPRF( 348)/WM( 33)/WM(  8)
      RC(IA, 349)=C2ND*RHO2*SPRF( 349)/WM( 31)/WM(  7)
      RC(IA, 350)=C2ND*RHO2*SPRF( 350)/WM( 34)/WM(  8)
      RC(IA, 351)=C2ND*RHO2*SPRF( 351)/WM( 31)/WM(  7)
      RC(IA, 352)=C2ND*RHO2*SPRF( 352)/WM( 36)/WM(  8)
      RC(IA, 353)=C2ND*RHO2*SPRF( 353)/WM( 17)/WM(  5)
      RC(IA, 354)=C1ST*RHO1*SPRF( 354)/WM( 34)
      RC(IA, 355)=C2ND*RHO2*SPRF( 355)/WM( 19)/WM(  7)
      RC(IA, 356)=C2ND*RHO2*SPRF( 356)/WM( 36)/WM(  5)
      RC(IA, 357)=C2ND*RHO2*TB009*SPRF( 357)/WM( 36)
      RC(IA, 358)=C3RD*RHO3*TB009*SPRF( 358)/WM( 32)/WM(  3)
      RC(IA, 359)=C2ND*RHO2*TB009*SPRF( 359)/WM( 36)
      RC(IA, 360)=C3RD*RHO3*TB009*SPRF( 360)/WM( 14)/WM( 12)
      RC(IA, 361)=C2ND*RHO2*SPRF( 361)/WM( 36)/WM(  2)
      RC(IA, 362)=C2ND*RHO2*SPRF( 362)/WM( 32)/WM(  7)
      RC(IA, 363)=C2ND*RHO2*SPRF( 363)/WM( 36)/WM(  9)
      RC(IA, 364)=C2ND*RHO2*SPRF( 364)/WM( 19)/WM( 10)
      RC(IA, 365)=C2ND*RHO2*SPRF( 365)/WM( 36)/WM(  3)
      RC(IA, 366)=C2ND*RHO2*SPRF( 366)/WM( 14)/WM( 29)
      RC(IA, 367)=C2ND*RHO2*SPRF( 367)/WM( 36)/WM(  3)
      RC(IA, 368)=C2ND*RHO2*SPRF( 368)/WM( 17)/WM(  6)
      RC(IA, 369)=C2ND*RHO2*SPRF( 369)/WM( 36)/WM(  5)
      RC(IA, 370)=C2ND*RHO2*SPRF( 370)/WM( 32)/WM(  6)
      RC(IA, 371)=C2ND*RHO2*SPRF( 371)/WM( 33)/WM(  2)
      RC(IA, 372)=C2ND*RHO2*SPRF( 372)/WM( 32)/WM(  7)
      RC(IA, 373)=C2ND*RHO2*SPRF( 373)/WM( 33)/WM(  4)
      RC(IA, 374)=C2ND*RHO2*SPRF( 374)/WM( 32)/WM(  5)
      RC(IA, 375)=C2ND*RHO2*SPRF( 375)/WM( 33)/WM(  3)
      RC(IA, 376)=C2ND*RHO2*SPRF( 376)/WM( 17)/WM(  6)
      RC(IA, 377)=C2ND*RHO2*SPRF( 377)/WM( 33)/WM(  3)
      RC(IA, 378)=C2ND*RHO2*SPRF( 378)/WM( 14)/WM( 29)
      RC(IA, 379)=C2ND*RHO2*SPRF( 379)/WM( 33)/WM(  7)
      RC(IA, 380)=C3RD*RHO3*SPRF( 380)/WM( 32)/WM(  5)/WM(  5)
      RC(IA, 381)=C2ND*RHO2*SPRF( 381)/WM( 33)/WM(  5)
      RC(IA, 382)=C2ND*RHO2*SPRF( 382)/WM( 32)/WM(  6)
      RC(IA, 383)=C2ND*RHO2*TB009*SPRF( 383)/WM( 33)
      RC(IA, 384)=C3RD*RHO3*TB009*SPRF( 384)/WM( 32)/WM(  3)
      RC(IA, 385)=C2ND*RHO2*SPRF( 385)/WM( 37)/WM(  4)
      RC(IA, 386)=C2ND*RHO2*SPRF( 386)/WM( 17)/WM(  9)
      RC(IA, 387)=C2ND*RHO2*SPRF( 387)/WM( 14)/WM( 22)
      RC(IA, 388)=C2ND*RHO2*SPRF( 388)/WM( 37)/WM(  3)
      RC(IA, 389)=C2ND*RHO2*SPRF( 389)/WM( 37)/WM(  4)
      RC(IA, 390)=C2ND*RHO2*SPRF( 390)/WM( 27)/WM( 14)
      RC(IA, 391)=C2ND*RHO2*SPRF( 391)/WM( 38)/WM(  3)
      RC(IA, 392)=C1ST*RHO1*SPRF( 392)/WM( 37)
      RC(IA, 393)=C2ND*RHO2*SPRF( 393)/WM( 38)/WM(  7)
      RC(IA, 394)=C2ND*RHO2*SPRF( 394)/WM( 37)/WM(  2)
      RC(IA, 395)=C2ND*RHO2*SPRF( 395)/WM( 37)/WM(  5)
      RC(IA, 396)=C2ND*RHO2*SPRF( 396)/WM( 38)/WM(  6)
      RC(IA, 397)=C2ND*RHO2*SPRF( 397)/WM( 38)/WM(  2)
      RC(IA, 398)=C2ND*RHO2*SPRF( 398)/WM( 26)/WM( 11)
      RC(IA, 399)=C2ND*RHO2*SPRF( 399)/WM( 37)/WM(  3)
      RC(IA, 400)=C1ST*RHO1*SPRF( 400)/WM( 39)
      RC(IA, 401)=C2ND*RHO2*SPRF( 401)/WM( 39)/WM(  3)
      RC(IA, 402)=C2ND*RHO2*SPRF( 402)/WM( 37)/WM(  1)
      RC(IA, 403)=C2ND*RHO2*SPRF( 403)/WM( 39)/WM(  2)
      RC(IA, 404)=C2ND*RHO2*SPRF( 404)/WM( 37)/WM(  7)
      RC(IA, 405)=C2ND*RHO2*SPRF( 405)/WM( 39)/WM( 14)
      RC(IA, 406)=C2ND*RHO2*SPRF( 406)/WM( 37)/WM( 13)
      RC(IA, 407)=C2ND*RHO2*SPRF( 407)/WM( 22)/WM( 14)
      RC(IA, 408)=C1ST*RHO1*SPRF( 408)/WM( 39)
      RC(IA, 409)=C2ND*RHO2*SPRF( 409)/WM( 39)/WM(  5)
      RC(IA, 410)=C2ND*RHO2*SPRF( 410)/WM( 37)/WM(  6)
      RC(IA, 411)=C2ND*RHO2*SPRF( 411)/WM( 38)/WM( 11)
      RC(IA, 412)=C2ND*RHO2*SPRF( 412)/WM( 37)/WM(  9)
      RC(IA, 413)=C2ND*RHO2*SPRF( 413)/WM( 38)/WM(  7)
      RC(IA, 414)=C3RD*RHO3*SPRF( 414)/WM(  5)/WM(  9)/WM( 23)
      RC(IA, 415)=C2ND*RHO2*SPRF( 415)/WM( 37)/WM(  2)
      RC(IA, 416)=C3RD*RHO3*SPRF( 416)/WM( 14)/WM( 11)/WM(  9)
      RC(IA, 417)=C2ND*RHO2*SPRF( 417)/WM( 40)/WM(  4)
      RC(IA, 418)=C2ND*RHO2*SPRF( 418)/WM( 19)/WM( 11)
      RC(IA, 419)=C2ND*RHO2*SPRF( 419)/WM( 40)/WM(  5)
      RC(IA, 420)=C2ND*RHO2*SPRF( 420)/WM( 39)/WM(  6)
      RC(IA, 421)=C2ND*RHO2*SPRF( 421)/WM( 40)/WM(  4)
      RC(IA, 422)=C3RD*RHO3*SPRF( 422)/WM( 26)/WM( 14)/WM(  3)
      RC(IA, 423)=C2ND*RHO2*SPRF( 423)/WM( 40)/WM(  3)
      RC(IA, 424)=C2ND*RHO2*SPRF( 424)/WM( 39)/WM(  1)
      RC(IA, 425)=C2ND*RHO2*SPRF( 425)/WM( 39)/WM(  3)
      RC(IA, 426)=C1ST*RHO1*SPRF( 426)/WM( 40)
      RC(IA, 427)=C2ND*RHO2*SPRF( 427)/WM( 39)/WM(  7)
      RC(IA, 428)=C2ND*RHO2*SPRF( 428)/WM( 40)/WM(  2)
      RC(IA, 429)=C2ND*RHO2*SPRF( 429)/WM( 39)/WM(  7)
      RC(IA, 430)=C3RD*RHO3*SPRF( 430)/WM(  5)/WM( 23)/WM( 12)
      RC(IA, 431)=C2ND*RHO2*SPRF( 431)/WM( 23)/WM( 14)
      RC(IA, 432)=C1ST*RHO1*SPRF( 432)/WM( 40)
      RC(IA, 433)=C2ND*RHO2*SPRF( 433)/WM( 40)/WM(  3)
      RC(IA, 434)=C2ND*RHO2*SPRF( 434)/WM( 17)/WM( 14)
      RC(IA, 435)=C2ND*RHO2*SPRF( 435)/WM( 14)/WM( 23)
      RC(IA, 436)=C2ND*RHO2*SPRF( 436)/WM( 39)/WM(  3)
      RC(IA, 437)=C1ST*RHO1*SPRF( 437)/WM( 41)
      RC(IA, 438)=C2ND*RHO2*SPRF( 438)/WM( 14)/WM( 19)
      RC(IA, 439)=C2ND*RHO2*SPRF( 439)/WM( 41)/WM(  2)
      RC(IA, 440)=C2ND*RHO2*SPRF( 440)/WM( 42)/WM(  7)
      RC(IA, 441)=C2ND*RHO2*SPRF( 441)/WM( 41)/WM(  2)
      RC(IA, 442)=C2ND*RHO2*SPRF( 442)/WM( 43)/WM(  7)
      RC(IA, 443)=C2ND*RHO2*SPRF( 443)/WM( 41)/WM(  3)
      RC(IA, 444)=C2ND*RHO2*SPRF( 444)/WM( 42)/WM(  1)
      RC(IA, 445)=C2ND*RHO2*SPRF( 445)/WM( 41)/WM(  3)
      RC(IA, 446)=C2ND*RHO2*SPRF( 446)/WM( 43)/WM(  1)
      RC(IA, 447)=C2ND*RHO2*SPRF( 447)/WM( 41)/WM(  4)
      RC(IA, 448)=C2ND*RHO2*SPRF( 448)/WM( 42)/WM(  5)
      RC(IA, 449)=C2ND*RHO2*SPRF( 449)/WM( 41)/WM(  4)
      RC(IA, 450)=C2ND*RHO2*SPRF( 450)/WM( 43)/WM(  5)
      RC(IA, 451)=C2ND*RHO2*SPRF( 451)/WM( 41)/WM(  5)
      RC(IA, 452)=C2ND*RHO2*SPRF( 452)/WM( 43)/WM(  6)
      RC(IA, 453)=C2ND*RHO2*SPRF( 453)/WM( 41)/WM(  5)
      RC(IA, 454)=C2ND*RHO2*SPRF( 454)/WM( 42)/WM(  6)
      RC(IA, 455)=C2ND*RHO2*SPRF( 455)/WM( 41)/WM(  7)
      RC(IA, 456)=C2ND*RHO2*SPRF( 456)/WM( 42)/WM(  8)
      RC(IA, 457)=C2ND*RHO2*SPRF( 457)/WM( 41)/WM(  7)
      RC(IA, 458)=C2ND*RHO2*SPRF( 458)/WM( 43)/WM(  8)
      RC(IA, 459)=C2ND*RHO2*SPRF( 459)*FC0(J,I)/WM( 42)/WM( 41)
      RC(IA, 460)=C2ND*RHO2*SPRF( 460)*FC0(J,I)/WM( 43)/WM( 41)
      RC(IA, 461)=C2ND*RHO2*SPRF( 461)/WM( 40)/WM(  3)
      RC(IA, 462)=C1ST*RHO1*SPRF( 462)/WM( 42)
      RC(IA, 463)=C2ND*RHO2*SPRF( 463)/WM( 42)/WM(  2)
      RC(IA, 464)=C2ND*RHO2*SPRF( 464)/WM( 40)/WM(  7)
      RC(IA, 465)=C1ST*RHO1*SPRF( 465)/WM( 43)
      RC(IA, 466)=C2ND*RHO2*SPRF( 466)/WM( 14)/WM( 17)
      RC(IA, 467)=C2ND*RHO2*SPRF( 467)/WM(  3)/WM( 40)
      RC(IA, 468)=C1ST*RHO1*SPRF( 468)/WM( 43)
      RC(IA, 469)=C2ND*RHO2*SPRF( 469)/WM( 43)/WM(  2)
      RC(IA, 470)=C2ND*RHO2*SPRF( 470)/WM( 40)/WM(  7)
      RC(IA, 471)=C1ST*RHO1*SPRF( 471)/WM( 44)
      RC(IA, 472)=C3RD*RHO3*SPRF( 472)/WM( 22)/WM( 23)/WM(  3)
      RC(IA, 473)=C1ST*RHO1*SPRF( 473)/WM( 44)
      RC(IA, 474)=C2ND*RHO2*SPRF( 474)/WM( 23)/WM( 23)
      RC(IA, 475)=C2ND*RHO2*SPRF( 475)/WM( 23)/WM( 23)
      RC(IA, 476)=C1ST*RHO1*SPRF( 476)/WM( 44)
      RC(IA, 477)=C2ND*RHO2*SPRF( 477)/WM( 44)/WM(  3)
      RC(IA, 478)=C2ND*RHO2*SPRF( 478)/WM( 23)/WM( 17)
      RC(IA, 479)=C2ND*RHO2*SPRF( 479)/WM( 44)/WM(  3)
      RC(IA, 480)=C3RD*RHO3*SPRF( 480)/WM(  1)/WM( 22)/WM( 23)
      RC(IA, 481)=C2ND*RHO2*SPRF( 481)/WM( 44)/WM(  5)
      RC(IA, 482)=C3RD*RHO3*SPRF( 482)/WM( 45)/WM(  3)/WM( 39)
      RC(IA, 483)=C2ND*RHO2*SPRF( 483)/WM( 44)/WM( 14)
      RC(IA, 484)=C3RD*RHO3*SPRF( 484)/WM( 13)/WM( 22)/WM( 23)
      RC(IA, 485)=C2ND*RHO2*SPRF( 485)/WM( 38)/WM( 14)
      RC(IA, 486)=C1ST*RHO1*SPRF( 486)/WM( 44)
      RC(IA, 487)=C1ST*RHO1*SPRF( 487)/WM( 46)
      RC(IA, 488)=C2ND*RHO2*SPRF( 488)/WM( 40)/WM( 22)
      RC(IA, 489)=C1ST*RHO1*SPRF( 489)/WM( 46)
      RC(IA, 490)=C2ND*RHO2*SPRF( 490)/WM( 37)/WM( 17)
      RC(IA, 491)=C1ST*RHO1*SPRF( 491)/WM( 46)
      RC(IA, 492)=C2ND*RHO2*SPRF( 492)/WM( 39)/WM( 23)
      RC(IA, 493)=C2ND*RHO2*SPRF( 493)/WM( 46)/WM(  2)
      RC(IA, 494)=C3RD*RHO3*SPRF( 494)/WM( 22)/WM( 39)/WM(  7)
      RC(IA, 495)=C2ND*RHO2*SPRF( 495)/WM( 46)/WM(  2)
      RC(IA, 496)=C3RD*RHO3*SPRF( 496)/WM( 23)/WM( 37)/WM(  7)
      RC(IA, 497)=C2ND*RHO2*SPRF( 497)/WM( 46)/WM(  7)
      RC(IA, 498)=C3RD*RHO3*SPRF( 498)/WM( 22)/WM( 39)/WM(  8)
      RC(IA, 499)=C2ND*RHO2*SPRF( 499)/WM( 46)/WM(  7)
      RC(IA, 500)=C3RD*RHO3*SPRF( 500)/WM( 23)/WM( 37)/WM(  8)
      RC(IA, 501)=C1ST*RHO1*SPRF( 501)/WM( 47)
      RC(IA, 503)=C2ND*RHO2*SPRF( 503)/WM( 47)/WM(  3)
      RC(IA, 505)=C2ND*RHO2*SPRF( 505)/WM( 47)/WM(  3)
      RC(IA, 507)=C2ND*RHO2*SPRF( 507)/WM( 47)/WM(  3)
      RC(IA, 509)=C2ND*RHO2*SPRF( 509)/WM( 47)/WM(  3)
      RC(IA, 510)=C3RD*RHO3*SPRF( 510)/WM(  1)/WM( 49)/WM( 19)
      RC(IA, 511)=C2ND*RHO2*SPRF( 511)/WM( 47)/WM(  5)
      RC(IA, 512)=C3RD*RHO3*SPRF( 512)/WM(  6)/WM( 49)/WM( 19)
      RC(IA, 513)=C2ND*RHO2*SPRF( 513)/WM( 47)/WM(  5)
      RC(IA, 515)=C2ND*RHO2*SPRF( 515)/WM( 47)/WM(  5)
      RC(IA, 517)=C2ND*RHO2*SPRF( 517)/WM( 47)/WM(  5)
      RC(IA, 519)=C1ST*RHO1*SPRF( 519)/WM( 49)
      RC(IA, 520)=C2ND*RHO2*SPRF( 520)/WM( 39)/WM( 19)
      RC(IA, 521)=C1ST*RHO1*SPRF( 521)/WM( 49)
      RC(IA, 522)=C2ND*RHO2*SPRF( 522)/WM( 40)/WM( 17)
      RC(IA, 523)=C2ND*RHO2*SPRF( 523)/WM( 49)/WM(  5)
      RC(IA, 524)=C3RD*RHO3*SPRF( 524)/WM(  6)/WM( 40)/WM( 23)
      RC(IA, 525)=C2ND*RHO2*SPRF( 525)/WM( 49)/WM(  3)
      RC(IA, 526)=C3RD*RHO3*SPRF( 526)/WM(  1)/WM( 40)/WM( 23)
      RC(IA, 527)=C2ND*RHO2*SPRF( 527)/WM( 49)/WM(  3)
      RC(IA, 528)=C3RD*RHO3*SPRF( 528)/WM( 17)/WM( 17)/WM( 14)
      RC(IA, 529)=C2ND*RHO2*SPRF( 529)/WM( 49)/WM(  3)
      RC(IA, 530)=C2ND*RHO2*SPRF( 530)/WM( 40)/WM( 19)
      RC(IA, 531)=C2ND*RHO2*SPRF( 531)/WM( 49)/WM(  3)
      RC(IA, 532)=C3RD*RHO3*SPRF( 532)/WM(  1)/WM( 17)/WM( 39)
      RC(IA, 533)=C2ND*RHO2*SPRF( 533)/WM( 49)/WM(  3)
      RC(IA, 534)=C3RD*RHO3*SPRF( 534)/WM(  1)/WM( 44)/WM( 14)
      RC(IA, 535)=C1ST*RHO1*SPRF( 535)/WM( 48)
      RC(IA, 536)=C2ND*RHO2*SPRF( 536)/WM( 39)/WM( 14)
      RC(IA, 537)=C2ND*RHO2*SPRF( 537)/WM( 48)/WM(  3)
      RC(IA, 538)=C2ND*RHO2*SPRF( 538)/WM( 17)/WM( 19)
      RC(IA, 539)=C2ND*RHO2*SPRF( 539)/WM( 48)/WM(  3)
      RC(IA, 540)=C2ND*RHO2*SPRF( 540)/WM( 40)/WM( 14)
      RC(IA, 541)=C2ND*RHO2*SPRF( 541)/WM( 48)/WM(  3)
      RC(IA, 542)=C3RD*RHO3*SPRF( 542)/WM(  1)/WM( 23)/WM( 17)
      RC(IA, 543)=C2ND*RHO2*SPRF( 543)/WM( 48)/WM(  5)
      RC(IA, 544)=C3RD*RHO3*SPRF( 544)/WM(  6)/WM( 44)/WM(  3)
C-------------------Source Terms in Soot Equations----------------------
C  Source terms 1,3,4,5 go into SMF equations & 2,6 go into SND equation
C--------C2H2 (22), O2 (02), OH (05), H (03), O (04), H2 (01), CO (09)
      DPST=0.0
      IF(STND(J,I).GE.1.0) 
     1   DPST=(6.0*STMF(J,I)/(PI*RSOOT*STND(J,I)))**0.333333
      ASST=PI*DPST*DPST*RSTR*RHO(J,I)*STND(J,I)
      SOURNU=C1ST*RHO1*AKST1*DEXP(-EAST1/TKD)/WM( 22)
      SOUR(J,I,1)=SOURNU
     1           +C1ST*RHO1*ASST*AKST2*DEXP(-EAST2/TKD)/WM( 22)
      SOUR(J,I,2)=2.0*SOURNU*AVNO/CMIN/WMSTR
      SOUR(J,I,3)=C1ST*RHO1*ASST*AKST3
     1           *(TKD**ALST3)*DEXP(-EAST3/TKD)/WM( 02)
      SOUR(J,I,4)=0.13*C1ST*ASST*AKST4
     1           *(TKD**ALST4)*DEXP(-EAST4/TKD)/WM( 05)/TBODY
      SOUR(J,I,5)=0.5*C1ST*ASST*AKST5
     1           *(TKD**ALST5)*DEXP(-EAST5/TKD)/WM( 04)/TBODY
      SOUR(J,I,6)=2.0*C1ST*RSTR*CAGL*DSQRT(DPST*6.0*SIGB*TKD/RSOOT)
     1           *RHO(J,I)*RHO(J,I)*STND(J,I)
C     SOUR(J,I,1)=0.0
C     SOUR(J,I,2)=0.0
C     SOUR(J,I,3)=0.0
C     SOUR(J,I,4)=0.0
C     SOUR(J,I,5)=0.0
C     SOUR(J,I,6)=0.0
  203 CONTINUE
      ENDDO
  202 CONTINUE
CCCCCCCCCCCCCCCCCCC     END REACTION RATES    CCCCCCCCCCCCCCCCCCCCCCCCCC
      RETURN
      END
C***********************************************************************
      SUBROUTINE SOLVESP(ISOR1,RELXSP,TOLRSP,SPNORM)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LJ2=LJ*2,LSP=52,LRX=544)
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FA0(LJ,LI),FA1(LJ,LI),FA2(LJ,LI),
     *  FA3(LJ,LI),FA4(LJ,LI),FA5(LJ,LI),FA6(LJ,LI),FA7(LJ,LI),
     *  FA8(LJ,LI),FA9(LJ,LI),FAA(LJ,LI),FAB(LJ,LI),FAC(LJ,LI),
     *  FAD(LJ,LI),FAE(LJ,LI),FAF(LJ,LI),FAG(LJ,LI),FAH(LJ,LI),
     *  FAI(LJ,LI),FAJ(LJ,LI),FB0(LJ,LI),FB1(LJ,LI),FB2(LJ,LI),
     *  FB3(LJ,LI),FB4(LJ,LI),FB5(LJ,LI),FB6(LJ,LI),FB7(LJ,LI),
     *  FB8(LJ,LI),FB9(LJ,LI),FBA(LJ,LI),FBB(LJ,LI),FBC(LJ,LI),
     *  FBD(LJ,LI),FBE(LJ,LI),FBF(LJ,LI),FBG(LJ,LI),FBH(LJ,LI),
     *  FBI(LJ,LI),FBJ(LJ,LI),FC0(LJ,LI),FC1(LJ,LI),FC2(LJ,LI),
     *  FC3(LJ,LI),FC4(LJ,LI),FC5(LJ,LI),FC6(LJ,LI),FC7(LJ,LI),
     *  FC8(LJ,LI),FC9(LJ,LI),FCA(LJ,LI),FCB(LJ,LI),
     *  U(LJ,LI),V(LJ,LI),W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),
     *  AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
      COMMON/DUMMY/
     *  EA0(LJ,LI),WA0(LJ,LI),PA0(LJ,LI),QA0(LJ,LI),EA1(LJ,LI),
     *  WA1(LJ,LI),PA1(LJ,LI),QA1(LJ,LI),EA2(LJ,LI),WA2(LJ,LI),
     *  PA2(LJ,LI),QA2(LJ,LI),EA3(LJ,LI),WA3(LJ,LI),PA3(LJ,LI),
     *  QA3(LJ,LI),EA4(LJ,LI),WA4(LJ,LI),PA4(LJ,LI),QA4(LJ,LI),
     *  EA5(LJ,LI),WA5(LJ,LI),PA5(LJ,LI),QA5(LJ,LI),EA6(LJ,LI),
     *  WA6(LJ,LI),PA6(LJ,LI),QA6(LJ,LI),EA7(LJ,LI),WA7(LJ,LI),
     *  PA7(LJ,LI),QA7(LJ,LI),EA8(LJ,LI),WA8(LJ,LI),PA8(LJ,LI),
     *  QA8(LJ,LI),EA9(LJ,LI),WA9(LJ,LI),PA9(LJ,LI),QA9(LJ,LI),
     *  EAA(LJ,LI),WAA(LJ,LI),PAA(LJ,LI),QAA(LJ,LI),EAB(LJ,LI),
     *  WAB(LJ,LI),PAB(LJ,LI),QAB(LJ,LI),EAC(LJ,LI),WAC(LJ,LI),
     *  PAC(LJ,LI),QAC(LJ,LI),EAD(LJ,LI),WAD(LJ,LI),PAD(LJ,LI)
      COMMON/DUM0/
     *  QAD(LJ,LI),EAE(LJ,LI),WAE(LJ,LI),PAE(LJ,LI),QAE(LJ,LI),
     *  EAF(LJ,LI),WAF(LJ,LI),PAF(LJ,LI),QAF(LJ,LI),EAG(LJ,LI),
     *  WAG(LJ,LI),PAG(LJ,LI),QAG(LJ,LI),EAH(LJ,LI),WAH(LJ,LI),
     *  PAH(LJ,LI),QAH(LJ,LI),EAI(LJ,LI),WAI(LJ,LI),PAI(LJ,LI),
     *  QAI(LJ,LI),EAJ(LJ,LI),WAJ(LJ,LI),PAJ(LJ,LI),QAJ(LJ,LI),
     *  EB0(LJ,LI),WB0(LJ,LI),PB0(LJ,LI),QB0(LJ,LI),EB1(LJ,LI),
     *  WB1(LJ,LI),PB1(LJ,LI),QB1(LJ,LI),EB2(LJ,LI),WB2(LJ,LI),
     *  PB2(LJ,LI),QB2(LJ,LI),EB3(LJ,LI),WB3(LJ,LI),PB3(LJ,LI),
     *  QB3(LJ,LI),EB4(LJ,LI),WB4(LJ,LI),PB4(LJ,LI),QB4(LJ,LI),
     *  EB5(LJ,LI),WB5(LJ,LI),PB5(LJ,LI),QB5(LJ,LI),EB6(LJ,LI),
     *  WB6(LJ,LI),PB6(LJ,LI),QB6(LJ,LI),EB7(LJ,LI),WB7(LJ,LI),
     *  PB7(LJ,LI),QB7(LJ,LI),EB8(LJ,LI),WB8(LJ,LI),PB8(LJ,LI),
     *  QB8(LJ,LI),EB9(LJ,LI),WB9(LJ,LI),PB9(LJ,LI),QB9(LJ,LI),
     *  EBA(LJ,LI),WBA(LJ,LI),PBA(LJ,LI),QBA(LJ,LI),EBB(LJ,LI),
     *  WBB(LJ,LI),PBB(LJ,LI),QBB(LJ,LI),EBC(LJ,LI),WBC(LJ,LI),
     *  PBC(LJ,LI),QBC(LJ,LI),EBD(LJ,LI),WBD(LJ,LI),PBD(LJ,LI),
     *  QBD(LJ,LI),EBE(LJ,LI),WBE(LJ,LI),PBE(LJ,LI),QBE(LJ,LI),
     *  EBF(LJ,LI),WBF(LJ,LI),PBF(LJ,LI),QBF(LJ,LI),EBG(LJ,LI),
     *  WBG(LJ,LI),PBG(LJ,LI),QBG(LJ,LI),EBH(LJ,LI),WBH(LJ,LI),
     *  PBH(LJ,LI),QBH(LJ,LI),EBI(LJ,LI),WBI(LJ,LI),PBI(LJ,LI),
     *  QBI(LJ,LI),EBJ(LJ,LI),WBJ(LJ,LI),PBJ(LJ,LI),QBJ(LJ,LI),
     *  EC0(LJ,LI),WC0(LJ,LI),PC0(LJ,LI),QC0(LJ,LI),EC1(LJ,LI),
     *  WC1(LJ,LI),PC1(LJ,LI),QC1(LJ,LI),EC2(LJ,LI),WC2(LJ,LI),
     *  PC2(LJ,LI),QC2(LJ,LI),EC3(LJ,LI),WC3(LJ,LI),PC3(LJ,LI),
     *  QC3(LJ,LI),EC4(LJ,LI),WC4(LJ,LI),PC4(LJ,LI),QC4(LJ,LI),
     *  EC5(LJ,LI),WC5(LJ,LI),PC5(LJ,LI),QC5(LJ,LI),EC6(LJ,LI),
     *  WC6(LJ,LI),PC6(LJ,LI),QC6(LJ,LI),EC7(LJ,LI),WC7(LJ,LI),
     *  PC7(LJ,LI),QC7(LJ,LI),EC8(LJ,LI),WC8(LJ,LI),PC8(LJ,LI),
     *  QC8(LJ,LI),EC9(LJ,LI),WC9(LJ,LI),PC9(LJ,LI),QC9(LJ,LI),
     *  ECA(LJ,LI),WCA(LJ,LI),PCA(LJ,LI),QCA(LJ,LI)
      COMMON/DUMX/ ZE(LJ,LI),ZW(LJ,LI),ZN(LJ,LI),ZS(LJ,LI)
      COMMON/DUMY/            SA0(LJ,LI),SA1(LJ,LI),SA2(LJ,LI),
     *  SA3(LJ,LI),SA4(LJ,LI),SA5(LJ,LI),SA6(LJ,LI),SA7(LJ,LI),
     *  SA8(LJ,LI),SA9(LJ,LI),SAA(LJ,LI),SAB(LJ,LI),SAC(LJ,LI),
     *  SAD(LJ,LI),SAE(LJ,LI),SAF(LJ,LI),SAG(LJ,LI),SAH(LJ,LI),
     *  SAI(LJ,LI),SAJ(LJ,LI),SB0(LJ,LI),SB1(LJ,LI),SB2(LJ,LI),
     *  SB3(LJ,LI),SB4(LJ,LI),SB5(LJ,LI),SB6(LJ,LI),SB7(LJ,LI),
     *  SB8(LJ,LI),SB9(LJ,LI),SBA(LJ,LI),SBB(LJ,LI),SBC(LJ,LI),
     *  SBD(LJ,LI),SBE(LJ,LI),SBF(LJ,LI),SBG(LJ,LI),SBH(LJ,LI),
     *  SBI(LJ,LI),SBJ(LJ,LI),SC0(LJ,LI),SC1(LJ,LI),SC2(LJ,LI),
     *  SC3(LJ,LI),SC4(LJ,LI),SC5(LJ,LI),SC6(LJ,LI),SC7(LJ,LI),
     *  SC8(LJ,LI),SC9(LJ,LI),SCA(LJ,LI),SCB(LJ,LI),
     *  FPZ(LJ,LI),EMU(LJ,LI),EKT(LJ,LI)
      COMMON/DUMS/ STE(LJ,LI),STW(LJ,LI),STN(LJ,LI),STS(LJ,LI),
     1  SST1(LJ,LI),SST2(LJ,LI),SOUR(LJ,LI,6)
      COMMON/DUMR/ RC(LJ,LI,LRX)
      SOR1=RELXSP
      SOR2=(1.0-SOR1)
      DO 300 ISORP=1,ISOR1
      RSD=0.0
      DO 304 J=2,LJ-1
      DO 305 I=2,LI-1
      IF(ISKIP(J,I).GE.3) GO TO 305
C----------------------  POINT RELAXATION SCHEME  ----------------------
      A0002=SA4(J,I)*SA3(J,I)*RC(J,I,   2)
      A0004=SA4(J,I)*SA2(J,I)*RC(J,I,   4)
      A0006=SA5(J,I)*SA2(J,I)*RC(J,I,   6)
      A0008=SA4(J,I)*SA4(J,I)*RC(J,I,   8)
      A0009=SA2(J,I)*SA2(J,I)*RC(J,I,   9)
      A0012=SA5(J,I)*RC(J,I,  12)
      A0013=SA3(J,I)*SA3(J,I)*RC(J,I,  13)
      A0016=SA4(J,I)*RC(J,I,  16)
      A0018=SA6(J,I)*RC(J,I,  18)
      A0020=SA6(J,I)*RC(J,I,  20)
      A0022=SA4(J,I)*SA4(J,I)*RC(J,I,  22)
      A0023=SA6(J,I)*SA2(J,I)*RC(J,I,  23)
      A0026=SA5(J,I)*SA3(J,I)*RC(J,I,  26)
      A0027=SA6(J,I)*SA3(J,I)*RC(J,I,  27)
      A0029=SA6(J,I)*SA4(J,I)*RC(J,I,  29)
      A0032=SA7(J,I)*RC(J,I,  32)
      A0033=SA6(J,I)*SA6(J,I)*RC(J,I,  33)
      A0035=SA7(J,I)*SA2(J,I)*RC(J,I,  35)
      A0038=SA5(J,I)*SA4(J,I)*RC(J,I,  38)
      A0040=SA5(J,I)*SA6(J,I)*RC(J,I,  40)
      A0042=SA6(J,I)*SA4(J,I)*RC(J,I,  42)
      A0043=SA8(J,I)*SA4(J,I)*RC(J,I,  43)
      A0045=SA8(J,I)*SA6(J,I)*RC(J,I,  45)
      A0048=SA9(J,I)*SA3(J,I)*RC(J,I,  48)
      A0049=SAA(J,I)*RC(J,I,  49)
      A0051=SAA(J,I)*SA2(J,I)*RC(J,I,  51)
      A0054=SA8(J,I)*SA4(J,I)*RC(J,I,  54)
      A0055=SAA(J,I)*SA3(J,I)*RC(J,I,  55)
      A0058=SA8(J,I)*SA5(J,I)*RC(J,I,  58)
      A0060=SA8(J,I)*SA6(J,I)*RC(J,I,  60)
      A0061=SAA(J,I)*SAD(J,I)*RC(J,I,  61)
      A0064=SAB(J,I)*RC(J,I,  64)
      A0065=SAB(J,I)*SA2(J,I)*RC(J,I,  65)
      A0068=SAA(J,I)*SA4(J,I)*RC(J,I,  68)
      A0070=SAA(J,I)*SA5(J,I)*RC(J,I,  70)
      A0072=SAA(J,I)*SA6(J,I)*RC(J,I,  72)
      A0074=SAA(J,I)*SA7(J,I)*RC(J,I,  74)
      A0075=SAC(J,I)*SA2(J,I)*RC(J,I,  75)
      A0078=SA5(J,I)*SAD(J,I)*RC(J,I,  78)
      A0080=SAD(J,I)*SA4(J,I)*RC(J,I,  80)
      A0082=SAD(J,I)*SA6(J,I)*RC(J,I,  82)
      A0084=SAD(J,I)*SA7(J,I)*RC(J,I,  84)
      A0085=SAD(J,I)*SA2(J,I)*RC(J,I,  85)
      A0087=SAD(J,I)*SA2(J,I)*RC(J,I,  87)
      A0090=SAF(J,I)*SA5(J,I)*RC(J,I,  90)
      A0091=SAD(J,I)*SA3(J,I)*RC(J,I,  91)
      A0093=SAD(J,I)*SAE(J,I)*RC(J,I,  93)
      A0095=SAD(J,I)*SA6(J,I)*RC(J,I,  95)
      A0098=SAB(J,I)*SA4(J,I)*RC(J,I,  98)
      A0100=SAH(J,I)*SA3(J,I)*RC(J,I, 100)
      A0101=SAD(J,I)*SAD(J,I)*RC(J,I, 101)
      A0103=SAD(J,I)*SAD(J,I)*RC(J,I, 103)
      A0106=SAC(J,I)*RC(J,I, 106)
      A0108=SAJ(J,I)*RC(J,I, 108)
      A0109=SAF(J,I)*SA4(J,I)*RC(J,I, 109)
      A0112=SA8(J,I)*SA4(J,I)*SA2(J,I)*RC(J,I, 112)
      A0113=SAF(J,I)*SA9(J,I)*RC(J,I, 113)
      A0115=SAF(J,I)*RC(J,I, 115)
      A0117=SAE(J,I)*SA2(J,I)*RC(J,I, 117)
      A0119=SAE(J,I)*SA4(J,I)*RC(J,I, 119)
      A0122=SB0(J,I)*SA5(J,I)*RC(J,I, 122)
      A0123=SAE(J,I)*SA3(J,I)*RC(J,I, 123)
      A0125=SAE(J,I)*SA3(J,I)*RC(J,I, 125)
      A0127=SAE(J,I)*SA1(J,I)*RC(J,I, 127)
      A0130=SA8(J,I)*SA4(J,I)*SA2(J,I)*RC(J,I, 130)
      A0131=SAE(J,I)*SAE(J,I)*RC(J,I, 131)
      A0133=SB0(J,I)*SA3(J,I)*RC(J,I, 133)
      A0136=SAA(J,I)*SA3(J,I)*RC(J,I, 136)
      A0137=SB0(J,I)*SA5(J,I)*RC(J,I, 137)
      A0139=SB0(J,I)*SA9(J,I)*RC(J,I, 139)
      A0141=SAH(J,I)*SA2(J,I)*RC(J,I, 141)
      A0144=SAF(J,I)*SA5(J,I)*RC(J,I, 144)
      A0146=SAB(J,I)*SA5(J,I)*RC(J,I, 146)
      A0148=SA4(J,I)*SAB(J,I)*RC(J,I, 148)
      A0150=SAB(J,I)*SA6(J,I)*RC(J,I, 150)
      A0151=SAH(J,I)*RC(J,I, 151)
      A0153=SAJ(J,I)*SA2(J,I)*RC(J,I, 153)
      A0156=SAI(J,I)*SA4(J,I)*RC(J,I, 156)
      A0158=SAI(J,I)*SA5(J,I)*RC(J,I, 158)
      A0159=SAJ(J,I)*SAD(J,I)*RC(J,I, 159)
      A0161=SAJ(J,I)*RC(J,I, 161)
      A0164=SAI(J,I)*SA7(J,I)*RC(J,I, 164)
      A0165=SAI(J,I)*SA2(J,I)*RC(J,I, 165)
      A0168=SAG(J,I)*SA4(J,I)*RC(J,I, 168)
      A0170=SAD(J,I)*SAB(J,I)*RC(J,I, 170)
      A0172=SAG(J,I)*SA6(J,I)*RC(J,I, 172)
      A0173=SAI(J,I)*RC(J,I, 173)
      A0175=SAG(J,I)*SA2(J,I)*RC(J,I, 175)
      A0178=SB2(J,I)*SA5(J,I)*RC(J,I, 178)
      A0180=SAD(J,I)*SAA(J,I)*RC(J,I, 180)
      A0181=SAG(J,I)*SA3(J,I)*RC(J,I, 181)
      A0184=SB2(J,I)*SAI(J,I)*RC(J,I, 184)
      A0186=SB2(J,I)*SA6(J,I)*RC(J,I, 186)
      A0187=SAG(J,I)*SA6(J,I)*RC(J,I, 187)
      A0190=SAD(J,I)*SA8(J,I)*SA7(J,I)*RC(J,I, 190)
      A0191=SAG(J,I)*RC(J,I, 191)
      A0193=SAG(J,I)*RC(J,I, 193)
      A0195=SB2(J,I)*SA2(J,I)*RC(J,I, 195)
      A0197=SB2(J,I)*RC(J,I, 197)
      A0200=SAB(J,I)*SAA(J,I)*RC(J,I, 200)
      A0202=SB3(J,I)*SA3(J,I)*RC(J,I, 202)
      A0204=SB1(J,I)*SA6(J,I)*RC(J,I, 204)
      A0205=SB1(J,I)*SA3(J,I)*RC(J,I, 205)
      A0208=SAE(J,I)*SA8(J,I)*RC(J,I, 208)
      A0210=SAB(J,I)*SA8(J,I)*RC(J,I, 210)
      A0211=SB1(J,I)*SA4(J,I)*RC(J,I, 211)
      A0214=SB7(J,I)*SA5(J,I)*RC(J,I, 214)
      A0216=SAD(J,I)*SA8(J,I)*RC(J,I, 216)
      A0218=SAE(J,I)*SA9(J,I)*RC(J,I, 218)
      A0220=SB6(J,I)*SA4(J,I)*RC(J,I, 220)
      A0221=SB5(J,I)*SAD(J,I)*RC(J,I, 221)
      A0224=SAF(J,I)*SA8(J,I)*RC(J,I, 224)
      A0225=SB6(J,I)*SA4(J,I)*RC(J,I, 225)
      A0227=SB6(J,I)*SA3(J,I)*RC(J,I, 227)
      A0230=SA8(J,I)*SA8(J,I)*SA4(J,I)*RC(J,I, 230)
      A0232=SA9(J,I)*SA8(J,I)*SA2(J,I)*RC(J,I, 232)
      A0233=SB7(J,I)*SA4(J,I)*RC(J,I, 233)
      A0236=SA8(J,I)*SB0(J,I)*RC(J,I, 236)
      A0238=SB6(J,I)*SA3(J,I)*RC(J,I, 238)
      A0240=SB0(J,I)*SA9(J,I)*RC(J,I, 240)
      A0242=SAA(J,I)*SA8(J,I)*RC(J,I, 242)
      A0243=SB8(J,I)*SA2(J,I)*RC(J,I, 243)
      A0246=SAD(J,I)*SA4(J,I)*RC(J,I, 246)
      A0248=SAB(J,I)*SA5(J,I)*RC(J,I, 248)
      A0250=SAB(J,I)*SA6(J,I)*RC(J,I, 250)
      A0251=SB8(J,I)*RC(J,I, 251)
      A0254=SB8(J,I)*RC(J,I, 254)
      A0256=SB8(J,I)*SA8(J,I)*RC(J,I, 256)
      A0258=SB8(J,I)*SA5(J,I)*RC(J,I, 258)
      A0260=SAH(J,I)*SA5(J,I)*RC(J,I, 260)
      A0261=SB9(J,I)*SA2(J,I)*RC(J,I, 261)
      A0263=SB9(J,I)*SA2(J,I)*RC(J,I, 263)
      A0266=SB8(J,I)*SA4(J,I)*RC(J,I, 266)
      A0268=SB8(J,I)*SA7(J,I)*RC(J,I, 268)
      A0270=SB8(J,I)*SA6(J,I)*RC(J,I, 270)
      A0271=SB9(J,I)*RC(J,I, 271)
      A0273=SB3(J,I)*RC(J,I, 273)
      A0276=SAD(J,I)*SAA(J,I)*RC(J,I, 276)
      A0277=SB3(J,I)*SA2(J,I)*RC(J,I, 277)
      A0280=SAB(J,I)*SAA(J,I)*RC(J,I, 280)
      A0282=SB5(J,I)*SA5(J,I)*RC(J,I, 282)
      A0284=SAB(J,I)*SA8(J,I)*SA4(J,I)*RC(J,I, 284)
      A0285=SB3(J,I)*SAD(J,I)*RC(J,I, 285)
      A0287=SB3(J,I)*SA6(J,I)*RC(J,I, 287)
      A0289=SB3(J,I)*SA6(J,I)*RC(J,I, 289)
      A0291=SB3(J,I)*RC(J,I, 291)
      A0293=SBB(J,I)*RC(J,I, 293)
      A0295=SBE(J,I)*RC(J,I, 295)
      A0298=SBE(J,I)*SA5(J,I)*RC(J,I, 298)
      A0300=SB3(J,I)*SA5(J,I)*RC(J,I, 300)
      A0302=SBE(J,I)*SA4(J,I)*RC(J,I, 302)
      A0304=SB3(J,I)*SA4(J,I)*RC(J,I, 304)
      A0305=SBB(J,I)*SA2(J,I)*RC(J,I, 305)
      A0307=SBB(J,I)*SA2(J,I)*RC(J,I, 307)
      A0309=SBB(J,I)*SAD(J,I)*RC(J,I, 309)
      A0311=SBB(J,I)*SAD(J,I)*RC(J,I, 311)
      A0314=SBE(J,I)*SA7(J,I)*RC(J,I, 314)
      A0316=SB3(J,I)*SA7(J,I)*RC(J,I, 316)
      A0318=SBE(J,I)*SA6(J,I)*RC(J,I, 318)
      A0319=SBA(J,I)*RC(J,I, 319)
      A0321=SBA(J,I)*RC(J,I, 321)
      A0324=SBD(J,I)*SA5(J,I)*RC(J,I, 324)
      A0326=SBC(J,I)*SA5(J,I)*RC(J,I, 326)
      A0328=SBF(J,I)*SA5(J,I)*RC(J,I, 328)
      A0329=SBA(J,I)*SA2(J,I)*RC(J,I, 329)
      A0331=SBA(J,I)*SA2(J,I)*RC(J,I, 331)
      A0333=SBA(J,I)*SA2(J,I)*RC(J,I, 333)
      A0336=SBD(J,I)*SA4(J,I)*RC(J,I, 336)
      A0338=SBC(J,I)*SA4(J,I)*RC(J,I, 338)
      A0340=SBF(J,I)*SA4(J,I)*RC(J,I, 340)
      A0341=SBA(J,I)*SAD(J,I)*RC(J,I, 341)
      A0343=SBA(J,I)*SAD(J,I)*RC(J,I, 343)
      A0345=SBA(J,I)*SAD(J,I)*RC(J,I, 345)
      A0348=SBC(J,I)*SA7(J,I)*RC(J,I, 348)
      A0350=SBD(J,I)*SA7(J,I)*RC(J,I, 350)
      A0352=SBF(J,I)*SA7(J,I)*RC(J,I, 352)
      A0354=SBD(J,I)*RC(J,I, 354)
      A0355=SAI(J,I)*SA6(J,I)*RC(J,I, 355)
      A0357=SBF(J,I)*RC(J,I, 357)
      A0359=SBF(J,I)*RC(J,I, 359)
      A0362=SBB(J,I)*SA6(J,I)*RC(J,I, 362)
      A0364=SAI(J,I)*SA9(J,I)*RC(J,I, 364)
      A0366=SAD(J,I)*SB8(J,I)*RC(J,I, 366)
      A0368=SAG(J,I)*SA5(J,I)*RC(J,I, 368)
      A0370=SBB(J,I)*SA5(J,I)*RC(J,I, 370)
      A0372=SBB(J,I)*SA6(J,I)*RC(J,I, 372)
      A0374=SBB(J,I)*SA4(J,I)*RC(J,I, 374)
      A0376=SAG(J,I)*SA5(J,I)*RC(J,I, 376)
      A0378=SAD(J,I)*SB8(J,I)*RC(J,I, 378)
      A0379=SBC(J,I)*SA6(J,I)*RC(J,I, 379)
      A0382=SBB(J,I)*SA5(J,I)*RC(J,I, 382)
      A0383=SBC(J,I)*RC(J,I, 383)
      A0386=SAG(J,I)*SA8(J,I)*RC(J,I, 386)
      A0387=SAD(J,I)*SB1(J,I)*RC(J,I, 387)
      A0390=SB6(J,I)*SAD(J,I)*RC(J,I, 390)
      A0392=SBG(J,I)*RC(J,I, 392)
      A0393=SBH(J,I)*SA6(J,I)*RC(J,I, 393)
      A0396=SBH(J,I)*SA5(J,I)*RC(J,I, 396)
      A0398=SB5(J,I)*SAA(J,I)*RC(J,I, 398)
      A0400=SBI(J,I)*RC(J,I, 400)
      A0401=SBI(J,I)*SA2(J,I)*RC(J,I, 401)
      A0404=SBG(J,I)*SA6(J,I)*RC(J,I, 404)
      A0405=SBI(J,I)*SAD(J,I)*RC(J,I, 405)
      A0408=SBI(J,I)*RC(J,I, 408)
      A0410=SBG(J,I)*SA5(J,I)*RC(J,I, 410)
      A0411=SBH(J,I)*SAA(J,I)*RC(J,I, 411)
      A0413=SBH(J,I)*SA6(J,I)*RC(J,I, 413)
      A0416=SAD(J,I)*SAA(J,I)*SA8(J,I)*RC(J,I, 416)
      A0418=SAI(J,I)*SAA(J,I)*RC(J,I, 418)
      A0420=SBI(J,I)*SA5(J,I)*RC(J,I, 420)
      A0421=SBJ(J,I)*SA3(J,I)*RC(J,I, 421)
      A0423=SBJ(J,I)*SA2(J,I)*RC(J,I, 423)
      A0426=SBJ(J,I)*RC(J,I, 426)
      A0427=SBI(J,I)*SA6(J,I)*RC(J,I, 427)
      A0429=SBI(J,I)*SA6(J,I)*RC(J,I, 429)
      A0432=SBJ(J,I)*RC(J,I, 432)
      A0434=SAG(J,I)*SAD(J,I)*RC(J,I, 434)
      A0435=SAD(J,I)*SB2(J,I)*RC(J,I, 435)
      A0437=SC0(J,I)*RC(J,I, 437)
      A0440=SC1(J,I)*SA6(J,I)*RC(J,I, 440)
      A0442=SC2(J,I)*SA6(J,I)*RC(J,I, 442)
      A0443=SC0(J,I)*SA2(J,I)*RC(J,I, 443)
      A0445=SC0(J,I)*SA2(J,I)*RC(J,I, 445)
      A0448=SC1(J,I)*SA4(J,I)*RC(J,I, 448)
      A0450=SC2(J,I)*SA4(J,I)*RC(J,I, 450)
      A0452=SC2(J,I)*SA5(J,I)*RC(J,I, 452)
      A0454=SC1(J,I)*SA5(J,I)*RC(J,I, 454)
      A0456=SC1(J,I)*SA7(J,I)*RC(J,I, 456)
      A0458=SC2(J,I)*SA7(J,I)*RC(J,I, 458)
      A0460=SC2(J,I)*RC(J,I, 460)
      A0462=SC1(J,I)*RC(J,I, 462)
      A0464=SBJ(J,I)*SA6(J,I)*RC(J,I, 464)
      A0465=SC2(J,I)*RC(J,I, 465)
      A0468=SC2(J,I)*RC(J,I, 468)
      A0470=SBJ(J,I)*SA6(J,I)*RC(J,I, 470)
      A0471=SC3(J,I)*RC(J,I, 471)
      A0473=SC3(J,I)*RC(J,I, 473)
      A0476=SC3(J,I)*RC(J,I, 476)
      A0478=SB2(J,I)*SAG(J,I)*RC(J,I, 478)
      A0479=SC3(J,I)*SA2(J,I)*RC(J,I, 479)
      A0481=SC3(J,I)*SA4(J,I)*RC(J,I, 481)
      A0483=SC3(J,I)*SAD(J,I)*RC(J,I, 483)
      A0486=SC3(J,I)*RC(J,I, 486)
      A0487=SC5(J,I)*RC(J,I, 487)
      A0489=SC5(J,I)*RC(J,I, 489)
      A0491=SC5(J,I)*RC(J,I, 491)
      A0494=SB1(J,I)*SBI(J,I)*SA6(J,I)*RC(J,I, 494)
      A0496=SB2(J,I)*SBG(J,I)*SA6(J,I)*RC(J,I, 496)
      A0498=SB1(J,I)*SBI(J,I)*SA7(J,I)*RC(J,I, 498)
      A0500=SB2(J,I)*SBG(J,I)*SA7(J,I)*RC(J,I, 500)
      A0501=SC6(J,I)*RC(J,I, 501)
      A0502=0.0
      A0503=SC6(J,I)*SA2(J,I)*RC(J,I, 503)
      A0504=0.0
      A0505=SC6(J,I)*SA2(J,I)*RC(J,I, 505)
      A0506=0.0
      A0507=SC6(J,I)*SA2(J,I)*RC(J,I, 507)
      A0508=0.0
      A0509=SC6(J,I)*SA2(J,I)*RC(J,I, 509)
      A0512=SA5(J,I)*SC8(J,I)*SAI(J,I)*RC(J,I, 512)
      A0513=SC6(J,I)*SA4(J,I)*RC(J,I, 513)
      A0514=0.0
      A0515=SC6(J,I)*SA4(J,I)*RC(J,I, 515)
      A0516=0.0
      A0517=SC6(J,I)*SA4(J,I)*RC(J,I, 517)
      A0518=0.0
      A0519=SC8(J,I)*RC(J,I, 519)
      A0521=SC8(J,I)*RC(J,I, 521)
      A0524=SA5(J,I)*SBJ(J,I)*SB2(J,I)*RC(J,I, 524)
      A0525=SC8(J,I)*SA2(J,I)*RC(J,I, 525)
      A0528=SAG(J,I)*SAG(J,I)*SAD(J,I)*RC(J,I, 528)
      A0530=SBJ(J,I)*SAI(J,I)*RC(J,I, 530)
      A0531=SC8(J,I)*SA2(J,I)*RC(J,I, 531)
      A0533=SC8(J,I)*SA2(J,I)*RC(J,I, 533)
      A0535=SC7(J,I)*RC(J,I, 535)
      A0538=SAG(J,I)*SAI(J,I)*RC(J,I, 538)
      A0540=SBJ(J,I)*SAD(J,I)*RC(J,I, 540)
      A0541=SC7(J,I)*SA2(J,I)*RC(J,I, 541)
      A0543=SC7(J,I)*SA4(J,I)*RC(J,I, 543)
      A0544=SA5(J,I)*SC3(J,I)*SA2(J,I)*RC(J,I, 544)
C----------------    SOLVING   1 ST SPECIES EQUATION    ----------------
      A0003=SA3(J,I)*RC(J,I,   3)
      A0005=SA4(J,I)*RC(J,I,   5)
      A0010=RC(J,I,  10)
      A0024=SA1(J,I)*RC(J,I,  24)
      A0036=SA6(J,I)*RC(J,I,  36)
      A0052=SA8(J,I)*RC(J,I,  52)
      A0066=SAA(J,I)*RC(J,I,  66)
      A0076=SAD(J,I)*RC(J,I,  76)
      A0086=SAE(J,I)*RC(J,I,  86)
      A0088=SAF(J,I)*RC(J,I,  88)
      A0102=SAG(J,I)*RC(J,I, 102)
      A0118=SB0(J,I)*RC(J,I, 118)
      A0126=SA8(J,I)*RC(J,I, 126)
      A0128=SA9(J,I)*RC(J,I, 128)
      A0142=SAB(J,I)*RC(J,I, 142)
      A0154=SAI(J,I)*RC(J,I, 154)
      A0166=SAG(J,I)*RC(J,I, 166)
      A0176=SB2(J,I)*RC(J,I, 176)
      A0194=SB1(J,I)*RC(J,I, 194)
      A0196=SB1(J,I)*RC(J,I, 196)
      A0244=SAB(J,I)*RC(J,I, 244)
      A0262=SB8(J,I)*RC(J,I, 262)
      A0264=SAH(J,I)*RC(J,I, 264)
      A0278=SB5(J,I)*RC(J,I, 278)
      A0306=SBE(J,I)*RC(J,I, 306)
      A0308=SB3(J,I)*RC(J,I, 308)
      A0330=SBD(J,I)*RC(J,I, 330)
      A0332=SBC(J,I)*RC(J,I, 332)
      A0334=SBF(J,I)*RC(J,I, 334)
      A0402=SBG(J,I)*RC(J,I, 402)
      A0424=SBI(J,I)*RC(J,I, 424)
      A0444=SC1(J,I)*RC(J,I, 444)
      A0446=SC2(J,I)*RC(J,I, 446)
      A0480=SB1(J,I)*SB2(J,I)*RC(J,I, 480)
      A0510=SC8(J,I)*SAI(J,I)*RC(J,I, 510)
      A0526=SBJ(J,I)*SB2(J,I)*RC(J,I, 526)
      A0532=SAG(J,I)*SBI(J,I)*RC(J,I, 532)
      A0534=SC3(J,I)*SAD(J,I)*RC(J,I, 534)
      A0542=SB2(J,I)*SAG(J,I)*RC(J,I, 542)
      RHS1A=+A0004+A0006+A0009+A0023+A0035+A0051+A0065+A0075+A0085+A0087
     1      +A0101+A0117+A0125+A0127+A0141+A0153+A0165+A0175+A0193+A0195
     2      +A0243+A0261+A0263+A0277+A0305+A0307+A0329+A0331+A0333+A0401
     3      +A0423+A0443+A0445+A0479+A0509+A0525+A0531+A0533+A0541
     4      +SOUR(J,I,1)*SB1(J,I)
      RHS1B=+1.000*A0503+1.000*A0505+1.000*A0507
      RHS1=RHO(J,I)*FA0(J,I)+DT*WM(  1)*(RHS1A+RHS1B)
      RHS2=-DT*WM(  1)*(
     1     +A0003+A0005+A0010+A0024+A0036+A0052+A0066+A0076+A0086+A0088
     2     +A0102+A0118+A0126+A0128+A0142+A0154+A0166+A0176+A0194+A0196
     3     +A0244+A0262+A0264+A0278+A0306+A0308+A0330+A0332+A0334+A0402
     4     +A0424+A0444+A0446+A0480+A0510+A0526+A0532+A0534+A0542
     *      )
      PPP=RHO(J,I)-(EA0(J,I)+WA0(J,I)+PA0(J,I)+QA0(J,I))
      FNEW=SOR2*SA0(J,I)+(RHS1-PA0(J,I)*SA0(J+1,I)-QA0(J,I)*SA0(J-1,I)
     1   -EA0(J,I)*SA0(J,I+1)-WA0(J,I)*SA0(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SA0(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SA0(J,I)=FNEW
      A0003=FNEW*A0003
      A0005=FNEW*A0005
      A0010=FNEW*A0010
      A0036=FNEW*A0036
      A0052=FNEW*A0052
      A0066=FNEW*A0066
      A0076=FNEW*A0076
      A0086=FNEW*A0086
      A0088=FNEW*A0088
      A0102=FNEW*A0102
      A0118=FNEW*A0118
      A0126=FNEW*A0126
      A0128=FNEW*A0128
      A0142=FNEW*A0142
      A0154=FNEW*A0154
      A0166=FNEW*A0166
      A0176=FNEW*A0176
      A0194=FNEW*A0194
      A0196=FNEW*A0196
      A0244=FNEW*A0244
      A0262=FNEW*A0262
      A0264=FNEW*A0264
      A0278=FNEW*A0278
      A0306=FNEW*A0306
      A0308=FNEW*A0308
      A0330=FNEW*A0330
      A0332=FNEW*A0332
      A0334=FNEW*A0334
      A0402=FNEW*A0402
      A0424=FNEW*A0424
      A0444=FNEW*A0444
      A0446=FNEW*A0446
      A0480=FNEW*A0480
      A0510=FNEW*A0510
      A0526=FNEW*A0526
      A0532=FNEW*A0532
      A0534=FNEW*A0534
      A0542=FNEW*A0542
C----------------    SOLVING   2 ND SPECIES EQUATION    ----------------
      A0001=SA2(J,I)*RC(J,I,   1)
      A0014=RC(J,I,  14)
      A0019=SA2(J,I)*RC(J,I,  19)
      A0024=SA0(J,I)*RC(J,I,  24)
      A0028=SA4(J,I)*RC(J,I,  28)
      A0030=SA5(J,I)*RC(J,I,  30)
      A0034=SA7(J,I)*RC(J,I,  34)
      A0047=SA8(J,I)*RC(J,I,  47)
      A0059=SAA(J,I)*RC(J,I,  59)
      A0071=SAB(J,I)*RC(J,I,  71)
      A0081=SAC(J,I)*RC(J,I,  81)
      A0097=SAD(J,I)*RC(J,I,  97)
      A0099=SAD(J,I)*RC(J,I,  99)
      A0111=SAF(J,I)*RC(J,I, 111)
      A0127=SAE(J,I)*RC(J,I, 127)
      A0129=SAE(J,I)*RC(J,I, 129)
      A0135=SB0(J,I)*RC(J,I, 135)
      A0149=SAH(J,I)*RC(J,I, 149)
      A0171=SAI(J,I)*RC(J,I, 171)
      A0185=SAG(J,I)*RC(J,I, 185)
      A0199=SB2(J,I)*RC(J,I, 199)
      A0201=SB2(J,I)*RC(J,I, 201)
      A0203=SB2(J,I)*RC(J,I, 203)
      A0209=SB1(J,I)*RC(J,I, 209)
      A0229=SB6(J,I)*RC(J,I, 229)
      A0231=SB6(J,I)*RC(J,I, 231)
      A0237=SB7(J,I)*RC(J,I, 237)
      A0239=SB7(J,I)*RC(J,I, 239)
      A0241=SB7(J,I)*RC(J,I, 241)
      A0249=SB8(J,I)*RC(J,I, 249)
      A0269=SB9(J,I)*RC(J,I, 269)
      A0283=SB3(J,I)*RC(J,I, 283)
      A0290=SBB(J,I)*RC(J,I, 290)
      A0317=SBB(J,I)*RC(J,I, 317)
      A0361=SBF(J,I)*RC(J,I, 361)
      A0371=SBC(J,I)*RC(J,I, 371)
      A0394=SBG(J,I)*RC(J,I, 394)
      A0397=SBH(J,I)*RC(J,I, 397)
      A0403=SBI(J,I)*RC(J,I, 403)
      A0415=SBG(J,I)*RC(J,I, 415)
      A0428=SBJ(J,I)*RC(J,I, 428)
      A0439=SC0(J,I)*RC(J,I, 439)
      A0441=SC0(J,I)*RC(J,I, 441)
      A0463=SC1(J,I)*RC(J,I, 463)
      A0469=SC2(J,I)*RC(J,I, 469)
      A0493=SC5(J,I)*RC(J,I, 493)
      A0495=SC5(J,I)*RC(J,I, 495)
      RHS1A=+A0002+A0013+A0020+A0023+A0027+A0029+A0033+A0048+A0060+A0072
     1      +A0082+A0098+A0100+A0112+A0128+A0130+A0136+A0150+A0172+A0186
     2      +A0200+A0202+A0204+A0210+A0230+A0232+A0238+A0240+A0242+A0250
     3      +A0270+A0284+A0289+A0318+A0362+A0372+A0393+A0398+A0404+A0416
     4      +A0427+A0440+A0442+A0464+A0470+A0494+A0496
      RHS1B=0.0
      RHS1=RHO(J,I)*FA1(J,I)+DT*WM(  2)*(RHS1A+RHS1B)
      RHS2=-DT*WM(  2)*(
     1     +A0001+A0014+A0019+A0024+A0028+A0030+A0034+A0047+A0059+A0071
     2     +A0081+A0097+A0099+A0111+A0127+A0129+A0135+A0149+A0171+A0185
     3     +A0199+A0201+A0203+A0209+A0229+A0231+A0237+A0239+A0241+A0249
     4     +A0269+A0283+A0290+A0317+A0361+A0371+A0394+A0397+A0403+A0415
     5     +A0428+A0439+A0441+A0463+A0469+A0493+A0495
     *     +0.5*SOUR(J,I,3))
      PPP=RHO(J,I)-(EA1(J,I)+WA1(J,I)+PA1(J,I)+QA1(J,I))
      FNEW=SOR2*SA1(J,I)+(RHS1-PA1(J,I)*SA1(J+1,I)-QA1(J,I)*SA1(J-1,I)
     1   -EA1(J,I)*SA1(J,I+1)-WA1(J,I)*SA1(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SA1(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SA1(J,I)=FNEW
      A0014=FNEW*A0014
      A0024=FNEW*A0024
      A0028=FNEW*A0028
      A0030=FNEW*A0030
      A0034=FNEW*A0034
      A0047=FNEW*A0047
      A0059=FNEW*A0059
      A0071=FNEW*A0071
      A0081=FNEW*A0081
      A0097=FNEW*A0097
      A0099=FNEW*A0099
      A0111=FNEW*A0111
      A0127=FNEW*A0127
      A0129=FNEW*A0129
      A0135=FNEW*A0135
      A0149=FNEW*A0149
      A0171=FNEW*A0171
      A0185=FNEW*A0185
      A0199=FNEW*A0199
      A0201=FNEW*A0201
      A0203=FNEW*A0203
      A0209=FNEW*A0209
      A0229=FNEW*A0229
      A0231=FNEW*A0231
      A0237=FNEW*A0237
      A0239=FNEW*A0239
      A0241=FNEW*A0241
      A0249=FNEW*A0249
      A0269=FNEW*A0269
      A0283=FNEW*A0283
      A0290=FNEW*A0290
      A0317=FNEW*A0317
      A0361=FNEW*A0361
      A0371=FNEW*A0371
      A0394=FNEW*A0394
      A0397=FNEW*A0397
      A0403=FNEW*A0403
      A0415=FNEW*A0415
      A0428=FNEW*A0428
      A0439=FNEW*A0439
      A0441=FNEW*A0441
      A0463=FNEW*A0463
      A0469=FNEW*A0469
      A0493=FNEW*A0493
      A0495=FNEW*A0495
C----------------    SOLVING   3 RD SPECIES EQUATION    ----------------
      A0001=SA1(J,I)*RC(J,I,   1)
      A0004=SA4(J,I)*RC(J,I,   4)
      A0006=SA5(J,I)*RC(J,I,   6)
      A0009=SA2(J,I)*RC(J,I,   9)
      A0011=SA4(J,I)*RC(J,I,  11)
      A0015=SA3(J,I)*RC(J,I,  15)
      A0019=SA1(J,I)*RC(J,I,  19)
      A0021=SA6(J,I)*RC(J,I,  21)
      A0023=SA6(J,I)*RC(J,I,  23)
      A0025=SA6(J,I)*RC(J,I,  25)
      A0035=SA7(J,I)*RC(J,I,  35)
      A0037=SA7(J,I)*RC(J,I,  37)
      A0044=SA9(J,I)*RC(J,I,  44)
      A0050=SA8(J,I)*RC(J,I,  50)
      A0051=SAA(J,I)*RC(J,I,  51)
      A0056=SA9(J,I)*RC(J,I,  56)
      A0063=SAA(J,I)*RC(J,I,  63)
      A0065=SAB(J,I)*RC(J,I,  65)
      A0075=SAC(J,I)*RC(J,I,  75)
      A0085=SAD(J,I)*RC(J,I,  85)
      A0087=SAD(J,I)*RC(J,I,  87)
      A0092=SAB(J,I)*RC(J,I,  92)
      A0094=SAG(J,I)*RC(J,I,  94)
      A0104=SAI(J,I)*RC(J,I, 104)
      A0105=SAD(J,I)*RC(J,I, 105)
      A0110=SAB(J,I)*RC(J,I, 110)
      A0112=SA8(J,I)*SA4(J,I)*RC(J,I, 112)
      A0117=SAE(J,I)*RC(J,I, 117)
      A0120=SAB(J,I)*RC(J,I, 120)
      A0124=SA8(J,I)*SA2(J,I)*RC(J,I, 124)
      A0130=SA8(J,I)*SA4(J,I)*RC(J,I, 130)
      A0132=SB1(J,I)*SA2(J,I)*RC(J,I, 132)
      A0134=SA8(J,I)*RC(J,I, 134)
      A0138=SAB(J,I)*RC(J,I, 138)
      A0141=SAH(J,I)*RC(J,I, 141)
      A0143=SAH(J,I)*RC(J,I, 143)
      A0152=SAB(J,I)*RC(J,I, 152)
      A0153=SAJ(J,I)*RC(J,I, 153)
      A0162=SAI(J,I)*RC(J,I, 162)
      A0165=SAI(J,I)*RC(J,I, 165)
      A0174=SAG(J,I)*RC(J,I, 174)
      A0175=SAG(J,I)*RC(J,I, 175)
      A0182=SB3(J,I)*RC(J,I, 182)
      A0192=SB2(J,I)*RC(J,I, 192)
      A0195=SB2(J,I)*RC(J,I, 195)
      A0198=SB1(J,I)*RC(J,I, 198)
      A0206=SB6(J,I)*RC(J,I, 206)
      A0212=SB5(J,I)*RC(J,I, 212)
      A0215=SB5(J,I)*RC(J,I, 215)
      A0223=SB6(J,I)*RC(J,I, 223)
      A0226=SAA(J,I)*SA8(J,I)*RC(J,I, 226)
      A0228=SA8(J,I)*SA8(J,I)*RC(J,I, 228)
      A0232=SA9(J,I)*SA8(J,I)*RC(J,I, 232)
      A0234=SB6(J,I)*RC(J,I, 234)
      A0243=SB8(J,I)*RC(J,I, 243)
      A0245=SB8(J,I)*RC(J,I, 245)
      A0252=SAB(J,I)*RC(J,I, 252)
      A0261=SB9(J,I)*RC(J,I, 261)
      A0263=SB9(J,I)*RC(J,I, 263)
      A0274=SB5(J,I)*RC(J,I, 274)
      A0275=SB3(J,I)*RC(J,I, 275)
      A0277=SB3(J,I)*RC(J,I, 277)
      A0286=SAI(J,I)*SA8(J,I)*RC(J,I, 286)
      A0305=SBB(J,I)*RC(J,I, 305)
      A0307=SBB(J,I)*RC(J,I, 307)
      A0329=SBA(J,I)*RC(J,I, 329)
      A0331=SBA(J,I)*RC(J,I, 331)
      A0333=SBA(J,I)*RC(J,I, 333)
      A0358=SBB(J,I)*RC(J,I, 358)
      A0365=SBF(J,I)*RC(J,I, 365)
      A0367=SBF(J,I)*RC(J,I, 367)
      A0375=SBC(J,I)*RC(J,I, 375)
      A0377=SBC(J,I)*RC(J,I, 377)
      A0384=SBB(J,I)*RC(J,I, 384)
      A0388=SBG(J,I)*RC(J,I, 388)
      A0391=SBH(J,I)*RC(J,I, 391)
      A0399=SBG(J,I)*RC(J,I, 399)
      A0401=SBI(J,I)*RC(J,I, 401)
      A0422=SB5(J,I)*SAD(J,I)*RC(J,I, 422)
      A0423=SBJ(J,I)*RC(J,I, 423)
      A0425=SBI(J,I)*RC(J,I, 425)
      A0433=SBJ(J,I)*RC(J,I, 433)
      A0436=SBI(J,I)*RC(J,I, 436)
      A0443=SC0(J,I)*RC(J,I, 443)
      A0445=SC0(J,I)*RC(J,I, 445)
      A0461=SBJ(J,I)*RC(J,I, 461)
      A0467=SBJ(J,I)*RC(J,I, 467)
      A0472=SB1(J,I)*SB2(J,I)*RC(J,I, 472)
      A0477=SC3(J,I)*RC(J,I, 477)
      A0479=SC3(J,I)*RC(J,I, 479)
      A0482=SC4(J,I)*SBI(J,I)*RC(J,I, 482)
      A0503=SC6(J,I)*RC(J,I, 503)
      A0505=SC6(J,I)*RC(J,I, 505)
      A0507=SC6(J,I)*RC(J,I, 507)
      A0509=SC6(J,I)*RC(J,I, 509)
      A0525=SC8(J,I)*RC(J,I, 525)
      A0527=SC8(J,I)*RC(J,I, 527)
      A0529=SC8(J,I)*RC(J,I, 529)
      A0531=SC8(J,I)*RC(J,I, 531)
      A0533=SC8(J,I)*RC(J,I, 533)
      A0537=SC7(J,I)*RC(J,I, 537)
      A0539=SC7(J,I)*RC(J,I, 539)
      A0541=SC7(J,I)*RC(J,I, 541)
      A0544=SA5(J,I)*SC3(J,I)*RC(J,I, 544)
      RHS1A=+A0002+A0003+A0005+A0010+A0010+A0012+A0016+A0020+A0022+A0024
     1      +A0026+A0036+A0038+A0043+A0049+A0052+A0055+A0064+A0066+A0076
     2      +A0086+A0088+A0091+A0093+A0103+A0106+A0109+A0111+A0118+A0119
     3      +A0123+A0123+A0129+A0131+A0131+A0133+A0137+A0142+A0144+A0151
     4      +A0154+A0161+A0166+A0173+A0176+A0181+A0191+A0196+A0197+A0205
     5      +A0211+A0216+A0224+A0225+A0227+A0231+A0233+A0244+A0246+A0251
     6      +A0262+A0264+A0273+A0276+A0278+A0285+A0306+A0308+A0330+A0332
     7      +A0334+A0357+A0366+A0368+A0376+A0378+A0383+A0387+A0392+A0400
     8      +A0402+A0421+A0424+A0426+A0434+A0435+A0444+A0446+A0462+A0468
     9      +A0471+A0478+A0480+A0481+A0504+A0506+A0508+A0510+A0526+A0528
     *      +A0530+A0532+A0534+A0538+A0540+A0542+A0543
     *      +SOUR(J,I,4)*SA4(J,I)
      RHS1B=0.0
      RHS1=RHO(J,I)*FA2(J,I)+DT*WM(  3)*(RHS1A+RHS1B)
      RHS2=-DT*WM(  3)*(
     1     +A0001+A0004+A0006+A0009+A0009+A0011+A0015+A0019+A0021+A0023
     2     +A0025+A0035+A0037+A0044+A0050+A0051+A0056+A0063+A0065+A0075
     3     +A0085+A0087+A0092+A0094+A0104+A0105+A0110+A0112+A0117+A0120
     4     +A0124+A0124+A0130+A0132+A0132+A0134+A0138+A0141+A0143+A0152
     5     +A0153+A0162+A0165+A0174+A0175+A0182+A0192+A0195+A0198+A0206
     6     +A0212+A0215+A0223+A0226+A0228+A0232+A0234+A0243+A0245+A0252
     7     +A0261+A0263+A0274+A0275+A0277+A0286+A0305+A0307+A0329+A0331
     8     +A0333+A0358+A0365+A0367+A0375+A0377+A0384+A0388+A0391+A0399
     9     +A0401+A0422+A0423+A0425+A0433+A0436+A0443+A0445+A0461+A0467
     *     +A0472+A0477+A0479+A0482+A0503+A0505+A0507+A0509+A0525+A0527
     *     +A0529+A0531+A0533+A0537+A0539+A0541+A0544
     *      )
      PPP=RHO(J,I)-(EA2(J,I)+WA2(J,I)+PA2(J,I)+QA2(J,I))
      FNEW=SOR2*SA2(J,I)+(RHS1-PA2(J,I)*SA2(J+1,I)-QA2(J,I)*SA2(J-1,I)
     1   -EA2(J,I)*SA2(J,I+1)-WA2(J,I)*SA2(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SA2(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SA2(J,I)=FNEW
      A0001=FNEW*A0001
      A0004=FNEW*A0004
      A0006=FNEW*A0006
      A0019=FNEW*A0019
      A0021=FNEW*A0021
      A0025=FNEW*A0025
      A0035=FNEW*A0035
      A0037=FNEW*A0037
      A0044=FNEW*A0044
      A0051=FNEW*A0051
      A0056=FNEW*A0056
      A0065=FNEW*A0065
      A0092=FNEW*A0092
      A0094=FNEW*A0094
      A0104=FNEW*A0104
      A0105=FNEW*A0105
      A0110=FNEW*A0110
      A0120=FNEW*A0120
      A0124=FNEW*A0124
      A0132=FNEW*A0132
      A0134=FNEW*A0134
      A0138=FNEW*A0138
      A0141=FNEW*A0141
      A0143=FNEW*A0143
      A0153=FNEW*A0153
      A0165=FNEW*A0165
      A0182=FNEW*A0182
      A0192=FNEW*A0192
      A0195=FNEW*A0195
      A0206=FNEW*A0206
      A0212=FNEW*A0212
      A0215=FNEW*A0215
      A0223=FNEW*A0223
      A0226=FNEW*A0226
      A0228=FNEW*A0228
      A0234=FNEW*A0234
      A0243=FNEW*A0243
      A0245=FNEW*A0245
      A0261=FNEW*A0261
      A0263=FNEW*A0263
      A0274=FNEW*A0274
      A0275=FNEW*A0275
      A0307=FNEW*A0307
      A0365=FNEW*A0365
      A0367=FNEW*A0367
      A0375=FNEW*A0375
      A0377=FNEW*A0377
      A0388=FNEW*A0388
      A0391=FNEW*A0391
      A0401=FNEW*A0401
      A0422=FNEW*A0422
      A0423=FNEW*A0423
      A0433=FNEW*A0433
      A0436=FNEW*A0436
      A0477=FNEW*A0477
      A0479=FNEW*A0479
      A0482=FNEW*A0482
      A0503=FNEW*A0503
      A0505=FNEW*A0505
      A0507=FNEW*A0507
      A0509=FNEW*A0509
      A0525=FNEW*A0525
      A0527=FNEW*A0527
      A0529=FNEW*A0529
      A0531=FNEW*A0531
      A0533=FNEW*A0533
      A0537=FNEW*A0537
      A0539=FNEW*A0539
      A0541=FNEW*A0541
      A0544=FNEW*A0544
C----------------    SOLVING   4 TH SPECIES EQUATION    ----------------
      A0002=SA4(J,I)*RC(J,I,   2)
      A0003=SA0(J,I)*RC(J,I,   3)
      A0007=SA5(J,I)*RC(J,I,   7)
      A0013=SA3(J,I)*RC(J,I,  13)
      A0015=SA2(J,I)*RC(J,I,  15)
      A0017=SA4(J,I)*RC(J,I,  17)
      A0026=SA5(J,I)*RC(J,I,  26)
      A0027=SA6(J,I)*RC(J,I,  27)
      A0041=SA7(J,I)*RC(J,I,  41)
      A0048=SA9(J,I)*RC(J,I,  48)
      A0053=SAA(J,I)*RC(J,I,  53)
      A0055=SAA(J,I)*RC(J,I,  55)
      A0067=SAB(J,I)*RC(J,I,  67)
      A0079=SAC(J,I)*RC(J,I,  79)
      A0091=SAD(J,I)*RC(J,I,  91)
      A0100=SAH(J,I)*RC(J,I, 100)
      A0123=SAE(J,I)*RC(J,I, 123)
      A0125=SAE(J,I)*RC(J,I, 125)
      A0133=SB0(J,I)*RC(J,I, 133)
      A0136=SAA(J,I)*RC(J,I, 136)
      A0147=SAH(J,I)*RC(J,I, 147)
      A0155=SAJ(J,I)*RC(J,I, 155)
      A0167=SAI(J,I)*RC(J,I, 167)
      A0169=SAI(J,I)*RC(J,I, 169)
      A0179=SAG(J,I)*RC(J,I, 179)
      A0181=SAG(J,I)*RC(J,I, 181)
      A0202=SB3(J,I)*RC(J,I, 202)
      A0205=SB1(J,I)*RC(J,I, 205)
      A0207=SB1(J,I)*RC(J,I, 207)
      A0217=SB5(J,I)*RC(J,I, 217)
      A0219=SB5(J,I)*RC(J,I, 219)
      A0227=SB6(J,I)*RC(J,I, 227)
      A0235=SB7(J,I)*RC(J,I, 235)
      A0238=SB6(J,I)*RC(J,I, 238)
      A0265=SB9(J,I)*RC(J,I, 265)
      A0279=SB3(J,I)*RC(J,I, 279)
      A0301=SBB(J,I)*RC(J,I, 301)
      A0303=SBB(J,I)*RC(J,I, 303)
      A0335=SBA(J,I)*RC(J,I, 335)
      A0337=SBA(J,I)*RC(J,I, 337)
      A0339=SBA(J,I)*RC(J,I, 339)
      A0373=SBC(J,I)*RC(J,I, 373)
      A0385=SBG(J,I)*RC(J,I, 385)
      A0389=SBG(J,I)*RC(J,I, 389)
      A0417=SBJ(J,I)*RC(J,I, 417)
      A0421=SBJ(J,I)*RC(J,I, 421)
      A0447=SC0(J,I)*RC(J,I, 447)
      A0449=SC0(J,I)*RC(J,I, 449)
      RHS1A=+A0001+A0004+A0008+A0014+A0014+A0016+A0018+A0025+A0028+A0042
     1      +A0047+A0054+A0056+A0068+A0080+A0092+A0099+A0124+A0126+A0134
     2      +A0135+A0148+A0156+A0168+A0170+A0180+A0182+A0201+A0206+A0208
     3      +A0218+A0220+A0228+A0236+A0237+A0266+A0280+A0302+A0304+A0336
     4      +A0338+A0340+A0374+A0386+A0390+A0418+A0422+A0448+A0450
      RHS1B=0.0
      RHS1=RHO(J,I)*FA3(J,I)+DT*WM(  4)*(RHS1A+RHS1B)
      RHS2=-DT*WM(  4)*(
     1     +A0002+A0003+A0007+A0013+A0013+A0015+A0017+A0026+A0027+A0041
     2     +A0048+A0053+A0055+A0067+A0079+A0091+A0100+A0123+A0125+A0133
     3     +A0136+A0147+A0155+A0167+A0169+A0179+A0181+A0202+A0205+A0207
     4     +A0217+A0219+A0227+A0235+A0238+A0265+A0279+A0301+A0303+A0335
     5     +A0337+A0339+A0373+A0385+A0389+A0417+A0421+A0447+A0449
     *     +SOUR(J,I,5))
      PPP=RHO(J,I)-(EA3(J,I)+WA3(J,I)+PA3(J,I)+QA3(J,I))
      FNEW=SOR2*SA3(J,I)+(RHS1-PA3(J,I)*SA3(J+1,I)-QA3(J,I)*SA3(J-1,I)
     1   -EA3(J,I)*SA3(J,I+1)-WA3(J,I)*SA3(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SA3(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SA3(J,I)=FNEW
      A0003=FNEW*A0003
      A0007=FNEW*A0007
      A0015=FNEW*A0015
      A0027=FNEW*A0027
      A0041=FNEW*A0041
      A0048=FNEW*A0048
      A0053=FNEW*A0053
      A0055=FNEW*A0055
      A0067=FNEW*A0067
      A0079=FNEW*A0079
      A0091=FNEW*A0091
      A0100=FNEW*A0100
      A0123=FNEW*A0123
      A0125=FNEW*A0125
      A0133=FNEW*A0133
      A0147=FNEW*A0147
      A0155=FNEW*A0155
      A0167=FNEW*A0167
      A0169=FNEW*A0169
      A0179=FNEW*A0179
      A0202=FNEW*A0202
      A0207=FNEW*A0207
      A0217=FNEW*A0217
      A0219=FNEW*A0219
      A0227=FNEW*A0227
      A0235=FNEW*A0235
      A0265=FNEW*A0265
      A0279=FNEW*A0279
      A0301=FNEW*A0301
      A0303=FNEW*A0303
      A0335=FNEW*A0335
      A0337=FNEW*A0337
      A0339=FNEW*A0339
      A0373=FNEW*A0373
      A0385=FNEW*A0385
      A0389=FNEW*A0389
      A0417=FNEW*A0417
      A0421=FNEW*A0421
      A0447=FNEW*A0447
      A0449=FNEW*A0449
C----------------    SOLVING   5 TH SPECIES EQUATION    ----------------
      A0002=SA3(J,I)*RC(J,I,   2)
      A0004=SA2(J,I)*RC(J,I,   4)
      A0005=SA0(J,I)*RC(J,I,   5)
      A0008=SA4(J,I)*RC(J,I,   8)
      A0011=SA2(J,I)*RC(J,I,  11)
      A0016=RC(J,I,  16)
      A0017=SA3(J,I)*RC(J,I,  17)
      A0022=SA4(J,I)*RC(J,I,  22)
      A0028=SA1(J,I)*RC(J,I,  28)
      A0029=SA6(J,I)*RC(J,I,  29)
      A0031=SA4(J,I)*RC(J,I,  31)
      A0038=SA5(J,I)*RC(J,I,  38)
      A0039=SA7(J,I)*RC(J,I,  39)
      A0042=SA6(J,I)*RC(J,I,  42)
      A0043=SA8(J,I)*RC(J,I,  43)
      A0046=SA9(J,I)*RC(J,I,  46)
      A0054=SA8(J,I)*RC(J,I,  54)
      A0057=SAA(J,I)*RC(J,I,  57)
      A0068=SAA(J,I)*RC(J,I,  68)
      A0069=SAB(J,I)*RC(J,I,  69)
      A0077=SAC(J,I)*RC(J,I,  77)
      A0080=SAD(J,I)*RC(J,I,  80)
      A0089=SAD(J,I)*RC(J,I,  89)
      A0096=SAH(J,I)*RC(J,I,  96)
      A0098=SAB(J,I)*RC(J,I,  98)
      A0109=SAF(J,I)*RC(J,I, 109)
      A0112=SA8(J,I)*SA2(J,I)*RC(J,I, 112)
      A0119=SAE(J,I)*RC(J,I, 119)
      A0121=SAE(J,I)*RC(J,I, 121)
      A0130=SA8(J,I)*SA2(J,I)*RC(J,I, 130)
      A0145=SAH(J,I)*RC(J,I, 145)
      A0148=SAB(J,I)*RC(J,I, 148)
      A0156=SAI(J,I)*RC(J,I, 156)
      A0157=SAJ(J,I)*RC(J,I, 157)
      A0168=SAG(J,I)*RC(J,I, 168)
      A0177=SAG(J,I)*RC(J,I, 177)
      A0188=SB4(J,I)*RC(J,I, 188)
      A0211=SB1(J,I)*RC(J,I, 211)
      A0213=SB1(J,I)*RC(J,I, 213)
      A0220=SB6(J,I)*RC(J,I, 220)
      A0225=SB6(J,I)*RC(J,I, 225)
      A0230=SA8(J,I)*SA8(J,I)*RC(J,I, 230)
      A0233=SB7(J,I)*RC(J,I, 233)
      A0246=SAD(J,I)*RC(J,I, 246)
      A0247=SB8(J,I)*RC(J,I, 247)
      A0255=SB5(J,I)*RC(J,I, 255)
      A0257=SB9(J,I)*RC(J,I, 257)
      A0259=SB9(J,I)*RC(J,I, 259)
      A0266=SB8(J,I)*RC(J,I, 266)
      A0272=SAD(J,I)*RC(J,I, 272)
      A0281=SB3(J,I)*RC(J,I, 281)
      A0284=SAB(J,I)*SA8(J,I)*RC(J,I, 284)
      A0288=SAB(J,I)*SAA(J,I)*RC(J,I, 288)
      A0297=SBB(J,I)*RC(J,I, 297)
      A0299=SBB(J,I)*RC(J,I, 299)
      A0302=SBE(J,I)*RC(J,I, 302)
      A0304=SB3(J,I)*RC(J,I, 304)
      A0323=SBA(J,I)*RC(J,I, 323)
      A0325=SBA(J,I)*RC(J,I, 325)
      A0327=SBA(J,I)*RC(J,I, 327)
      A0336=SBD(J,I)*RC(J,I, 336)
      A0338=SBC(J,I)*RC(J,I, 338)
      A0340=SBF(J,I)*RC(J,I, 340)
      A0353=SAG(J,I)*RC(J,I, 353)
      A0356=SBF(J,I)*RC(J,I, 356)
      A0369=SBF(J,I)*RC(J,I, 369)
      A0374=SBB(J,I)*RC(J,I, 374)
      A0380=SBB(J,I)*SA4(J,I)*RC(J,I, 380)
      A0381=SBC(J,I)*RC(J,I, 381)
      A0395=SBG(J,I)*RC(J,I, 395)
      A0409=SBI(J,I)*RC(J,I, 409)
      A0414=SA8(J,I)*SB2(J,I)*RC(J,I, 414)
      A0419=SBJ(J,I)*RC(J,I, 419)
      A0430=SB2(J,I)*SAB(J,I)*RC(J,I, 430)
      A0448=SC1(J,I)*RC(J,I, 448)
      A0450=SC2(J,I)*RC(J,I, 450)
      A0451=SC0(J,I)*RC(J,I, 451)
      A0453=SC0(J,I)*RC(J,I, 453)
      A0481=SC3(J,I)*RC(J,I, 481)
      A0511=SC6(J,I)*RC(J,I, 511)
      A0513=SC6(J,I)*RC(J,I, 513)
      A0515=SC6(J,I)*RC(J,I, 515)
      A0517=SC6(J,I)*RC(J,I, 517)
      A0523=SC8(J,I)*RC(J,I, 523)
      A0543=SC7(J,I)*RC(J,I, 543)
      RHS1A=+A0001+A0003+A0006+A0007+A0007+A0012+A0015+A0018+A0021+A0021
     1      +A0027+A0030+A0032+A0032+A0037+A0040+A0041+A0044+A0045+A0053
     2      +A0058+A0067+A0070+A0078+A0079+A0090+A0095+A0097+A0110+A0111
     3      +A0120+A0122+A0129+A0146+A0147+A0155+A0158+A0167+A0178+A0187
     4      +A0212+A0214+A0219+A0226+A0229+A0234+A0245+A0248+A0256+A0258
     5      +A0260+A0265+A0271+A0282+A0283+A0287+A0298+A0300+A0301+A0303
     6      +A0324+A0326+A0328+A0335+A0337+A0339+A0354+A0355+A0370+A0373
     7      +A0379+A0379+A0382+A0396+A0410+A0413+A0420+A0429+A0447+A0449
     8      +A0452+A0454+A0482+A0512+A0514+A0516+A0518+A0524+A0544
      RHS1B=0.0
      RHS1=RHO(J,I)*FA4(J,I)+DT*WM(  5)*(RHS1A+RHS1B)
      RHS2=-DT*WM(  5)*(
     1     +A0002+A0004+A0005+A0008+A0008+A0011+A0016+A0017+A0022+A0022
     2     +A0028+A0029+A0031+A0031+A0038+A0039+A0042+A0043+A0046+A0054
     3     +A0057+A0068+A0069+A0077+A0080+A0089+A0096+A0098+A0109+A0112
     4     +A0119+A0121+A0130+A0145+A0148+A0156+A0157+A0168+A0177+A0188
     5     +A0211+A0213+A0220+A0225+A0230+A0233+A0246+A0247+A0255+A0257
     6     +A0259+A0266+A0272+A0281+A0284+A0288+A0297+A0299+A0302+A0304
     7     +A0323+A0325+A0327+A0336+A0338+A0340+A0353+A0356+A0369+A0374
     8     +A0380+A0380+A0381+A0395+A0409+A0414+A0419+A0430+A0448+A0450
     9     +A0451+A0453+A0481+A0511+A0513+A0515+A0517+A0523+A0543
     *     +SOUR(J,I,4)
     *      )
      PPP=RHO(J,I)-(EA4(J,I)+WA4(J,I)+PA4(J,I)+QA4(J,I))
      FNEW=SOR2*SA4(J,I)+(RHS1-PA4(J,I)*SA4(J+1,I)-QA4(J,I)*SA4(J-1,I)
     1   -EA4(J,I)*SA4(J,I+1)-WA4(J,I)*SA4(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SA4(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SA4(J,I)=FNEW
      A0005=FNEW*A0005
      A0008=FNEW*A0008
      A0011=FNEW*A0011
      A0017=FNEW*A0017
      A0022=FNEW*A0022
      A0028=FNEW*A0028
      A0029=FNEW*A0029
      A0031=FNEW*A0031
      A0039=FNEW*A0039
      A0046=FNEW*A0046
      A0057=FNEW*A0057
      A0069=FNEW*A0069
      A0077=FNEW*A0077
      A0080=FNEW*A0080
      A0089=FNEW*A0089
      A0096=FNEW*A0096
      A0109=FNEW*A0109
      A0119=FNEW*A0119
      A0121=FNEW*A0121
      A0145=FNEW*A0145
      A0157=FNEW*A0157
      A0177=FNEW*A0177
      A0188=FNEW*A0188
      A0213=FNEW*A0213
      A0220=FNEW*A0220
      A0225=FNEW*A0225
      A0233=FNEW*A0233
      A0247=FNEW*A0247
      A0255=FNEW*A0255
      A0257=FNEW*A0257
      A0259=FNEW*A0259
      A0281=FNEW*A0281
      A0288=FNEW*A0288
      A0297=FNEW*A0297
      A0299=FNEW*A0299
      A0302=FNEW*A0302
      A0323=FNEW*A0323
      A0325=FNEW*A0325
      A0327=FNEW*A0327
      A0336=FNEW*A0336
      A0338=FNEW*A0338
      A0340=FNEW*A0340
      A0356=FNEW*A0356
      A0369=FNEW*A0369
      A0380=FNEW*A0380
      A0381=FNEW*A0381
      A0395=FNEW*A0395
      A0409=FNEW*A0409
      A0414=FNEW*A0414
      A0419=FNEW*A0419
      A0430=FNEW*A0430
      A0448=FNEW*A0448
      A0450=FNEW*A0450
      A0451=FNEW*A0451
      A0453=FNEW*A0453
      A0481=FNEW*A0481
      A0511=FNEW*A0511
      A0513=FNEW*A0513
      A0515=FNEW*A0515
      A0517=FNEW*A0517
      A0523=FNEW*A0523
      A0543=FNEW*A0543
C----------------    SOLVING   6 TH SPECIES EQUATION    ----------------
      A0006=SA2(J,I)*RC(J,I,   6)
      A0007=SA3(J,I)*RC(J,I,   7)
      A0012=RC(J,I,  12)
      A0026=SA3(J,I)*RC(J,I,  26)
      A0030=SA1(J,I)*RC(J,I,  30)
      A0038=SA4(J,I)*RC(J,I,  38)
      A0040=SA6(J,I)*RC(J,I,  40)
      A0058=SA8(J,I)*RC(J,I,  58)
      A0070=SAA(J,I)*RC(J,I,  70)
      A0078=SAD(J,I)*RC(J,I,  78)
      A0090=SAF(J,I)*RC(J,I,  90)
      A0122=SB0(J,I)*RC(J,I, 122)
      A0137=SB0(J,I)*RC(J,I, 137)
      A0144=SAF(J,I)*RC(J,I, 144)
      A0146=SAB(J,I)*RC(J,I, 146)
      A0158=SAI(J,I)*RC(J,I, 158)
      A0178=SB2(J,I)*RC(J,I, 178)
      A0214=SB7(J,I)*RC(J,I, 214)
      A0248=SAB(J,I)*RC(J,I, 248)
      A0258=SB8(J,I)*RC(J,I, 258)
      A0260=SAH(J,I)*RC(J,I, 260)
      A0282=SB5(J,I)*RC(J,I, 282)
      A0298=SBE(J,I)*RC(J,I, 298)
      A0300=SB3(J,I)*RC(J,I, 300)
      A0322=SAG(J,I)*RC(J,I, 322)
      A0324=SBD(J,I)*RC(J,I, 324)
      A0326=SBC(J,I)*RC(J,I, 326)
      A0328=SBF(J,I)*RC(J,I, 328)
      A0368=SAG(J,I)*RC(J,I, 368)
      A0370=SBB(J,I)*RC(J,I, 370)
      A0376=SAG(J,I)*RC(J,I, 376)
      A0382=SBB(J,I)*RC(J,I, 382)
      A0396=SBH(J,I)*RC(J,I, 396)
      A0410=SBG(J,I)*RC(J,I, 410)
      A0420=SBI(J,I)*RC(J,I, 420)
      A0452=SC2(J,I)*RC(J,I, 452)
      A0454=SC1(J,I)*RC(J,I, 454)
      A0512=SC8(J,I)*SAI(J,I)*RC(J,I, 512)
      A0524=SBJ(J,I)*SB2(J,I)*RC(J,I, 524)
      A0544=SC3(J,I)*SA2(J,I)*RC(J,I, 544)
      RHS1A=+A0005+A0008+A0011+A0025+A0029+A0037+A0039+A0057+A0069+A0077
     1      +A0089+A0121+A0138+A0143+A0145+A0157+A0177+A0213+A0247+A0257
     2      +A0259+A0281+A0297+A0299+A0321+A0323+A0325+A0327+A0367+A0369
     3      +A0375+A0381+A0395+A0409+A0419+A0451+A0453+A0511+A0523+A0543
      RHS1B=+1.000*A0513+1.000*A0515+1.000*A0517
      RHS1=RHO(J,I)*FA5(J,I)+DT*WM(  6)*(RHS1A+RHS1B)
      RHS2=-DT*WM(  6)*(
     1     +A0006+A0007+A0012+A0026+A0030+A0038+A0040+A0058+A0070+A0078
     2     +A0090+A0122+A0137+A0144+A0146+A0158+A0178+A0214+A0248+A0258
     3     +A0260+A0282+A0298+A0300+A0322+A0324+A0326+A0328+A0368+A0370
     4     +A0376+A0382+A0396+A0410+A0420+A0452+A0454+A0512+A0524+A0544
     *      )
      PPP=RHO(J,I)-(EA5(J,I)+WA5(J,I)+PA5(J,I)+QA5(J,I))
      FNEW=SOR2*SA5(J,I)+(RHS1-PA5(J,I)*SA5(J+1,I)-QA5(J,I)*SA5(J-1,I)
     1   -EA5(J,I)*SA5(J,I+1)-WA5(J,I)*SA5(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SA5(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SA5(J,I)=FNEW
      A0026=FNEW*A0026
      A0030=FNEW*A0030
      A0038=FNEW*A0038
      A0078=FNEW*A0078
      A0090=FNEW*A0090
      A0122=FNEW*A0122
      A0137=FNEW*A0137
      A0178=FNEW*A0178
      A0214=FNEW*A0214
      A0282=FNEW*A0282
      A0298=FNEW*A0298
      A0324=FNEW*A0324
      A0326=FNEW*A0326
      A0328=FNEW*A0328
      A0396=FNEW*A0396
      A0452=FNEW*A0452
      A0454=FNEW*A0454
      A0544=FNEW*A0544
C----------------    SOLVING   7 TH SPECIES EQUATION    ----------------
      A0018=RC(J,I,  18)
      A0020=RC(J,I,  20)
      A0021=SA2(J,I)*RC(J,I,  21)
      A0023=SA2(J,I)*RC(J,I,  23)
      A0025=SA2(J,I)*RC(J,I,  25)
      A0027=SA3(J,I)*RC(J,I,  27)
      A0029=SA4(J,I)*RC(J,I,  29)
      A0033=SA6(J,I)*RC(J,I,  33)
      A0036=SA0(J,I)*RC(J,I,  36)
      A0040=SA5(J,I)*RC(J,I,  40)
      A0042=SA4(J,I)*RC(J,I,  42)
      A0045=SA8(J,I)*RC(J,I,  45)
      A0060=SA8(J,I)*RC(J,I,  60)
      A0072=SAA(J,I)*RC(J,I,  72)
      A0073=SAB(J,I)*RC(J,I,  73)
      A0082=SAD(J,I)*RC(J,I,  82)
      A0083=SAC(J,I)*RC(J,I,  83)
      A0095=SAD(J,I)*RC(J,I,  95)
      A0150=SAB(J,I)*RC(J,I, 150)
      A0163=SAJ(J,I)*RC(J,I, 163)
      A0172=SAG(J,I)*RC(J,I, 172)
      A0186=SB2(J,I)*RC(J,I, 186)
      A0187=SAG(J,I)*RC(J,I, 187)
      A0189=SB4(J,I)*RC(J,I, 189)
      A0204=SB1(J,I)*RC(J,I, 204)
      A0250=SAB(J,I)*RC(J,I, 250)
      A0267=SB9(J,I)*RC(J,I, 267)
      A0270=SB8(J,I)*RC(J,I, 270)
      A0287=SB3(J,I)*RC(J,I, 287)
      A0289=SB3(J,I)*RC(J,I, 289)
      A0313=SBB(J,I)*RC(J,I, 313)
      A0315=SBB(J,I)*RC(J,I, 315)
      A0318=SBE(J,I)*RC(J,I, 318)
      A0347=SBA(J,I)*RC(J,I, 347)
      A0349=SBA(J,I)*RC(J,I, 349)
      A0351=SBA(J,I)*RC(J,I, 351)
      A0355=SAI(J,I)*RC(J,I, 355)
      A0362=SBB(J,I)*RC(J,I, 362)
      A0372=SBB(J,I)*RC(J,I, 372)
      A0379=SBC(J,I)*RC(J,I, 379)
      A0393=SBH(J,I)*RC(J,I, 393)
      A0404=SBG(J,I)*RC(J,I, 404)
      A0413=SBH(J,I)*RC(J,I, 413)
      A0427=SBI(J,I)*RC(J,I, 427)
      A0429=SBI(J,I)*RC(J,I, 429)
      A0440=SC1(J,I)*RC(J,I, 440)
      A0442=SC2(J,I)*RC(J,I, 442)
      A0455=SC0(J,I)*RC(J,I, 455)
      A0457=SC0(J,I)*RC(J,I, 457)
      A0464=SBJ(J,I)*RC(J,I, 464)
      A0470=SBJ(J,I)*RC(J,I, 470)
      A0494=SB1(J,I)*SBI(J,I)*RC(J,I, 494)
      A0496=SB2(J,I)*SBG(J,I)*RC(J,I, 496)
      A0497=SC5(J,I)*RC(J,I, 497)
      A0499=SC5(J,I)*RC(J,I, 499)
      RHS1A=+A0017+A0019+A0022+A0024+A0026+A0028+A0030+A0034+A0034+A0035
     1      +A0039+A0041+A0046+A0059+A0071+A0074+A0081+A0084+A0096+A0149
     2      +A0164+A0171+A0185+A0188+A0190+A0203+A0249+A0268+A0269+A0288
     3      +A0290+A0314+A0316+A0317+A0348+A0350+A0352+A0356+A0361+A0371
     4      +A0380+A0394+A0403+A0414+A0428+A0430+A0439+A0441+A0456+A0458
     5      +A0463+A0469+A0493+A0495+A0498+A0500
      RHS1B=0.0
      RHS1=RHO(J,I)*FA6(J,I)+DT*WM(  7)*(RHS1A+RHS1B)
      RHS2=-DT*WM(  7)*(
     1     +A0018+A0020+A0021+A0023+A0025+A0027+A0029+A0033+A0033+A0036
     2     +A0040+A0042+A0045+A0060+A0072+A0073+A0082+A0083+A0095+A0150
     3     +A0163+A0172+A0186+A0187+A0189+A0204+A0250+A0267+A0270+A0287
     4     +A0289+A0313+A0315+A0318+A0347+A0349+A0351+A0355+A0362+A0372
     5     +A0379+A0393+A0404+A0413+A0427+A0429+A0440+A0442+A0455+A0457
     6     +A0464+A0470+A0494+A0496+A0497+A0499
     *      )
      PPP=RHO(J,I)-(EA6(J,I)+WA6(J,I)+PA6(J,I)+QA6(J,I))
      FNEW=SOR2*SA6(J,I)+(RHS1-PA6(J,I)*SA6(J+1,I)-QA6(J,I)*SA6(J-1,I)
     1   -EA6(J,I)*SA6(J,I+1)-WA6(J,I)*SA6(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SA6(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SA6(J,I)=FNEW
      A0033=FNEW*A0033
      A0036=FNEW*A0036
      A0040=FNEW*A0040
      A0042=FNEW*A0042
      A0073=FNEW*A0073
      A0082=FNEW*A0082
      A0083=FNEW*A0083
      A0163=FNEW*A0163
      A0186=FNEW*A0186
      A0189=FNEW*A0189
      A0267=FNEW*A0267
      A0287=FNEW*A0287
      A0313=FNEW*A0313
      A0315=FNEW*A0315
      A0318=FNEW*A0318
      A0347=FNEW*A0347
      A0349=FNEW*A0349
      A0351=FNEW*A0351
      A0379=FNEW*A0379
      A0393=FNEW*A0393
      A0413=FNEW*A0413
      A0429=FNEW*A0429
      A0440=FNEW*A0440
      A0442=FNEW*A0442
      A0455=FNEW*A0455
      A0457=FNEW*A0457
      A0497=FNEW*A0497
      A0499=FNEW*A0499
C----------------    SOLVING   8 TH SPECIES EQUATION    ----------------
      A0032=RC(J,I,  32)
      A0034=SA1(J,I)*RC(J,I,  34)
      A0035=SA2(J,I)*RC(J,I,  35)
      A0037=SA2(J,I)*RC(J,I,  37)
      A0039=SA4(J,I)*RC(J,I,  39)
      A0041=SA3(J,I)*RC(J,I,  41)
      A0074=SAA(J,I)*RC(J,I,  74)
      A0084=SAD(J,I)*RC(J,I,  84)
      A0164=SAI(J,I)*RC(J,I, 164)
      A0190=SAD(J,I)*SA8(J,I)*RC(J,I, 190)
      A0268=SB8(J,I)*RC(J,I, 268)
      A0314=SBE(J,I)*RC(J,I, 314)
      A0316=SB3(J,I)*RC(J,I, 316)
      A0348=SBC(J,I)*RC(J,I, 348)
      A0350=SBD(J,I)*RC(J,I, 350)
      A0352=SBF(J,I)*RC(J,I, 352)
      A0456=SC1(J,I)*RC(J,I, 456)
      A0458=SC2(J,I)*RC(J,I, 458)
      A0498=SB1(J,I)*SBI(J,I)*RC(J,I, 498)
      A0500=SB2(J,I)*SBG(J,I)*RC(J,I, 500)
      RHS1A=+A0031+A0033+A0036+A0038+A0040+A0042+A0073+A0083+A0163+A0189
     1      +A0267+A0313+A0315+A0347+A0349+A0351+A0455+A0457+A0497+A0499
      RHS1B=0.0
      RHS1=RHO(J,I)*FA7(J,I)+DT*WM(  8)*(RHS1A+RHS1B)
      RHS2=-DT*WM(  8)*(
     1     +A0032+A0034+A0035+A0037+A0039+A0041+A0074+A0084+A0164+A0190
     2     +A0268+A0314+A0316+A0348+A0350+A0352+A0456+A0458+A0498+A0500
     *      )
      PPP=RHO(J,I)-(EA7(J,I)+WA7(J,I)+PA7(J,I)+QA7(J,I))
      FNEW=SOR2*SA7(J,I)+(RHS1-PA7(J,I)*SA7(J+1,I)-QA7(J,I)*SA7(J-1,I)
     1   -EA7(J,I)*SA7(J,I+1)-WA7(J,I)*SA7(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SA7(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SA7(J,I)=FNEW
      A0084=FNEW*A0084
      A0314=FNEW*A0314
      A0348=FNEW*A0348
      A0350=FNEW*A0350
      A0352=FNEW*A0352
      A0456=FNEW*A0456
      A0458=FNEW*A0458
C----------------    SOLVING   9 TH SPECIES EQUATION    ----------------
      A0043=SA4(J,I)*RC(J,I,  43)
      A0045=SA6(J,I)*RC(J,I,  45)
      A0047=SA1(J,I)*RC(J,I,  47)
      A0050=SA2(J,I)*RC(J,I,  50)
      A0052=SA0(J,I)*RC(J,I,  52)
      A0054=SA4(J,I)*RC(J,I,  54)
      A0058=SA5(J,I)*RC(J,I,  58)
      A0060=SA6(J,I)*RC(J,I,  60)
      A0062=SAC(J,I)*RC(J,I,  62)
      A0112=SA4(J,I)*SA2(J,I)*RC(J,I, 112)
      A0114=SAB(J,I)*RC(J,I, 114)
      A0124=SA2(J,I)*SA2(J,I)*RC(J,I, 124)
      A0126=SA0(J,I)*RC(J,I, 126)
      A0130=SA4(J,I)*SA2(J,I)*RC(J,I, 130)
      A0134=SA2(J,I)*RC(J,I, 134)
      A0140=SAA(J,I)*RC(J,I, 140)
      A0190=SAD(J,I)*SA7(J,I)*RC(J,I, 190)
      A0208=SAE(J,I)*RC(J,I, 208)
      A0210=SAB(J,I)*RC(J,I, 210)
      A0216=SAD(J,I)*RC(J,I, 216)
      A0222=SAI(J,I)*RC(J,I, 222)
      A0224=SAF(J,I)*RC(J,I, 224)
      A0226=SAA(J,I)*SA2(J,I)*RC(J,I, 226)
      A0228=SA8(J,I)*SA2(J,I)*RC(J,I, 228)
      A0230=SA8(J,I)*SA4(J,I)*RC(J,I, 230)
      A0232=SA9(J,I)*SA2(J,I)*RC(J,I, 232)
      A0236=SB0(J,I)*RC(J,I, 236)
      A0242=SAA(J,I)*RC(J,I, 242)
      A0256=SB8(J,I)*RC(J,I, 256)
      A0284=SAB(J,I)*SA4(J,I)*RC(J,I, 284)
      A0286=SAI(J,I)*SA2(J,I)*RC(J,I, 286)
      A0292=SAD(J,I)*RC(J,I, 292)
      A0296=SAD(J,I)*RC(J,I, 296)
      A0363=SBF(J,I)*RC(J,I, 363)
      A0386=SAG(J,I)*RC(J,I, 386)
      A0412=SBG(J,I)*RC(J,I, 412)
      A0414=SA4(J,I)*SB2(J,I)*RC(J,I, 414)
      A0416=SAD(J,I)*SAA(J,I)*RC(J,I, 416)
      RHS1A=+A0044+A0046+A0048+A0049+A0051+A0053+A0057+A0059+A0061+A0111
     1      +A0113+A0123+A0125+A0129+A0133+A0139+A0189+A0207+A0209+A0215
     2      +A0221+A0223+A0225+A0227+A0227+A0229+A0229+A0231+A0235+A0241
     3      +A0255+A0283+A0285+A0291+A0295+A0364+A0385+A0411+A0413+A0415
     *      +0.5*SOUR(J,I,3)*SA1(J,I)+SOUR(J,I,4)*SA4(J,I)
     *      +SOUR(J,I,5)*SA3(J,I)
      RHS1B=0.0
      RHS1=RHO(J,I)*FA8(J,I)+DT*WM(  9)*(RHS1A+RHS1B)
      RHS2=-DT*WM(  9)*(
     1     +A0043+A0045+A0047+A0050+A0052+A0054+A0058+A0060+A0062+A0112
     2     +A0114+A0124+A0126+A0130+A0134+A0140+A0190+A0208+A0210+A0216
     3     +A0222+A0224+A0226+A0228+A0228+A0230+A0230+A0232+A0236+A0242
     4     +A0256+A0284+A0286+A0292+A0296+A0363+A0386+A0412+A0414+A0416
     *      )
      PPP=RHO(J,I)-(EA8(J,I)+WA8(J,I)+PA8(J,I)+QA8(J,I))
      FNEW=SOR2*SA8(J,I)+(RHS1-PA8(J,I)*SA8(J+1,I)-QA8(J,I)*SA8(J-1,I)
     1   -EA8(J,I)*SA8(J,I+1)-WA8(J,I)*SA8(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SA8(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SA8(J,I)=FNEW
      A0043=FNEW*A0043
      A0045=FNEW*A0045
      A0047=FNEW*A0047
      A0050=FNEW*A0050
      A0052=FNEW*A0052
      A0054=FNEW*A0054
      A0058=FNEW*A0058
      A0060=FNEW*A0060
      A0062=FNEW*A0062
      A0112=FNEW*A0112
      A0114=FNEW*A0114
      A0124=FNEW*A0124
      A0126=FNEW*A0126
      A0130=FNEW*A0130
      A0134=FNEW*A0134
      A0140=FNEW*A0140
      A0222=FNEW*A0222
      A0228=FNEW*A0228
      A0230=FNEW*A0230
      A0256=FNEW*A0256
      A0286=FNEW*A0286
      A0363=FNEW*A0363
      A0412=FNEW*A0412
C----------------    SOLVING  10 TH SPECIES EQUATION    ----------------
      A0044=SA2(J,I)*RC(J,I,  44)
      A0046=SA4(J,I)*RC(J,I,  46)
      A0048=SA3(J,I)*RC(J,I,  48)
      A0056=SA2(J,I)*RC(J,I,  56)
      A0113=SAF(J,I)*RC(J,I, 113)
      A0128=SA0(J,I)*RC(J,I, 128)
      A0139=SB0(J,I)*RC(J,I, 139)
      A0218=SAE(J,I)*RC(J,I, 218)
      A0232=SA8(J,I)*SA2(J,I)*RC(J,I, 232)
      A0240=SB0(J,I)*RC(J,I, 240)
      A0364=SAI(J,I)*RC(J,I, 364)
      RHS1A=+A0043+A0045+A0047+A0055+A0114+A0127+A0140+A0217+A0231+A0239
     1      +A0363
      RHS1B=0.0
      RHS1=RHO(J,I)*FA9(J,I)+DT*WM( 10)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 10)*(
     1     +A0044+A0046+A0048+A0056+A0113+A0128+A0139+A0218+A0232+A0240
     2     +A0364
     *      )
      PPP=RHO(J,I)-(EA9(J,I)+WA9(J,I)+PA9(J,I)+QA9(J,I))
      FNEW=SOR2*SA9(J,I)+(RHS1-PA9(J,I)*SA9(J+1,I)-QA9(J,I)*SA9(J-1,I)
     1   -EA9(J,I)*SA9(J,I+1)-WA9(J,I)*SA9(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SA9(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SA9(J,I)=FNEW
      A0056=FNEW*A0056
      A0113=FNEW*A0113
      A0128=FNEW*A0128
      A0139=FNEW*A0139
      A0232=FNEW*A0232
C----------------    SOLVING  11 TH SPECIES EQUATION    ----------------
      A0049=RC(J,I,  49)
      A0051=SA2(J,I)*RC(J,I,  51)
      A0053=SA3(J,I)*RC(J,I,  53)
      A0055=SA3(J,I)*RC(J,I,  55)
      A0057=SA4(J,I)*RC(J,I,  57)
      A0059=SA1(J,I)*RC(J,I,  59)
      A0061=SAD(J,I)*RC(J,I,  61)
      A0063=SA2(J,I)*RC(J,I,  63)
      A0066=SA0(J,I)*RC(J,I,  66)
      A0068=SA4(J,I)*RC(J,I,  68)
      A0070=SA5(J,I)*RC(J,I,  70)
      A0072=SA6(J,I)*RC(J,I,  72)
      A0074=SA7(J,I)*RC(J,I,  74)
      A0136=SA3(J,I)*RC(J,I, 136)
      A0140=SA8(J,I)*RC(J,I, 140)
      A0180=SAD(J,I)*RC(J,I, 180)
      A0200=SAB(J,I)*RC(J,I, 200)
      A0226=SA8(J,I)*SA2(J,I)*RC(J,I, 226)
      A0242=SA8(J,I)*RC(J,I, 242)
      A0276=SAD(J,I)*RC(J,I, 276)
      A0280=SAB(J,I)*RC(J,I, 280)
      A0288=SAB(J,I)*SA4(J,I)*RC(J,I, 288)
      A0294=SAD(J,I)*RC(J,I, 294)
      A0398=SB5(J,I)*RC(J,I, 398)
      A0411=SBH(J,I)*RC(J,I, 411)
      A0416=SAD(J,I)*SA8(J,I)*RC(J,I, 416)
      A0418=SAI(J,I)*RC(J,I, 418)
      RHS1A=+A0050+A0052+A0054+A0056+A0058+A0060+A0062+A0064+A0065+A0067
     1      +A0069+A0071+A0073+A0135+A0139+A0179+A0199+A0225+A0241+A0275
     2      +A0279+A0287+A0293+A0397+A0412+A0415+A0417
      RHS1B=0.0
      RHS1=RHO(J,I)*FAA(J,I)+DT*WM( 11)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 11)*(
     1     +A0049+A0051+A0053+A0055+A0057+A0059+A0061+A0063+A0066+A0068
     2     +A0070+A0072+A0074+A0136+A0140+A0180+A0200+A0226+A0242+A0276
     3     +A0280+A0288+A0294+A0398+A0411+A0416+A0418
     *      )
      PPP=RHO(J,I)-(EAA(J,I)+WAA(J,I)+PAA(J,I)+QAA(J,I))
      FNEW=SOR2*SAA(J,I)+(RHS1-PAA(J,I)*SAA(J+1,I)-QAA(J,I)*SAA(J-1,I)
     1   -EAA(J,I)*SAA(J,I+1)-WAA(J,I)*SAA(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SAA(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SAA(J,I)=FNEW
      A0061=FNEW*A0061
      A0063=FNEW*A0063
      A0066=FNEW*A0066
      A0068=FNEW*A0068
      A0070=FNEW*A0070
      A0072=FNEW*A0072
      A0074=FNEW*A0074
      A0136=FNEW*A0136
      A0140=FNEW*A0140
      A0226=FNEW*A0226
      A0242=FNEW*A0242
      A0411=FNEW*A0411
C----------------    SOLVING  12 TH SPECIES EQUATION    ----------------
      A0064=RC(J,I,  64)
      A0065=SA2(J,I)*RC(J,I,  65)
      A0067=SA3(J,I)*RC(J,I,  67)
      A0069=SA4(J,I)*RC(J,I,  69)
      A0071=SA1(J,I)*RC(J,I,  71)
      A0073=SA6(J,I)*RC(J,I,  73)
      A0092=SA2(J,I)*RC(J,I,  92)
      A0098=SA4(J,I)*RC(J,I,  98)
      A0110=SA2(J,I)*RC(J,I, 110)
      A0114=SA8(J,I)*RC(J,I, 114)
      A0120=SA2(J,I)*RC(J,I, 120)
      A0138=SA2(J,I)*RC(J,I, 138)
      A0142=SA0(J,I)*RC(J,I, 142)
      A0146=SA5(J,I)*RC(J,I, 146)
      A0148=SA4(J,I)*RC(J,I, 148)
      A0150=SA6(J,I)*RC(J,I, 150)
      A0152=SA2(J,I)*RC(J,I, 152)
      A0170=SAD(J,I)*RC(J,I, 170)
      A0200=SAA(J,I)*RC(J,I, 200)
      A0210=SA8(J,I)*RC(J,I, 210)
      A0244=SA0(J,I)*RC(J,I, 244)
      A0248=SA5(J,I)*RC(J,I, 248)
      A0250=SA6(J,I)*RC(J,I, 250)
      A0252=SA2(J,I)*RC(J,I, 252)
      A0280=SAA(J,I)*RC(J,I, 280)
      A0284=SA8(J,I)*SA4(J,I)*RC(J,I, 284)
      A0288=SAA(J,I)*SA4(J,I)*RC(J,I, 288)
      A0360=SAD(J,I)*RC(J,I, 360)
      A0430=SA4(J,I)*SB2(J,I)*RC(J,I, 430)
      RHS1A=+A0063+A0066+A0068+A0070+A0072+A0074+A0091+A0097+A0109+A0113
     1      +A0119+A0137+A0141+A0145+A0147+A0149+A0151+A0169+A0199+A0209
     2      +A0243+A0247+A0249+A0251+A0279+A0283+A0287+A0359+A0429
      RHS1B=0.0
      RHS1=RHO(J,I)*FAB(J,I)+DT*WM( 12)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 12)*(
     1     +A0064+A0065+A0067+A0069+A0071+A0073+A0092+A0098+A0110+A0114
     2     +A0120+A0138+A0142+A0146+A0148+A0150+A0152+A0170+A0200+A0210
     3     +A0244+A0248+A0250+A0252+A0280+A0284+A0288+A0360+A0430
     *      )
      PPP=RHO(J,I)-(EAB(J,I)+WAB(J,I)+PAB(J,I)+QAB(J,I))
      FNEW=SOR2*SAB(J,I)+(RHS1-PAB(J,I)*SAB(J+1,I)-QAB(J,I)*SAB(J-1,I)
     1   -EAB(J,I)*SAB(J,I+1)-WAB(J,I)*SAB(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SAB(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SAB(J,I)=FNEW
      A0092=FNEW*A0092
      A0098=FNEW*A0098
      A0110=FNEW*A0110
      A0114=FNEW*A0114
      A0120=FNEW*A0120
      A0138=FNEW*A0138
      A0142=FNEW*A0142
      A0146=FNEW*A0146
      A0148=FNEW*A0148
      A0150=FNEW*A0150
      A0152=FNEW*A0152
      A0200=FNEW*A0200
      A0210=FNEW*A0210
      A0244=FNEW*A0244
      A0248=FNEW*A0248
      A0250=FNEW*A0250
      A0252=FNEW*A0252
      A0280=FNEW*A0280
      A0284=FNEW*A0284
      A0288=FNEW*A0288
C----------------    SOLVING  13 TH SPECIES EQUATION    ----------------
      A0062=SA8(J,I)*RC(J,I,  62)
      A0075=SA2(J,I)*RC(J,I,  75)
      A0077=SA4(J,I)*RC(J,I,  77)
      A0079=SA3(J,I)*RC(J,I,  79)
      A0081=SA1(J,I)*RC(J,I,  81)
      A0083=SA6(J,I)*RC(J,I,  83)
      A0106=RC(J,I, 106)
      A0160=SAI(J,I)*RC(J,I, 160)
      A0310=SBE(J,I)*RC(J,I, 310)
      A0312=SB3(J,I)*RC(J,I, 312)
      A0342=SBD(J,I)*RC(J,I, 342)
      A0344=SBC(J,I)*RC(J,I, 344)
      A0346=SBF(J,I)*RC(J,I, 346)
      A0406=SBG(J,I)*RC(J,I, 406)
      A0484=SB1(J,I)*SB2(J,I)*RC(J,I, 484)
      RHS1A=+A0061+A0076+A0078+A0080+A0082+A0084+A0105+A0159+A0309+A0311
     1      +A0341+A0343+A0345+A0405+A0483
      RHS1B=0.0
      RHS1=RHO(J,I)*FAC(J,I)+DT*WM( 13)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 13)*(
     1     +A0062+A0075+A0077+A0079+A0081+A0083+A0106+A0160+A0310+A0312
     2     +A0342+A0344+A0346+A0406+A0484
     *      )
      PPP=RHO(J,I)-(EAC(J,I)+WAC(J,I)+PAC(J,I)+QAC(J,I))
      FNEW=SOR2*SAC(J,I)+(RHS1-PAC(J,I)*SAC(J+1,I)-QAC(J,I)*SAC(J-1,I)
     1   -EAC(J,I)*SAC(J,I+1)-WAC(J,I)*SAC(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SAC(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SAC(J,I)=FNEW
      A0062=FNEW*A0062
      A0075=FNEW*A0075
      A0077=FNEW*A0077
      A0079=FNEW*A0079
      A0081=FNEW*A0081
      A0083=FNEW*A0083
      A0106=FNEW*A0106
      A0160=FNEW*A0160
      A0310=FNEW*A0310
      A0312=FNEW*A0312
      A0342=FNEW*A0342
      A0344=FNEW*A0344
      A0346=FNEW*A0346
      A0406=FNEW*A0406
      A0484=FNEW*A0484
C----------------    SOLVING  14 TH SPECIES EQUATION    ----------------
      A0061=SAA(J,I)*RC(J,I,  61)
      A0076=SA0(J,I)*RC(J,I,  76)
      A0078=SA5(J,I)*RC(J,I,  78)
      A0080=SA4(J,I)*RC(J,I,  80)
      A0082=SA6(J,I)*RC(J,I,  82)
      A0084=SA7(J,I)*RC(J,I,  84)
      A0085=SA2(J,I)*RC(J,I,  85)
      A0087=SA2(J,I)*RC(J,I,  87)
      A0089=SA4(J,I)*RC(J,I,  89)
      A0091=SA3(J,I)*RC(J,I,  91)
      A0093=SAE(J,I)*RC(J,I,  93)
      A0095=SA6(J,I)*RC(J,I,  95)
      A0097=SA1(J,I)*RC(J,I,  97)
      A0099=SA1(J,I)*RC(J,I,  99)
      A0101=SAD(J,I)*RC(J,I, 101)
      A0103=SAD(J,I)*RC(J,I, 103)
      A0105=SA2(J,I)*RC(J,I, 105)
      A0107=SAD(J,I)*RC(J,I, 107)
      A0159=SAJ(J,I)*RC(J,I, 159)
      A0170=SAB(J,I)*RC(J,I, 170)
      A0180=SAA(J,I)*RC(J,I, 180)
      A0190=SA8(J,I)*SA7(J,I)*RC(J,I, 190)
      A0216=SA8(J,I)*RC(J,I, 216)
      A0221=SB5(J,I)*RC(J,I, 221)
      A0246=SA4(J,I)*RC(J,I, 246)
      A0272=SA4(J,I)*RC(J,I, 272)
      A0276=SAA(J,I)*RC(J,I, 276)
      A0285=SB3(J,I)*RC(J,I, 285)
      A0292=SA8(J,I)*RC(J,I, 292)
      A0294=SAA(J,I)*RC(J,I, 294)
      A0296=SA8(J,I)*RC(J,I, 296)
      A0309=SBB(J,I)*RC(J,I, 309)
      A0311=SBB(J,I)*RC(J,I, 311)
      A0320=SB8(J,I)*RC(J,I, 320)
      A0341=SBA(J,I)*RC(J,I, 341)
      A0343=SBA(J,I)*RC(J,I, 343)
      A0345=SBA(J,I)*RC(J,I, 345)
      A0360=SAB(J,I)*RC(J,I, 360)
      A0366=SB8(J,I)*RC(J,I, 366)
      A0378=SB8(J,I)*RC(J,I, 378)
      A0387=SB1(J,I)*RC(J,I, 387)
      A0390=SB6(J,I)*RC(J,I, 390)
      A0405=SBI(J,I)*RC(J,I, 405)
      A0407=SB1(J,I)*RC(J,I, 407)
      A0416=SAA(J,I)*SA8(J,I)*RC(J,I, 416)
      A0422=SB5(J,I)*SA2(J,I)*RC(J,I, 422)
      A0431=SB2(J,I)*RC(J,I, 431)
      A0434=SAG(J,I)*RC(J,I, 434)
      A0435=SB2(J,I)*RC(J,I, 435)
      A0438=SAI(J,I)*RC(J,I, 438)
      A0466=SAG(J,I)*RC(J,I, 466)
      A0483=SC3(J,I)*RC(J,I, 483)
      A0485=SBH(J,I)*RC(J,I, 485)
      A0528=SAG(J,I)*SAG(J,I)*RC(J,I, 528)
      A0534=SA0(J,I)*SC3(J,I)*RC(J,I, 534)
      A0536=SBI(J,I)*RC(J,I, 536)
      A0540=SBJ(J,I)*RC(J,I, 540)
      RHS1A=+A0062+A0075+A0077+A0079+A0081+A0083+A0086+A0088+A0090+A0092
     1      +A0094+A0096+A0098+A0100+A0102+A0102+A0104+A0104+A0106+A0108
     2      +A0108+A0160+A0169+A0179+A0189+A0215+A0222+A0245+A0271+A0275
     3      +A0286+A0291+A0293+A0295+A0310+A0312+A0319+A0342+A0344+A0346
     4      +A0359+A0365+A0377+A0388+A0389+A0406+A0408+A0415+A0421+A0432
     5      +A0433+A0436+A0437+A0465+A0484+A0486+A0527+A0533+A0535+A0539
      RHS1B=+1.000*A0501+1.000*A0505+1.000*A0507+1.000*A0513+1.000*A0517
     *
      RHS1=RHO(J,I)*FAD(J,I)+DT*WM( 14)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 14)*(
     1     +A0061+A0076+A0078+A0080+A0082+A0084+A0085+A0087+A0089+A0091
     2     +A0093+A0095+A0097+A0099+A0101+A0101+A0103+A0103+A0105+A0107
     3     +A0107+A0159+A0170+A0180+A0190+A0216+A0221+A0246+A0272+A0276
     4     +A0285+A0292+A0294+A0296+A0309+A0311+A0320+A0341+A0343+A0345
     5     +A0360+A0366+A0378+A0387+A0390+A0405+A0407+A0416+A0422+A0431
     6     +A0434+A0435+A0438+A0466+A0483+A0485+A0528+A0534+A0536+A0540
     *      )
      PPP=RHO(J,I)-(EAD(J,I)+WAD(J,I)+PAD(J,I)+QAD(J,I))
      FNEW=SOR2*SAD(J,I)+(RHS1-PAD(J,I)*SAD(J+1,I)-QAD(J,I)*SAD(J-1,I)
     1   -EAD(J,I)*SAD(J,I+1)-WAD(J,I)*SAD(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SAD(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SAD(J,I)=FNEW
      A0085=FNEW*A0085
      A0087=FNEW*A0087
      A0089=FNEW*A0089
      A0095=FNEW*A0095
      A0099=FNEW*A0099
      A0101=FNEW*A0101
      A0103=FNEW*A0103
      A0107=FNEW*A0107
      A0159=FNEW*A0159
      A0170=FNEW*A0170
      A0180=FNEW*A0180
      A0190=FNEW*A0190
      A0216=FNEW*A0216
      A0221=FNEW*A0221
      A0246=FNEW*A0246
      A0272=FNEW*A0272
      A0276=FNEW*A0276
      A0285=FNEW*A0285
      A0292=FNEW*A0292
      A0294=FNEW*A0294
      A0296=FNEW*A0296
      A0311=FNEW*A0311
      A0360=FNEW*A0360
      A0405=FNEW*A0405
      A0416=FNEW*A0416
      A0483=FNEW*A0483
C----------------    SOLVING  15 TH SPECIES EQUATION    ----------------
      A0086=SA0(J,I)*RC(J,I,  86)
      A0093=SAD(J,I)*RC(J,I,  93)
      A0116=RC(J,I, 116)
      A0117=SA2(J,I)*RC(J,I, 117)
      A0119=SA4(J,I)*RC(J,I, 119)
      A0121=SA4(J,I)*RC(J,I, 121)
      A0123=SA3(J,I)*RC(J,I, 123)
      A0125=SA3(J,I)*RC(J,I, 125)
      A0127=SA1(J,I)*RC(J,I, 127)
      A0129=SA1(J,I)*RC(J,I, 129)
      A0131=SAE(J,I)*RC(J,I, 131)
      A0208=SA8(J,I)*RC(J,I, 208)
      A0218=SA9(J,I)*RC(J,I, 218)
      RHS1A=+A0085+A0094+A0115+A0118+A0120+A0122+A0124+A0126+A0128+A0130
     1      +A0132+A0132+A0207+A0217
      RHS1B=0.0
      RHS1=RHO(J,I)*FAE(J,I)+DT*WM( 15)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 15)*(
     1     +A0086+A0093+A0116+A0117+A0119+A0121+A0123+A0125+A0127+A0129
     2     +A0131+A0131+A0208+A0218
     *      )
      PPP=RHO(J,I)-(EAE(J,I)+WAE(J,I)+PAE(J,I)+QAE(J,I))
      FNEW=SOR2*SAE(J,I)+(RHS1-PAE(J,I)*SAE(J+1,I)-QAE(J,I)*SAE(J-1,I)
     1   -EAE(J,I)*SAE(J,I+1)-WAE(J,I)*SAE(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SAE(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SAE(J,I)=FNEW
      A0093=FNEW*A0093
      A0116=FNEW*A0116
      A0117=FNEW*A0117
      A0121=FNEW*A0121
      A0131=FNEW*A0131
      A0208=FNEW*A0208
      A0218=FNEW*A0218
C----------------    SOLVING  16 TH SPECIES EQUATION    ----------------
      A0088=SA0(J,I)*RC(J,I,  88)
      A0090=SA5(J,I)*RC(J,I,  90)
      A0109=SA4(J,I)*RC(J,I, 109)
      A0111=SA1(J,I)*RC(J,I, 111)
      A0113=SA9(J,I)*RC(J,I, 113)
      A0115=RC(J,I, 115)
      A0144=SA5(J,I)*RC(J,I, 144)
      A0224=SA8(J,I)*RC(J,I, 224)
      RHS1A=+A0087+A0089+A0110+A0112+A0114+A0116+A0143+A0223
      RHS1B=0.0
      RHS1=RHO(J,I)*FAF(J,I)+DT*WM( 16)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 16)*(
     1     +A0088+A0090+A0109+A0111+A0113+A0115+A0144+A0224
     *      )
      PPP=RHO(J,I)-(EAF(J,I)+WAF(J,I)+PAF(J,I)+QAF(J,I))
      FNEW=SOR2*SAF(J,I)+(RHS1-PAF(J,I)*SAF(J+1,I)-QAF(J,I)*SAF(J-1,I)
     1   -EAF(J,I)*SAF(J,I+1)-WAF(J,I)*SAF(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SAF(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SAF(J,I)=FNEW
      A0144=FNEW*A0144
      A0224=FNEW*A0224
C----------------    SOLVING  17 TH SPECIES EQUATION    ----------------
      A0094=SA2(J,I)*RC(J,I,  94)
      A0102=SA0(J,I)*RC(J,I, 102)
      A0166=SA0(J,I)*RC(J,I, 166)
      A0168=SA4(J,I)*RC(J,I, 168)
      A0172=SA6(J,I)*RC(J,I, 172)
      A0174=SA2(J,I)*RC(J,I, 174)
      A0175=SA2(J,I)*RC(J,I, 175)
      A0177=SA4(J,I)*RC(J,I, 177)
      A0179=SA3(J,I)*RC(J,I, 179)
      A0181=SA3(J,I)*RC(J,I, 181)
      A0183=SAG(J,I)*RC(J,I, 183)
      A0185=SA1(J,I)*RC(J,I, 185)
      A0187=SA6(J,I)*RC(J,I, 187)
      A0191=RC(J,I, 191)
      A0193=RC(J,I, 193)
      A0322=SA5(J,I)*RC(J,I, 322)
      A0353=SA4(J,I)*RC(J,I, 353)
      A0368=SA5(J,I)*RC(J,I, 368)
      A0376=SA5(J,I)*RC(J,I, 376)
      A0386=SA8(J,I)*RC(J,I, 386)
      A0434=SAD(J,I)*RC(J,I, 434)
      A0466=SAD(J,I)*RC(J,I, 466)
      A0478=SB2(J,I)*RC(J,I, 478)
      A0490=SBG(J,I)*RC(J,I, 490)
      A0522=SBJ(J,I)*RC(J,I, 522)
      A0528=SAG(J,I)*SAD(J,I)*RC(J,I, 528)
      A0532=SA0(J,I)*SBI(J,I)*RC(J,I, 532)
      A0538=SAI(J,I)*RC(J,I, 538)
      A0542=SA0(J,I)*SB2(J,I)*RC(J,I, 542)
      RHS1A=+A0093+A0101+A0165+A0167+A0171+A0173+A0176+A0178+A0180+A0182
     1      +A0184+A0184+A0186+A0188+A0192+A0194+A0321+A0354+A0367+A0375
     2      +A0385+A0433+A0465+A0477+A0489+A0521+A0527+A0527+A0531+A0537
     3      +A0541
      RHS1B=+1.000*A0501+1.000*A0501+1.000*A0503+1.000*A0505+1.000*A0507
     *      +1.000*A0507+1.000*A0507+1.000*A0513+1.000*A0515+1.000*A0517
     *      +1.000*A0517+1.000*A0517
      RHS1=RHO(J,I)*FAG(J,I)+DT*WM( 17)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 17)*(
     1     +A0094+A0102+A0166+A0168+A0172+A0174+A0175+A0177+A0179+A0181
     2     +A0183+A0183+A0185+A0187+A0191+A0193+A0322+A0353+A0368+A0376
     3     +A0386+A0434+A0466+A0478+A0490+A0522+A0528+A0528+A0532+A0538
     4     +A0542
     *      )
      PPP=RHO(J,I)-(EAG(J,I)+WAG(J,I)+PAG(J,I)+QAG(J,I))
      FNEW=SOR2*SAG(J,I)+(RHS1-PAG(J,I)*SAG(J+1,I)-QAG(J,I)*SAG(J-1,I)
     1   -EAG(J,I)*SAG(J,I+1)-WAG(J,I)*SAG(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SAG(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SAG(J,I)=FNEW
      A0166=FNEW*A0166
      A0168=FNEW*A0168
      A0172=FNEW*A0172
      A0174=FNEW*A0174
      A0175=FNEW*A0175
      A0177=FNEW*A0177
      A0181=FNEW*A0181
      A0183=FNEW*A0183
      A0185=FNEW*A0185
      A0187=FNEW*A0187
      A0191=FNEW*A0191
      A0193=FNEW*A0193
      A0322=FNEW*A0322
      A0353=FNEW*A0353
      A0368=FNEW*A0368
      A0376=FNEW*A0376
      A0386=FNEW*A0386
      A0434=FNEW*A0434
      A0466=FNEW*A0466
      A0528=FNEW*A0528
C----------------    SOLVING  18 TH SPECIES EQUATION    ----------------
      A0096=SA4(J,I)*RC(J,I,  96)
      A0100=SA3(J,I)*RC(J,I, 100)
      A0141=SA2(J,I)*RC(J,I, 141)
      A0143=SA2(J,I)*RC(J,I, 143)
      A0145=SA4(J,I)*RC(J,I, 145)
      A0147=SA3(J,I)*RC(J,I, 147)
      A0149=SA1(J,I)*RC(J,I, 149)
      A0151=RC(J,I, 151)
      A0253=RC(J,I, 253)
      A0260=SA5(J,I)*RC(J,I, 260)
      A0264=SA0(J,I)*RC(J,I, 264)
      RHS1A=+A0095+A0099+A0142+A0144+A0146+A0148+A0150+A0152+A0254+A0259
     1      +A0263
      RHS1B=0.0
      RHS1=RHO(J,I)*FAH(J,I)+DT*WM( 18)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 18)*(
     1     +A0096+A0100+A0141+A0143+A0145+A0147+A0149+A0151+A0253+A0260
     2     +A0264
     *      )
      PPP=RHO(J,I)-(EAH(J,I)+WAH(J,I)+PAH(J,I)+QAH(J,I))
      FNEW=SOR2*SAH(J,I)+(RHS1-PAH(J,I)*SAH(J+1,I)-QAH(J,I)*SAH(J-1,I)
     1   -EAH(J,I)*SAH(J,I+1)-WAH(J,I)*SAH(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SAH(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SAH(J,I)=FNEW
      A0253=FNEW*A0253
      A0260=FNEW*A0260
      A0264=FNEW*A0264
C----------------    SOLVING  19 TH SPECIES EQUATION    ----------------
      A0104=SA2(J,I)*RC(J,I, 104)
      A0154=SA0(J,I)*RC(J,I, 154)
      A0156=SA4(J,I)*RC(J,I, 156)
      A0158=SA5(J,I)*RC(J,I, 158)
      A0160=SAC(J,I)*RC(J,I, 160)
      A0162=SA2(J,I)*RC(J,I, 162)
      A0164=SA7(J,I)*RC(J,I, 164)
      A0165=SA2(J,I)*RC(J,I, 165)
      A0167=SA3(J,I)*RC(J,I, 167)
      A0169=SA3(J,I)*RC(J,I, 169)
      A0171=SA1(J,I)*RC(J,I, 171)
      A0173=RC(J,I, 173)
      A0184=SB2(J,I)*RC(J,I, 184)
      A0222=SA8(J,I)*RC(J,I, 222)
      A0286=SA8(J,I)*SA2(J,I)*RC(J,I, 286)
      A0355=SA6(J,I)*RC(J,I, 355)
      A0364=SA9(J,I)*RC(J,I, 364)
      A0418=SAA(J,I)*RC(J,I, 418)
      A0438=SAD(J,I)*RC(J,I, 438)
      A0510=SA0(J,I)*SC8(J,I)*RC(J,I, 510)
      A0512=SA5(J,I)*SC8(J,I)*RC(J,I, 512)
      A0520=SBI(J,I)*RC(J,I, 520)
      A0530=SBJ(J,I)*RC(J,I, 530)
      A0538=SAG(J,I)*RC(J,I, 538)
      RHS1A=+A0103+A0153+A0155+A0157+A0159+A0161+A0163+A0166+A0168+A0170
     1      +A0172+A0174+A0183+A0221+A0285+A0356+A0363+A0417+A0437+A0509
     2      +A0511+A0519+A0529+A0537
      RHS1B=+1.000*A0501+1.000*A0503+1.000*A0515
      RHS1=RHO(J,I)*FAI(J,I)+DT*WM( 19)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 19)*(
     1     +A0104+A0154+A0156+A0158+A0160+A0162+A0164+A0165+A0167+A0169
     2     +A0171+A0173+A0184+A0222+A0286+A0355+A0364+A0418+A0438+A0510
     3     +A0512+A0520+A0530+A0538
     *      )
      PPP=RHO(J,I)-(EAI(J,I)+WAI(J,I)+PAI(J,I)+QAI(J,I))
      FNEW=SOR2*SAI(J,I)+(RHS1-PAI(J,I)*SAI(J+1,I)-QAI(J,I)*SAI(J-1,I)
     1   -EAI(J,I)*SAI(J,I+1)-WAI(J,I)*SAI(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SAI(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SAI(J,I)=FNEW
      A0154=FNEW*A0154
      A0156=FNEW*A0156
      A0158=FNEW*A0158
      A0160=FNEW*A0160
      A0162=FNEW*A0162
      A0164=FNEW*A0164
      A0222=FNEW*A0222
      A0286=FNEW*A0286
      A0355=FNEW*A0355
      A0364=FNEW*A0364
      A0418=FNEW*A0418
      A0438=FNEW*A0438
      A0510=FNEW*A0510
      A0512=FNEW*A0512
      A0538=FNEW*A0538
C----------------    SOLVING  20 TH SPECIES EQUATION    ----------------
      A0108=RC(J,I, 108)
      A0153=SA2(J,I)*RC(J,I, 153)
      A0155=SA3(J,I)*RC(J,I, 155)
      A0157=SA4(J,I)*RC(J,I, 157)
      A0159=SAD(J,I)*RC(J,I, 159)
      A0161=RC(J,I, 161)
      A0163=SA6(J,I)*RC(J,I, 163)
      RHS1A=+A0107+A0154+A0156+A0158+A0160+A0162+A0164
      RHS1B=0.0
      RHS1=RHO(J,I)*FAJ(J,I)+DT*WM( 20)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 20)*(
     1     +A0108+A0153+A0155+A0157+A0159+A0161+A0163
     *      )
      PPP=RHO(J,I)-(EAJ(J,I)+WAJ(J,I)+PAJ(J,I)+QAJ(J,I))
      FNEW=SOR2*SAJ(J,I)+(RHS1-PAJ(J,I)*SAJ(J+1,I)-QAJ(J,I)*SAJ(J-1,I)
     1   -EAJ(J,I)*SAJ(J,I+1)-WAJ(J,I)*SAJ(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SAJ(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SAJ(J,I)=FNEW
C----------------    SOLVING  21 ST SPECIES EQUATION    ----------------
      A0118=SA0(J,I)*RC(J,I, 118)
      A0122=SA5(J,I)*RC(J,I, 122)
      A0133=SA3(J,I)*RC(J,I, 133)
      A0135=SA1(J,I)*RC(J,I, 135)
      A0137=SA5(J,I)*RC(J,I, 137)
      A0139=SA9(J,I)*RC(J,I, 139)
      A0236=SA8(J,I)*RC(J,I, 236)
      A0240=SA9(J,I)*RC(J,I, 240)
      RHS1A=+A0117+A0121+A0134+A0136+A0138+A0140+A0235+A0239
      RHS1B=0.0
      RHS1=RHO(J,I)*FB0(J,I)+DT*WM( 21)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 21)*(
     1     +A0118+A0122+A0133+A0135+A0137+A0139+A0236+A0240
     *      )
      PPP=RHO(J,I)-(EB0(J,I)+WB0(J,I)+PB0(J,I)+QB0(J,I))
      FNEW=SOR2*SB0(J,I)+(RHS1-PB0(J,I)*SB0(J+1,I)-QB0(J,I)*SB0(J-1,I)
     1   -EB0(J,I)*SB0(J,I+1)-WB0(J,I)*SB0(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SB0(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SB0(J,I)=FNEW
      A0236=FNEW*A0236
      A0240=FNEW*A0240
C----------------    SOLVING  22 ND SPECIES EQUATION    ----------------
      A0132=SA2(J,I)*SA2(J,I)*RC(J,I, 132)
      A0194=SA0(J,I)*RC(J,I, 194)
      A0196=SA0(J,I)*RC(J,I, 196)
      A0198=SA2(J,I)*RC(J,I, 198)
      A0204=SA6(J,I)*RC(J,I, 204)
      A0205=SA3(J,I)*RC(J,I, 205)
      A0207=SA3(J,I)*RC(J,I, 207)
      A0209=SA1(J,I)*RC(J,I, 209)
      A0211=SA4(J,I)*RC(J,I, 211)
      A0213=SA4(J,I)*RC(J,I, 213)
      A0387=SAD(J,I)*RC(J,I, 387)
      A0407=SAD(J,I)*RC(J,I, 407)
      A0472=SB2(J,I)*SA2(J,I)*RC(J,I, 472)
      A0480=SA0(J,I)*SB2(J,I)*RC(J,I, 480)
      A0484=SAC(J,I)*SB2(J,I)*RC(J,I, 484)
      A0488=SBJ(J,I)*RC(J,I, 488)
      A0494=SBI(J,I)*SA6(J,I)*RC(J,I, 494)
      A0498=SBI(J,I)*SA7(J,I)*RC(J,I, 498)
      RHS1A=+A0131+A0193+A0195+A0197+A0203+A0206+A0208+A0210+A0212+A0214
     1      +A0388+A0408+A0471+A0479+A0483+A0487+A0493+A0497
      RHS1B=0.0
      RHS1=RHO(J,I)*FB1(J,I)+DT*WM( 22)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 22)*(
     1     +A0132+A0194+A0196+A0198+A0204+A0205+A0207+A0209+A0211+A0213
     2     +A0387+A0407+A0472+A0480+A0484+A0488+A0494+A0498
     *     +SOUR(J,I,1)
     *      )
      PPP=RHO(J,I)-(EB1(J,I)+WB1(J,I)+PB1(J,I)+QB1(J,I))
      FNEW=SOR2*SB1(J,I)+(RHS1-PB1(J,I)*SB1(J+1,I)-QB1(J,I)*SB1(J-1,I)
     1   -EB1(J,I)*SB1(J,I+1)-WB1(J,I)*SB1(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SB1(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SB1(J,I)=FNEW
      A0196=FNEW*A0196
      A0198=FNEW*A0198
      A0204=FNEW*A0204
      A0205=FNEW*A0205
      A0211=FNEW*A0211
      A0213=FNEW*A0213
      A0387=FNEW*A0387
      A0407=FNEW*A0407
C----------------    SOLVING  23 RD SPECIES EQUATION    ----------------
      A0176=SA0(J,I)*RC(J,I, 176)
      A0178=SA5(J,I)*RC(J,I, 178)
      A0184=SAI(J,I)*RC(J,I, 184)
      A0186=SA6(J,I)*RC(J,I, 186)
      A0192=SA2(J,I)*RC(J,I, 192)
      A0195=SA2(J,I)*RC(J,I, 195)
      A0197=RC(J,I, 197)
      A0199=SA1(J,I)*RC(J,I, 199)
      A0201=SA1(J,I)*RC(J,I, 201)
      A0203=SA1(J,I)*RC(J,I, 203)
      A0414=SA4(J,I)*SA8(J,I)*RC(J,I, 414)
      A0430=SA4(J,I)*SAB(J,I)*RC(J,I, 430)
      A0431=SAD(J,I)*RC(J,I, 431)
      A0435=SAD(J,I)*RC(J,I, 435)
      A0472=SB1(J,I)*SA2(J,I)*RC(J,I, 472)
      A0474=SB2(J,I)*RC(J,I, 474)
      A0475=SB2(J,I)*RC(J,I, 475)
      A0478=SAG(J,I)*RC(J,I, 478)
      A0480=SA0(J,I)*SB1(J,I)*RC(J,I, 480)
      A0484=SAC(J,I)*SB1(J,I)*RC(J,I, 484)
      A0492=SBI(J,I)*RC(J,I, 492)
      A0496=SBG(J,I)*SA6(J,I)*RC(J,I, 496)
      A0500=SBG(J,I)*SA7(J,I)*RC(J,I, 500)
      A0524=SA5(J,I)*SBJ(J,I)*RC(J,I, 524)
      A0526=SA0(J,I)*SBJ(J,I)*RC(J,I, 526)
      A0542=SA0(J,I)*SAG(J,I)*RC(J,I, 542)
      RHS1A=+A0175+A0177+A0183+A0185+A0191+A0196+A0198+A0200+A0202+A0204
     1      +A0413+A0429+A0432+A0436+A0471+A0473+A0473+A0476+A0476+A0477
     2      +A0479+A0483+A0491+A0495+A0499+A0523+A0525+A0541
      RHS1B=0.0
      RHS1=RHO(J,I)*FB2(J,I)+DT*WM( 23)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 23)*(
     1     +A0176+A0178+A0184+A0186+A0192+A0195+A0197+A0199+A0201+A0203
     2     +A0414+A0430+A0431+A0435+A0472+A0474+A0474+A0475+A0475+A0478
     3     +A0480+A0484+A0492+A0496+A0500+A0524+A0526+A0542
     *      )
      PPP=RHO(J,I)-(EB2(J,I)+WB2(J,I)+PB2(J,I)+QB2(J,I))
      FNEW=SOR2*SB2(J,I)+(RHS1-PB2(J,I)*SB2(J+1,I)-QB2(J,I)*SB2(J-1,I)
     1   -EB2(J,I)*SB2(J,I+1)-WB2(J,I)*SB2(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SB2(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SB2(J,I)=FNEW
      A0201=FNEW*A0201
      A0414=FNEW*A0414
      A0430=FNEW*A0430
      A0431=FNEW*A0431
      A0435=FNEW*A0435
      A0472=FNEW*A0472
      A0474=FNEW*A0474
      A0475=FNEW*A0475
      A0478=FNEW*A0478
      A0480=FNEW*A0480
      A0484=FNEW*A0484
      A0542=FNEW*A0542
C----------------    SOLVING  24 TH SPECIES EQUATION    ----------------
      A0182=SA2(J,I)*RC(J,I, 182)
      A0202=SA3(J,I)*RC(J,I, 202)
      A0273=RC(J,I, 273)
      A0275=SA2(J,I)*RC(J,I, 275)
      A0277=SA2(J,I)*RC(J,I, 277)
      A0279=SA3(J,I)*RC(J,I, 279)
      A0281=SA4(J,I)*RC(J,I, 281)
      A0283=SA1(J,I)*RC(J,I, 283)
      A0285=SAD(J,I)*RC(J,I, 285)
      A0287=SA6(J,I)*RC(J,I, 287)
      A0289=SA6(J,I)*RC(J,I, 289)
      A0291=RC(J,I, 291)
      A0300=SA5(J,I)*RC(J,I, 300)
      A0304=SA4(J,I)*RC(J,I, 304)
      A0308=SA0(J,I)*RC(J,I, 308)
      A0312=SAC(J,I)*RC(J,I, 312)
      A0316=SA7(J,I)*RC(J,I, 316)
      RHS1A=+A0181+A0201+A0274+A0276+A0278+A0280+A0282+A0284+A0286+A0288
     1      +A0290+A0292+A0299+A0303+A0307+A0311+A0315
      RHS1B=0.0
      RHS1=RHO(J,I)*FB3(J,I)+DT*WM( 24)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 24)*(
     1     +A0182+A0202+A0273+A0275+A0277+A0279+A0281+A0283+A0285+A0287
     2     +A0289+A0291+A0300+A0304+A0308+A0312+A0316
     *      )
      PPP=RHO(J,I)-(EB3(J,I)+WB3(J,I)+PB3(J,I)+QB3(J,I))
      FNEW=SOR2*SB3(J,I)+(RHS1-PB3(J,I)*SB3(J+1,I)-QB3(J,I)*SB3(J-1,I)
     1   -EB3(J,I)*SB3(J,I+1)-WB3(J,I)*SB3(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SB3(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SB3(J,I)=FNEW
      A0273=FNEW*A0273
      A0277=FNEW*A0277
      A0281=FNEW*A0281
      A0289=FNEW*A0289
      A0300=FNEW*A0300
      A0304=FNEW*A0304
      A0308=FNEW*A0308
      A0312=FNEW*A0312
      A0316=FNEW*A0316
C----------------    SOLVING  25 TH SPECIES EQUATION    ----------------
      A0188=SA4(J,I)*RC(J,I, 188)
      A0189=SA6(J,I)*RC(J,I, 189)
      RHS1A=+A0187+A0190
      RHS1B=0.0
      RHS1=RHO(J,I)*FB4(J,I)+DT*WM( 25)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 25)*(
     1     +A0188+A0189
     *      )
      PPP=RHO(J,I)-(EB4(J,I)+WB4(J,I)+PB4(J,I)+QB4(J,I))
      FNEW=SOR2*SB4(J,I)+(RHS1-PB4(J,I)*SB4(J+1,I)-QB4(J,I)*SB4(J-1,I)
     1   -EB4(J,I)*SB4(J,I+1)-WB4(J,I)*SB4(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SB4(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SB4(J,I)=FNEW
C----------------    SOLVING  26 TH SPECIES EQUATION    ----------------
      A0212=SA2(J,I)*RC(J,I, 212)
      A0215=SA2(J,I)*RC(J,I, 215)
      A0217=SA3(J,I)*RC(J,I, 217)
      A0219=SA3(J,I)*RC(J,I, 219)
      A0221=SAD(J,I)*RC(J,I, 221)
      A0255=SA4(J,I)*RC(J,I, 255)
      A0274=SA2(J,I)*RC(J,I, 274)
      A0278=SA0(J,I)*RC(J,I, 278)
      A0282=SA5(J,I)*RC(J,I, 282)
      A0398=SAA(J,I)*RC(J,I, 398)
      A0422=SAD(J,I)*SA2(J,I)*RC(J,I, 422)
      RHS1A=+A0211+A0216+A0218+A0220+A0222+A0256+A0273+A0277+A0281+A0397
     1      +A0421
      RHS1B=0.0
      RHS1=RHO(J,I)*FB5(J,I)+DT*WM( 26)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 26)*(
     1     +A0212+A0215+A0217+A0219+A0221+A0255+A0274+A0278+A0282+A0398
     2     +A0422
     *      )
      PPP=RHO(J,I)-(EB5(J,I)+WB5(J,I)+PB5(J,I)+QB5(J,I))
      FNEW=SOR2*SB5(J,I)+(RHS1-PB5(J,I)*SB5(J+1,I)-QB5(J,I)*SB5(J-1,I)
     1   -EB5(J,I)*SB5(J,I+1)-WB5(J,I)*SB5(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SB5(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SB5(J,I)=FNEW
      A0219=FNEW*A0219
      A0255=FNEW*A0255
      A0398=FNEW*A0398
      A0422=FNEW*A0422
C----------------    SOLVING  27 TH SPECIES EQUATION    ----------------
      A0206=SA2(J,I)*RC(J,I, 206)
      A0220=SA4(J,I)*RC(J,I, 220)
      A0223=SA2(J,I)*RC(J,I, 223)
      A0225=SA4(J,I)*RC(J,I, 225)
      A0227=SA3(J,I)*RC(J,I, 227)
      A0229=SA1(J,I)*RC(J,I, 229)
      A0231=SA1(J,I)*RC(J,I, 231)
      A0234=SA2(J,I)*RC(J,I, 234)
      A0238=SA3(J,I)*RC(J,I, 238)
      A0390=SAD(J,I)*RC(J,I, 390)
      RHS1A=+A0205+A0219+A0224+A0226+A0228+A0230+A0232+A0233+A0237+A0389
      RHS1B=0.0
      RHS1=RHO(J,I)*FB6(J,I)+DT*WM( 27)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 27)*(
     1     +A0206+A0220+A0223+A0225+A0227+A0229+A0231+A0234+A0238+A0390
     *      )
      PPP=RHO(J,I)-(EB6(J,I)+WB6(J,I)+PB6(J,I)+QB6(J,I))
      FNEW=SOR2*SB6(J,I)+(RHS1-PB6(J,I)*SB6(J+1,I)-QB6(J,I)*SB6(J-1,I)
     1   -EB6(J,I)*SB6(J,I+1)-WB6(J,I)*SB6(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SB6(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SB6(J,I)=FNEW
      A0234=FNEW*A0234
      A0238=FNEW*A0238
      A0390=FNEW*A0390
C----------------    SOLVING  28 TH SPECIES EQUATION    ----------------
      A0214=SA5(J,I)*RC(J,I, 214)
      A0233=SA4(J,I)*RC(J,I, 233)
      A0235=SA3(J,I)*RC(J,I, 235)
      A0237=SA1(J,I)*RC(J,I, 237)
      A0239=SA1(J,I)*RC(J,I, 239)
      A0241=SA1(J,I)*RC(J,I, 241)
      RHS1A=+A0213+A0234+A0236+A0238+A0240+A0242
      RHS1B=0.0
      RHS1=RHO(J,I)*FB7(J,I)+DT*WM( 28)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 28)*(
     1     +A0214+A0233+A0235+A0237+A0239+A0241
     *      )
      PPP=RHO(J,I)-(EB7(J,I)+WB7(J,I)+PB7(J,I)+QB7(J,I))
      FNEW=SOR2*SB7(J,I)+(RHS1-PB7(J,I)*SB7(J+1,I)-QB7(J,I)*SB7(J-1,I)
     1   -EB7(J,I)*SB7(J,I+1)-WB7(J,I)*SB7(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SB7(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SB7(J,I)=FNEW
C----------------    SOLVING  29 TH SPECIES EQUATION    ----------------
      A0243=SA2(J,I)*RC(J,I, 243)
      A0245=SA2(J,I)*RC(J,I, 245)
      A0247=SA4(J,I)*RC(J,I, 247)
      A0249=SA1(J,I)*RC(J,I, 249)
      A0251=RC(J,I, 251)
      A0254=RC(J,I, 254)
      A0256=SA8(J,I)*RC(J,I, 256)
      A0258=SA5(J,I)*RC(J,I, 258)
      A0262=SA0(J,I)*RC(J,I, 262)
      A0266=SA4(J,I)*RC(J,I, 266)
      A0268=SA7(J,I)*RC(J,I, 268)
      A0270=SA6(J,I)*RC(J,I, 270)
      A0320=SAD(J,I)*RC(J,I, 320)
      A0366=SAD(J,I)*RC(J,I, 366)
      A0378=SAD(J,I)*RC(J,I, 378)
      RHS1A=+A0244+A0246+A0248+A0250+A0252+A0253+A0255+A0257+A0261+A0265
     1      +A0267+A0269+A0319+A0365+A0377
      RHS1B=0.0
      RHS1=RHO(J,I)*FB8(J,I)+DT*WM( 29)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 29)*(
     1     +A0243+A0245+A0247+A0249+A0251+A0254+A0256+A0258+A0262+A0266
     2     +A0268+A0270+A0320+A0366+A0378
     *      )
      PPP=RHO(J,I)-(EB8(J,I)+WB8(J,I)+PB8(J,I)+QB8(J,I))
      FNEW=SOR2*SB8(J,I)+(RHS1-PB8(J,I)*SB8(J+1,I)-QB8(J,I)*SB8(J-1,I)
     1   -EB8(J,I)*SB8(J,I+1)-WB8(J,I)*SB8(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SB8(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SB8(J,I)=FNEW
      A0258=FNEW*A0258
      A0262=FNEW*A0262
      A0266=FNEW*A0266
      A0268=FNEW*A0268
      A0270=FNEW*A0270
      A0320=FNEW*A0320
      A0366=FNEW*A0366
      A0378=FNEW*A0378
C----------------    SOLVING  30 TH SPECIES EQUATION    ----------------
      A0257=SA4(J,I)*RC(J,I, 257)
      A0259=SA4(J,I)*RC(J,I, 259)
      A0261=SA2(J,I)*RC(J,I, 261)
      A0263=SA2(J,I)*RC(J,I, 263)
      A0265=SA3(J,I)*RC(J,I, 265)
      A0267=SA6(J,I)*RC(J,I, 267)
      A0269=SA1(J,I)*RC(J,I, 269)
      A0271=RC(J,I, 271)
      RHS1A=+A0258+A0260+A0262+A0264+A0266+A0268+A0270+A0272
      RHS1B=0.0
      RHS1=RHO(J,I)*FB9(J,I)+DT*WM( 30)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 30)*(
     1     +A0257+A0259+A0261+A0263+A0265+A0267+A0269+A0271
     *      )
      PPP=RHO(J,I)-(EB9(J,I)+WB9(J,I)+PB9(J,I)+QB9(J,I))
      FNEW=SOR2*SB9(J,I)+(RHS1-PB9(J,I)*SB9(J+1,I)-QB9(J,I)*SB9(J-1,I)
     1   -EB9(J,I)*SB9(J,I+1)-WB9(J,I)*SB9(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SB9(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SB9(J,I)=FNEW
C----------------    SOLVING  31 TH SPECIES EQUATION    ----------------
      A0319=RC(J,I, 319)
      A0321=RC(J,I, 321)
      A0323=SA4(J,I)*RC(J,I, 323)
      A0325=SA4(J,I)*RC(J,I, 325)
      A0327=SA4(J,I)*RC(J,I, 327)
      A0329=SA2(J,I)*RC(J,I, 329)
      A0331=SA2(J,I)*RC(J,I, 331)
      A0333=SA2(J,I)*RC(J,I, 333)
      A0335=SA3(J,I)*RC(J,I, 335)
      A0337=SA3(J,I)*RC(J,I, 337)
      A0339=SA3(J,I)*RC(J,I, 339)
      A0341=SAD(J,I)*RC(J,I, 341)
      A0343=SAD(J,I)*RC(J,I, 343)
      A0345=SAD(J,I)*RC(J,I, 345)
      A0347=SA6(J,I)*RC(J,I, 347)
      A0349=SA6(J,I)*RC(J,I, 349)
      A0351=SA6(J,I)*RC(J,I, 351)
      RHS1A=+A0320+A0322+A0324+A0326+A0328+A0330+A0332+A0334+A0336+A0338
     1      +A0340+A0342+A0344+A0346+A0348+A0350+A0352
      RHS1B=0.0
      RHS1=RHO(J,I)*FBA(J,I)+DT*WM( 31)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 31)*(
     1     +A0319+A0321+A0323+A0325+A0327+A0329+A0331+A0333+A0335+A0337
     2     +A0339+A0341+A0343+A0345+A0347+A0349+A0351
     *      )
      PPP=RHO(J,I)-(EBA(J,I)+WBA(J,I)+PBA(J,I)+QBA(J,I))
      FNEW=SOR2*SBA(J,I)+(RHS1-PBA(J,I)*SBA(J+1,I)-QBA(J,I)*SBA(J-1,I)
     1   -EBA(J,I)*SBA(J,I+1)-WBA(J,I)*SBA(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SBA(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SBA(J,I)=FNEW
      A0323=FNEW*A0323
      A0325=FNEW*A0325
      A0327=FNEW*A0327
      A0329=FNEW*A0329
      A0331=FNEW*A0331
      A0333=FNEW*A0333
      A0335=FNEW*A0335
      A0337=FNEW*A0337
      A0339=FNEW*A0339
      A0341=FNEW*A0341
      A0343=FNEW*A0343
      A0345=FNEW*A0345
      A0347=FNEW*A0347
      A0349=FNEW*A0349
      A0351=FNEW*A0351
C----------------    SOLVING  32 TH SPECIES EQUATION    ----------------
      A0290=SA1(J,I)*RC(J,I, 290)
      A0293=RC(J,I, 293)
      A0297=SA4(J,I)*RC(J,I, 297)
      A0299=SA4(J,I)*RC(J,I, 299)
      A0301=SA3(J,I)*RC(J,I, 301)
      A0303=SA3(J,I)*RC(J,I, 303)
      A0305=SA2(J,I)*RC(J,I, 305)
      A0307=SA2(J,I)*RC(J,I, 307)
      A0309=SAD(J,I)*RC(J,I, 309)
      A0311=SAD(J,I)*RC(J,I, 311)
      A0313=SA6(J,I)*RC(J,I, 313)
      A0315=SA6(J,I)*RC(J,I, 315)
      A0317=SA1(J,I)*RC(J,I, 317)
      A0358=SA2(J,I)*RC(J,I, 358)
      A0362=SA6(J,I)*RC(J,I, 362)
      A0370=SA5(J,I)*RC(J,I, 370)
      A0372=SA6(J,I)*RC(J,I, 372)
      A0374=SA4(J,I)*RC(J,I, 374)
      A0380=SA4(J,I)*SA4(J,I)*RC(J,I, 380)
      A0382=SA5(J,I)*RC(J,I, 382)
      A0384=SA2(J,I)*RC(J,I, 384)
      RHS1A=+A0289+A0294+A0298+A0300+A0302+A0304+A0306+A0308+A0310+A0312
     1      +A0314+A0316+A0318+A0357+A0361+A0369+A0371+A0373+A0379+A0381
     2      +A0383
      RHS1B=0.0
      RHS1=RHO(J,I)*FBB(J,I)+DT*WM( 32)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 32)*(
     1     +A0290+A0293+A0297+A0299+A0301+A0303+A0305+A0307+A0309+A0311
     2     +A0313+A0315+A0317+A0358+A0362+A0370+A0372+A0374+A0380+A0382
     3     +A0384
     *      )
      PPP=RHO(J,I)-(EBB(J,I)+WBB(J,I)+PBB(J,I)+QBB(J,I))
      FNEW=SOR2*SBB(J,I)+(RHS1-PBB(J,I)*SBB(J+1,I)-QBB(J,I)*SBB(J-1,I)
     1   -EBB(J,I)*SBB(J,I+1)-WBB(J,I)*SBB(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SBB(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SBB(J,I)=FNEW
      A0297=FNEW*A0297
      A0301=FNEW*A0301
      A0305=FNEW*A0305
      A0309=FNEW*A0309
      A0313=FNEW*A0313
      A0317=FNEW*A0317
      A0358=FNEW*A0358
      A0362=FNEW*A0362
      A0370=FNEW*A0370
      A0372=FNEW*A0372
      A0374=FNEW*A0374
      A0380=FNEW*A0380
      A0382=FNEW*A0382
      A0384=FNEW*A0384
C----------------    SOLVING  33 TH SPECIES EQUATION    ----------------
      A0326=SA5(J,I)*RC(J,I, 326)
      A0332=SA0(J,I)*RC(J,I, 332)
      A0338=SA4(J,I)*RC(J,I, 338)
      A0344=SAC(J,I)*RC(J,I, 344)
      A0348=SA7(J,I)*RC(J,I, 348)
      A0371=SA1(J,I)*RC(J,I, 371)
      A0373=SA3(J,I)*RC(J,I, 373)
      A0375=SA2(J,I)*RC(J,I, 375)
      A0377=SA2(J,I)*RC(J,I, 377)
      A0379=SA6(J,I)*RC(J,I, 379)
      A0381=SA4(J,I)*RC(J,I, 381)
      A0383=RC(J,I, 383)
      RHS1A=+A0325+A0331+A0337+A0343+A0347+A0372+A0374+A0376+A0378+A0380
     1      +A0382+A0384
      RHS1B=0.0
      RHS1=RHO(J,I)*FBC(J,I)+DT*WM( 33)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 33)*(
     1     +A0326+A0332+A0338+A0344+A0348+A0371+A0373+A0375+A0377+A0379
     2     +A0381+A0383
     *      )
      PPP=RHO(J,I)-(EBC(J,I)+WBC(J,I)+PBC(J,I)+QBC(J,I))
      FNEW=SOR2*SBC(J,I)+(RHS1-PBC(J,I)*SBC(J+1,I)-QBC(J,I)*SBC(J-1,I)
     1   -EBC(J,I)*SBC(J,I+1)-WBC(J,I)*SBC(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SBC(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SBC(J,I)=FNEW
C----------------    SOLVING  34 TH SPECIES EQUATION    ----------------
      A0324=SA5(J,I)*RC(J,I, 324)
      A0330=SA0(J,I)*RC(J,I, 330)
      A0336=SA4(J,I)*RC(J,I, 336)
      A0342=SAC(J,I)*RC(J,I, 342)
      A0350=SA7(J,I)*RC(J,I, 350)
      A0354=RC(J,I, 354)
      RHS1A=+A0323+A0329+A0335+A0341+A0349+A0353
      RHS1B=0.0
      RHS1=RHO(J,I)*FBD(J,I)+DT*WM( 34)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 34)*(
     1     +A0324+A0330+A0336+A0342+A0350+A0354
     *      )
      PPP=RHO(J,I)-(EBD(J,I)+WBD(J,I)+PBD(J,I)+QBD(J,I))
      FNEW=SOR2*SBD(J,I)+(RHS1-PBD(J,I)*SBD(J+1,I)-QBD(J,I)*SBD(J-1,I)
     1   -EBD(J,I)*SBD(J,I+1)-WBD(J,I)*SBD(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SBD(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SBD(J,I)=FNEW
C----------------    SOLVING  35 TH SPECIES EQUATION    ----------------
      A0295=RC(J,I, 295)
      A0298=SA5(J,I)*RC(J,I, 298)
      A0302=SA4(J,I)*RC(J,I, 302)
      A0306=SA0(J,I)*RC(J,I, 306)
      A0310=SAC(J,I)*RC(J,I, 310)
      A0314=SA7(J,I)*RC(J,I, 314)
      A0318=SA6(J,I)*RC(J,I, 318)
      RHS1A=+A0296+A0297+A0301+A0305+A0309+A0313+A0317
      RHS1B=0.0
      RHS1=RHO(J,I)*FBE(J,I)+DT*WM( 35)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 35)*(
     1     +A0295+A0298+A0302+A0306+A0310+A0314+A0318
     *      )
      PPP=RHO(J,I)-(EBE(J,I)+WBE(J,I)+PBE(J,I)+QBE(J,I))
      FNEW=SOR2*SBE(J,I)+(RHS1-PBE(J,I)*SBE(J+1,I)-QBE(J,I)*SBE(J-1,I)
     1   -EBE(J,I)*SBE(J,I+1)-WBE(J,I)*SBE(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SBE(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SBE(J,I)=FNEW
C----------------    SOLVING  36 TH SPECIES EQUATION    ----------------
      A0328=SA5(J,I)*RC(J,I, 328)
      A0334=SA0(J,I)*RC(J,I, 334)
      A0340=SA4(J,I)*RC(J,I, 340)
      A0346=SAC(J,I)*RC(J,I, 346)
      A0352=SA7(J,I)*RC(J,I, 352)
      A0356=SA4(J,I)*RC(J,I, 356)
      A0357=RC(J,I, 357)
      A0359=RC(J,I, 359)
      A0361=SA1(J,I)*RC(J,I, 361)
      A0363=SA8(J,I)*RC(J,I, 363)
      A0365=SA2(J,I)*RC(J,I, 365)
      A0367=SA2(J,I)*RC(J,I, 367)
      A0369=SA4(J,I)*RC(J,I, 369)
      RHS1A=+A0327+A0333+A0339+A0345+A0351+A0355+A0358+A0360+A0362+A0364
     1      +A0366+A0368+A0370
      RHS1B=0.0
      RHS1=RHO(J,I)*FBF(J,I)+DT*WM( 36)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 36)*(
     1     +A0328+A0334+A0340+A0346+A0352+A0356+A0357+A0359+A0361+A0363
     2     +A0365+A0367+A0369
     *      )
      PPP=RHO(J,I)-(EBF(J,I)+WBF(J,I)+PBF(J,I)+QBF(J,I))
      FNEW=SOR2*SBF(J,I)+(RHS1-PBF(J,I)*SBF(J+1,I)-QBF(J,I)*SBF(J-1,I)
     1   -EBF(J,I)*SBF(J,I+1)-WBF(J,I)*SBF(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SBF(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SBF(J,I)=FNEW
C----------------    SOLVING  37 TH SPECIES EQUATION    ----------------
      A0385=SA3(J,I)*RC(J,I, 385)
      A0388=SA2(J,I)*RC(J,I, 388)
      A0389=SA3(J,I)*RC(J,I, 389)
      A0392=RC(J,I, 392)
      A0394=SA1(J,I)*RC(J,I, 394)
      A0395=SA4(J,I)*RC(J,I, 395)
      A0399=SA2(J,I)*RC(J,I, 399)
      A0402=SA0(J,I)*RC(J,I, 402)
      A0404=SA6(J,I)*RC(J,I, 404)
      A0406=SAC(J,I)*RC(J,I, 406)
      A0410=SA5(J,I)*RC(J,I, 410)
      A0412=SA8(J,I)*RC(J,I, 412)
      A0415=SA1(J,I)*RC(J,I, 415)
      A0490=SAG(J,I)*RC(J,I, 490)
      A0496=SB2(J,I)*SA6(J,I)*RC(J,I, 496)
      A0500=SB2(J,I)*SA7(J,I)*RC(J,I, 500)
      RHS1A=+A0386+A0387+A0390+A0391+A0393+A0396+A0400+A0401+A0403+A0405
     1      +A0409+A0411+A0416+A0489+A0495+A0499
      RHS1B=0.0
      RHS1=RHO(J,I)*FBG(J,I)+DT*WM( 37)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 37)*(
     1     +A0385+A0388+A0389+A0392+A0394+A0395+A0399+A0402+A0404+A0406
     2     +A0410+A0412+A0415+A0490+A0496+A0500
     *      )
      PPP=RHO(J,I)-(EBG(J,I)+WBG(J,I)+PBG(J,I)+QBG(J,I))
      FNEW=SOR2*SBG(J,I)+(RHS1-PBG(J,I)*SBG(J+1,I)-QBG(J,I)*SBG(J-1,I)
     1   -EBG(J,I)*SBG(J,I+1)-WBG(J,I)*SBG(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SBG(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SBG(J,I)=FNEW
      A0392=FNEW*A0392
      A0394=FNEW*A0394
      A0395=FNEW*A0395
      A0399=FNEW*A0399
      A0402=FNEW*A0402
      A0404=FNEW*A0404
      A0406=FNEW*A0406
      A0410=FNEW*A0410
      A0412=FNEW*A0412
      A0490=FNEW*A0490
      A0496=FNEW*A0496
      A0500=FNEW*A0500
C----------------    SOLVING  38 TH SPECIES EQUATION    ----------------
      A0391=SA2(J,I)*RC(J,I, 391)
      A0393=SA6(J,I)*RC(J,I, 393)
      A0396=SA5(J,I)*RC(J,I, 396)
      A0397=SA1(J,I)*RC(J,I, 397)
      A0411=SAA(J,I)*RC(J,I, 411)
      A0413=SA6(J,I)*RC(J,I, 413)
      A0485=SAD(J,I)*RC(J,I, 485)
      RHS1A=+A0392+A0394+A0395+A0398+A0412+A0414+A0486
      RHS1B=0.0
      RHS1=RHO(J,I)*FBH(J,I)+DT*WM( 38)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 38)*(
     1     +A0391+A0393+A0396+A0397+A0411+A0413+A0485
     *      )
      PPP=RHO(J,I)-(EBH(J,I)+WBH(J,I)+PBH(J,I)+QBH(J,I))
      FNEW=SOR2*SBH(J,I)+(RHS1-PBH(J,I)*SBH(J+1,I)-QBH(J,I)*SBH(J-1,I)
     1   -EBH(J,I)*SBH(J,I+1)-WBH(J,I)*SBH(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SBH(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SBH(J,I)=FNEW
      A0485=FNEW*A0485
C----------------    SOLVING  39 TH SPECIES EQUATION    ----------------
      A0400=RC(J,I, 400)
      A0401=SA2(J,I)*RC(J,I, 401)
      A0403=SA1(J,I)*RC(J,I, 403)
      A0405=SAD(J,I)*RC(J,I, 405)
      A0408=RC(J,I, 408)
      A0409=SA4(J,I)*RC(J,I, 409)
      A0420=SA5(J,I)*RC(J,I, 420)
      A0424=SA0(J,I)*RC(J,I, 424)
      A0425=SA2(J,I)*RC(J,I, 425)
      A0427=SA6(J,I)*RC(J,I, 427)
      A0429=SA6(J,I)*RC(J,I, 429)
      A0436=SA2(J,I)*RC(J,I, 436)
      A0482=SC4(J,I)*SA2(J,I)*RC(J,I, 482)
      A0492=SB2(J,I)*RC(J,I, 492)
      A0494=SB1(J,I)*SA6(J,I)*RC(J,I, 494)
      A0498=SB1(J,I)*SA7(J,I)*RC(J,I, 498)
      A0520=SAI(J,I)*RC(J,I, 520)
      A0532=SA0(J,I)*SAG(J,I)*RC(J,I, 532)
      A0536=SAD(J,I)*RC(J,I, 536)
      RHS1A=+A0399+A0402+A0404+A0406+A0407+A0410+A0419+A0423+A0426+A0428
     1      +A0430+A0435+A0481+A0491+A0493+A0497+A0519+A0531+A0535
      RHS1B=0.0
      RHS1=RHO(J,I)*FBI(J,I)+DT*WM( 39)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 39)*(
     1     +A0400+A0401+A0403+A0405+A0408+A0409+A0420+A0424+A0425+A0427
     2     +A0429+A0436+A0482+A0492+A0494+A0498+A0520+A0532+A0536
     *      )
      PPP=RHO(J,I)-(EBI(J,I)+WBI(J,I)+PBI(J,I)+QBI(J,I))
      FNEW=SOR2*SBI(J,I)+(RHS1-PBI(J,I)*SBI(J+1,I)-QBI(J,I)*SBI(J-1,I)
     1   -EBI(J,I)*SBI(J,I+1)-WBI(J,I)*SBI(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SBI(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SBI(J,I)=FNEW
      A0420=FNEW*A0420
      A0424=FNEW*A0424
      A0425=FNEW*A0425
      A0427=FNEW*A0427
      A0482=FNEW*A0482
      A0492=FNEW*A0492
      A0494=FNEW*A0494
      A0498=FNEW*A0498
      A0520=FNEW*A0520
      A0532=FNEW*A0532
      A0536=FNEW*A0536
C----------------    SOLVING  40 TH SPECIES EQUATION    ----------------
      A0417=SA3(J,I)*RC(J,I, 417)
      A0419=SA4(J,I)*RC(J,I, 419)
      A0421=SA3(J,I)*RC(J,I, 421)
      A0423=SA2(J,I)*RC(J,I, 423)
      A0426=RC(J,I, 426)
      A0428=SA1(J,I)*RC(J,I, 428)
      A0432=RC(J,I, 432)
      A0433=SA2(J,I)*RC(J,I, 433)
      A0461=SA2(J,I)*RC(J,I, 461)
      A0464=SA6(J,I)*RC(J,I, 464)
      A0467=SA2(J,I)*RC(J,I, 467)
      A0470=SA6(J,I)*RC(J,I, 470)
      A0488=SB1(J,I)*RC(J,I, 488)
      A0522=SAG(J,I)*RC(J,I, 522)
      A0524=SA5(J,I)*SB2(J,I)*RC(J,I, 524)
      A0526=SA0(J,I)*SB2(J,I)*RC(J,I, 526)
      A0530=SAI(J,I)*RC(J,I, 530)
      A0540=SAD(J,I)*RC(J,I, 540)
      RHS1A=+A0418+A0420+A0422+A0424+A0425+A0427+A0431+A0434+A0462+A0463
     1      +A0468+A0469+A0487+A0521+A0523+A0525+A0529+A0539
      RHS1B=+1.000*A0503+1.000*A0515
      RHS1=RHO(J,I)*FBJ(J,I)+DT*WM( 40)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 40)*(
     1     +A0417+A0419+A0421+A0423+A0426+A0428+A0432+A0433+A0461+A0464
     2     +A0467+A0470+A0488+A0522+A0524+A0526+A0530+A0540
     *      )
      PPP=RHO(J,I)-(EBJ(J,I)+WBJ(J,I)+PBJ(J,I)+QBJ(J,I))
      FNEW=SOR2*SBJ(J,I)+(RHS1-PBJ(J,I)*SBJ(J+1,I)-QBJ(J,I)*SBJ(J-1,I)
     1   -EBJ(J,I)*SBJ(J,I+1)-WBJ(J,I)*SBJ(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SBJ(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SBJ(J,I)=FNEW
      A0461=FNEW*A0461
      A0464=FNEW*A0464
      A0467=FNEW*A0467
      A0470=FNEW*A0470
      A0488=FNEW*A0488
      A0522=FNEW*A0522
      A0524=FNEW*A0524
      A0526=FNEW*A0526
      A0530=FNEW*A0530
      A0540=FNEW*A0540
C----------------    SOLVING  41 ST SPECIES EQUATION    ----------------
      A0437=RC(J,I, 437)
      A0439=SA1(J,I)*RC(J,I, 439)
      A0441=SA1(J,I)*RC(J,I, 441)
      A0443=SA2(J,I)*RC(J,I, 443)
      A0445=SA2(J,I)*RC(J,I, 445)
      A0447=SA3(J,I)*RC(J,I, 447)
      A0449=SA3(J,I)*RC(J,I, 449)
      A0451=SA4(J,I)*RC(J,I, 451)
      A0453=SA4(J,I)*RC(J,I, 453)
      A0455=SA6(J,I)*RC(J,I, 455)
      A0457=SA6(J,I)*RC(J,I, 457)
      RHS1A=+A0438+A0440+A0442+A0444+A0446+A0448+A0450+A0452+A0454+A0456
     1      +A0458
      RHS1B=0.0
      RHS1=RHO(J,I)*FC0(J,I)+DT*WM( 41)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 41)*(
     1     +A0437+A0439+A0441+A0443+A0445+A0447+A0449+A0451+A0453+A0455
     2     +A0457
     *      )
      PPP=RHO(J,I)-(EC0(J,I)+WC0(J,I)+PC0(J,I)+QC0(J,I))
      FNEW=SOR2*SC0(J,I)+(RHS1-PC0(J,I)*SC0(J+1,I)-QC0(J,I)*SC0(J-1,I)
     1   -EC0(J,I)*SC0(J,I+1)-WC0(J,I)*SC0(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SC0(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SC0(J,I)=FNEW
      A0439=FNEW*A0439
      A0441=FNEW*A0441
      A0443=FNEW*A0443
      A0445=FNEW*A0445
      A0447=FNEW*A0447
      A0449=FNEW*A0449
      A0451=FNEW*A0451
      A0453=FNEW*A0453
      A0455=FNEW*A0455
      A0457=FNEW*A0457
C----------------    SOLVING  42 ND SPECIES EQUATION    ----------------
      A0440=SA6(J,I)*RC(J,I, 440)
      A0444=SA0(J,I)*RC(J,I, 444)
      A0448=SA4(J,I)*RC(J,I, 448)
      A0454=SA5(J,I)*RC(J,I, 454)
      A0456=SA7(J,I)*RC(J,I, 456)
      A0459=RC(J,I, 459)
      A0462=RC(J,I, 462)
      A0463=SA1(J,I)*RC(J,I, 463)
      RHS1A=+A0439+A0443+A0447+A0453+A0455+A0460+A0461+A0464
      RHS1B=0.0
      RHS1=RHO(J,I)*FC1(J,I)+DT*WM( 42)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 42)*(
     1     +A0440+A0444+A0448+A0454+A0456+A0459+A0462+A0463
     *      )
      PPP=RHO(J,I)-(EC1(J,I)+WC1(J,I)+PC1(J,I)+QC1(J,I))
      FNEW=SOR2*SC1(J,I)+(RHS1-PC1(J,I)*SC1(J+1,I)-QC1(J,I)*SC1(J-1,I)
     1   -EC1(J,I)*SC1(J,I+1)-WC1(J,I)*SC1(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SC1(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SC1(J,I)=FNEW
      A0459=FNEW*A0459
C----------------    SOLVING  43 RD SPECIES EQUATION    ----------------
      A0442=SA6(J,I)*RC(J,I, 442)
      A0446=SA0(J,I)*RC(J,I, 446)
      A0450=SA4(J,I)*RC(J,I, 450)
      A0452=SA5(J,I)*RC(J,I, 452)
      A0458=SA7(J,I)*RC(J,I, 458)
      A0460=RC(J,I, 460)
      A0465=RC(J,I, 465)
      A0468=RC(J,I, 468)
      A0469=SA1(J,I)*RC(J,I, 469)
      RHS1A=+A0441+A0445+A0449+A0451+A0457+A0459+A0466+A0467+A0470
      RHS1B=0.0
      RHS1=RHO(J,I)*FC2(J,I)+DT*WM( 43)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 43)*(
     1     +A0442+A0446+A0450+A0452+A0458+A0460+A0465+A0468+A0469
     *      )
      PPP=RHO(J,I)-(EC2(J,I)+WC2(J,I)+PC2(J,I)+QC2(J,I))
      FNEW=SOR2*SC2(J,I)+(RHS1-PC2(J,I)*SC2(J+1,I)-QC2(J,I)*SC2(J-1,I)
     1   -EC2(J,I)*SC2(J,I+1)-WC2(J,I)*SC2(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SC2(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SC2(J,I)=FNEW
C----------------    SOLVING  44 TH SPECIES EQUATION    ----------------
      A0471=RC(J,I, 471)
      A0473=RC(J,I, 473)
      A0476=RC(J,I, 476)
      A0477=SA2(J,I)*RC(J,I, 477)
      A0479=SA2(J,I)*RC(J,I, 479)
      A0481=SA4(J,I)*RC(J,I, 481)
      A0483=SAD(J,I)*RC(J,I, 483)
      A0486=RC(J,I, 486)
      A0534=SA0(J,I)*SAD(J,I)*RC(J,I, 534)
      A0544=SA5(J,I)*SA2(J,I)*RC(J,I, 544)
      RHS1A=+A0472+A0474+A0475+A0478+A0480+A0482+A0484+A0485+A0533+A0543
      RHS1B=0.0
      RHS1=RHO(J,I)*FC3(J,I)+DT*WM( 44)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 44)*(
     1     +A0471+A0473+A0476+A0477+A0479+A0481+A0483+A0486+A0534+A0544
     *      )
      PPP=RHO(J,I)-(EC3(J,I)+WC3(J,I)+PC3(J,I)+QC3(J,I))
      FNEW=SOR2*SC3(J,I)+(RHS1-PC3(J,I)*SC3(J+1,I)-QC3(J,I)*SC3(J-1,I)
     1   -EC3(J,I)*SC3(J,I+1)-WC3(J,I)*SC3(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SC3(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SC3(J,I)=FNEW
      A0481=FNEW*A0481
      A0534=FNEW*A0534
      A0544=FNEW*A0544
C----------------    SOLVING  45 TH SPECIES EQUATION    ----------------
      A0482=SA2(J,I)*SBI(J,I)*RC(J,I, 482)
      RHS1A=+A0481
      RHS1B=0.0
      RHS1=RHO(J,I)*FC4(J,I)+DT*WM( 45)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 45)*(
     1     +A0482
     *      )
      PPP=RHO(J,I)-(EC4(J,I)+WC4(J,I)+PC4(J,I)+QC4(J,I))
      FNEW=SOR2*SC4(J,I)+(RHS1-PC4(J,I)*SC4(J+1,I)-QC4(J,I)*SC4(J-1,I)
     1   -EC4(J,I)*SC4(J,I+1)-WC4(J,I)*SC4(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SC4(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SC4(J,I)=FNEW
C----------------    SOLVING  46 TH SPECIES EQUATION    ----------------
      A0487=RC(J,I, 487)
      A0489=RC(J,I, 489)
      A0491=RC(J,I, 491)
      A0493=SA1(J,I)*RC(J,I, 493)
      A0495=SA1(J,I)*RC(J,I, 495)
      A0497=SA6(J,I)*RC(J,I, 497)
      A0499=SA6(J,I)*RC(J,I, 499)
      RHS1A=+A0488+A0490+A0492+A0494+A0496+A0498+A0500
      RHS1B=0.0
      RHS1=RHO(J,I)*FC5(J,I)+DT*WM( 46)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 46)*(
     1     +A0487+A0489+A0491+A0493+A0495+A0497+A0499
     *      )
      PPP=RHO(J,I)-(EC5(J,I)+WC5(J,I)+PC5(J,I)+QC5(J,I))
      FNEW=SOR2*SC5(J,I)+(RHS1-PC5(J,I)*SC5(J+1,I)-QC5(J,I)*SC5(J-1,I)
     1   -EC5(J,I)*SC5(J,I+1)-WC5(J,I)*SC5(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SC5(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SC5(J,I)=FNEW
C----------------    SOLVING  47 TH SPECIES EQUATION    ----------------
      A0501=RC(J,I, 501)
      A0503=SA2(J,I)*RC(J,I, 503)
      A0505=SA2(J,I)*RC(J,I, 505)
      A0507=SA2(J,I)*RC(J,I, 507)
      A0509=SA2(J,I)*RC(J,I, 509)
      A0511=SA4(J,I)*RC(J,I, 511)
      A0513=SA4(J,I)*RC(J,I, 513)
      A0515=SA4(J,I)*RC(J,I, 515)
      A0517=SA4(J,I)*RC(J,I, 517)
      RHS1A=+A0502+A0504+A0506+A0508+A0510+A0512+A0514+A0516+A0518
      RHS1B=0.0
      RHS1=RHO(J,I)*FC6(J,I)+DT*WM( 47)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 47)*(
     1     +A0501+A0503+A0505+A0507+A0509+A0511+A0513+A0515+A0517
     *      )
      PPP=RHO(J,I)-(EC6(J,I)+WC6(J,I)+PC6(J,I)+QC6(J,I))
      FNEW=SOR2*SC6(J,I)+(RHS1-PC6(J,I)*SC6(J+1,I)-QC6(J,I)*SC6(J-1,I)
     1   -EC6(J,I)*SC6(J,I+1)-WC6(J,I)*SC6(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SC6(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SC6(J,I)=FNEW
      A0501=FNEW*A0501
      A0503=FNEW*A0503
      A0505=FNEW*A0505
      A0507=FNEW*A0507
      A0509=FNEW*A0509
      A0511=FNEW*A0511
      A0513=FNEW*A0513
      A0515=FNEW*A0515
      A0517=FNEW*A0517
C----------------    SOLVING  48 TH SPECIES EQUATION    ----------------
      A0535=RC(J,I, 535)
      A0537=SA2(J,I)*RC(J,I, 537)
      A0539=SA2(J,I)*RC(J,I, 539)
      A0541=SA2(J,I)*RC(J,I, 541)
      A0543=SA4(J,I)*RC(J,I, 543)
      RHS1A=+A0536+A0538+A0540+A0542+A0544
      RHS1B=+1.000*A0505+1.000*A0513
      RHS1=RHO(J,I)*FC7(J,I)+DT*WM( 48)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 48)*(
     1     +A0535+A0537+A0539+A0541+A0543
     *      )
      PPP=RHO(J,I)-(EC7(J,I)+WC7(J,I)+PC7(J,I)+QC7(J,I))
      FNEW=SOR2*SC7(J,I)+(RHS1-PC7(J,I)*SC7(J+1,I)-QC7(J,I)*SC7(J-1,I)
     1   -EC7(J,I)*SC7(J,I+1)-WC7(J,I)*SC7(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SC7(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SC7(J,I)=FNEW
C----------------    SOLVING  49 TH SPECIES EQUATION    ----------------
      A0510=SA0(J,I)*SAI(J,I)*RC(J,I, 510)
      A0512=SA5(J,I)*SAI(J,I)*RC(J,I, 512)
      A0519=RC(J,I, 519)
      A0521=RC(J,I, 521)
      A0523=SA4(J,I)*RC(J,I, 523)
      A0525=SA2(J,I)*RC(J,I, 525)
      A0527=SA2(J,I)*RC(J,I, 527)
      A0529=SA2(J,I)*RC(J,I, 529)
      A0531=SA2(J,I)*RC(J,I, 531)
      A0533=SA2(J,I)*RC(J,I, 533)
      RHS1A=+A0509+A0511+A0520+A0522+A0524+A0526+A0528+A0530+A0532+A0534
      RHS1B=0.0
      RHS1=RHO(J,I)*FC8(J,I)+DT*WM( 49)*(RHS1A+RHS1B)
      RHS2=-DT*WM( 49)*(
     1     +A0510+A0512+A0519+A0521+A0523+A0525+A0527+A0529+A0531+A0533
     *      )
      PPP=RHO(J,I)-(EC8(J,I)+WC8(J,I)+PC8(J,I)+QC8(J,I))
      FNEW=SOR2*SC8(J,I)+(RHS1-PC8(J,I)*SC8(J+1,I)-QC8(J,I)*SC8(J-1,I)
     1   -EC8(J,I)*SC8(J,I+1)-WC8(J,I)*SC8(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SC8(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SC8(J,I)=FNEW
C----------------    SOLVING  50 TH SPECIES EQUATION    ----------------
      RHS1B=0.0
      RHS1=RHO(J,I)*FC9(J,I)+DT*WM( 50)*(RHS1A+RHS1B)
      RHS2=0.0
      PPP=RHO(J,I)-(EC9(J,I)+WC9(J,I)+PC9(J,I)+QC9(J,I))
      FNEW=SOR2*SC9(J,I)+(RHS1-PC9(J,I)*SC9(J+1,I)-QC9(J,I)*SC9(J-1,I)
     1   -EC9(J,I)*SC9(J,I+1)-WC9(J,I)*SC9(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SC9(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SC9(J,I)=FNEW
C----------------    SOLVING  51 TH SPECIES EQUATION    ----------------
      RHS1B=0.0
      RHS1=RHO(J,I)*FCA(J,I)+DT*WM( 51)*(RHS1A+RHS1B)
      RHS2=0.0
      PPP=RHO(J,I)-(ECA(J,I)+WCA(J,I)+PCA(J,I)+QCA(J,I))
      FNEW=SOR2*SCA(J,I)+(RHS1-PCA(J,I)*SCA(J+1,I)-QCA(J,I)*SCA(J-1,I)
     1   -ECA(J,I)*SCA(J,I+1)-WCA(J,I)*SCA(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SCA(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SCA(J,I)=FNEW
C---------------    SOLVING  SOOT MASS FRACTION EQUATION    ------------
C------------ Molecular Weight of Soot is equal to that of C -----------
      RHS1=RHO(J,I)*STMF(J,I)
     1    +DT*(0.012011/WMSTR)*(2.0*SOUR(J,I,1)*SB1(J,I)
     2    -SOUR(J,I,3)*SA1(J,I)-SOUR(J,I,4)*SA4(J,I)
     3    -SOUR(J,I,5)*SA3(J,I))
      RHS2=0.0
      PPP=RHO(J,I)-(STE(J,I)+STW(J,I)+STN(J,I)+STS(J,I))
      FNEW=SOR2*SST1(J,I)+(RHS1-STN(J,I)*SST1(J+1,I)
     1   -STS(J,I)*SST1(J-1,I)
     2   -STE(J,I)*SST1(J,I+1)-STW(J,I)*SST1(J,I-1))*SOR1/(PPP-RHS2)
      RSD=RSD+DABS(SST1(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SST1(J,I)=FNEW
C---------------    SOLVING  SOOT NUMBER DENSITY EQUATION    -----------
      RHS1=RHO(J,I)*STND(J,I)
     1    +DT*SOUR(J,I,2)*SB1(J,I)
      RHS2=-DT*SOUR(J,I,6)
      PPP=RHO(J,I)-(STE(J,I)+STW(J,I)+STN(J,I)+STS(J,I))
      FNEW=SOR2*SST2(J,I)+(RHS1-STN(J,I)*SST2(J+1,I)
     1   -STS(J,I)*SST2(J-1,I)
     2   -STE(J,I)*SST2(J,I+1)-STW(J,I)*SST2(J,I-1))*SOR1/(PPP-RHS2)
      RSDY=RSDY+DABS(SST2(J,I)-FNEW)
      IF(FNEW.LT.0.0) FNEW=0.0
      SST2(J,I)=FNEW
  305 CONTINUE
  304 CONTINUE
      IF((RSD/SPNORM).LE.TOLRSP) GO TO 301
  300 CONTINUE
      ISOR1=ISORP-1
  301 CONTINUE
      ISOR1=ISORP
C------------------------   N2  N2  N2   -------------------------------
      DO 394 J=2,LJ-1
      DO 390 I=2,LI-1
      IF(ISKIP(J,I).GE.3) GO TO 390
C--------------------------SOLVING FOR N2------------------------------
      FTOT=+SA0(J,I)+SA1(J,I)+SA2(J,I)+SA3(J,I)+SA4(J,I)
     *     +SA5(J,I)+SA6(J,I)+SA7(J,I)+SA8(J,I)+SA9(J,I)
     *     +SAA(J,I)+SAB(J,I)+SAC(J,I)+SAD(J,I)+SAE(J,I)
     *     +SAF(J,I)+SAG(J,I)+SAH(J,I)+SAI(J,I)+SAJ(J,I)
     *     +SB0(J,I)+SB1(J,I)+SB2(J,I)+SB3(J,I)+SB4(J,I)
     *     +SB5(J,I)+SB6(J,I)+SB7(J,I)+SB8(J,I)+SB9(J,I)
     *     +SBA(J,I)+SBB(J,I)+SBC(J,I)+SBD(J,I)+SBE(J,I)
     *     +SBF(J,I)+SBG(J,I)+SBH(J,I)+SBI(J,I)+SBJ(J,I)
     *     +SC0(J,I)+SC1(J,I)+SC2(J,I)+SC3(J,I)+SC4(J,I)
     *     +SC5(J,I)+SC6(J,I)+SC7(J,I)+SC8(J,I)+SC9(J,I)
     *     +SCA(J,I)
     *     +SST1(J,I)
      SCB(J,I)=1.0-FTOT
      IF(FTOT.GT.1.0) THEN
            SA0(J,I)=SA0(J,I)/FTOT
            SA1(J,I)=SA1(J,I)/FTOT
            SA2(J,I)=SA2(J,I)/FTOT
            SA3(J,I)=SA3(J,I)/FTOT
            SA4(J,I)=SA4(J,I)/FTOT
            SA5(J,I)=SA5(J,I)/FTOT
            SA6(J,I)=SA6(J,I)/FTOT
            SA7(J,I)=SA7(J,I)/FTOT
            SA8(J,I)=SA8(J,I)/FTOT
            SA9(J,I)=SA9(J,I)/FTOT
            SAA(J,I)=SAA(J,I)/FTOT
            SAB(J,I)=SAB(J,I)/FTOT
            SAC(J,I)=SAC(J,I)/FTOT
            SAD(J,I)=SAD(J,I)/FTOT
            SAE(J,I)=SAE(J,I)/FTOT
            SAF(J,I)=SAF(J,I)/FTOT
            SAG(J,I)=SAG(J,I)/FTOT
            SAH(J,I)=SAH(J,I)/FTOT
            SAI(J,I)=SAI(J,I)/FTOT
            SAJ(J,I)=SAJ(J,I)/FTOT
            SB0(J,I)=SB0(J,I)/FTOT
            SB1(J,I)=SB1(J,I)/FTOT
            SB2(J,I)=SB2(J,I)/FTOT
            SB3(J,I)=SB3(J,I)/FTOT
            SB4(J,I)=SB4(J,I)/FTOT
            SB5(J,I)=SB5(J,I)/FTOT
            SB6(J,I)=SB6(J,I)/FTOT
            SB7(J,I)=SB7(J,I)/FTOT
            SB8(J,I)=SB8(J,I)/FTOT
            SB9(J,I)=SB9(J,I)/FTOT
            SBA(J,I)=SBA(J,I)/FTOT
            SBB(J,I)=SBB(J,I)/FTOT
            SBC(J,I)=SBC(J,I)/FTOT
            SBD(J,I)=SBD(J,I)/FTOT
            SBE(J,I)=SBE(J,I)/FTOT
            SBF(J,I)=SBF(J,I)/FTOT
            SBG(J,I)=SBG(J,I)/FTOT
            SBH(J,I)=SBH(J,I)/FTOT
            SBI(J,I)=SBI(J,I)/FTOT
            SBJ(J,I)=SBJ(J,I)/FTOT
            SC0(J,I)=SC0(J,I)/FTOT
            SC1(J,I)=SC1(J,I)/FTOT
            SC2(J,I)=SC2(J,I)/FTOT
            SC3(J,I)=SC3(J,I)/FTOT
            SC4(J,I)=SC4(J,I)/FTOT
            SC5(J,I)=SC5(J,I)/FTOT
            SC6(J,I)=SC6(J,I)/FTOT
            SC7(J,I)=SC7(J,I)/FTOT
            SC8(J,I)=SC8(J,I)/FTOT
            SC9(J,I)=SC9(J,I)/FTOT
            SCA(J,I)=SCA(J,I)/FTOT
            SST1(J,I)=SST1(J,I)/FTOT
            SCB(J,I)=0.0
            END IF
      FCB(J,I)=SCB(J,I)
  390 CONTINUE
  394 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE SOLVEH(ISOR2,RELXH,TOLRH,SIGH,SIGSP,HNORM)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LJ2=LJ*2,LSP=52,LRX=544,
     1          LPD=LE*(149-LSP-LSP-1))
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),
     1  W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDF(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/CB08/ TPOL1,TPOL2,POLSP(14,LSP),CISP(12,LSP)
      COMMON/CB09/ IBOT(LI),ITOP(LI),JLFT(LJ),JRGT(LJ),
     1       FBOT(8+LSP,LI),FTOP(8+LSP,LI),FLFT(8+LSP,LJ),FRGT(8+LSP,LJ)
      COMMON/CB10/ FBXP(4,LJ),FBXM(4,LJ),FBYP(4,LI),FBYM(4,LI)
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
      COMMON/REAC/ LREV,NRLIND,ITBEF(LRX),ALF(LRX),AKF(LRX),EXAF(LRX),
     1             ALLOW(LRX),AKLOW(LRX),EALOW(LRX),TROE(LRX,4)
      COMMON/HFORM/ HFO(LSP)
      COMMON/HRLS/ QDOT(LJ,LI)
      COMMON/DUM0/ RDHSP(LJ,LI,LSP),RHSZ(LJ,LI),HENTH(LJ,LI,LSP),
     1             DUM(LPD)
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
      COMMON/DUMD/ QDSP(LJ,LI,LSP)
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
      COMMON/DUMX/ ZE(LJ,LI),ZW(LJ,LI),ZN(LJ,LI),ZS(LJ,LI)
      COMMON/DUMY/ SSP(LE,LSP),FPZ(LJ,LI),EMU(LJ,LI),EKT(LJ,LI)
      COMMON/DUMR/ RC(LE,LRX)
      COMMON/DUMS/ STE(LJ,LI),STW(LJ,LI),STN(LJ,LI),STS(LJ,LI),
     1  SST1(LJ,LI),SST2(LJ,LI),SOUR(LE,6)
      DIMENSION RA(LRX),PRD(LSP),SRA(10)
      SYM=DFLOAT(ISYM)
      HCONVR=GASC/WMSTR/HSTR
      CRAD=ALSTR/RSTR/VSTR/VSTR/VSTR
      RCONV=VSTR/WMSTR/ALSTR*1.0D-06
      TKDMIN=TSTR-1.0
      DO 402 I=1,LI
      DO J=1,LJ
      TKD=TK(J,I)*TSTR
      TKD2=TKD*TKD
      TKD3=TKD*TKD2
      TKD4=TKD*TKD3
      TKD5=TKD*TKD4
      IPOLY=1
      IF(TKD.GT.TPOL1) IPOLY=8
      DO 404 ISP=1,LSP
      HSSP=(POLSP(IPOLY,ISP)*TKD+POLSP(IPOLY+1,ISP)*TKD2/2.0
     1     +POLSP(IPOLY+2,ISP)*TKD3/3.0+POLSP(IPOLY+3,ISP)*TKD4/4.0
     2     +POLSP(IPOLY+4,ISP)*TKD5/5.0)/WM(ISP)*HCONVR
      RDHSP(J,I,ISP)=(EMU(J,I)*VDF(J,I,ISP)+REI*TMU(J,I)/SIGSP
     1               -EKT(J,I))*HSSP
      HENTH(J,I,ISP)=HSSP*WMSTR*HSTR+POLSP(IPOLY+5,ISP)/WM(ISP)*GASC
  404 CONTINUE
      ENDDO
  402 CONTINUE
C-----------------   SOURCE TERM IN ENTHALPY EQUATION   
      DO 410 J=2,LJ-1
      DO I=2,LI-1
      IF(ISKIP(J,I).GE.10) GO TO 411
      IA=(I-1)*LJ+J
      RDHE=0.0
      RDHW=0.0
      RDHN=0.0
      RDHS=0.0
      RDHP=0.0
      TMOLE=0.0
      DO 412 ISP=1,LSP
      TMOLE=TMOLE+SSP(IA,ISP)/WM(ISP)
      RDHE=RDHE+(RDHSP(J,I+1,ISP)+RDHSP(J,I,ISP))
     1   *(SSP(IA+LJ,ISP)-SSP(IA,ISP))
      RDHW=RDHW+(RDHSP(J,I,ISP)+RDHSP(J,I-1,ISP))
     1   *(SSP(IA,ISP)-SSP(IA-LJ,ISP))
      RDHN=RDHN+(RDHSP(J+1,I,ISP)+RDHSP(J,I,ISP))
     1   *(SSP(IA+1,ISP)-SSP(IA,ISP))
      RDHS=RDHS+(RDHSP(J,I,ISP)+RDHSP(J-1,I,ISP))
     1   *(SSP(IA,ISP)-SSP(IA-1,ISP))
      RDHP=RDHP+(RDHSP(J,I,ISP))
     1   *(SSP(IA+1,ISP)-SSP(IA-1,ISP))
  412 CONTINUE
C---------------------------- RADIATION --------------------------------
C------------------------CH4, CO2, CO & H2O & SOOT----------------------
      TKD=TK(J,I)*TSTR
      TKDI=1000.0/TKD
      TKDIS=TKDI*TKDI
      PLANK1=6.6334+TKD*(-3.5686D-03+TKD*(1.6682D-08+TKD*(2.5611D-10
     1       -2.6558D-14*TKD)))
      PLANK2=1.8741D+01-TKDI*(1.2131D+02-2.7350D+02*TKDI
     1  +1.9405D+02*TKDIS-5.6310D+01*TKDIS*TKDI+5.81690*TKDIS*TKDIS)
      IF(TKD.LE.750.0) THEN
      PLANK3=4.7869+TKD*(-6.953D-02+TKD*(2.95775D-04+TKD*(-4.25732D-07
     1       +2.02894D-10*TKD)))
             ELSE
      PLANK3=1.009D+01+TKD*(-1.183D-02+TKD*(4.7753D-06+TKD*(-5.87209D-10
     1       -2.5334D-14*TKD)))
             END IF
      PLANK4=-0.23093-TKDI*(1.12390-9.41530*TKDI+2.99880*TKDIS
     1       -0.51382*TKDIS*TKDI+1.86840D-05*TKDIS*TKDIS)
C-----------------CH4, CO2, CO, and H2O
      PLANKC=(PLANK1*SSP(IA,13)/WM(13)+PLANK2*SSP(IA,10)/WM(10)
     1       +PLANK3*SSP(IA,09)/WM(09)+PLANK4*SSP(IA,06)/WM(06))/TMOLE
C
C----*****----ONLY 40% OF THE CALCULATED RADIATION INCLUDED----*****----
C
      QRAD=CRAD*0.4*(4.0*5.669D-08*PLANKC*(TKD**4.0-8.1D+09)
     1    +3.334D-04*STMF(J,I)*(RHO(J,I)/RSOOT)*TKD**5)
     
C-----------------------------------------------------------------------
      RHSZ(J,I)=DT*0.5*((RDHE/XXS(I+1)-RDHW/XXS(I))/XXC(I)
     1         +(RDHN/YYS(J+1)-RDHS/YYS(J))/YYC(J)
     2         +SYM*RDHP/Y(J,I)/YYC(J))-QRAD*DT
  411 CONTINUE
      ENDDO
  410 CONTINUE
C-----------------   RATES FOR ENTHALPY EQUATION   ---------------------
      DO 414 J=2,LJ-1
      DO I=2,LI-1
      IF(ISKIP(J,I).GE.10) GO TO 415
      IA=(I-1)*LJ+J
      RA(0001)=RC(IA,   1)*SSP(IA,  3)*SSP(IA,  2)
      RA(0002)=RC(IA,   2)*SSP(IA,  5)*SSP(IA,  4)
      RA(0003)=RC(IA,   3)*SSP(IA,  1)*SSP(IA,  4)
      RA(0004)=RC(IA,   4)*SSP(IA,  5)*SSP(IA,  3)
      RA(0005)=RC(IA,   5)*SSP(IA,  1)*SSP(IA,  5)
      RA(0006)=RC(IA,   6)*SSP(IA,  6)*SSP(IA,  3)
      RA(0007)=RC(IA,   7)*SSP(IA,  6)*SSP(IA,  4)
      RA(0008)=RC(IA,   8)*SSP(IA,  5)*SSP(IA,  5)
      RA(0009)=RC(IA,   9)*SSP(IA,  3)*SSP(IA,  3)
      RA(0010)=RC(IA,  10)*SSP(IA,  1)
      RA(0011)=RC(IA,  11)*SSP(IA,  3)*SSP(IA,  5)
      RA(0012)=RC(IA,  12)*SSP(IA,  6)
      RA(0013)=RC(IA,  13)*SSP(IA,  4)*SSP(IA,  4)
      RA(0014)=RC(IA,  14)*SSP(IA,  2)
      RA(0015)=RC(IA,  15)*SSP(IA,  3)*SSP(IA,  4)
      RA(0016)=RC(IA,  16)*SSP(IA,  5)
      RA(0017)=RC(IA,  17)*SSP(IA,  4)*SSP(IA,  5)
      RA(0018)=RC(IA,  18)*SSP(IA,  7)
      RA(0019)=RC(IA,  19)*SSP(IA,  3)*SSP(IA,  2)
      RA(0020)=RC(IA,  20)*SSP(IA,  7)
      RA(0021)=RC(IA,  21)*SSP(IA,  7)*SSP(IA,  3)
      RA(0022)=RC(IA,  22)*SSP(IA,  5)*SSP(IA,  5)
      RA(0023)=RC(IA,  23)*SSP(IA,  7)*SSP(IA,  3)
      RA(0024)=RC(IA,  24)*SSP(IA,  1)*SSP(IA,  2)
      RA(0025)=RC(IA,  25)*SSP(IA,  7)*SSP(IA,  3)
      RA(0026)=RC(IA,  26)*SSP(IA,  6)*SSP(IA,  4)
      RA(0027)=RC(IA,  27)*SSP(IA,  7)*SSP(IA,  4)
      RA(0028)=RC(IA,  28)*SSP(IA,  5)*SSP(IA,  2)
      RA(0029)=RC(IA,  29)*SSP(IA,  7)*SSP(IA,  5)
      RA(0030)=RC(IA,  30)*SSP(IA,  6)*SSP(IA,  2)
      RA(0031)=RC(IA,  31)*SSP(IA,  5)*SSP(IA,  5)
      RA(0032)=RC(IA,  32)*SSP(IA,  8)
      RA(0033)=RC(IA,  33)*SSP(IA,  7)*SSP(IA,  7)
      RA(0034)=RC(IA,  34)*SSP(IA,  8)*SSP(IA,  2)
      RA(0035)=RC(IA,  35)*SSP(IA,  8)*SSP(IA,  3)
      RA(0036)=RC(IA,  36)*SSP(IA,  7)*SSP(IA,  1)
      RA(0037)=RC(IA,  37)*SSP(IA,  8)*SSP(IA,  3)
      RA(0038)=RC(IA,  38)*SSP(IA,  6)*SSP(IA,  5)
      RA(0039)=RC(IA,  39)*SSP(IA,  8)*SSP(IA,  5)
      RA(0040)=RC(IA,  40)*SSP(IA,  6)*SSP(IA,  7)
      RA(0041)=RC(IA,  41)*SSP(IA,  8)*SSP(IA,  4)
      RA(0042)=RC(IA,  42)*SSP(IA,  7)*SSP(IA,  5)
      RA(0043)=RC(IA,  43)*SSP(IA,  9)*SSP(IA,  5)
      RA(0044)=RC(IA,  44)*SSP(IA, 10)*SSP(IA,  3)
      RA(0045)=RC(IA,  45)*SSP(IA,  9)*SSP(IA,  7)
      RA(0046)=RC(IA,  46)*SSP(IA, 10)*SSP(IA,  5)
      RA(0047)=RC(IA,  47)*SSP(IA,  9)*SSP(IA,  2)
      RA(0048)=RC(IA,  48)*SSP(IA, 10)*SSP(IA,  4)
      RA(0049)=RC(IA,  49)*SSP(IA, 11)
      RA(0050)=RC(IA,  50)*SSP(IA,  9)*SSP(IA,  3)
      RA(0051)=RC(IA,  51)*SSP(IA, 11)*SSP(IA,  3)
      RA(0052)=RC(IA,  52)*SSP(IA,  9)*SSP(IA,  1)
      RA(0053)=RC(IA,  53)*SSP(IA, 11)*SSP(IA,  4)
      RA(0054)=RC(IA,  54)*SSP(IA,  9)*SSP(IA,  5)
      RA(0055)=RC(IA,  55)*SSP(IA, 11)*SSP(IA,  4)
      RA(0056)=RC(IA,  56)*SSP(IA, 10)*SSP(IA,  3)
      RA(0057)=RC(IA,  57)*SSP(IA, 11)*SSP(IA,  5)
      RA(0058)=RC(IA,  58)*SSP(IA,  9)*SSP(IA,  6)
      RA(0059)=RC(IA,  59)*SSP(IA, 11)*SSP(IA,  2)
      RA(0060)=RC(IA,  60)*SSP(IA,  9)*SSP(IA,  7)
      RA(0061)=RC(IA,  61)*SSP(IA, 11)*SSP(IA, 14)
      RA(0062)=RC(IA,  62)*SSP(IA,  9)*SSP(IA, 13)
      RA(0063)=RC(IA,  63)*SSP(IA,  3)*SSP(IA, 11)
      RA(0064)=RC(IA,  64)*SSP(IA, 12)
      RA(0065)=RC(IA,  65)*SSP(IA, 12)*SSP(IA,  3)
      RA(0066)=RC(IA,  66)*SSP(IA, 11)*SSP(IA,  1)
      RA(0067)=RC(IA,  67)*SSP(IA, 12)*SSP(IA,  4)
      RA(0068)=RC(IA,  68)*SSP(IA, 11)*SSP(IA,  5)
      RA(0069)=RC(IA,  69)*SSP(IA, 12)*SSP(IA,  5)
      RA(0070)=RC(IA,  70)*SSP(IA, 11)*SSP(IA,  6)
      RA(0071)=RC(IA,  71)*SSP(IA, 12)*SSP(IA,  2)
      RA(0072)=RC(IA,  72)*SSP(IA, 11)*SSP(IA,  7)
      RA(0073)=RC(IA,  73)*SSP(IA, 12)*SSP(IA,  7)
      RA(0074)=RC(IA,  74)*SSP(IA, 11)*SSP(IA,  8)
      RA(0075)=RC(IA,  75)*SSP(IA, 13)*SSP(IA,  3)
      RA(0076)=RC(IA,  76)*SSP(IA,  1)*SSP(IA, 14)
      RA(0077)=RC(IA,  77)*SSP(IA, 13)*SSP(IA,  5)
      RA(0078)=RC(IA,  78)*SSP(IA,  6)*SSP(IA, 14)
      RA(0079)=RC(IA,  79)*SSP(IA, 13)*SSP(IA,  4)
      RA(0080)=RC(IA,  80)*SSP(IA, 14)*SSP(IA,  5)
      RA(0081)=RC(IA,  81)*SSP(IA, 13)*SSP(IA,  2)
      RA(0082)=RC(IA,  82)*SSP(IA, 14)*SSP(IA,  7)
      RA(0083)=RC(IA,  83)*SSP(IA, 13)*SSP(IA,  7)
      RA(0084)=RC(IA,  84)*SSP(IA, 14)*SSP(IA,  8)
      RA(0085)=RC(IA,  85)*SSP(IA, 14)*SSP(IA,  3)
      RA(0086)=RC(IA,  86)*SSP(IA, 15)*SSP(IA,  1)
      RA(0087)=RC(IA,  87)*SSP(IA, 14)*SSP(IA,  3)
      RA(0088)=RC(IA,  88)*SSP(IA, 16)*SSP(IA,  1)
      RA(0089)=RC(IA,  89)*SSP(IA, 14)*SSP(IA,  5)
      RA(0090)=RC(IA,  90)*SSP(IA, 16)*SSP(IA,  6)
      RA(0091)=RC(IA,  91)*SSP(IA, 14)*SSP(IA,  4)
      RA(0092)=RC(IA,  92)*SSP(IA, 12)*SSP(IA,  3)
      RA(0093)=RC(IA,  93)*SSP(IA, 14)*SSP(IA, 15)
      RA(0094)=RC(IA,  94)*SSP(IA, 17)*SSP(IA,  3)
      RA(0095)=RC(IA,  95)*SSP(IA, 14)*SSP(IA,  7)
      RA(0096)=RC(IA,  96)*SSP(IA, 18)*SSP(IA,  5)
      RA(0097)=RC(IA,  97)*SSP(IA, 14)*SSP(IA,  2)
      RA(0098)=RC(IA,  98)*SSP(IA, 12)*SSP(IA,  5)
      RA(0099)=RC(IA,  99)*SSP(IA, 14)*SSP(IA,  2)
      RA(0100)=RC(IA, 100)*SSP(IA, 18)*SSP(IA,  4)
      RA(0101)=RC(IA, 101)*SSP(IA, 14)*SSP(IA, 14)
      RA(0102)=RC(IA, 102)*SSP(IA, 17)*SSP(IA,  1)
      RA(0103)=RC(IA, 103)*SSP(IA, 14)*SSP(IA, 14)
      RA(0104)=RC(IA, 104)*SSP(IA, 19)*SSP(IA,  3)
      RA(0105)=RC(IA, 105)*SSP(IA,  3)*SSP(IA, 14)
      RA(0106)=RC(IA, 106)*SSP(IA, 13)
      RA(0107)=RC(IA, 107)*SSP(IA, 14)*SSP(IA, 14)
      RA(0108)=RC(IA, 108)*SSP(IA, 20)
      RA(0109)=RC(IA, 109)*SSP(IA, 16)*SSP(IA,  5)
      RA(0110)=RC(IA, 110)*SSP(IA, 12)*SSP(IA,  3)
      RA(0111)=RC(IA, 111)*SSP(IA, 16)*SSP(IA,  2)
      RA(0112)=RC(IA, 112)*SSP(IA,  9)*SSP(IA,  5)*SSP(IA,  3)
      RA(0113)=RC(IA, 113)*SSP(IA, 16)*SSP(IA, 10)
      RA(0114)=RC(IA, 114)*SSP(IA,  9)*SSP(IA, 12)
      RA(0115)=RC(IA, 115)*SSP(IA, 16)
      RA(0116)=RC(IA, 116)*SSP(IA, 15)
      RA(0117)=RC(IA, 117)*SSP(IA, 15)*SSP(IA,  3)
      RA(0118)=RC(IA, 118)*SSP(IA, 21)*SSP(IA,  1)
      RA(0119)=RC(IA, 119)*SSP(IA, 15)*SSP(IA,  5)
      RA(0120)=RC(IA, 120)*SSP(IA, 12)*SSP(IA,  3)
      RA(0121)=RC(IA, 121)*SSP(IA, 15)*SSP(IA,  5)
      RA(0122)=RC(IA, 122)*SSP(IA, 21)*SSP(IA,  6)
      RA(0123)=RC(IA, 123)*SSP(IA, 15)*SSP(IA,  4)
      RA(0124)=RC(IA, 124)*SSP(IA,  9)*SSP(IA,  3)*SSP(IA,  3)
      RA(0125)=RC(IA, 125)*SSP(IA, 15)*SSP(IA,  4)
      RA(0126)=RC(IA, 126)*SSP(IA,  9)*SSP(IA,  1)
      RA(0127)=RC(IA, 127)*SSP(IA, 15)*SSP(IA,  2)
      RA(0128)=RC(IA, 128)*SSP(IA, 10)*SSP(IA,  1)
      RA(0129)=RC(IA, 129)*SSP(IA, 15)*SSP(IA,  2)
      RA(0130)=RC(IA, 130)*SSP(IA,  9)*SSP(IA,  5)*SSP(IA,  3)
      RA(0131)=RC(IA, 131)*SSP(IA, 15)*SSP(IA, 15)
      RA(0132)=RC(IA, 132)*SSP(IA, 22)*SSP(IA,  3)*SSP(IA,  3)
      RA(0133)=RC(IA, 133)*SSP(IA, 21)*SSP(IA,  4)
      RA(0134)=RC(IA, 134)*SSP(IA,  9)*SSP(IA,  3)
      RA(0135)=RC(IA, 135)*SSP(IA, 21)*SSP(IA,  2)
      RA(0136)=RC(IA, 136)*SSP(IA, 11)*SSP(IA,  4)
      RA(0137)=RC(IA, 137)*SSP(IA, 21)*SSP(IA,  6)
      RA(0138)=RC(IA, 138)*SSP(IA, 12)*SSP(IA,  3)
      RA(0139)=RC(IA, 139)*SSP(IA, 21)*SSP(IA, 10)
      RA(0140)=RC(IA, 140)*SSP(IA, 11)*SSP(IA,  9)
      RA(0141)=RC(IA, 141)*SSP(IA, 18)*SSP(IA,  3)
      RA(0142)=RC(IA, 142)*SSP(IA, 12)*SSP(IA,  1)
      RA(0143)=RC(IA, 143)*SSP(IA, 18)*SSP(IA,  3)
      RA(0144)=RC(IA, 144)*SSP(IA, 16)*SSP(IA,  6)
      RA(0145)=RC(IA, 145)*SSP(IA, 18)*SSP(IA,  5)
      RA(0146)=RC(IA, 146)*SSP(IA, 12)*SSP(IA,  6)
      RA(0147)=RC(IA, 147)*SSP(IA, 18)*SSP(IA,  4)
      RA(0148)=RC(IA, 148)*SSP(IA,  5)*SSP(IA, 12)
      RA(0149)=RC(IA, 149)*SSP(IA, 18)*SSP(IA,  2)
      RA(0150)=RC(IA, 150)*SSP(IA, 12)*SSP(IA,  7)
      RA(0151)=RC(IA, 151)*SSP(IA, 18)
      RA(0152)=RC(IA, 152)*SSP(IA, 12)*SSP(IA,  3)
      RA(0153)=RC(IA, 153)*SSP(IA, 20)*SSP(IA,  3)
      RA(0154)=RC(IA, 154)*SSP(IA, 19)*SSP(IA,  1)
      RA(0155)=RC(IA, 155)*SSP(IA, 20)*SSP(IA,  4)
      RA(0156)=RC(IA, 156)*SSP(IA, 19)*SSP(IA,  5)
      RA(0157)=RC(IA, 157)*SSP(IA, 20)*SSP(IA,  5)
      RA(0158)=RC(IA, 158)*SSP(IA, 19)*SSP(IA,  6)
      RA(0159)=RC(IA, 159)*SSP(IA, 20)*SSP(IA, 14)
      RA(0160)=RC(IA, 160)*SSP(IA, 19)*SSP(IA, 13)
      RA(0161)=RC(IA, 161)*SSP(IA, 20)
      RA(0162)=RC(IA, 162)*SSP(IA, 19)*SSP(IA,  3)
      RA(0163)=RC(IA, 163)*SSP(IA, 20)*SSP(IA,  7)
      RA(0164)=RC(IA, 164)*SSP(IA, 19)*SSP(IA,  8)
      RA(0165)=RC(IA, 165)*SSP(IA, 19)*SSP(IA,  3)
      RA(0166)=RC(IA, 166)*SSP(IA, 17)*SSP(IA,  1)
      RA(0167)=RC(IA, 167)*SSP(IA, 19)*SSP(IA,  4)
      RA(0168)=RC(IA, 168)*SSP(IA, 17)*SSP(IA,  5)
      RA(0169)=RC(IA, 169)*SSP(IA, 19)*SSP(IA,  4)
      RA(0170)=RC(IA, 170)*SSP(IA, 14)*SSP(IA, 12)
      RA(0171)=RC(IA, 171)*SSP(IA, 19)*SSP(IA,  2)
      RA(0172)=RC(IA, 172)*SSP(IA, 17)*SSP(IA,  7)
      RA(0173)=RC(IA, 173)*SSP(IA, 19)
      RA(0174)=RC(IA, 174)*SSP(IA, 17)*SSP(IA,  3)
      RA(0175)=RC(IA, 175)*SSP(IA, 17)*SSP(IA,  3)
      RA(0176)=RC(IA, 176)*SSP(IA, 23)*SSP(IA,  1)
      RA(0177)=RC(IA, 177)*SSP(IA, 17)*SSP(IA,  5)
      RA(0178)=RC(IA, 178)*SSP(IA, 23)*SSP(IA,  6)
      RA(0179)=RC(IA, 179)*SSP(IA, 17)*SSP(IA,  4)
      RA(0180)=RC(IA, 180)*SSP(IA, 14)*SSP(IA, 11)
      RA(0181)=RC(IA, 181)*SSP(IA, 17)*SSP(IA,  4)
      RA(0182)=RC(IA, 182)*SSP(IA, 24)*SSP(IA,  3)
      RA(0183)=RC(IA, 183)*SSP(IA, 17)*SSP(IA, 17)
      RA(0184)=RC(IA, 184)*SSP(IA, 23)*SSP(IA, 19)
      RA(0185)=RC(IA, 185)*SSP(IA, 17)*SSP(IA,  2)
      RA(0186)=RC(IA, 186)*SSP(IA, 23)*SSP(IA,  7)
      RA(0187)=RC(IA, 187)*SSP(IA, 17)*SSP(IA,  7)
      RA(0188)=RC(IA, 188)*SSP(IA, 25)*SSP(IA,  5)
      RA(0189)=RC(IA, 189)*SSP(IA, 25)*SSP(IA,  7)
      RA(0190)=RC(IA, 190)*SSP(IA, 14)*SSP(IA,  9)*SSP(IA,  8)
      RA(0191)=RC(IA, 191)*SSP(IA, 17)
      RA(0192)=RC(IA, 192)*SSP(IA, 23)*SSP(IA,  3)
      RA(0193)=RC(IA, 193)*SSP(IA, 17)
      RA(0194)=RC(IA, 194)*SSP(IA, 22)*SSP(IA,  1)
      RA(0195)=RC(IA, 195)*SSP(IA, 23)*SSP(IA,  3)
      RA(0196)=RC(IA, 196)*SSP(IA, 22)*SSP(IA,  1)
      RA(0197)=RC(IA, 197)*SSP(IA, 23)
      RA(0198)=RC(IA, 198)*SSP(IA, 22)*SSP(IA,  3)
      RA(0199)=RC(IA, 199)*SSP(IA, 23)*SSP(IA,  2)
      RA(0200)=RC(IA, 200)*SSP(IA, 12)*SSP(IA, 11)
      RA(0201)=RC(IA, 201)*SSP(IA, 23)*SSP(IA,  2)
      RA(0202)=RC(IA, 202)*SSP(IA, 24)*SSP(IA,  4)
      RA(0203)=RC(IA, 203)*SSP(IA, 23)*SSP(IA,  2)
      RA(0204)=RC(IA, 204)*SSP(IA, 22)*SSP(IA,  7)
      RA(0205)=RC(IA, 205)*SSP(IA, 22)*SSP(IA,  4)
      RA(0206)=RC(IA, 206)*SSP(IA, 27)*SSP(IA,  3)
      RA(0207)=RC(IA, 207)*SSP(IA, 22)*SSP(IA,  4)
      RA(0208)=RC(IA, 208)*SSP(IA, 15)*SSP(IA,  9)
      RA(0209)=RC(IA, 209)*SSP(IA, 22)*SSP(IA,  2)
      RA(0210)=RC(IA, 210)*SSP(IA, 12)*SSP(IA,  9)
      RA(0211)=RC(IA, 211)*SSP(IA, 22)*SSP(IA,  5)
      RA(0212)=RC(IA, 212)*SSP(IA, 26)*SSP(IA,  3)
      RA(0213)=RC(IA, 213)*SSP(IA, 22)*SSP(IA,  5)
      RA(0214)=RC(IA, 214)*SSP(IA, 28)*SSP(IA,  6)
      RA(0215)=RC(IA, 215)*SSP(IA, 26)*SSP(IA,  3)
      RA(0216)=RC(IA, 216)*SSP(IA, 14)*SSP(IA,  9)
      RA(0217)=RC(IA, 217)*SSP(IA, 26)*SSP(IA,  4)
      RA(0218)=RC(IA, 218)*SSP(IA, 15)*SSP(IA, 10)
      RA(0219)=RC(IA, 219)*SSP(IA, 26)*SSP(IA,  4)
      RA(0220)=RC(IA, 220)*SSP(IA, 27)*SSP(IA,  5)
      RA(0221)=RC(IA, 221)*SSP(IA, 26)*SSP(IA, 14)
      RA(0222)=RC(IA, 222)*SSP(IA, 19)*SSP(IA,  9)
      RA(0223)=RC(IA, 223)*SSP(IA, 27)*SSP(IA,  3)
      RA(0224)=RC(IA, 224)*SSP(IA, 16)*SSP(IA,  9)
      RA(0225)=RC(IA, 225)*SSP(IA, 27)*SSP(IA,  5)
      RA(0226)=RC(IA, 226)*SSP(IA, 11)*SSP(IA,  9)*SSP(IA,  3)
      RA(0227)=RC(IA, 227)*SSP(IA, 27)*SSP(IA,  4)
      RA(0228)=RC(IA, 228)*SSP(IA,  9)*SSP(IA,  9)*SSP(IA,  3)
      RA(0229)=RC(IA, 229)*SSP(IA, 27)*SSP(IA,  2)
      RA(0230)=RC(IA, 230)*SSP(IA,  9)*SSP(IA,  9)*SSP(IA,  5)
      RA(0231)=RC(IA, 231)*SSP(IA, 27)*SSP(IA,  2)
      RA(0232)=RC(IA, 232)*SSP(IA, 10)*SSP(IA,  9)*SSP(IA,  3)
      RA(0233)=RC(IA, 233)*SSP(IA, 28)*SSP(IA,  5)
      RA(0234)=RC(IA, 234)*SSP(IA, 27)*SSP(IA,  3)
      RA(0235)=RC(IA, 235)*SSP(IA, 28)*SSP(IA,  4)
      RA(0236)=RC(IA, 236)*SSP(IA,  9)*SSP(IA, 21)
      RA(0237)=RC(IA, 237)*SSP(IA, 28)*SSP(IA,  2)
      RA(0238)=RC(IA, 238)*SSP(IA, 27)*SSP(IA,  4)
      RA(0239)=RC(IA, 239)*SSP(IA, 28)*SSP(IA,  2)
      RA(0240)=RC(IA, 240)*SSP(IA, 21)*SSP(IA, 10)
      RA(0241)=RC(IA, 241)*SSP(IA, 28)*SSP(IA,  2)
      RA(0242)=RC(IA, 242)*SSP(IA, 11)*SSP(IA,  9)
      RA(0243)=RC(IA, 243)*SSP(IA, 29)*SSP(IA,  3)
      RA(0244)=RC(IA, 244)*SSP(IA, 12)*SSP(IA,  1)
      RA(0245)=RC(IA, 245)*SSP(IA, 29)*SSP(IA,  3)
      RA(0246)=RC(IA, 246)*SSP(IA, 14)*SSP(IA,  5)
      RA(0247)=RC(IA, 247)*SSP(IA, 29)*SSP(IA,  5)
      RA(0248)=RC(IA, 248)*SSP(IA, 12)*SSP(IA,  6)
      RA(0249)=RC(IA, 249)*SSP(IA, 29)*SSP(IA,  2)
      RA(0250)=RC(IA, 250)*SSP(IA, 12)*SSP(IA,  7)
      RA(0251)=RC(IA, 251)*SSP(IA, 29)
      RA(0252)=RC(IA, 252)*SSP(IA, 12)*SSP(IA,  3)
      RA(0253)=RC(IA, 253)*SSP(IA, 18)
      RA(0254)=RC(IA, 254)*SSP(IA, 29)
      RA(0255)=RC(IA, 255)*SSP(IA, 26)*SSP(IA,  5)
      RA(0256)=RC(IA, 256)*SSP(IA, 29)*SSP(IA,  9)
      RA(0257)=RC(IA, 257)*SSP(IA, 30)*SSP(IA,  5)
      RA(0258)=RC(IA, 258)*SSP(IA, 29)*SSP(IA,  6)
      RA(0259)=RC(IA, 259)*SSP(IA, 30)*SSP(IA,  5)
      RA(0260)=RC(IA, 260)*SSP(IA, 18)*SSP(IA,  6)
      RA(0261)=RC(IA, 261)*SSP(IA, 30)*SSP(IA,  3)
      RA(0262)=RC(IA, 262)*SSP(IA, 29)*SSP(IA,  1)
      RA(0263)=RC(IA, 263)*SSP(IA, 30)*SSP(IA,  3)
      RA(0264)=RC(IA, 264)*SSP(IA, 18)*SSP(IA,  1)
      RA(0265)=RC(IA, 265)*SSP(IA, 30)*SSP(IA,  4)
      RA(0266)=RC(IA, 266)*SSP(IA, 29)*SSP(IA,  5)
      RA(0267)=RC(IA, 267)*SSP(IA, 30)*SSP(IA,  7)
      RA(0268)=RC(IA, 268)*SSP(IA, 29)*SSP(IA,  8)
      RA(0269)=RC(IA, 269)*SSP(IA, 30)*SSP(IA,  2)
      RA(0270)=RC(IA, 270)*SSP(IA, 29)*SSP(IA,  7)
      RA(0271)=RC(IA, 271)*SSP(IA, 30)
      RA(0272)=RC(IA, 272)*SSP(IA, 14)*SSP(IA,  5)
      RA(0273)=RC(IA, 273)*SSP(IA, 24)
      RA(0274)=RC(IA, 274)*SSP(IA, 26)*SSP(IA,  3)
      RA(0275)=RC(IA, 275)*SSP(IA, 24)*SSP(IA,  3)
      RA(0276)=RC(IA, 276)*SSP(IA, 14)*SSP(IA, 11)
      RA(0277)=RC(IA, 277)*SSP(IA, 24)*SSP(IA,  3)
      RA(0278)=RC(IA, 278)*SSP(IA, 26)*SSP(IA,  1)
      RA(0279)=RC(IA, 279)*SSP(IA, 24)*SSP(IA,  4)
      RA(0280)=RC(IA, 280)*SSP(IA, 12)*SSP(IA, 11)
      RA(0281)=RC(IA, 281)*SSP(IA, 24)*SSP(IA,  5)
      RA(0282)=RC(IA, 282)*SSP(IA, 26)*SSP(IA,  6)
      RA(0283)=RC(IA, 283)*SSP(IA, 24)*SSP(IA,  2)
      RA(0284)=RC(IA, 284)*SSP(IA, 12)*SSP(IA,  9)*SSP(IA,  5)
      RA(0285)=RC(IA, 285)*SSP(IA, 24)*SSP(IA, 14)
      RA(0286)=RC(IA, 286)*SSP(IA, 19)*SSP(IA,  9)*SSP(IA,  3)
      RA(0287)=RC(IA, 287)*SSP(IA, 24)*SSP(IA,  7)
      RA(0288)=RC(IA, 288)*SSP(IA, 12)*SSP(IA, 11)*SSP(IA,  5)
      RA(0289)=RC(IA, 289)*SSP(IA, 24)*SSP(IA,  7)
      RA(0290)=RC(IA, 290)*SSP(IA, 32)*SSP(IA,  2)
      RA(0291)=RC(IA, 291)*SSP(IA, 24)
      RA(0292)=RC(IA, 292)*SSP(IA, 14)*SSP(IA,  9)
      RA(0293)=RC(IA, 293)*SSP(IA, 32)
      RA(0294)=RC(IA, 294)*SSP(IA, 14)*SSP(IA, 11)
      RA(0295)=RC(IA, 295)*SSP(IA, 35)
      RA(0296)=RC(IA, 296)*SSP(IA, 14)*SSP(IA,  9)
      RA(0297)=RC(IA, 297)*SSP(IA, 32)*SSP(IA,  5)
      RA(0298)=RC(IA, 298)*SSP(IA, 35)*SSP(IA,  6)
      RA(0299)=RC(IA, 299)*SSP(IA, 32)*SSP(IA,  5)
      RA(0300)=RC(IA, 300)*SSP(IA, 24)*SSP(IA,  6)
      RA(0301)=RC(IA, 301)*SSP(IA, 32)*SSP(IA,  4)
      RA(0302)=RC(IA, 302)*SSP(IA, 35)*SSP(IA,  5)
      RA(0303)=RC(IA, 303)*SSP(IA, 32)*SSP(IA,  4)
      RA(0304)=RC(IA, 304)*SSP(IA, 24)*SSP(IA,  5)
      RA(0305)=RC(IA, 305)*SSP(IA, 32)*SSP(IA,  3)
      RA(0306)=RC(IA, 306)*SSP(IA, 35)*SSP(IA,  1)
      RA(0307)=RC(IA, 307)*SSP(IA, 32)*SSP(IA,  3)
      RA(0308)=RC(IA, 308)*SSP(IA, 24)*SSP(IA,  1)
      RA(0309)=RC(IA, 309)*SSP(IA, 32)*SSP(IA, 14)
      RA(0310)=RC(IA, 310)*SSP(IA, 35)*SSP(IA, 13)
      RA(0311)=RC(IA, 311)*SSP(IA, 32)*SSP(IA, 14)
      RA(0312)=RC(IA, 312)*SSP(IA, 24)*SSP(IA, 13)
      RA(0313)=RC(IA, 313)*SSP(IA, 32)*SSP(IA,  7)
      RA(0314)=RC(IA, 314)*SSP(IA, 35)*SSP(IA,  8)
      RA(0315)=RC(IA, 315)*SSP(IA, 32)*SSP(IA,  7)
      RA(0316)=RC(IA, 316)*SSP(IA, 24)*SSP(IA,  8)
      RA(0317)=RC(IA, 317)*SSP(IA, 32)*SSP(IA,  2)
      RA(0318)=RC(IA, 318)*SSP(IA, 35)*SSP(IA,  7)
      RA(0319)=RC(IA, 319)*SSP(IA, 31)
      RA(0320)=RC(IA, 320)*SSP(IA, 14)*SSP(IA, 29)
      RA(0321)=RC(IA, 321)*SSP(IA, 31)
      RA(0322)=RC(IA, 322)*SSP(IA, 17)*SSP(IA,  6)
      RA(0323)=RC(IA, 323)*SSP(IA, 31)*SSP(IA,  5)
      RA(0324)=RC(IA, 324)*SSP(IA, 34)*SSP(IA,  6)
      RA(0325)=RC(IA, 325)*SSP(IA, 31)*SSP(IA,  5)
      RA(0326)=RC(IA, 326)*SSP(IA, 33)*SSP(IA,  6)
      RA(0327)=RC(IA, 327)*SSP(IA, 31)*SSP(IA,  5)
      RA(0328)=RC(IA, 328)*SSP(IA, 36)*SSP(IA,  6)
      RA(0329)=RC(IA, 329)*SSP(IA, 31)*SSP(IA,  3)
      RA(0330)=RC(IA, 330)*SSP(IA, 34)*SSP(IA,  1)
      RA(0331)=RC(IA, 331)*SSP(IA, 31)*SSP(IA,  3)
      RA(0332)=RC(IA, 332)*SSP(IA, 33)*SSP(IA,  1)
      RA(0333)=RC(IA, 333)*SSP(IA, 31)*SSP(IA,  3)
      RA(0334)=RC(IA, 334)*SSP(IA, 36)*SSP(IA,  1)
      RA(0335)=RC(IA, 335)*SSP(IA, 31)*SSP(IA,  4)
      RA(0336)=RC(IA, 336)*SSP(IA, 34)*SSP(IA,  5)
      RA(0337)=RC(IA, 337)*SSP(IA, 31)*SSP(IA,  4)
      RA(0338)=RC(IA, 338)*SSP(IA, 33)*SSP(IA,  5)
      RA(0339)=RC(IA, 339)*SSP(IA, 31)*SSP(IA,  4)
      RA(0340)=RC(IA, 340)*SSP(IA, 36)*SSP(IA,  5)
      RA(0341)=RC(IA, 341)*SSP(IA, 31)*SSP(IA, 14)
      RA(0342)=RC(IA, 342)*SSP(IA, 34)*SSP(IA, 13)
      RA(0343)=RC(IA, 343)*SSP(IA, 31)*SSP(IA, 14)
      RA(0344)=RC(IA, 344)*SSP(IA, 33)*SSP(IA, 13)
      RA(0345)=RC(IA, 345)*SSP(IA, 31)*SSP(IA, 14)
      RA(0346)=RC(IA, 346)*SSP(IA, 36)*SSP(IA, 13)
      RA(0347)=RC(IA, 347)*SSP(IA, 31)*SSP(IA,  7)
      RA(0348)=RC(IA, 348)*SSP(IA, 33)*SSP(IA,  8)
      RA(0349)=RC(IA, 349)*SSP(IA, 31)*SSP(IA,  7)
      RA(0350)=RC(IA, 350)*SSP(IA, 34)*SSP(IA,  8)
      RA(0351)=RC(IA, 351)*SSP(IA, 31)*SSP(IA,  7)
      RA(0352)=RC(IA, 352)*SSP(IA, 36)*SSP(IA,  8)
      RA(0353)=RC(IA, 353)*SSP(IA, 17)*SSP(IA,  5)
      RA(0354)=RC(IA, 354)*SSP(IA, 34)
      RA(0355)=RC(IA, 355)*SSP(IA, 19)*SSP(IA,  7)
      RA(0356)=RC(IA, 356)*SSP(IA, 36)*SSP(IA,  5)
      RA(0357)=RC(IA, 357)*SSP(IA, 36)
      RA(0358)=RC(IA, 358)*SSP(IA, 32)*SSP(IA,  3)
      RA(0359)=RC(IA, 359)*SSP(IA, 36)
      RA(0360)=RC(IA, 360)*SSP(IA, 14)*SSP(IA, 12)
      RA(0361)=RC(IA, 361)*SSP(IA, 36)*SSP(IA,  2)
      RA(0362)=RC(IA, 362)*SSP(IA, 32)*SSP(IA,  7)
      RA(0363)=RC(IA, 363)*SSP(IA, 36)*SSP(IA,  9)
      RA(0364)=RC(IA, 364)*SSP(IA, 19)*SSP(IA, 10)
      RA(0365)=RC(IA, 365)*SSP(IA, 36)*SSP(IA,  3)
      RA(0366)=RC(IA, 366)*SSP(IA, 14)*SSP(IA, 29)
      RA(0367)=RC(IA, 367)*SSP(IA, 36)*SSP(IA,  3)
      RA(0368)=RC(IA, 368)*SSP(IA, 17)*SSP(IA,  6)
      RA(0369)=RC(IA, 369)*SSP(IA, 36)*SSP(IA,  5)
      RA(0370)=RC(IA, 370)*SSP(IA, 32)*SSP(IA,  6)
      RA(0371)=RC(IA, 371)*SSP(IA, 33)*SSP(IA,  2)
      RA(0372)=RC(IA, 372)*SSP(IA, 32)*SSP(IA,  7)
      RA(0373)=RC(IA, 373)*SSP(IA, 33)*SSP(IA,  4)
      RA(0374)=RC(IA, 374)*SSP(IA, 32)*SSP(IA,  5)
      RA(0375)=RC(IA, 375)*SSP(IA, 33)*SSP(IA,  3)
      RA(0376)=RC(IA, 376)*SSP(IA, 17)*SSP(IA,  6)
      RA(0377)=RC(IA, 377)*SSP(IA, 33)*SSP(IA,  3)
      RA(0378)=RC(IA, 378)*SSP(IA, 14)*SSP(IA, 29)
      RA(0379)=RC(IA, 379)*SSP(IA, 33)*SSP(IA,  7)
      RA(0380)=RC(IA, 380)*SSP(IA, 32)*SSP(IA,  5)*SSP(IA,  5)
      RA(0381)=RC(IA, 381)*SSP(IA, 33)*SSP(IA,  5)
      RA(0382)=RC(IA, 382)*SSP(IA, 32)*SSP(IA,  6)
      RA(0383)=RC(IA, 383)*SSP(IA, 33)
      RA(0384)=RC(IA, 384)*SSP(IA, 32)*SSP(IA,  3)
      RA(0385)=RC(IA, 385)*SSP(IA, 37)*SSP(IA,  4)
      RA(0386)=RC(IA, 386)*SSP(IA, 17)*SSP(IA,  9)
      RA(0387)=RC(IA, 387)*SSP(IA, 14)*SSP(IA, 22)
      RA(0388)=RC(IA, 388)*SSP(IA, 37)*SSP(IA,  3)
      RA(0389)=RC(IA, 389)*SSP(IA, 37)*SSP(IA,  4)
      RA(0390)=RC(IA, 390)*SSP(IA, 27)*SSP(IA, 14)
      RA(0391)=RC(IA, 391)*SSP(IA, 38)*SSP(IA,  3)
      RA(0392)=RC(IA, 392)*SSP(IA, 37)
      RA(0393)=RC(IA, 393)*SSP(IA, 38)*SSP(IA,  7)
      RA(0394)=RC(IA, 394)*SSP(IA, 37)*SSP(IA,  2)
      RA(0395)=RC(IA, 395)*SSP(IA, 37)*SSP(IA,  5)
      RA(0396)=RC(IA, 396)*SSP(IA, 38)*SSP(IA,  6)
      RA(0397)=RC(IA, 397)*SSP(IA, 38)*SSP(IA,  2)
      RA(0398)=RC(IA, 398)*SSP(IA, 26)*SSP(IA, 11)
      RA(0399)=RC(IA, 399)*SSP(IA, 37)*SSP(IA,  3)
      RA(0400)=RC(IA, 400)*SSP(IA, 39)
      RA(0401)=RC(IA, 401)*SSP(IA, 39)*SSP(IA,  3)
      RA(0402)=RC(IA, 402)*SSP(IA, 37)*SSP(IA,  1)
      RA(0403)=RC(IA, 403)*SSP(IA, 39)*SSP(IA,  2)
      RA(0404)=RC(IA, 404)*SSP(IA, 37)*SSP(IA,  7)
      RA(0405)=RC(IA, 405)*SSP(IA, 39)*SSP(IA, 14)
      RA(0406)=RC(IA, 406)*SSP(IA, 37)*SSP(IA, 13)
      RA(0407)=RC(IA, 407)*SSP(IA, 22)*SSP(IA, 14)
      RA(0408)=RC(IA, 408)*SSP(IA, 39)
      RA(0409)=RC(IA, 409)*SSP(IA, 39)*SSP(IA,  5)
      RA(0410)=RC(IA, 410)*SSP(IA, 37)*SSP(IA,  6)
      RA(0411)=RC(IA, 411)*SSP(IA, 38)*SSP(IA, 11)
      RA(0412)=RC(IA, 412)*SSP(IA, 37)*SSP(IA,  9)
      RA(0413)=RC(IA, 413)*SSP(IA, 38)*SSP(IA,  7)
      RA(0414)=RC(IA, 414)*SSP(IA,  5)*SSP(IA,  9)*SSP(IA, 23)
      RA(0415)=RC(IA, 415)*SSP(IA, 37)*SSP(IA,  2)
      RA(0416)=RC(IA, 416)*SSP(IA, 14)*SSP(IA, 11)*SSP(IA,  9)
      RA(0417)=RC(IA, 417)*SSP(IA, 40)*SSP(IA,  4)
      RA(0418)=RC(IA, 418)*SSP(IA, 19)*SSP(IA, 11)
      RA(0419)=RC(IA, 419)*SSP(IA, 40)*SSP(IA,  5)
      RA(0420)=RC(IA, 420)*SSP(IA, 39)*SSP(IA,  6)
      RA(0421)=RC(IA, 421)*SSP(IA, 40)*SSP(IA,  4)
      RA(0422)=RC(IA, 422)*SSP(IA, 26)*SSP(IA, 14)*SSP(IA,  3)
      RA(0423)=RC(IA, 423)*SSP(IA, 40)*SSP(IA,  3)
      RA(0424)=RC(IA, 424)*SSP(IA, 39)*SSP(IA,  1)
      RA(0425)=RC(IA, 425)*SSP(IA, 39)*SSP(IA,  3)
      RA(0426)=RC(IA, 426)*SSP(IA, 40)
      RA(0427)=RC(IA, 427)*SSP(IA, 39)*SSP(IA,  7)
      RA(0428)=RC(IA, 428)*SSP(IA, 40)*SSP(IA,  2)
      RA(0429)=RC(IA, 429)*SSP(IA, 39)*SSP(IA,  7)
      RA(0430)=RC(IA, 430)*SSP(IA,  5)*SSP(IA, 23)*SSP(IA, 12)
      RA(0431)=RC(IA, 431)*SSP(IA, 23)*SSP(IA, 14)
      RA(0432)=RC(IA, 432)*SSP(IA, 40)
      RA(0433)=RC(IA, 433)*SSP(IA, 40)*SSP(IA,  3)
      RA(0434)=RC(IA, 434)*SSP(IA, 17)*SSP(IA, 14)
      RA(0435)=RC(IA, 435)*SSP(IA, 14)*SSP(IA, 23)
      RA(0436)=RC(IA, 436)*SSP(IA, 39)*SSP(IA,  3)
      RA(0437)=RC(IA, 437)*SSP(IA, 41)
      RA(0438)=RC(IA, 438)*SSP(IA, 14)*SSP(IA, 19)
      RA(0439)=RC(IA, 439)*SSP(IA, 41)*SSP(IA,  2)
      RA(0440)=RC(IA, 440)*SSP(IA, 42)*SSP(IA,  7)
      RA(0441)=RC(IA, 441)*SSP(IA, 41)*SSP(IA,  2)
      RA(0442)=RC(IA, 442)*SSP(IA, 43)*SSP(IA,  7)
      RA(0443)=RC(IA, 443)*SSP(IA, 41)*SSP(IA,  3)
      RA(0444)=RC(IA, 444)*SSP(IA, 42)*SSP(IA,  1)
      RA(0445)=RC(IA, 445)*SSP(IA, 41)*SSP(IA,  3)
      RA(0446)=RC(IA, 446)*SSP(IA, 43)*SSP(IA,  1)
      RA(0447)=RC(IA, 447)*SSP(IA, 41)*SSP(IA,  4)
      RA(0448)=RC(IA, 448)*SSP(IA, 42)*SSP(IA,  5)
      RA(0449)=RC(IA, 449)*SSP(IA, 41)*SSP(IA,  4)
      RA(0450)=RC(IA, 450)*SSP(IA, 43)*SSP(IA,  5)
      RA(0451)=RC(IA, 451)*SSP(IA, 41)*SSP(IA,  5)
      RA(0452)=RC(IA, 452)*SSP(IA, 43)*SSP(IA,  6)
      RA(0453)=RC(IA, 453)*SSP(IA, 41)*SSP(IA,  5)
      RA(0454)=RC(IA, 454)*SSP(IA, 42)*SSP(IA,  6)
      RA(0455)=RC(IA, 455)*SSP(IA, 41)*SSP(IA,  7)
      RA(0456)=RC(IA, 456)*SSP(IA, 42)*SSP(IA,  8)
      RA(0457)=RC(IA, 457)*SSP(IA, 41)*SSP(IA,  7)
      RA(0458)=RC(IA, 458)*SSP(IA, 43)*SSP(IA,  8)
      RA(0459)=RC(IA, 459)*SSP(IA, 42)
      RA(0460)=RC(IA, 460)*SSP(IA, 43)
      RA(0461)=RC(IA, 461)*SSP(IA, 40)*SSP(IA,  3)
      RA(0462)=RC(IA, 462)*SSP(IA, 42)
      RA(0463)=RC(IA, 463)*SSP(IA, 42)*SSP(IA,  2)
      RA(0464)=RC(IA, 464)*SSP(IA, 40)*SSP(IA,  7)
      RA(0465)=RC(IA, 465)*SSP(IA, 43)
      RA(0466)=RC(IA, 466)*SSP(IA, 14)*SSP(IA, 17)
      RA(0467)=RC(IA, 467)*SSP(IA,  3)*SSP(IA, 40)
      RA(0468)=RC(IA, 468)*SSP(IA, 43)
      RA(0469)=RC(IA, 469)*SSP(IA, 43)*SSP(IA,  2)
      RA(0470)=RC(IA, 470)*SSP(IA, 40)*SSP(IA,  7)
      RA(0471)=RC(IA, 471)*SSP(IA, 44)
      RA(0472)=RC(IA, 472)*SSP(IA, 22)*SSP(IA, 23)*SSP(IA,  3)
      RA(0473)=RC(IA, 473)*SSP(IA, 44)
      RA(0474)=RC(IA, 474)*SSP(IA, 23)*SSP(IA, 23)
      RA(0475)=RC(IA, 475)*SSP(IA, 23)*SSP(IA, 23)
      RA(0476)=RC(IA, 476)*SSP(IA, 44)
      RA(0477)=RC(IA, 477)*SSP(IA, 44)*SSP(IA,  3)
      RA(0478)=RC(IA, 478)*SSP(IA, 23)*SSP(IA, 17)
      RA(0479)=RC(IA, 479)*SSP(IA, 44)*SSP(IA,  3)
      RA(0480)=RC(IA, 480)*SSP(IA,  1)*SSP(IA, 22)*SSP(IA, 23)
      RA(0481)=RC(IA, 481)*SSP(IA, 44)*SSP(IA,  5)
      RA(0482)=RC(IA, 482)*SSP(IA, 45)*SSP(IA,  3)*SSP(IA, 39)
      RA(0483)=RC(IA, 483)*SSP(IA, 44)*SSP(IA, 14)
      RA(0484)=RC(IA, 484)*SSP(IA, 13)*SSP(IA, 22)*SSP(IA, 23)
      RA(0485)=RC(IA, 485)*SSP(IA, 38)*SSP(IA, 14)
      RA(0486)=RC(IA, 486)*SSP(IA, 44)
      RA(0487)=RC(IA, 487)*SSP(IA, 46)
      RA(0488)=RC(IA, 488)*SSP(IA, 40)*SSP(IA, 22)
      RA(0489)=RC(IA, 489)*SSP(IA, 46)
      RA(0490)=RC(IA, 490)*SSP(IA, 37)*SSP(IA, 17)
      RA(0491)=RC(IA, 491)*SSP(IA, 46)
      RA(0492)=RC(IA, 492)*SSP(IA, 39)*SSP(IA, 23)
      RA(0493)=RC(IA, 493)*SSP(IA, 46)*SSP(IA,  2)
      RA(0494)=RC(IA, 494)*SSP(IA, 22)*SSP(IA, 39)*SSP(IA,  7)
      RA(0495)=RC(IA, 495)*SSP(IA, 46)*SSP(IA,  2)
      RA(0496)=RC(IA, 496)*SSP(IA, 23)*SSP(IA, 37)*SSP(IA,  7)
      RA(0497)=RC(IA, 497)*SSP(IA, 46)*SSP(IA,  7)
      RA(0498)=RC(IA, 498)*SSP(IA, 22)*SSP(IA, 39)*SSP(IA,  8)
      RA(0499)=RC(IA, 499)*SSP(IA, 46)*SSP(IA,  7)
      RA(0500)=RC(IA, 500)*SSP(IA, 23)*SSP(IA, 37)*SSP(IA,  8)
      RA(0501)=RC(IA, 501)*SSP(IA, 47)
      RA(0502)=0.0
      RA(0503)=RC(IA, 503)*SSP(IA, 47)*SSP(IA,  3)
      RA(0504)=0.0
      RA(0505)=RC(IA, 505)*SSP(IA, 47)*SSP(IA,  3)
      RA(0506)=0.0
      RA(0507)=RC(IA, 507)*SSP(IA, 47)*SSP(IA,  3)
      RA(0508)=0.0
      RA(0509)=RC(IA, 509)*SSP(IA, 47)*SSP(IA,  3)
      RA(0510)=RC(IA, 510)*SSP(IA,  1)*SSP(IA, 49)*SSP(IA, 19)
      RA(0511)=RC(IA, 511)*SSP(IA, 47)*SSP(IA,  5)
      RA(0512)=RC(IA, 512)*SSP(IA,  6)*SSP(IA, 49)*SSP(IA, 19)
      RA(0513)=RC(IA, 513)*SSP(IA, 47)*SSP(IA,  5)
      RA(0514)=0.0
      RA(0515)=RC(IA, 515)*SSP(IA, 47)*SSP(IA,  5)
      RA(0516)=0.0
      RA(0517)=RC(IA, 517)*SSP(IA, 47)*SSP(IA,  5)
      RA(0518)=0.0
      RA(0519)=RC(IA, 519)*SSP(IA, 49)
      RA(0520)=RC(IA, 520)*SSP(IA, 39)*SSP(IA, 19)
      RA(0521)=RC(IA, 521)*SSP(IA, 49)
      RA(0522)=RC(IA, 522)*SSP(IA, 40)*SSP(IA, 17)
      RA(0523)=RC(IA, 523)*SSP(IA, 49)*SSP(IA,  5)
      RA(0524)=RC(IA, 524)*SSP(IA,  6)*SSP(IA, 40)*SSP(IA, 23)
      RA(0525)=RC(IA, 525)*SSP(IA, 49)*SSP(IA,  3)
      RA(0526)=RC(IA, 526)*SSP(IA,  1)*SSP(IA, 40)*SSP(IA, 23)
      RA(0527)=RC(IA, 527)*SSP(IA, 49)*SSP(IA,  3)
      RA(0528)=RC(IA, 528)*SSP(IA, 17)*SSP(IA, 17)*SSP(IA, 14)
      RA(0529)=RC(IA, 529)*SSP(IA, 49)*SSP(IA,  3)
      RA(0530)=RC(IA, 530)*SSP(IA, 40)*SSP(IA, 19)
      RA(0531)=RC(IA, 531)*SSP(IA, 49)*SSP(IA,  3)
      RA(0532)=RC(IA, 532)*SSP(IA,  1)*SSP(IA, 17)*SSP(IA, 39)
      RA(0533)=RC(IA, 533)*SSP(IA, 49)*SSP(IA,  3)
      RA(0534)=RC(IA, 534)*SSP(IA,  1)*SSP(IA, 44)*SSP(IA, 14)
      RA(0535)=RC(IA, 535)*SSP(IA, 48)
      RA(0536)=RC(IA, 536)*SSP(IA, 39)*SSP(IA, 14)
      RA(0537)=RC(IA, 537)*SSP(IA, 48)*SSP(IA,  3)
      RA(0538)=RC(IA, 538)*SSP(IA, 17)*SSP(IA, 19)
      RA(0539)=RC(IA, 539)*SSP(IA, 48)*SSP(IA,  3)
      RA(0540)=RC(IA, 540)*SSP(IA, 40)*SSP(IA, 14)
      RA(0541)=RC(IA, 541)*SSP(IA, 48)*SSP(IA,  3)
      RA(0542)=RC(IA, 542)*SSP(IA,  1)*SSP(IA, 23)*SSP(IA, 17)
      RA(0543)=RC(IA, 543)*SSP(IA, 48)*SSP(IA,  5)
      RA(0544)=RC(IA, 544)*SSP(IA,  6)*SSP(IA, 44)*SSP(IA,  3)
      SRA(01)=SOUR(IA,1)*SSP(IA, 22)
C     SRA(02)=SOUR(IA,2)
C     SRA(03)=SOUR(IA,3)*0.5*SSP(IA, 02)
C     SRA(04)=SOUR(IA,4)*SSP(IA, 05)
C     SRA(05)=SOUR(IA,5)*SSP(IA, 04)
C     SRA(06)=SOUR(IA,6)
      SRA(01)=0.0
      SRA(02)=0.0
      SRA(03)=0.0
      SRA(04)=0.0
      SRA(05)=0.0
      SRA(06)=0.0
C-------------------------   PRODUCTION RATES  -------------------------
      PRD(  1)=WM(  1)*HFO(  1)*(
     1  +RA(0004)+RA(0006)+RA(0009)+RA(0023)+RA(0035)+RA(0051)+RA(0065)
     2  +RA(0075)+RA(0085)+RA(0087)+RA(0101)+RA(0117)+RA(0125)+RA(0127)
     3  +RA(0141)+RA(0153)+RA(0165)+RA(0175)+RA(0193)+RA(0195)+RA(0243)
     4  +RA(0261)+RA(0263)+RA(0277)+RA(0305)+RA(0307)+RA(0329)+RA(0331)
     5  +RA(0333)+RA(0401)+RA(0423)+RA(0443)+RA(0445)+RA(0479)+RA(0509)
     6  +RA(0525)+RA(0531)+RA(0533)+RA(0541)
     7  -RA(0003)-RA(0005)-RA(0010)-RA(0024)-RA(0036)-RA(0052)-RA(0066)
     8  -RA(0076)-RA(0086)-RA(0088)-RA(0102)-RA(0118)-RA(0126)-RA(0128)
     9  -RA(0142)-RA(0154)-RA(0166)-RA(0176)-RA(0194)-RA(0196)-RA(0244)
     *  -RA(0262)-RA(0264)-RA(0278)-RA(0306)-RA(0308)-RA(0330)-RA(0332)
     *  -RA(0334)-RA(0402)-RA(0424)-RA(0444)-RA(0446)-RA(0480)-RA(0510)
     *  -RA(0526)-RA(0532)-RA(0534)-RA(0542)
     *  +1.000*RA(0503)+1.000*RA(0505)+1.000*RA(0507)
     *  +SRA(01))
      PRD(  2)=WM(  2)*HFO(  2)*(
     1  +RA(0002)+RA(0013)+RA(0020)+RA(0023)+RA(0027)+RA(0029)+RA(0033)
     2  +RA(0048)+RA(0060)+RA(0072)+RA(0082)+RA(0098)+RA(0100)+RA(0112)
     3  +RA(0128)+RA(0130)+RA(0136)+RA(0150)+RA(0172)+RA(0186)+RA(0200)
     4  +RA(0202)+RA(0204)+RA(0210)+RA(0230)+RA(0232)+RA(0238)+RA(0240)
     5  +RA(0242)+RA(0250)+RA(0270)+RA(0284)+RA(0289)+RA(0318)+RA(0362)
     6  +RA(0372)+RA(0393)+RA(0398)+RA(0404)+RA(0416)+RA(0427)+RA(0440)
     7  +RA(0442)+RA(0464)+RA(0470)+RA(0494)+RA(0496)
     8  -RA(0001)-RA(0014)-RA(0019)-RA(0024)-RA(0028)-RA(0030)-RA(0034)
     9  -RA(0047)-RA(0059)-RA(0071)-RA(0081)-RA(0097)-RA(0099)-RA(0111)
     *  -RA(0127)-RA(0129)-RA(0135)-RA(0149)-RA(0171)-RA(0185)-RA(0199)
     *  -RA(0201)-RA(0203)-RA(0209)-RA(0229)-RA(0231)-RA(0237)-RA(0239)
     *  -RA(0241)-RA(0249)-RA(0269)-RA(0283)-RA(0290)-RA(0317)-RA(0361)
     *  -RA(0371)-RA(0394)-RA(0397)-RA(0403)-RA(0415)-RA(0428)-RA(0439)
     *  -RA(0441)-RA(0463)-RA(0469)-RA(0493)-RA(0495)
     *  -SRA(02))
      PRD(  3)=WM(  3)*HFO(  3)*(
     1  +RA(0002)+RA(0003)+RA(0005)+RA(0010)+RA(0010)+RA(0012)+RA(0016)
     2  +RA(0020)+RA(0022)+RA(0024)+RA(0026)+RA(0036)+RA(0038)+RA(0043)
     3  +RA(0049)+RA(0052)+RA(0055)+RA(0064)+RA(0066)+RA(0076)+RA(0086)
     4  +RA(0088)+RA(0091)+RA(0093)+RA(0103)+RA(0106)+RA(0109)+RA(0111)
     5  +RA(0118)+RA(0119)+RA(0123)+RA(0123)+RA(0129)+RA(0131)+RA(0131)
     6  +RA(0133)+RA(0137)+RA(0142)+RA(0144)+RA(0151)+RA(0154)+RA(0161)
     7  +RA(0166)+RA(0173)+RA(0176)+RA(0181)+RA(0191)+RA(0196)+RA(0197)
     8  +RA(0205)+RA(0211)+RA(0216)+RA(0224)+RA(0225)+RA(0227)+RA(0231)
     9  +RA(0233)+RA(0244)+RA(0246)+RA(0251)+RA(0262)+RA(0264)+RA(0273)
     *  +RA(0276)+RA(0278)+RA(0285)+RA(0306)+RA(0308)+RA(0330)+RA(0332)
     *  +RA(0334)+RA(0357)+RA(0366)+RA(0368)+RA(0376)+RA(0378)+RA(0383)
     *  +RA(0387)+RA(0392)+RA(0400)+RA(0402)+RA(0421)+RA(0424)+RA(0426)
     *  +RA(0434)+RA(0435)+RA(0444)+RA(0446)+RA(0462)+RA(0468)+RA(0471)
     *  +RA(0478)+RA(0480)+RA(0481)+RA(0504)+RA(0506)+RA(0508)+RA(0510)
     *  +RA(0526)+RA(0528)+RA(0530)+RA(0532)+RA(0534)+RA(0538)+RA(0540)
     *  +RA(0542)+RA(0543)
     *  -RA(0001)-RA(0004)-RA(0006)-RA(0009)-RA(0009)-RA(0011)-RA(0015)
     *  -RA(0019)-RA(0021)-RA(0023)-RA(0025)-RA(0035)-RA(0037)-RA(0044)
     *  -RA(0050)-RA(0051)-RA(0056)-RA(0063)-RA(0065)-RA(0075)-RA(0085)
     *  -RA(0087)-RA(0092)-RA(0094)-RA(0104)-RA(0105)-RA(0110)-RA(0112)
     *  -RA(0117)-RA(0120)-RA(0124)-RA(0124)-RA(0130)-RA(0132)-RA(0132)
     *  -RA(0134)-RA(0138)-RA(0141)-RA(0143)-RA(0152)-RA(0153)-RA(0162)
     *  -RA(0165)-RA(0174)-RA(0175)-RA(0182)-RA(0192)-RA(0195)-RA(0198)
     *  -RA(0206)-RA(0212)-RA(0215)-RA(0223)-RA(0226)-RA(0228)-RA(0232)
     *  -RA(0234)-RA(0243)-RA(0245)-RA(0252)-RA(0261)-RA(0263)-RA(0274)
     *  -RA(0275)-RA(0277)-RA(0286)-RA(0305)-RA(0307)-RA(0329)-RA(0331)
     *  -RA(0333)-RA(0358)-RA(0365)-RA(0367)-RA(0375)-RA(0377)-RA(0384)
     *  -RA(0388)-RA(0391)-RA(0399)-RA(0401)-RA(0422)-RA(0423)-RA(0425)
     *  -RA(0433)-RA(0436)-RA(0443)-RA(0445)-RA(0461)-RA(0467)-RA(0472)
     *  -RA(0477)-RA(0479)-RA(0482)-RA(0503)-RA(0505)-RA(0507)-RA(0509)
     *  -RA(0525)-RA(0527)-RA(0529)-RA(0531)-RA(0533)-RA(0537)-RA(0539)
     *  -RA(0541)-RA(0544)
     *  +SRA(04))
      PRD(  4)=WM(  4)*HFO(  4)*(
     1  +RA(0001)+RA(0004)+RA(0008)+RA(0014)+RA(0014)+RA(0016)+RA(0018)
     2  +RA(0025)+RA(0028)+RA(0042)+RA(0047)+RA(0054)+RA(0056)+RA(0068)
     3  +RA(0080)+RA(0092)+RA(0099)+RA(0124)+RA(0126)+RA(0134)+RA(0135)
     4  +RA(0148)+RA(0156)+RA(0168)+RA(0170)+RA(0180)+RA(0182)+RA(0201)
     5  +RA(0206)+RA(0208)+RA(0218)+RA(0220)+RA(0228)+RA(0236)+RA(0237)
     6  +RA(0266)+RA(0280)+RA(0302)+RA(0304)+RA(0336)+RA(0338)+RA(0340)
     7  +RA(0374)+RA(0386)+RA(0390)+RA(0418)+RA(0422)+RA(0448)+RA(0450)
     8  -RA(0002)-RA(0003)-RA(0007)-RA(0013)-RA(0013)-RA(0015)-RA(0017)
     9  -RA(0026)-RA(0027)-RA(0041)-RA(0048)-RA(0053)-RA(0055)-RA(0067)
     *  -RA(0079)-RA(0091)-RA(0100)-RA(0123)-RA(0125)-RA(0133)-RA(0136)
     *  -RA(0147)-RA(0155)-RA(0167)-RA(0169)-RA(0179)-RA(0181)-RA(0202)
     *  -RA(0205)-RA(0207)-RA(0217)-RA(0219)-RA(0227)-RA(0235)-RA(0238)
     *  -RA(0265)-RA(0279)-RA(0301)-RA(0303)-RA(0335)-RA(0337)-RA(0339)
     *  -RA(0373)-RA(0385)-RA(0389)-RA(0417)-RA(0421)-RA(0447)-RA(0449)
     *  -SRA(05))
      PRD(  5)=WM(  5)*HFO(  5)*(
     1  +RA(0001)+RA(0003)+RA(0006)+RA(0007)+RA(0007)+RA(0012)+RA(0015)
     2  +RA(0018)+RA(0021)+RA(0021)+RA(0027)+RA(0030)+RA(0032)+RA(0032)
     3  +RA(0037)+RA(0040)+RA(0041)+RA(0044)+RA(0045)+RA(0053)+RA(0058)
     4  +RA(0067)+RA(0070)+RA(0078)+RA(0079)+RA(0090)+RA(0095)+RA(0097)
     5  +RA(0110)+RA(0111)+RA(0120)+RA(0122)+RA(0129)+RA(0146)+RA(0147)
     6  +RA(0155)+RA(0158)+RA(0167)+RA(0178)+RA(0187)+RA(0212)+RA(0214)
     7  +RA(0219)+RA(0226)+RA(0229)+RA(0234)+RA(0245)+RA(0248)+RA(0256)
     8  +RA(0258)+RA(0260)+RA(0265)+RA(0271)+RA(0282)+RA(0283)+RA(0287)
     9  +RA(0298)+RA(0300)+RA(0301)+RA(0303)+RA(0324)+RA(0326)+RA(0328)
     *  +RA(0335)+RA(0337)+RA(0339)+RA(0354)+RA(0355)+RA(0370)+RA(0373)
     *  +RA(0379)+RA(0379)+RA(0382)+RA(0396)+RA(0410)+RA(0413)+RA(0420)
     *  +RA(0429)+RA(0447)+RA(0449)+RA(0452)+RA(0454)+RA(0482)+RA(0512)
     *  +RA(0514)+RA(0516)+RA(0518)+RA(0524)+RA(0544)
     *  -RA(0002)-RA(0004)-RA(0005)-RA(0008)-RA(0008)-RA(0011)-RA(0016)
     *  -RA(0017)-RA(0022)-RA(0022)-RA(0028)-RA(0029)-RA(0031)-RA(0031)
     *  -RA(0038)-RA(0039)-RA(0042)-RA(0043)-RA(0046)-RA(0054)-RA(0057)
     *  -RA(0068)-RA(0069)-RA(0077)-RA(0080)-RA(0089)-RA(0096)-RA(0098)
     *  -RA(0109)-RA(0112)-RA(0119)-RA(0121)-RA(0130)-RA(0145)-RA(0148)
     *  -RA(0156)-RA(0157)-RA(0168)-RA(0177)-RA(0188)-RA(0211)-RA(0213)
     *  -RA(0220)-RA(0225)-RA(0230)-RA(0233)-RA(0246)-RA(0247)-RA(0255)
     *  -RA(0257)-RA(0259)-RA(0266)-RA(0272)-RA(0281)-RA(0284)-RA(0288)
     *  -RA(0297)-RA(0299)-RA(0302)-RA(0304)-RA(0323)-RA(0325)-RA(0327)
     *  -RA(0336)-RA(0338)-RA(0340)-RA(0353)-RA(0356)-RA(0369)-RA(0374)
     *  -RA(0380)-RA(0380)-RA(0381)-RA(0395)-RA(0409)-RA(0414)-RA(0419)
     *  -RA(0430)-RA(0448)-RA(0450)-RA(0451)-RA(0453)-RA(0481)-RA(0511)
     *  -RA(0513)-RA(0515)-RA(0517)-RA(0523)-RA(0543)
     *  -SRA(04))
      PRD(  6)=WM(  6)*HFO(  6)*(
     1  +RA(0005)+RA(0008)+RA(0011)+RA(0025)+RA(0029)+RA(0037)+RA(0039)
     2  +RA(0057)+RA(0069)+RA(0077)+RA(0089)+RA(0121)+RA(0138)+RA(0143)
     3  +RA(0145)+RA(0157)+RA(0177)+RA(0213)+RA(0247)+RA(0257)+RA(0259)
     4  +RA(0281)+RA(0297)+RA(0299)+RA(0321)+RA(0323)+RA(0325)+RA(0327)
     5  +RA(0367)+RA(0369)+RA(0375)+RA(0381)+RA(0395)+RA(0409)+RA(0419)
     6  +RA(0451)+RA(0453)+RA(0511)+RA(0523)+RA(0543)
     7  -RA(0006)-RA(0007)-RA(0012)-RA(0026)-RA(0030)-RA(0038)-RA(0040)
     8  -RA(0058)-RA(0070)-RA(0078)-RA(0090)-RA(0122)-RA(0137)-RA(0144)
     9  -RA(0146)-RA(0158)-RA(0178)-RA(0214)-RA(0248)-RA(0258)-RA(0260)
     *  -RA(0282)-RA(0298)-RA(0300)-RA(0322)-RA(0324)-RA(0326)-RA(0328)
     *  -RA(0368)-RA(0370)-RA(0376)-RA(0382)-RA(0396)-RA(0410)-RA(0420)
     *  -RA(0452)-RA(0454)-RA(0512)-RA(0524)-RA(0544)
     *  +1.000*RA(0513)+1.000*RA(0515)+1.000*RA(0517)
     *  )
      PRD(  7)=WM(  7)*HFO(  7)*(
     1  +RA(0017)+RA(0019)+RA(0022)+RA(0024)+RA(0026)+RA(0028)+RA(0030)
     2  +RA(0034)+RA(0034)+RA(0035)+RA(0039)+RA(0041)+RA(0046)+RA(0059)
     3  +RA(0071)+RA(0074)+RA(0081)+RA(0084)+RA(0096)+RA(0149)+RA(0164)
     4  +RA(0171)+RA(0185)+RA(0188)+RA(0190)+RA(0203)+RA(0249)+RA(0268)
     5  +RA(0269)+RA(0288)+RA(0290)+RA(0314)+RA(0316)+RA(0317)+RA(0348)
     6  +RA(0350)+RA(0352)+RA(0356)+RA(0361)+RA(0371)+RA(0380)+RA(0394)
     7  +RA(0403)+RA(0414)+RA(0428)+RA(0430)+RA(0439)+RA(0441)+RA(0456)
     8  +RA(0458)+RA(0463)+RA(0469)+RA(0493)+RA(0495)+RA(0498)+RA(0500)
     9  -RA(0018)-RA(0020)-RA(0021)-RA(0023)-RA(0025)-RA(0027)-RA(0029)
     *  -RA(0033)-RA(0033)-RA(0036)-RA(0040)-RA(0042)-RA(0045)-RA(0060)
     *  -RA(0072)-RA(0073)-RA(0082)-RA(0083)-RA(0095)-RA(0150)-RA(0163)
     *  -RA(0172)-RA(0186)-RA(0187)-RA(0189)-RA(0204)-RA(0250)-RA(0267)
     *  -RA(0270)-RA(0287)-RA(0289)-RA(0313)-RA(0315)-RA(0318)-RA(0347)
     *  -RA(0349)-RA(0351)-RA(0355)-RA(0362)-RA(0372)-RA(0379)-RA(0393)
     *  -RA(0404)-RA(0413)-RA(0427)-RA(0429)-RA(0440)-RA(0442)-RA(0455)
     *  -RA(0457)-RA(0464)-RA(0470)-RA(0494)-RA(0496)-RA(0497)-RA(0499)
     *  )
      PRD(  8)=WM(  8)*HFO(  8)*(
     1  +RA(0031)+RA(0033)+RA(0036)+RA(0038)+RA(0040)+RA(0042)+RA(0073)
     2  +RA(0083)+RA(0163)+RA(0189)+RA(0267)+RA(0313)+RA(0315)+RA(0347)
     3  +RA(0349)+RA(0351)+RA(0455)+RA(0457)+RA(0497)+RA(0499)
     4  -RA(0032)-RA(0034)-RA(0035)-RA(0037)-RA(0039)-RA(0041)-RA(0074)
     5  -RA(0084)-RA(0164)-RA(0190)-RA(0268)-RA(0314)-RA(0316)-RA(0348)
     6  -RA(0350)-RA(0352)-RA(0456)-RA(0458)-RA(0498)-RA(0500)
     *  )
      PRD(  9)=WM(  9)*HFO(  9)*(
     1  +RA(0044)+RA(0046)+RA(0048)+RA(0049)+RA(0051)+RA(0053)+RA(0057)
     2  +RA(0059)+RA(0061)+RA(0111)+RA(0113)+RA(0123)+RA(0125)+RA(0129)
     3  +RA(0133)+RA(0139)+RA(0189)+RA(0207)+RA(0209)+RA(0215)+RA(0221)
     4  +RA(0223)+RA(0225)+RA(0227)+RA(0227)+RA(0229)+RA(0229)+RA(0231)
     5  +RA(0235)+RA(0241)+RA(0255)+RA(0283)+RA(0285)+RA(0291)+RA(0295)
     6  +RA(0364)+RA(0385)+RA(0411)+RA(0413)+RA(0415)
     7  -RA(0043)-RA(0045)-RA(0047)-RA(0050)-RA(0052)-RA(0054)-RA(0058)
     8  -RA(0060)-RA(0062)-RA(0112)-RA(0114)-RA(0124)-RA(0126)-RA(0130)
     9  -RA(0134)-RA(0140)-RA(0190)-RA(0208)-RA(0210)-RA(0216)-RA(0222)
     *  -RA(0224)-RA(0226)-RA(0228)-RA(0228)-RA(0230)-RA(0230)-RA(0232)
     *  -RA(0236)-RA(0242)-RA(0256)-RA(0284)-RA(0286)-RA(0292)-RA(0296)
     *  -RA(0363)-RA(0386)-RA(0412)-RA(0414)-RA(0416)
     *  +SRA(03)+SRA(04)+SRA(05))
      PRD( 10)=WM( 10)*HFO( 10)*(
     1  +RA(0043)+RA(0045)+RA(0047)+RA(0055)+RA(0114)+RA(0127)+RA(0140)
     2  +RA(0217)+RA(0231)+RA(0239)+RA(0363)
     3  -RA(0044)-RA(0046)-RA(0048)-RA(0056)-RA(0113)-RA(0128)-RA(0139)
     4  -RA(0218)-RA(0232)-RA(0240)-RA(0364)
     *  )
      PRD( 11)=WM( 11)*HFO( 11)*(
     1  +RA(0050)+RA(0052)+RA(0054)+RA(0056)+RA(0058)+RA(0060)+RA(0062)
     2  +RA(0064)+RA(0065)+RA(0067)+RA(0069)+RA(0071)+RA(0073)+RA(0135)
     3  +RA(0139)+RA(0179)+RA(0199)+RA(0225)+RA(0241)+RA(0275)+RA(0279)
     4  +RA(0287)+RA(0293)+RA(0397)+RA(0412)+RA(0415)+RA(0417)
     5  -RA(0049)-RA(0051)-RA(0053)-RA(0055)-RA(0057)-RA(0059)-RA(0061)
     6  -RA(0063)-RA(0066)-RA(0068)-RA(0070)-RA(0072)-RA(0074)-RA(0136)
     7  -RA(0140)-RA(0180)-RA(0200)-RA(0226)-RA(0242)-RA(0276)-RA(0280)
     8  -RA(0288)-RA(0294)-RA(0398)-RA(0411)-RA(0416)-RA(0418)
     *  )
      PRD( 12)=WM( 12)*HFO( 12)*(
     1  +RA(0063)+RA(0066)+RA(0068)+RA(0070)+RA(0072)+RA(0074)+RA(0091)
     2  +RA(0097)+RA(0109)+RA(0113)+RA(0119)+RA(0137)+RA(0141)+RA(0145)
     3  +RA(0147)+RA(0149)+RA(0151)+RA(0169)+RA(0199)+RA(0209)+RA(0243)
     4  +RA(0247)+RA(0249)+RA(0251)+RA(0279)+RA(0283)+RA(0287)+RA(0359)
     5  +RA(0429)
     6  -RA(0064)-RA(0065)-RA(0067)-RA(0069)-RA(0071)-RA(0073)-RA(0092)
     7  -RA(0098)-RA(0110)-RA(0114)-RA(0120)-RA(0138)-RA(0142)-RA(0146)
     8  -RA(0148)-RA(0150)-RA(0152)-RA(0170)-RA(0200)-RA(0210)-RA(0244)
     9  -RA(0248)-RA(0250)-RA(0252)-RA(0280)-RA(0284)-RA(0288)-RA(0360)
     *  -RA(0430)
     *  )
      PRD( 13)=WM( 13)*HFO( 13)*(
     1  +RA(0061)+RA(0076)+RA(0078)+RA(0080)+RA(0082)+RA(0084)+RA(0105)
     2  +RA(0159)+RA(0309)+RA(0311)+RA(0341)+RA(0343)+RA(0345)+RA(0405)
     3  +RA(0483)
     4  -RA(0062)-RA(0075)-RA(0077)-RA(0079)-RA(0081)-RA(0083)-RA(0106)
     5  -RA(0160)-RA(0310)-RA(0312)-RA(0342)-RA(0344)-RA(0346)-RA(0406)
     6  -RA(0484)
     *  )
      PRD( 14)=WM( 14)*HFO( 14)*(
     1  +RA(0062)+RA(0075)+RA(0077)+RA(0079)+RA(0081)+RA(0083)+RA(0086)
     2  +RA(0088)+RA(0090)+RA(0092)+RA(0094)+RA(0096)+RA(0098)+RA(0100)
     3  +RA(0102)+RA(0102)+RA(0104)+RA(0104)+RA(0106)+RA(0108)+RA(0108)
     4  +RA(0160)+RA(0169)+RA(0179)+RA(0189)+RA(0215)+RA(0222)+RA(0245)
     5  +RA(0271)+RA(0275)+RA(0286)+RA(0291)+RA(0293)+RA(0295)+RA(0310)
     6  +RA(0312)+RA(0319)+RA(0342)+RA(0344)+RA(0346)+RA(0359)+RA(0365)
     7  +RA(0377)+RA(0388)+RA(0389)+RA(0406)+RA(0408)+RA(0415)+RA(0421)
     8  +RA(0432)+RA(0433)+RA(0436)+RA(0437)+RA(0465)+RA(0484)+RA(0486)
     9  +RA(0527)+RA(0533)+RA(0535)+RA(0539)
     *  -RA(0061)-RA(0076)-RA(0078)-RA(0080)-RA(0082)-RA(0084)-RA(0085)
     *  -RA(0087)-RA(0089)-RA(0091)-RA(0093)-RA(0095)-RA(0097)-RA(0099)
     *  -RA(0101)-RA(0101)-RA(0103)-RA(0103)-RA(0105)-RA(0107)-RA(0107)
     *  -RA(0159)-RA(0170)-RA(0180)-RA(0190)-RA(0216)-RA(0221)-RA(0246)
     *  -RA(0272)-RA(0276)-RA(0285)-RA(0292)-RA(0294)-RA(0296)-RA(0309)
     *  -RA(0311)-RA(0320)-RA(0341)-RA(0343)-RA(0345)-RA(0360)-RA(0366)
     *  -RA(0378)-RA(0387)-RA(0390)-RA(0405)-RA(0407)-RA(0416)-RA(0422)
     *  -RA(0431)-RA(0434)-RA(0435)-RA(0438)-RA(0466)-RA(0483)-RA(0485)
     *  -RA(0528)-RA(0534)-RA(0536)-RA(0540)
     *  +1.000*RA(0501)+1.000*RA(0505)+1.000*RA(0507)+1.000*RA(0513)
     *  +1.000*RA(0517)
     *  )
      PRD( 15)=WM( 15)*HFO( 15)*(
     1  +RA(0085)+RA(0094)+RA(0115)+RA(0118)+RA(0120)+RA(0122)+RA(0124)
     2  +RA(0126)+RA(0128)+RA(0130)+RA(0132)+RA(0132)+RA(0207)+RA(0217)
     3  -RA(0086)-RA(0093)-RA(0116)-RA(0117)-RA(0119)-RA(0121)-RA(0123)
     4  -RA(0125)-RA(0127)-RA(0129)-RA(0131)-RA(0131)-RA(0208)-RA(0218)
     *  )
      PRD( 16)=WM( 16)*HFO( 16)*(
     1  +RA(0087)+RA(0089)+RA(0110)+RA(0112)+RA(0114)+RA(0116)+RA(0143)
     2  +RA(0223)
     3  -RA(0088)-RA(0090)-RA(0109)-RA(0111)-RA(0113)-RA(0115)-RA(0144)
     4  -RA(0224)
     *  )
      PRD( 17)=WM( 17)*HFO( 17)*(
     1  +RA(0093)+RA(0101)+RA(0165)+RA(0167)+RA(0171)+RA(0173)+RA(0176)
     2  +RA(0178)+RA(0180)+RA(0182)+RA(0184)+RA(0184)+RA(0186)+RA(0188)
     3  +RA(0192)+RA(0194)+RA(0321)+RA(0354)+RA(0367)+RA(0375)+RA(0385)
     4  +RA(0433)+RA(0465)+RA(0477)+RA(0489)+RA(0521)+RA(0527)+RA(0527)
     5  +RA(0531)+RA(0537)+RA(0541)
     6  -RA(0094)-RA(0102)-RA(0166)-RA(0168)-RA(0172)-RA(0174)-RA(0175)
     7  -RA(0177)-RA(0179)-RA(0181)-RA(0183)-RA(0183)-RA(0185)-RA(0187)
     8  -RA(0191)-RA(0193)-RA(0322)-RA(0353)-RA(0368)-RA(0376)-RA(0386)
     9  -RA(0434)-RA(0466)-RA(0478)-RA(0490)-RA(0522)-RA(0528)-RA(0528)
     *  -RA(0532)-RA(0538)-RA(0542)
     *  +1.000*RA(0501)+1.000*RA(0501)+1.000*RA(0503)+1.000*RA(0505)
     *  +1.000*RA(0507)+1.000*RA(0507)+1.000*RA(0507)+1.000*RA(0513)
     *  +1.000*RA(0515)+1.000*RA(0517)+1.000*RA(0517)+1.000*RA(0517)
     *
     *  )
      PRD( 18)=WM( 18)*HFO( 18)*(
     1  +RA(0095)+RA(0099)+RA(0142)+RA(0144)+RA(0146)+RA(0148)+RA(0150)
     2  +RA(0152)+RA(0254)+RA(0259)+RA(0263)
     3  -RA(0096)-RA(0100)-RA(0141)-RA(0143)-RA(0145)-RA(0147)-RA(0149)
     4  -RA(0151)-RA(0253)-RA(0260)-RA(0264)
     *  )
      PRD( 19)=WM( 19)*HFO( 19)*(
     1  +RA(0103)+RA(0153)+RA(0155)+RA(0157)+RA(0159)+RA(0161)+RA(0163)
     2  +RA(0166)+RA(0168)+RA(0170)+RA(0172)+RA(0174)+RA(0183)+RA(0221)
     3  +RA(0285)+RA(0356)+RA(0363)+RA(0417)+RA(0437)+RA(0509)+RA(0511)
     4  +RA(0519)+RA(0529)+RA(0537)
     5  -RA(0104)-RA(0154)-RA(0156)-RA(0158)-RA(0160)-RA(0162)-RA(0164)
     6  -RA(0165)-RA(0167)-RA(0169)-RA(0171)-RA(0173)-RA(0184)-RA(0222)
     7  -RA(0286)-RA(0355)-RA(0364)-RA(0418)-RA(0438)-RA(0510)-RA(0512)
     8  -RA(0520)-RA(0530)-RA(0538)
     *  +1.000*RA(0501)+1.000*RA(0503)+1.000*RA(0515)
     *  )
      PRD( 20)=WM( 20)*HFO( 20)*(
     1  +RA(0107)+RA(0154)+RA(0156)+RA(0158)+RA(0160)+RA(0162)+RA(0164)
     2  -RA(0108)-RA(0153)-RA(0155)-RA(0157)-RA(0159)-RA(0161)-RA(0163)
     *  )
      PRD( 21)=WM( 21)*HFO( 21)*(
     1  +RA(0117)+RA(0121)+RA(0134)+RA(0136)+RA(0138)+RA(0140)+RA(0235)
     2  +RA(0239)
     3  -RA(0118)-RA(0122)-RA(0133)-RA(0135)-RA(0137)-RA(0139)-RA(0236)
     4  -RA(0240)
     *  )
      PRD( 22)=WM( 22)*HFO( 22)*(
     1  +RA(0131)+RA(0193)+RA(0195)+RA(0197)+RA(0203)+RA(0206)+RA(0208)
     2  +RA(0210)+RA(0212)+RA(0214)+RA(0388)+RA(0408)+RA(0471)+RA(0479)
     3  +RA(0483)+RA(0487)+RA(0493)+RA(0497)
     4  -RA(0132)-RA(0194)-RA(0196)-RA(0198)-RA(0204)-RA(0205)-RA(0207)
     5  -RA(0209)-RA(0211)-RA(0213)-RA(0387)-RA(0407)-RA(0472)-RA(0480)
     6  -RA(0484)-RA(0488)-RA(0494)-RA(0498)
     *  -SRA(01))
      PRD( 23)=WM( 23)*HFO( 23)*(
     1  +RA(0175)+RA(0177)+RA(0183)+RA(0185)+RA(0191)+RA(0196)+RA(0198)
     2  +RA(0200)+RA(0202)+RA(0204)+RA(0413)+RA(0429)+RA(0432)+RA(0436)
     3  +RA(0471)+RA(0473)+RA(0473)+RA(0476)+RA(0476)+RA(0477)+RA(0479)
     4  +RA(0483)+RA(0491)+RA(0495)+RA(0499)+RA(0523)+RA(0525)+RA(0541)
     5  -RA(0176)-RA(0178)-RA(0184)-RA(0186)-RA(0192)-RA(0195)-RA(0197)
     6  -RA(0199)-RA(0201)-RA(0203)-RA(0414)-RA(0430)-RA(0431)-RA(0435)
     7  -RA(0472)-RA(0474)-RA(0474)-RA(0475)-RA(0475)-RA(0478)-RA(0480)
     8  -RA(0484)-RA(0492)-RA(0496)-RA(0500)-RA(0524)-RA(0526)-RA(0542)
     *  )
      PRD( 24)=WM( 24)*HFO( 24)*(
     1  +RA(0181)+RA(0201)+RA(0274)+RA(0276)+RA(0278)+RA(0280)+RA(0282)
     2  +RA(0284)+RA(0286)+RA(0288)+RA(0290)+RA(0292)+RA(0299)+RA(0303)
     3  +RA(0307)+RA(0311)+RA(0315)
     4  -RA(0182)-RA(0202)-RA(0273)-RA(0275)-RA(0277)-RA(0279)-RA(0281)
     5  -RA(0283)-RA(0285)-RA(0287)-RA(0289)-RA(0291)-RA(0300)-RA(0304)
     6  -RA(0308)-RA(0312)-RA(0316)
     *  )
      PRD( 25)=WM( 25)*HFO( 25)*(
     1  +RA(0187)+RA(0190)
     2  -RA(0188)-RA(0189)
     *  )
      PRD( 26)=WM( 26)*HFO( 26)*(
     1  +RA(0211)+RA(0216)+RA(0218)+RA(0220)+RA(0222)+RA(0256)+RA(0273)
     2  +RA(0277)+RA(0281)+RA(0397)+RA(0421)
     3  -RA(0212)-RA(0215)-RA(0217)-RA(0219)-RA(0221)-RA(0255)-RA(0274)
     4  -RA(0278)-RA(0282)-RA(0398)-RA(0422)
     *  )
      PRD( 27)=WM( 27)*HFO( 27)*(
     1  +RA(0205)+RA(0219)+RA(0224)+RA(0226)+RA(0228)+RA(0230)+RA(0232)
     2  +RA(0233)+RA(0237)+RA(0389)
     3  -RA(0206)-RA(0220)-RA(0223)-RA(0225)-RA(0227)-RA(0229)-RA(0231)
     4  -RA(0234)-RA(0238)-RA(0390)
     *  )
      PRD( 28)=WM( 28)*HFO( 28)*(
     1  +RA(0213)+RA(0234)+RA(0236)+RA(0238)+RA(0240)+RA(0242)
     2  -RA(0214)-RA(0233)-RA(0235)-RA(0237)-RA(0239)-RA(0241)
     *  )
      PRD( 29)=WM( 29)*HFO( 29)*(
     1  +RA(0244)+RA(0246)+RA(0248)+RA(0250)+RA(0252)+RA(0253)+RA(0255)
     2  +RA(0257)+RA(0261)+RA(0265)+RA(0267)+RA(0269)+RA(0319)+RA(0365)
     3  +RA(0377)
     4  -RA(0243)-RA(0245)-RA(0247)-RA(0249)-RA(0251)-RA(0254)-RA(0256)
     5  -RA(0258)-RA(0262)-RA(0266)-RA(0268)-RA(0270)-RA(0320)-RA(0366)
     6  -RA(0378)
     *  )
      PRD( 30)=WM( 30)*HFO( 30)*(
     1  +RA(0258)+RA(0260)+RA(0262)+RA(0264)+RA(0266)+RA(0268)+RA(0270)
     2  +RA(0272)
     3  -RA(0257)-RA(0259)-RA(0261)-RA(0263)-RA(0265)-RA(0267)-RA(0269)
     4  -RA(0271)
     *  )
      PRD( 31)=WM( 31)*HFO( 31)*(
     1  +RA(0320)+RA(0322)+RA(0324)+RA(0326)+RA(0328)+RA(0330)+RA(0332)
     2  +RA(0334)+RA(0336)+RA(0338)+RA(0340)+RA(0342)+RA(0344)+RA(0346)
     3  +RA(0348)+RA(0350)+RA(0352)
     4  -RA(0319)-RA(0321)-RA(0323)-RA(0325)-RA(0327)-RA(0329)-RA(0331)
     5  -RA(0333)-RA(0335)-RA(0337)-RA(0339)-RA(0341)-RA(0343)-RA(0345)
     6  -RA(0347)-RA(0349)-RA(0351)
     *  )
      PRD( 32)=WM( 32)*HFO( 32)*(
     1  +RA(0289)+RA(0294)+RA(0298)+RA(0300)+RA(0302)+RA(0304)+RA(0306)
     2  +RA(0308)+RA(0310)+RA(0312)+RA(0314)+RA(0316)+RA(0318)+RA(0357)
     3  +RA(0361)+RA(0369)+RA(0371)+RA(0373)+RA(0379)+RA(0381)+RA(0383)
     4  -RA(0290)-RA(0293)-RA(0297)-RA(0299)-RA(0301)-RA(0303)-RA(0305)
     5  -RA(0307)-RA(0309)-RA(0311)-RA(0313)-RA(0315)-RA(0317)-RA(0358)
     6  -RA(0362)-RA(0370)-RA(0372)-RA(0374)-RA(0380)-RA(0382)-RA(0384)
     *  )
      PRD( 33)=WM( 33)*HFO( 33)*(
     1  +RA(0325)+RA(0331)+RA(0337)+RA(0343)+RA(0347)+RA(0372)+RA(0374)
     2  +RA(0376)+RA(0378)+RA(0380)+RA(0382)+RA(0384)
     3  -RA(0326)-RA(0332)-RA(0338)-RA(0344)-RA(0348)-RA(0371)-RA(0373)
     4  -RA(0375)-RA(0377)-RA(0379)-RA(0381)-RA(0383)
     *  )
      PRD( 34)=WM( 34)*HFO( 34)*(
     1  +RA(0323)+RA(0329)+RA(0335)+RA(0341)+RA(0349)+RA(0353)
     2  -RA(0324)-RA(0330)-RA(0336)-RA(0342)-RA(0350)-RA(0354)
     *  )
      PRD( 35)=WM( 35)*HFO( 35)*(
     1  +RA(0296)+RA(0297)+RA(0301)+RA(0305)+RA(0309)+RA(0313)+RA(0317)
     2  -RA(0295)-RA(0298)-RA(0302)-RA(0306)-RA(0310)-RA(0314)-RA(0318)
     *  )
      PRD( 36)=WM( 36)*HFO( 36)*(
     1  +RA(0327)+RA(0333)+RA(0339)+RA(0345)+RA(0351)+RA(0355)+RA(0358)
     2  +RA(0360)+RA(0362)+RA(0364)+RA(0366)+RA(0368)+RA(0370)
     3  -RA(0328)-RA(0334)-RA(0340)-RA(0346)-RA(0352)-RA(0356)-RA(0357)
     4  -RA(0359)-RA(0361)-RA(0363)-RA(0365)-RA(0367)-RA(0369)
     *  )
      PRD( 37)=WM( 37)*HFO( 37)*(
     1  +RA(0386)+RA(0387)+RA(0390)+RA(0391)+RA(0393)+RA(0396)+RA(0400)
     2  +RA(0401)+RA(0403)+RA(0405)+RA(0409)+RA(0411)+RA(0416)+RA(0489)
     3  +RA(0495)+RA(0499)
     4  -RA(0385)-RA(0388)-RA(0389)-RA(0392)-RA(0394)-RA(0395)-RA(0399)
     5  -RA(0402)-RA(0404)-RA(0406)-RA(0410)-RA(0412)-RA(0415)-RA(0490)
     6  -RA(0496)-RA(0500)
     *  )
      PRD( 38)=WM( 38)*HFO( 38)*(
     1  +RA(0392)+RA(0394)+RA(0395)+RA(0398)+RA(0412)+RA(0414)+RA(0486)
     2  -RA(0391)-RA(0393)-RA(0396)-RA(0397)-RA(0411)-RA(0413)-RA(0485)
     *  )
      PRD( 39)=WM( 39)*HFO( 39)*(
     1  +RA(0399)+RA(0402)+RA(0404)+RA(0406)+RA(0407)+RA(0410)+RA(0419)
     2  +RA(0423)+RA(0426)+RA(0428)+RA(0430)+RA(0435)+RA(0481)+RA(0491)
     3  +RA(0493)+RA(0497)+RA(0519)+RA(0531)+RA(0535)
     4  -RA(0400)-RA(0401)-RA(0403)-RA(0405)-RA(0408)-RA(0409)-RA(0420)
     5  -RA(0424)-RA(0425)-RA(0427)-RA(0429)-RA(0436)-RA(0482)-RA(0492)
     6  -RA(0494)-RA(0498)-RA(0520)-RA(0532)-RA(0536)
     *  )
      PRD( 40)=WM( 40)*HFO( 40)*(
     1  +RA(0418)+RA(0420)+RA(0422)+RA(0424)+RA(0425)+RA(0427)+RA(0431)
     2  +RA(0434)+RA(0462)+RA(0463)+RA(0468)+RA(0469)+RA(0487)+RA(0521)
     3  +RA(0523)+RA(0525)+RA(0529)+RA(0539)
     4  -RA(0417)-RA(0419)-RA(0421)-RA(0423)-RA(0426)-RA(0428)-RA(0432)
     5  -RA(0433)-RA(0461)-RA(0464)-RA(0467)-RA(0470)-RA(0488)-RA(0522)
     6  -RA(0524)-RA(0526)-RA(0530)-RA(0540)
     *  +1.000*RA(0503)+1.000*RA(0515)
     *  )
      PRD( 41)=WM( 41)*HFO( 41)*(
     1  +RA(0438)+RA(0440)+RA(0442)+RA(0444)+RA(0446)+RA(0448)+RA(0450)
     2  +RA(0452)+RA(0454)+RA(0456)+RA(0458)
     3  -RA(0437)-RA(0439)-RA(0441)-RA(0443)-RA(0445)-RA(0447)-RA(0449)
     4  -RA(0451)-RA(0453)-RA(0455)-RA(0457)
     *  )
      PRD( 42)=WM( 42)*HFO( 42)*(
     1  +RA(0439)+RA(0443)+RA(0447)+RA(0453)+RA(0455)+RA(0460)+RA(0461)
     2  +RA(0464)
     3  -RA(0440)-RA(0444)-RA(0448)-RA(0454)-RA(0456)-RA(0459)-RA(0462)
     4  -RA(0463)
     *  )
      PRD( 43)=WM( 43)*HFO( 43)*(
     1  +RA(0441)+RA(0445)+RA(0449)+RA(0451)+RA(0457)+RA(0459)+RA(0466)
     2  +RA(0467)+RA(0470)
     3  -RA(0442)-RA(0446)-RA(0450)-RA(0452)-RA(0458)-RA(0460)-RA(0465)
     4  -RA(0468)-RA(0469)
     *  )
      PRD( 44)=WM( 44)*HFO( 44)*(
     1  +RA(0472)+RA(0474)+RA(0475)+RA(0478)+RA(0480)+RA(0482)+RA(0484)
     2  +RA(0485)+RA(0533)+RA(0543)
     3  -RA(0471)-RA(0473)-RA(0476)-RA(0477)-RA(0479)-RA(0481)-RA(0483)
     4  -RA(0486)-RA(0534)-RA(0544)
     *  )
      PRD( 45)=WM( 45)*HFO( 45)*(
     1  +RA(0481)
     2  -RA(0482)
     *  )
      PRD( 46)=WM( 46)*HFO( 46)*(
     1  +RA(0488)+RA(0490)+RA(0492)+RA(0494)+RA(0496)+RA(0498)+RA(0500)
     2  -RA(0487)-RA(0489)-RA(0491)-RA(0493)-RA(0495)-RA(0497)-RA(0499)
     *  )
      PRD( 47)=WM( 47)*HFO( 47)*(
     1  +RA(0502)+RA(0504)+RA(0506)+RA(0508)+RA(0510)+RA(0512)+RA(0514)
     2  +RA(0516)+RA(0518)
     3  -RA(0501)-RA(0503)-RA(0505)-RA(0507)-RA(0509)-RA(0511)-RA(0513)
     4  -RA(0515)-RA(0517)
     *  )
      PRD( 48)=WM( 48)*HFO( 48)*(
     1  +RA(0536)+RA(0538)+RA(0540)+RA(0542)+RA(0544)
     2  -RA(0535)-RA(0537)-RA(0539)-RA(0541)-RA(0543)
     *  +1.000*RA(0505)+1.000*RA(0513)
     *  )
      PRD( 49)=WM( 49)*HFO( 49)*(
     1  +RA(0509)+RA(0511)+RA(0520)+RA(0522)+RA(0524)+RA(0526)+RA(0528)
     2  +RA(0530)+RA(0532)+RA(0534)
     3  -RA(0510)-RA(0512)-RA(0519)-RA(0521)-RA(0523)-RA(0525)-RA(0527)
     4  -RA(0529)-RA(0531)-RA(0533)
     *  )
      PRD(50)=0.0
      PRD(51)=0.0
      PRD(52)=0.0
      PROD=0.0
      HRR=0.0
      DO 416 ISP=1,LSP
      PROD=PROD+PRD(ISP)
      HRR=HRR+PRD(ISP)*HENTH(J,I,ISP)/HFO(ISP)
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
      QDSP(J,I,ISP)=-RCONV*PRD(ISP)/HFO(ISP)
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
  416 CONTINUE
      RHSZ(J,I)=RHSZ(J,I)+RHO(J,I)*HT(J,I)-DT*PROD
      QDOT(J,I)=-HRR*RCONV
  415 CONTINUE
      ENDDO
  414 CONTINUE
C--------------------    SOLVE ENTHALPY EQUATION    --------------------
      SOR1=RELXH
      SOR2=(1.0-SOR1)
      DO ISORP=1,ISOR2
      RSDZ=0.0
C----------------------   POINT RELAXATION SCHEME  ----------------------
      DO 420 J=2,LJ-1
      DO 422 I=2,LI-1
      IF(ISKIP(J,I).GE.3) GO TO 422
      PPP=RHO(J,I)-(ZE(J,I)+ZW(J,I)+ZN(J,I)+ZS(J,I))
      FNEW=SOR2*FPZ(J,I)+(RHSZ(J,I)-ZN(J,I)*FPZ(J+1,I)
     1   -ZS(J,I)*FPZ(J-1,I)-ZE(J,I)*FPZ(J,I+1)
     2   -ZW(J,I)*FPZ(J,I-1))*SOR1/PPP
      RSDZ=RSDZ+DABS(FPZ(J,I)-FNEW)
      FPZ(J,I)=FNEW
  422 CONTINUE
  420 CONTINUE
      IF((RSDZ/HNORM).LE.TOLRH) GO TO 401
      ENDDO
      ISOR2=ISORP-1
  401 CONTINUE
      ISOR2=ISORP
C--------------------              OVER             --------------------
C--------------------    CALCULATING TEMPERATURE    --------------------
      IOVER=1
      DO I=2,LI-1
      DO J=2,LJ-1
      IF(ISKIP(J,I).GE.3) GO TO 511
      IA=(I-1)*LJ+J
      ENTH=(FPZ(J,I)
     1    -0.5*(U(J,I)*U(J,I)+V(J,I)*V(J,I)+W(J,I)*W(J,I)))*HSTR
      DO ISP=1,LSP
      ENTH=ENTH+SSP(IA,ISP)*HFO(ISP)*HSTR
      ENDDO
      IPL=1
      TPOLY=TPOL1
      TKD=TK(J,I)*TSTR
  514 CONTINUE
      IF(TKD.GT.TPOLY) TKD=TPOLY/2.0
      A1=0.0
      A2=0.0
      A3=0.0
      A4=0.0
      A5=0.0
      A6=0.0
      DO ISP=1,LSP
      A1=A1+POLSP(IPL,ISP)*SSP(IA,ISP)/WM(ISP)
      A2=A2+POLSP(IPL+1,ISP)*SSP(IA,ISP)/WM(ISP)/2.0
      A3=A3+POLSP(IPL+2,ISP)*SSP(IA,ISP)/WM(ISP)/3.0
      A4=A4+POLSP(IPL+3,ISP)*SSP(IA,ISP)/WM(ISP)/4.0
      A5=A5+POLSP(IPL+4,ISP)*SSP(IA,ISP)/WM(ISP)/5.0
      A6=A6+POLSP(IPL+5,ISP)*SSP(IA,ISP)/WM(ISP)
      ENDDO
      FUNT=((A1+(A2+(A3+(A4+A5*TPOLY)*TPOLY)*TPOLY)*TPOLY)*TPOLY
     1    +A6)*GASC/WMSTR-ENTH
      IF(FUNT.GE.1.0D-07) GO TO 516
      IF(IPL.GE.8) THEN
                   IF(IOVER.LE.0) WRITE(*,712) I,HT(J,I)
                   IOVER=1
                   GO TO 516
                   END IF
      IPL=8
      IF(TKD.LE.TPOLY) TKD=TPOL2/2.0
      TPOLY=TPOL2
      GO TO 514
  516 CONTINUE
      DO INEWT=1,100
      DTKD=((A1+(A2+(A3+(A4+A5*TKD)*TKD)*TKD)*TKD)*TKD+A6-ENTH*WMSTR
     1   /GASC)/(A1+(2.0*A2+(3.0*A3+(4.0*A4+5.0*A5*TKD)*TKD)*TKD)*TKD)
      IF(DABS(DTKD).GT.0.001.AND.INEWT.EQ.100) WRITE(*,9427) 
     1   I,J,A1,A2,A3,A4,A5,A6,TKD,ENTH,WMSTR,GASC
      
      
 9427 FORMAT('DEBUG(9427):', 5X,I4, 2X,I4, 4x, 
     1   g0,2x,g0,2x,g0,2x,g0,2x,g0,2x,g0,2x,g0,2x,g0,2x,g0,2x,g0,2x)
      IF(DABS(DTKD).LE.0.001) GO TO 522
      TKD=TKD-DTKD
      ENDDO
      IF(IOVER.LE.1) WRITE(*,714) I, HT(J,I)
      IOVER=2
  712 FORMAT(5X,'TEMPERATURE IS GOING OUT-OF-LIMIT',I4,2X,D12.5)
  714 FORMAT(5X,'UNABLE TO FIND THE ROOT',I4,2X,D12.5)
  522 CONTINUE
      IF(TKD.GT.TPOL2) TKD=TPOL2
      IF(TKD.LT.TKDMIN) TKD=TKDMIN
      TK(J,I)=TKD/TSTR
  523 CONTINUE
  511 CONTINUE
      ENDDO
      ENDDO
      RETURN
      END
C***********************************************************************
      SUBROUTINE KESOLV(ISOR,RELXKE,TOLRKE,ISCHEM,SIGK,SIGE,RESDK,RESDE)
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
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
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
  110     CONTINUE
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
  112     CONTINUE
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
  200 CONTINUE
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
C***********************************************************************
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
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/FINJ/NFINJ,IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30)
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
            DO 530 I=IBODYM+1,IBODYP-1
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
            DO 602 J=JFINJM,JFINJP
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
C***********************************************************************
      SUBROUTINE PERTRB(IANOIS,IBNOIS,JANOIS,JBNOIS,IRAND)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)        
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),
     1  W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/NOISE/ INOISE,IVXTYP,FNOISE,ANOISE(200),VXTYP(5),TTIME
      PI2=2.0*3.141592654
      ICENT=(IANOIS+IBNOIS)/2
      JCENT=(JANOIS+JBNOIS)/2
      FT=TTIME*FNOISE
      FTIME=PI2*(FT-DFLOAT(INT(FT)))
      RADO=DSQRT((X(JBNOIS,IBNOIS)-X(JCENT,ICENT))**2
     1    +(Y(JBNOIS,IBNOIS)-Y(JCENT,ICENT))**2)
      DO 20 I=IANOIS,IBNOIS
      DO 22 J=JANOIS,JBNOIS
C      ARAND=RAN(IRAND)
C     ARAND=RANF()
C     UFLUC=ANOISE*(ARAND-0.5)*2.0
C     U(J,I)=U(J,I)+UFLUC
       IF(I.EQ.ICENT.AND.J.EQ.JCENT) GO TO 22
       YA=Y(J,I)-Y(JCENT,ICENT)
       XA=X(J,I)-X(JCENT,ICENT)
C      RAD=DSQRT(XA*XA+YA*YA)
      RAD=RADO
       R1=DCOS(FTIME*RADO/RAD)/RAD
       R2=DSIN(FTIME*RADO/RAD)/RAD
       U(J,I)=U(J,I)-ANOISE(1)*(YA*R1+XA*R2)
       V(J,I)=V(J,I)-ANOISE(1)*(XA*R1-YA*R2)
   22 CONTINUE
   20 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE FIRE(INSP,IGKEEP)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)                                   
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),
     1  W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      IFIRE=0
      DO 12 I=2,LI-1
      DO 12 J=2,LJ-1
      IF(FSP(J,I,INSP).GE.0.6.AND.FSP(J,I,INSP).LE.0.95) THEN
             IFIRE=IFIRE+1
             ISKIP(J,I)=-9
             END IF
   12 CONTINUE
      IF(IFIRE.EQ.0) WRITE(*,812) 
  812 FORMAT(5X,'FIRE IS NOT INTRODUCED INTO THE FLOW')
      IF(IFIRE.GE.1) WRITE(*,814) IFIRE,ITR,ITR+IGKEEP
  814 FORMAT(5X,50('&'),
     1  /5X,'&&',46X,'&&',
     2  /5X,'&&',5X,'FIRE IS INITIATED AT ',I4,' GRID POINTS',4X,'&&',
     3  /5X,'&&',8X,'BETWEEN ITERATIONS ',I3,' AND ',I3,8X,'&&',
     4  /5X,'&&',46X,'&&'/5X,50('&')/)
      RETURN
      END
C***********************************************************************
      SUBROUTINE PRTRACK(NOPR,NOINJ,NEWINJ,DTPR,PRMOVE,
     1                   PDIA,PDEN,PTHR,PVEL)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LSP=52,LNPR=200,LNMX=2000)              
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),
     1  W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1 HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDFSP(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/TRACK/ IPICK(101),JPICK(101),XPR(LNPR),YPR(LNPR),
     1              PATHX(LNMX,LNPR),PATHY(LNMX,LNPR),
     2              PUVEL(LNMX,LNPR),PVVEL(LNMX,LNPR),
     3              GUVEL(LNMX,LNPR),GVVEL(LNMX,LNPR)
      DIMENSION IAP(LNMX),JAP(LNMX)
C-----------------------------------------------------------------------
C    ITRACK - n th trace,
C    DTPR - non-dimensional time between n-th & (n-1)-th traces
C-----------------------------------------------------------------------
      PDLMT=1.0D-07/ALSTR
      LI3=LI-2
      LJ3=LJ-2
      NINT=10
      DTINT=DTPR/DFLOAT(NINT)
      NOINJ=NOINJ+NEWINJ
      IF(NOINJ.GT.LNMX) NOINJ=LNMX
      IF(NOINJ.LE.0) RETURN
      XMIN=X(1,1)+0.5*XXS(2)
      XMAX=X(1,LI)-0.5*XXS(LI)
      YMIN=Y(1,1)+0.5*YYS(2)
      YMAX=Y(LJ,1)-0.5*YYS(LJ)
      AXMAX=100.0/X(1,LI)
      AYMAX=100.0/(Y(LJ,1)-Y(1,1))
      IF(NEWINJ.EQ.1) THEN
          DO 10 N=1,NOPR
          DO 11 ITRACK=NOINJ,2,-1
          PATHX(ITRACK,N)=PATHX(ITRACK-1,N)
          PATHY(ITRACK,N)=PATHY(ITRACK-1,N)
          PUVEL(ITRACK,N)=PUVEL(ITRACK-1,N)
          PVVEL(ITRACK,N)=PVVEL(ITRACK-1,N)
   11     CONTINUE
          PATHX(1,N)=XPR(N)+PRMOVE
          PATHY(1,N)=YPR(N)
   10     CONTINUE
          END IF
      DO 100 N=1,NOPR
      DO 12 ITRACK=1,NOINJ
      XX=PATHX(ITRACK,N)
      YY=PATHY(ITRACK,N)
      IAP(ITRACK)=0
      JAP(ITRACK)=0
      IF(XX.LT.-999.0) GO TO 12
      IA1=INT((XX-X(1,1))*AXMAX)+1
      IA2=IPICK(IA1)
      IA3=IPICK(IA1+1)
      DO 14 I=IA2,IA3
      IF(X(1,I).GT.XX) GO TO 15
   14 CONTINUE
   15 CONTINUE
      IAA=I-1
      JA1=INT((YY-Y(1,1))*AYMAX)+1
      JA2=JPICK(JA1)
      JA3=JPICK(JA1+1)
      DO 16 J=JA2,JA3
      IF(Y(J,1).GT.YY) GO TO 17
   16 CONTINUE
   17 CONTINUE
      JAA=J-1
      IAP(ITRACK)=IAA
      JAP(ITRACK)=JAA
   12 CONTINUE
      DO 22 ITRACK=1,NOINJ
      IA=IAP(ITRACK)
      JA=JAP(ITRACK)
      IB=IA
      JB=JA
      IF(IA.LE.0.OR.JA.LE.0) GO TO 22
      XX=PATHX(ITRACK,N)
      YY=PATHY(ITRACK,N)
      IF((XX-X(JB,IB)).GE.0.5*XXS(IB+1)) IA=IA+1
      IF((YY-Y(JB,IB)).GE.0.5*YYS(JB+1)) JA=JA+1
C---------------------- NINE-POINT INTERPOLATION -----------------------
      IF(IA.GE.3.AND.IA.LE.LI3.AND.JA.GE.3.AND.JA.LE.LJ3) THEN
           XU1=XXC(IB)
           XU2=XXC(IB+1)
           YU1=YYS(JA)
           YU2=YYS(JA+1)
           XUP=XX-0.5*(X(1,IB+1)+X(1,IB))
           YUP=YY-Y(JA,1)
           XU12=XU1*XU2
           YU12=YU1*YU2
           U1=U(JA,IB+1)
           U2=((XU1*XU1*U(JA,IB+2)-XU2*XU2*U(JA,IB))/(XU1+XU2)
     1       -(XU1-XU2)*U1)/XU12*XUP
           U3=((XU1*U(JA,IB+2)+XU2*U(JA,IB))/(XU1+XU2)-U1)/XU12
     1       *XUP*XUP
           U4=((YU1*YU1*U(JA+1,IB+1)-YU2*YU2*U(JA-1,IB+1))/(YU1+YU2)
     1       -(YU1-YU2)*U1)/YU12*YUP
           U5=((YU1*U(JA+1,IB+1)+YU2*U(JA-1,IB+1))/(YU1+YU2)-U1)/YU12
     1       *YUP*YUP
           U6=(U(JA+1,IB+2)+U(JA-1,IB)-U(JA-1,IB+2)-U(JA+1,IB))
     1       /((XU1+XU2)*(YU1+YU2))*XUP*YUP
           UGS=U1+U2+U3+U4+U5+U6
           YV1=YYC(JB)
           YV2=YYC(JB+1)
           XV1=XXS(IA)
           XV2=XXS(IA+1)
           YVP=YY-0.5*(Y(JB+1,1)+Y(JB,1))
           XVP=XX-X(1,IA)
           XV12=XV1*XV2
           YV12=YV1*YV2
           V1=V(JB+1,IA)
           V2=((XV1*XV1*V(JB+1,IA+1)-XV2*XV2*V(JB+1,IA-1))/(XV1+XV2)
     1       -(XV1-XV2)*V1)/XV12*XVP
           V3=((XV1*V(JB+1,IA+1)+XV2*V(JB+1,IA-1))/(XV1+XV2)-V1)/XV12
     1       *XVP*XVP
           V4=((YV1*YV1*V(JB+2,IA)-YV2*YV2*V(JB,IA))/(YV1+YV2)
     1       -(YV1-YV2)*V1)/YV12*YVP
           V5=((YV1*V(JB+2,IA)+YV2*V(JB,IA))/(YV1+YV2)-V1)/YV12
     1       *YVP*YVP
           V6=(V(JB+2,IA+1)+V(JB,IA-1)-V(JB,IA+1)-V(JB+2,IA-1))
     1       /((XV1+XV2)*(YV1+YV2))*XVP*YVP
           VGS=V1+V2+V3+V4+V5+V6
           IF(PDIA.GT.0.0) THEN
           XSP=XX-X(1,IA)
           YSP=YY-Y(JA,1)
           XS1P2=XXS(IA)+XXS(IA+1)
           XS1M2=XXS(IA)-XXS(IA+1)
           XS11=XXS(IA)*XXS(IA)
           XS22=XXS(IA+1)*XXS(IA+1)
           XS12=XXS(IA)*XXS(IA+1)
           YS1P2=YYS(JA)+YYS(JA+1)
           YS1M2=YYS(JA)-YYS(JA+1)
           YS11=YYS(JA)*YYS(JA)
           YS22=YYS(JA+1)*YYS(JA+1)
           YS12=YYS(JA)*YYS(JA+1)
           RGS=RHO(JA,IA)+((XS11*RHO(JA,IA+1)-XS22*RHO(JA,IA-1))/XS1P2
     1        -XS1M2*RHO(JA,IA))/XS12*XSP
     2        +((XXS(IA)*RHO(JA,IA+1)+XXS(IA+1)*RHO(JA,IA-1))/XS1P2
     3        -RHO(JA,IA))/XS12*XSP*XSP
     4        +((YS11*RHO(JA+1,IA)-YS22*RHO(JA-1,IA))/YS1P2
     5        -YS1M2*RHO(JA,IA))/YS12*YSP
     6        +((YYS(JA)*RHO(JA+1,IA)+YYS(JA+1)*RHO(JA-1,IA))/YS1P2
     7        -RHO(JA,IA))/YS12*YSP*YSP
     8        +(RHO(JA+1,IA+1)+RHO(JA-1,IA-1)-RHO(JA-1,IA+1)
     9        -RHO(JA+1,IA-1))/(XS1P2*YS1P2)*XSP*YSP
           TGS=TK(JA,IA)+((XS11*TK(JA,IA+1)-XS22*TK(JA,IA-1))/XS1P2
     1        -XS1M2*TK(JA,IA))/XS12*XSP
     2        +((XXS(IA)*TK(JA,IA+1)+XXS(IA+1)*TK(JA,IA-1))/XS1P2
     3        -TK(JA,IA))/XS12*XSP*XSP
     4        +((YS11*TK(JA+1,IA)-YS22*TK(JA-1,IA))/YS1P2
     5        -YS1M2*TK(JA,IA))/YS12*YSP
     6        +((YYS(JA)*TK(JA+1,IA)+YYS(JA+1)*TK(JA-1,IA))/YS1P2
     7        -TK(JA,IA))/YS12*YSP*YSP
     8        +(TK(JA+1,IA+1)+TK(JA-1,IA-1)-TK(JA-1,IA+1)
     9        -TK(JA+1,IA-1))/(XS1P2*YS1P2)*XSP*YSP
          VMGS=VMU(JA,IA)+((XS11*VMU(JA,IA+1)-XS22*VMU(JA,IA-1))/XS1P2
     1        -XS1M2*VMU(JA,IA))/XS12*XSP
     2        +((XXS(IA)*VMU(JA,IA+1)+XXS(IA+1)*VMU(JA,IA-1))/XS1P2
     3        -VMU(JA,IA))/XS12*XSP*XSP
     4        +((YS11*VMU(JA+1,IA)-YS22*VMU(JA-1,IA))/YS1P2
     5        -YS1M2*VMU(JA,IA))/YS12*YSP
     6        +((YYS(JA)*VMU(JA+1,IA)+YYS(JA+1)*VMU(JA-1,IA))/YS1P2
     7        -VMU(JA,IA))/YS12*YSP*YSP
     8        +(VMU(JA+1,IA+1)+VMU(JA-1,IA-1)-VMU(JA-1,IA+1)
     9        -VMU(JA+1,IA-1))/(XS1P2*YS1P2)*XSP*YSP
          VTGS=VTC(JA,IA)+((XS11*VTC(JA,IA+1)-XS22*VTC(JA,IA-1))/XS1P2
     1        -XS1M2*VTC(JA,IA))/XS12*XSP
     2        +((XXS(IA)*VTC(JA,IA+1)+XXS(IA+1)*VTC(JA,IA-1))/XS1P2
     3        -VTC(JA,IA))/XS12*XSP*XSP
     4        +((YS11*VTC(JA+1,IA)-YS22*VTC(JA-1,IA))/YS1P2
     5        -YS1M2*VTC(JA,IA))/YS12*YSP
     6        +((YYS(JA)*VTC(JA+1,IA)+YYS(JA+1)*VTC(JA-1,IA))/YS1P2
     7        -VTC(JA,IA))/YS12*YSP*YSP
     8        +(VTC(JA+1,IA+1)+VTC(JA-1,IA-1)-VTC(JA-1,IA+1)
     9        -VTC(JA+1,IA-1))/(XS1P2*YS1P2)*XSP*YSP
               END IF 
           ELSE
           AL1=(XX-0.5*(X(JA,IA)+X(JA,IA-1)))/XXC(IA)
           UA=AL1*U(JA,IA+1)+(1.0-AL1)*U(JA,IA)
           UB=AL1*U(JA+1,IA+1)+(1.0-AL1)*U(JA+1,IA)  
           BL1=(YY-Y(JA,IA))/YYS(JA+1)
           UGS=BL1*UB+(1.0-BL1)*UA
           AL1=(YY-0.5*(Y(JA,IA)+Y(JA-1,IA)))/YYC(JA)
           VA=AL1*V(JA+1,IA)+(1.0-AL1)*V(JA,IA)
           VB=AL1*V(JA+1,IA+1)+(1.0-AL1)*V(JA,IA+1)  
           BL1=(XX-X(JA,IA))/XXS(IA+1)
           VGS=BL1*VB+(1.0-BL1)*VA
           IF(PDIA.GT.0.0) THEN 
           AL1=(XX-X(JA,IA))/XXS(IA+1)
           BL1=(YY-Y(JA,IA))/YYS(JA+1)
           RGS=BL1*(AL1*RHO(JA+1,IA+1)+(1.0-AL1)*RHO(JA+1,IA))
     1       +(1.0-BL1)*(AL1*RHO(JA,IA+1)+(1.0-AL1)*RHO(JA,IA))
           TGS=BL1*(AL1*TK(JA+1,IA+1)+(1.0-AL1)*TK(JA+1,IA))
     1       +(1.0-BL1)*(AL1*TK(JA,IA+1)+(1.0-AL1)*TK(JA,IA))
           VMGS=BL1*(AL1*VMU(JA+1,IA+1)+(1.0-AL1)*VMU(JA+1,IA))
     1       +(1.0-BL1)*(AL1*VMU(JA,IA+1)+(1.0-AL1)*VMU(JA,IA))
           VTGS=BL1*(AL1*VTC(JA+1,IA+1)+(1.0-AL1)*VTC(JA+1,IA))
     1       +(1.0-BL1)*(AL1*VTC(JA,IA+1)+(1.0-AL1)*VTC(JA,IA))
               END IF 
           END IF
      IF(NEWINJ.EQ.1.AND.ITRACK.EQ.1) THEN
                      UPR=PVEL*UGS
                      VPR=PVEL*VGS
                      ELSE 
                      UPR=PUVEL(ITRACK,N)
                      VPR=PVVEL(ITRACK,N)
                      END IF 
      UTP=0.0
      VTP=0.0
      IF(PDIA.GT.0.0) THEN
            IAA=IA+1
            JAA=JA+1
            AL1=(XX-X(JA,IA))/XXS(IAA)
            BL1=(YY-Y(JA,IA))/YYS(JAA)
            DTBTX=(BL1*(TK(JAA,IAA)-TK(JAA,IA))
     1           +(1.0-BL1)*(TK(JA,IAA)-TK(JA,IA)))/TGS/XXS(IAA)
            DTBTY=(AL1*(TK(JAA,IAA)-TK(JA,IAA))
     1           +(1.0-AL1)*(TK(JAA,IA)-TK(JA,IA)))/TGS/YYS(JAA)
            IF(PDIA.LE.PDLMT) THEN
                UPR=UGS-REI*0.5385226*VMGS*DTBTX/RGS
                VPR=VGS-REI*0.5385226*VMGS*DTBTY/RGS
                ELSE
                CSD=0.0
                FGRX=0.0
                FGRY=0.0
                CFTP=0.0
                IF(PTHR.GT.0.0) THEN
                    RBAR=0.0
                    DO 30 ISP=1,LSP
                    RBAR=RBAR+0.25*(FSP(JA,IA,ISP)+FSP(JA,IA+1,ISP)
     1              +FSP(JA+1,IA,ISP)+FSP(JA+1,IA+1,ISP))/WM(ISP)/WMSTR
   30               CONTINUE
                    AKN=REI*2.0*VMGS/(0.491*RGS*
     1               DSQRT(8.0*GASC*RBAR*TGS*TSTR/3.1415927)/VSTR)/PDIA
                    CFTP=-REI*REI*42.12*VMGS*VMGS*(VTGS/PTHR+2.18*AKN)/
     1               (RGS*PDEN*PDIA*PDIA*(1.0+3.42*AKN)
     2               *(1.0+2.0*VTGS/PTHR+4.36*AKN))
                    END IF
                FGRX=0.0
                FGRY=0.0
                IF(IGRAV.EQ.1.OR.IGRAV.EQ.2) FGRX=BETA4*(RGS/PDEN-1.0)
                IF(IGRAV.EQ.3.OR.IGRAV.EQ.4) FGRY=BETA4*(RGS/PDEN-1.0)
                REP=RGS*DSQRT((UGS-UPR)**2+(VGS-VPR)**2)*PDIA/VMGS/REI
                CSD=0.0
                IF(REP.GE.1.0D-20) CSD=
     1            18.0*(1.0+0.1667*REP**0.667)*RGS/PDIA/PDEN/REP
                CSDX=CSD*DABS(UGS-UPR)
                CSDY=CSD*DABS(VGS-VPR)
                UPR1=UPR
                VPR1=VPR
                DO 32 IINT=1,NINT
                UPR=(UPR+DTINT*(CSDX*UGS+FGRX+CFTP*DTBTX))
     1              /(1.0+DTINT*CSDX)
                VPR=(VPR+DTINT*(CSDY*VGS+FGRY+CFTP*DTBTY))
     1              /(1.0+DTINT*CSDY)
   32           CONTINUE
                END IF
            ELSE
	    UPR=UGS
	    VPR=VGS
	    END IF
      GUVEL(ITRACK,N)=TGS*TSTR
      GVVEL(ITRACK,N)=VGS
      PUVEL(ITRACK,N)=UPR
      PVVEL(ITRACK,N)=VPR
      XB=XX+DTPR*UPR
      YB=YY+DTPR*VPR
      PATHX(ITRACK,N)=XB
      PATHY(ITRACK,N)=YB
      IF((XB.GE.XMAX).OR.(YB.GE.YMAX).OR.(XB.LE.XMIN).OR.
     1   (YB.LE.YMIN)) PATHX(ITRACK,N)=-999999.0
   22 CONTINUE
  100 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE PRNTDK(MP,INFUEL)                      
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)                                   
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LE,LSP-1),FSN2(LE),
     6  U(LJ,LI),V(LJ,LI),W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),
     7  AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/HRLS/ QDOT(LJ,LI)
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/SOOT/ RSOOT,STMF(LE),STND(LE),STDF(LE)
C---add 200000 (for swirl), 400000 (for soot), 100000 (for bodies)
      IFLAME=200000+LSP*100+INFUEL
      IFLAME=IFLAME+400000
      IF(NBODY.GE.1) IFLAME=IFLAME+100000
      INERT=0
      WRITE(MP,100) IFLAME
      WRITE(MP,102) LI,LJ,ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,INERT,ITR,
     1              DT
      WRITE(MP,104) ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      IF(NBODY.GE.1) 
     1 WRITE(MP,103) NBODY,(IBDM(N),JBDM(N),IBDP(N),JBDP(N),N=1,NBODY)
      WRITE(MP,104) X,Y,RHO,FSP,U,V,W,P,QDOT,TK
      WRITE(MP,104) STMF,STND
      IF(IFLOW.EQ.2) WRITE(MP,104) AK,EPS
  100 FORMAT(I6,2('-'),'FLAME DATA FROM UNICOND-HEPT-SD(',
     1    ' 52 SPECIES & 544 REACTIONS OF SANDIEGO MECH) DATA',2('-'))
  102 FORMAT(10I6,F12.10)
  103 FORMAT(21I6)
  104 FORMAT(8(1PE14.7,1X))
      RETURN
      END
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
C***********************************************************************
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
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
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
      DO 1320 J=1,LJ-1
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
      DO 10 J=1,LJ
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
   10 CONTINUE
C----------------------------WRITING THE DATA---------------------------
      WRITE(13,102) ITR,TTIME
      WRITE(13,105)
  102 FORMAT(1X,I6,F15.10)
      DO 1340 I=1,LI
      DO 1340 J=1,LJ
C---I,J,XRG(J,I),YRG(J,I),
C---    RHO(J,I),U(J,I),V(J,I),T(J,I),QDOT(J,I),
C---    MU(J,I),TC(J,I),CP(J,I),
C---    FSP(J,I,ISP),DSP(J,I,ISP),QSP(J,I,ISP),ENTH(J,I,ISP)
      WRITE(13,104) I,J,XRG(J,I),YRG(J,I),RHORG(J,I),
     *          URG(J,I),VRG(J,I),
     *          TKRG(J,I),QRG(J,I),VMURG(J,I),VTCRG(J,I),CP(J,I),
     *         (FSPRG(J,I,ISP),VDFRG(J,I,ISP),QSPRG(J,I,ISP),
     *          ENTH(J,I,ISP),ISP=1,LSP)
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
C--NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW-NEW--
C***********************************************************************
      SUBROUTINE READFF(MP,IREAD,TTIME0)                      
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52) 
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),
     1  U(LJ,LI),V(LJ,LI),W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),
     2  AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDFSP(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/CB08/ TPOL1,TPOL2,POLSP(14,LSP),CISP(12,LSP)
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
      COMMON/HFORM/ HFO(LSP)
      COMMON/DUM0/ DATA2(LE*149)
      DIMENSION DAT1(20),ITEMP(600),JTEMP(200)
      READ(MP,102) IFLAME
      IRFL=IFLAME
                           NSOOT=0
      IF(IFLAME.GT.400000) THEN
                           NSOOT=1
                           IFLAME=IFLAME-400000
                           END IF
                           NSWRL=0
      IF(IFLAME.GT.200000) THEN
                           NSWRL=1
                           IFLAME=IFLAME-200000
                           END IF
      ISWRLA=0
      IF(NSWRL.EQ.0) READ(MP,101) ILA,JLA,ISYMA,IFLOWA,
     1    ITHRMA,ICHEMA,IGRAVA,INERTA,ITRA,DTA
      IF(NSWRL.EQ.1) READ(MP,102) ILA,JLA,ISYMA,IFLOWA,ISWRLA,
     1    ITHRMA,ICHEMA,IGRAVA,INERTA,ITRA,DTA
      READ(MP,104) DAT1
      IF(IFLAME.GE.100000) THEN
          IFLAME=IFLAME-100000
          READ(MP,103) NN,(IB1,JB1,IB2,JB2,N=1,NN)
          END IF
      WRITE(*,109) IRFL,ILA,JLA,ISYMA,IFLOWA,ISWRLA,
     1    ITHRMA,ICHEMA,IGRAVA,INERTA,ITRA,DTA*DAT1(1)/DAT1(2)
  109 FORMAT(6X,'IT CONTAINS--',I6,2X,2I4,7I2,I6,F10.7)
      WRITE(16,110) ILA,JLA,DAT1(2)
  110 FORMAT(1X,'READ DATA, IL=',I6,',  JL=',I6,',  VREF=',F7.2/)
C     TTIME0=DTA*DFLOAT(ITRA-1)*DAT1(1)/DAT1(2)
      TTIME0=0.0
      ILJL=ILA*JLA
      LLSP=IFLAME/100
      IDATA2=ILJL*(LLSP+7)+NSWRL*ILJL
      READ(MP,104) (DATA2(I),I=1,IDATA2)
      IF(NSOOT.EQ.1) THEN
           READ(MP,104) (DATA2(I),I=IDATA2+1,IDATA2+ILJL*2)
           IDATA2=IDATA2+ILJL*2
           END IF
      IF(IFLOWA.EQ.2) 
     1   READ(MP,104) (DATA2(I),I=IDATA2+1,IDATA2+ILJL*2)
  101 FORMAT(9I6,F12.10)
  102 FORMAT(10I6,F12.10)
  103 FORMAT(21I6)
  104 FORMAT(8(1PE14.7,1X))
      NI=ILJL*NSWRL
      NJ=ILJL*2*NSOOT
C--------------------------ONE TO ONE TRANSFER--------------------------
      IF(IREAD.EQ.2) THEN
          IF((LI.NE.ILA).OR.(LJ.NE.JLA).OR.(LSP.NE.LLSP)) THEN
            WRITE(*,922) ILA,JLA,LLSP,IFLAME
  922       FORMAT(5X,'ILA = ',I3,',  JLA = ',I3,',  LLSP = ',I3,
     1      '; but IREAD = ',I2/5X,'****RESETTING IREAD = 1****')
            IREAD=1
            GO TO 920
            END IF
          DO 910 I=1,LI
          DO 910 J=1,LJ
          IA=(I-1)*LJ+J
          RHO(J,I)=DATA2(IA+ILJL*2)
          TOT=0.0
          DO 912 ISP=1,LLSP-1
          FSP(J,I,ISP)=DATA2(IA+ILJL*2+ILJL*ISP)
          TOT=TOT+FSP(J,I,ISP)
  912     CONTINUE
          FSP(J,I,LSP)=1.0-TOT
          U(J,I)=DATA2(IA+ILJL*2+ILJL*LLSP)
          V(J,I)=DATA2(IA+ILJL*3+ILJL*LLSP)
          W(J,I)=0.0
          IF(NSWRL.EQ.1) W(J,I)=DATA2(IA+ILJL*4+ILJL*LLSP)
          P(J,I)=DATA2(IA+ILJL*4+NI+ILJL*LLSP)
          HT(J,I)=DATA2(IA+ILJL*5+NI+ILJL*LLSP)
          TK(J,I)=DATA2(IA+ILJL*6+NI+ILJL*LLSP)
          STMF(J,I)=0.0
          STND(J,I)=0.0
          IF(NSOOT.EQ.1) THEN
            STMF(J,I)=DATA2(IA+ILJL*7+NI+ILJL*LLSP)
            STND(J,I)=DATA2(IA+ILJL*8+NI+ILJL*LLSP)
          END IF
          AK(J,I)=0.0
          EPS(J,I)=0.0
          IF(IFLOWA.EQ.2) THEN
            AK(J,I)=DATA2(IA+ILJL*7+NI+NJ+ILJL*LLSP)
            EPS(J,I)=DATA2(IA+ILJL*8+NI+NJ+ILJL*LLSP)
          END IF
  910     CONTINUE
          GO TO 397
          END IF
  920     CONTINUE
C----------------CONVERT READ X & Y INTO THE PRESENT MESH SYSTEM--------
      GCONV=DAT1(1)/ALSTR
      DO 120 I=1,ILJL+ILJL
      DATA2(I)=DATA2(I)*GCONV
      IF(DATA2(I).GT.1.0D-15) DATA2(I)=DATA2(I)-1.0D-15
  120 CONTINUE
C--------------------------------OVER-----------------------------------
C-------------- CONVERT READ GLOBAL -->> FINITE-RATE CHEMISTRY ---------
      IF(LLSP.EQ.5) THEN
      DO 132 I=1,ILJL
      FFU=DATA2(ILJL*3+I)
      FO2=DATA2(ILJL*4+I)
      FCO2=DATA2(ILJL*5+I)
      FH2O=DATA2(ILJL*6+I)
      FU=DATA2(ILJL*7+I)
      FV=DATA2(ILJL*8+I)
      FW=0.0
      IF(NSWRL.EQ.1) FW=DATA2(ILJL*9+I)
      FP=DATA2(ILJL*9+NI+I)
      FH=DATA2(ILJL*10+NI+I)
      FT=DATA2(ILJL*11+NI+I)
      FSTMF=0.0
      FSTND=0.0
      IF(NSOOT.EQ.1) THEN
         FSTMF=DATA2(ILJL*12+NI+I)
         FSTND=DATA2(ILJL*13+NI+I)
         END IF
      FAK=0.0
      FEPS=0.0
      IF(IFLOWA.EQ.2) THEN
         FAK=DATA2(ILJL*12+NI+NJ+I)
         FEPS=DATA2(ILJL*13+NI+NJ+I)
         END IF
      DO 133 ISP=1,LSP-1
      II=ILJL*(ISP+2)+I
      DATA2(II)=0.0
  133 CONTINUE      
      RADADD=0.0001*FH2O
      ILJLU=ILJL*(LSP+2)
      DATA2(ILJLU+ILJL*9+I)=FEPS
      DATA2(ILJLU+ILJL*8+I)=FAK
      DATA2(ILJLU+ILJL*7+I)=FSTND
      DATA2(ILJLU+ILJL*6+I)=FSTMF
      DATA2(ILJLU+ILJL*5+I)=FT
      DATA2(ILJLU+ILJL*4+I)=FH
      DATA2(ILJLU+ILJL*3+I)=FP
      DATA2(ILJLU+ILJL*2+I)=FW
      DATA2(ILJLU+ILJL+I)=FV
      DATA2(ILJLU+I)=FU
      DATA2(ILJL*08+I)=FH2O-3.0*RADADD
      DATA2(ILJL*07+I)=RADADD
      DATA2(ILJL*06+I)=RADADD
      DATA2(ILJL*05+I)=RADADD
      DATA2(ILJL*12+I)=FCO2
      DATA2(ILJL*04+I)=FO2
      IF(IFLAME.EQ.511) DATA2(ILJL*49+I)=FFU
      IF(IFLAME.EQ.510) DATA2(ILJL*11+I)=FFU
      IF(IFLAME.EQ.509) DATA2(ILJL*14+I)=FFU
      IF(IFLAME.EQ.508) DATA2(ILJL*42+I)=FFU
      IF(IFLAME.EQ.507) DATA2(ILJL*43+I)=FFU
      IF(IFLAME.EQ.506) DATA2(ILJL*22+I)=FFU
      IF(IFLAME.EQ.505) DATA2(ILJL*19+I)=FFU
      IF(IFLAME.EQ.504) DATA2(ILJL*24+I)=FFU
      IF(IFLAME.EQ.503) DATA2(ILJL*32+I)=FFU
      IF(IFLAME.EQ.502) DATA2(ILJL*15+I)=FFU
      IF(IFLAME.EQ.501) DATA2(ILJL*03+I)=FFU
  132 CONTINUE
      END IF
C--------------------------------OVER-----------------------------------
      JB=1
      DO 210 J=1,LJ
      YNEW=Y(J,1)
      JA=JB-1
      IF(JA.LT.1) JA=1
      DO 212 JB=JA,JLA
      YOLD=DATA2(ILJL+JB)
      JTEMP(J)=JB
      IF(YOLD.GT.YNEW) GO TO 214
  212 CONTINUE
  214 CONTINUE
  210 CONTINUE
      IB=1
      DO 220 I=1,LI
      XNEW=X(1,I)
      IA=IB-1
      IF(IA.LT.1) IA=1
      DO 222 IB=IA,ILA
      IIB=(IB-1)*JLA+1
      IIA=IIB-JLA
      XOLD=DATA2(IIB)
      ITEMP(I)=IIB
      IF(XOLD.GT.XNEW) GO TO 224
  222 CONTINUE
  224 CONTINUE
  220 CONTINUE
      AF=DAT1(2)/VSTR
      AFS=AF*AF
C--------------------------INTERPOLATE SCALARS--------------------------
      DO 310 I=1,LI
      XNEW=X(1,I)
      IIB=ITEMP(I)
      IIA=IIB-JLA
      XOLD=DATA2(IIB)
      XOLDM=DATA2(IIA)
      AX1=(XNEW-XOLDM)/(XOLD-XOLDM)
      AX2=(XOLD-XNEW)/(XOLD-XOLDM)
C     IF(AX2.LT.0.0) THEN
                     AX1=1.0
                     AX2=0.0
C                    END IF
      DO 312 J=1,LJ
      YNEW=Y(J,1)
      JJB=JTEMP(J)
      JJA=JJB-1
      YOLD=DATA2(ILJL+JJB)
      YOLDM=DATA2(ILJL+JJA)
      AY1=(YNEW-YOLDM)/(YOLD-YOLDM)
      AY2=(YOLD-YNEW)/(YOLD-YOLDM)
C     IF(AY2.LT.0.0) THEN
                     AY1=1.0
                     AY2=0.0
C                    END IF
      IXY=IIA+JJA-1+ILJL*2
      RHO(J,I)=AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1        +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1)
      FTOT=0.0
      DO 314 ISP=1,LSP-1
      IXY=IIA+JJA-1+ILJL*(ISP+2)
      FSP(J,I,ISP)=AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1            +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1)
      FTOT=FTOT+FSP(J,I,ISP)
  314 CONTINUE
      IXY=IIA+JJA-1+ILJL*(LSP+4)
      W(J,I)=AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1      +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1)
      IXY=IIA+JJA-1+ILJL*(LSP+5)
      P(J,I)=(AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1       +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1))*AFS
      IXY=IIA+JJA-1+ILJL*(LSP+6)
      HT(J,I)=(AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1        +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1))*AFS
      IXY=IIA+JJA-1+ILJL*(LSP+7)
      TK(J,I)=AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1       +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1)
      IXY=IIA+JJA-1+ILJL*(LSP+8)
      STMF(J,I)=AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1       +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1)
      IXY=IIA+JJA-1+ILJL*(LSP+9)
      STND(J,I)=AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1        +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1)
      FTOT=FTOT+STMF(J,I)
      FSP(J,I,LSP)=1.0-FTOT
      IF(IFLOWA.NE.2) GO TO 312
      IXY=IIA+JJA-1+ILJL*(LSP+10)
      AK(J,I)=AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1       +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1)
      IXY=IIA+JJA-1+ILJL*(LSP+11)
      EPS(J,I)=AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1        +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1)
  312 CONTINUE
  310 CONTINUE
C-----------------------------STAGGERED GRID----------------------------
C-----------------------------INTERPOLATE U-----------------------------
      IB=2
      DO 320 I=2,LI
      XNEW=0.5*(X(1,I-1)+X(1,I))
      IA=IB-1
      IF(IA.LT.2) IA=2
      DO 322 IB=IA,ILA
      IIB=(IB-1)*JLA+1
      IIA=IIB-JLA
      XOLD=0.5*(DATA2(IIA)+DATA2(IIB))
      IF(XOLD.GT.XNEW) GO TO 324
  322 CONTINUE
  324 CONTINUE
      IIM=IIA-JLA
      IF(IIM.LT.1) IIM=1
      XOLDM=0.5*(DATA2(IIM)+DATA2(IIA))
      AX1=(XNEW-XOLDM)/(XOLD-XOLDM)
      AX2=(XOLD-XNEW)/(XOLD-XOLDM)
      IF(AX2.LT.0.0) THEN
                     AX1=1.0
                     AX2=0.0
                     END IF
      DO 326 J=1,LJ
      YNEW=Y(J,1)
      JJB=JTEMP(J)
      JJA=JJB-1
      YOLD=DATA2(ILJL+JJB)
      YOLDM=DATA2(ILJL+JJA)
      AY1=(YNEW-YOLDM)/(YOLD-YOLDM)
      AY2=(YOLD-YNEW)/(YOLD-YOLDM)
      IF(AY2.LT.0.0) THEN
                     AY1=1.0
                     AY2=0.0
                     END IF
      IXY=IIA+JJA-1+ILJL*(LSP+2)
      U(J,I)=(AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1       +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1))*AF
  326 CONTINUE
  320 CONTINUE
C-----------------------------INTERPOLATE V-----------------------------
      JTEMP(1)=1
      JB=2
      DO 330 J=2,LJ
      YNEW=0.5*(Y(J-1,1)+Y(J,1))
      JA=JB-1
      IF(JA.LT.2) JA=2
      DO 332 JB=JA,JLA
      YOLD=0.5*(DATA2(ILJL+JB-1)+DATA2(ILJL+JB))
      JTEMP(J)=JB
      IF(YOLD.GT.YNEW) GO TO 334
  332 CONTINUE
  334 CONTINUE
  330 CONTINUE
      DO 336 I=1,LI
      XNEW=X(1,I)
      IIB=ITEMP(I)
      IIA=IIB-JLA
      XOLD=DATA2(IIB)
      XOLDM=DATA2(IIA)
      AX1=(XNEW-XOLDM)/(XOLD-XOLDM)
      AX2=(XOLD-XNEW)/(XOLD-XOLDM)
      IF(AX2.LT.0.0) THEN
                     AX1=1.0
                     AX2=0.0
                     END IF
      DO 338 J=2,LJ
      YNEW=Y(J,1)
      JJB=JTEMP(J)
      JJA=JJB-1
      JJM=JJA-1
      IF(JJM.LT.1) JJM=1
      YOLD=0.5*(DATA2(ILJL+JJA)+DATA2(ILJL+JJB))
      YOLDM=0.5*(DATA2(ILJL+JJM)+DATA2(ILJL+JJA))
      AY1=(YNEW-YOLDM)/(YOLD-YOLDM)
      AY2=(YOLD-YNEW)/(YOLD-YOLDM)
      IF(AY2.LT.0.0) THEN
                     AY1=1.0
                     AY2=0.0
                     END IF
      IXY=IIA+JJA-1+ILJL*(LSP+3)
      V(J,I)=(AX2*AY2*DATA2(IXY)+AX1*AY2*DATA2(IXY+JLA)
     1       +AX2*AY1*DATA2(IXY+1)+AX1*AY1*DATA2(IXY+JLA+1))*AF
  338 CONTINUE
  336 CONTINUE
C------------------CALCULATE ENTHALPY USING NEW VARIABLES----------------
  397 CONTINUE
      DO 398 I=1,LI
      DO 398 J=1,LJ
      TKD=TK(J,I)*TSTR
      TKD2=TKD*TKD
      TKD3=TKD*TKD2
      TKD4=TKD*TKD3
      TKD5=TKD*TKD4
      IPOLY=1
      IF(TKD.GT.TPOL1) IPOLY=8
      HNEW=0.0
      DO 399 ISP=1,LSP
      HNEW=HNEW+FSP(J,I,ISP)*((POLSP(IPOLY,ISP)*TKD
     1    +POLSP(IPOLY+1,ISP)*TKD2/2.0+POLSP(IPOLY+2,ISP)*TKD3/3.0
     1    +POLSP(IPOLY+3,ISP)*TKD4/4.0+POLSP(IPOLY+4,ISP)*TKD5/5.0
     2    +POLSP(IPOLY+5,ISP))*GASC/WMSTR/HSTR/WM(ISP)-HFO(ISP))
  399 CONTINUE
      HT(J,I)=HNEW
     1       +0.5*(U(J,I)*U(J,I)+V(J,I)*V(J,I)+W(J,I)*W(J,I))
  398 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE PRNTFF(INSP,ITRANS)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)                                   
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),
     *  U(LJ,LI),V(LJ,LI),W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),
     *  AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDFSP(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/CB08/ TPOL1,TPOL2,POLSP(14,LSP),CISP(12,LSP)
      COMMON/HRLS/ QDOT(LJ,LI)
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
      DIMENSION IP(16),A(16)
      CHARACTER CAPT1(LSP+16)*18                                             
      DATA CAPT1 
     1 /'     PRESSURE     ','      DENSITY     ','   TEMPERATURE    ',    
     2  '   U VELOCITY     ','    V VELOCITY    ','  SWIRL VELOCITY  ',
     3  '  KINETIC ENERGY  ','  KE DISSIPATION  ',' N2 MASS FRACTION ',
     *  'SOOT MASS FRACTION','SOOT NUMBE DENSITY',
     *  'LAMINAR VISCOSITY ','TURBULT VISCOSITY ','HEAT CONDUCTIVITY ',
     *  'CH4-DIFFUSION CO. ','O2-DIFFUSION CO.  ','   HRR (J/CC/S)   ',
     *  'H2                ','O2                ','H                 ',
     *  'O                 ','OH                ','H2O               ',
     *  'HO2               ','H2O2              ','CO                ',
     *  'CO2               ','HCO               ','CH2O              ',
     *  'CH4               ','CH3               ','T-CH2             ',
     *  'S-CH2             ','C2H4              ','CH3O              ',
     *  'C2H5              ','C2H6              ','CH                ',
     *  'C2H2              ','C2H3              ','CH2CHO            ',
     *  'C2H4O             ','CH2CO             ','HCCO              ',
     *  'C2H               ','CH2OH             ','CH3OH             ',
     *  'C2H5OH            ','CH3CHO            ','CH3CHOH           ',
     *  'CH2CH2OH          ','CH3CO             ','CH3CH2O           ',
     *  'C3H4              ','C3H3              ','C3H5              ',
     *  'C3H6              ','C3H8              ','I-C3H7            ',
     *  'N-C3H7            ','C4H6              ','CHO               ',
     *  'C5H8              ','C7H16             ','C4H8              ',
     *  'C5H10             ','AR                ','HE                '/
      KA=LI
      KB=LJ
      IF(ITRANS.EQ.1) THEN
                      KA=LJ
                      KB=LI
                      END IF
      IVAL=8
      IP(1)=1                                                           
      DO 10 I=2,IVAL-1                                                      
      IP(I)=INT(DFLOAT(KA*(I-1))/DFLOAT(IVAL-1))                             
      IF(IP(I).EQ.IP(I-1)) IP(I)=IP(I)+1
      IF(IP(I).GT.LI) IP(I)=KA
   10 CONTINUE
      IP(IVAL)=KA
      WRITE(16,810)                                                      
  810 FORMAT(/20X,'GEOMETRICAL DATA IN millimeters'/)                          
      WRITE(16,812)                                                      
  812 FORMAT(2X,'NO.',4X,'X-Dis.',5X,'DX',7X,'R-outer ',     
     1       5X,'NO.',3X,'Y-Dis.',5X,'DY',7X,'Velocity')     
      AA=1000.0*ALSTR
      IJ=LI
      IF(LJ.GT.LI) IJ=LJ
      DO 12 I=1,IJ
      DXX=0.0
      DYY=0.0
      IF((I.LE.LI).AND.(I.LE.LJ)) THEN
          IF(I.GE.2) DXX=(X(1,I)-X(1,I-1))*AA
          IF(I.GE.2) DYY=(Y(I,1)-Y(I-1,1))*AA
          WRITE(16,814) I,X(1,I)*AA,DXX,
     1    Y(LJ,I)*AA,I,Y(I,1)*AA,DYY,U(I,1)*VSTR
          END IF
      IF((I.LE.LI).AND.(I.GT.LJ)) THEN
          DXX=(X(1,I)-X(1,I-1))*AA
          WRITE(16,814) I,X(1,I)*AA,DXX,Y(LJ,I)*AA
          END IF
      IF((I.GT.LI).AND.(I.LE.LJ)) THEN
          DYY=(Y(I,1)-Y(I-1,1))*AA
          WRITE(16,815) I,Y(I,1)*AA,DYY,U(I,1)*VSTR
          END IF
  814 FORMAT(2X,I3,F10.3,1X,F10.7,F10.5,5X,I3,F10.5,1X,F10.7,1X,F10.6)                    
  815 FORMAT(42X,I3,F10.5,2X,F10.7,1X,F10.6)                    
   12 CONTINUE                                                          
      DO 20 II=1,LSP+16                                                     
      IF(II.LE.17) THEN
                   WRITE(16,818) ITR,CAPT1(II)
                   ELSE
                   WRITE(16,819) ITR,II-17,CAPT1(II)
                   END IF
  818 FORMAT(/1X,14('@'),'  ITR=',I7,7X,A18,10X,14('@')/)                            
  819 FORMAT(/1X,14('@'),'  ITR=',I7,3X,I3,'-',A18,10X,14('@')/)                            
      IF(ITRANS.EQ.0) WRITE(16,830) (IP(I),I=1,IVAL)                                       
      IF(ITRANS.EQ.1) WRITE(16,831) (IP(I),I=1,IVAL)                                       
  830 FORMAT(1X,'(J)   ',16('I=',I3,4X))                         
  831 FORMAT(1X,'(I)   ',16('J=',I3,4X))                         
      DO 22 K=1,KB                                                      
      DO 24 IA=1,IVAL
      IB1=IP(IA)                                                        
      J=K*(1-ITRANS)+IB1*ITRANS
      I=K*ITRANS+IB1*(1-ITRANS)
      IF(II.EQ.1) A(IA)=P(J,I)*PSTR                       
      IF(II.EQ.2) A(IA)=RHO(J,I)*RSTR                                 
      IF(II.EQ.3) A(IA)=TK(J,I)*TSTR                       
      IF(II.EQ.4) A(IA)=U(J,I)*VSTR                       
      IF(II.EQ.5) A(IA)=V(J,I)*VSTR                          
      IF(II.EQ.6) A(IA)=W(J,I)*VSTR                          
      IF(II.EQ.7) A(IA)=AK(J,I)*AKSTR                          
      IF(II.EQ.8) A(IA)=EPS(J,I)*EPSTR                          
      IF(II.EQ.9) A(IA)=FSP(J,I,LSP)
      IF(II.EQ.10) A(IA)=STMF(J,I)                     
      IF(II.EQ.11) A(IA)=STND(J,I)                       
      IF(II.EQ.12) A(IA)=VMU(J,I)*AMSTR                       
      IF(II.EQ.13) A(IA)=TMU(J,I)*AMSTR                        
      IF(II.EQ.14) A(IA)=VTC(J,I)*ACSTR                       
      IF(II.EQ.15) A(IA)=VDFSP(J,I,INSP)*DSTR                   
      IF(II.EQ.16) A(IA)=VDFSP(J,I,2)*DSTR                   
      IF(II.EQ.17) A(IA)=QDOT(J,I)                   
      IF(II.GE.18) A(IA)=FSP(J,I,II-17)                     
   24 CONTINUE                                                          
      IF(II.EQ.1) WRITE(16,824) K,(A(IA),IA=1,IVAL)                                
      IF(II.EQ.2) WRITE(16,821) K,(A(IA),IA=1,IVAL)                                
      IF(II.EQ.3) WRITE(16,822) K,(A(IA),IA=1,IVAL)                                
      IF((II.GE.4).AND.(II.LE.6)) WRITE(16,823) K,(A(IA),IA=1,IVAL)                                
      IF((II.GE.7).AND.(II.LE.8)) WRITE(16,826) K,(A(I),IA=1,IVAL)                                
      IF(II.EQ.9) WRITE(16,825) K,(A(IA),IA=1,IVAL)                                
      IF(II.GE.10) WRITE(16,826) K,(A(IA),IA=1,IVAL)                                
  821 FORMAT(1X,I3,16F9.5)                                      
  822 FORMAT(1X,I3,16F9.2)                                      
  823 FORMAT(1X,I3,16F9.5)                                      
  824 FORMAT(1X,I3,1P16D9.2)                                      
  825 FORMAT(1X,I3,16F9.5)                                      
  826 FORMAT(1X,I3,1P16D9.2)                                      
   22 CONTINUE                                                          
   20 CONTINUE                                                          
      RETURN                                                            
      END                                                               
C***********************************************************************
      SUBROUTINE PRESET(ICODE,RTIN,RTOT,ALENG,
     1                  XBDM,YBDM,XBDP,YBDP,XFIM,YFIM,XFIP,YFIP,
     2                  NIREGN,IREGN,XREGN,NJREGN,JREGN,YREGN)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)                                   
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),
     1  W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/FINJ/NFINJ,IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30)
      DIMENSION XBDM(10),YBDM(10),XBDP(10),YBDP(10)
      DIMENSION XFIM(10),YFIM(10),XFIP(10),YFIP(10)
      DIMENSION IREGN(10),XREGN(10),JREGN(10),YREGN(10)
      DIMENSION AAI(LI),AAJ(LJ)
C-----------------------------------------------------------------------
      BETAX=0.0
      BETAY=0.0
      BTMIN=0.001
      BTMAX=1.55
C----------------   UNIFORM MESH IN Y-DIRECTION   ----------------------
      IF(NJREGN.EQ.0) THEN
                      DO 120 J=1,LJ
                      AAJ(J)=DFLOAT(J-1)/DFLOAT(LJ-1)
  120                 CONTINUE
                      GO TO 129
                      END IF
C----------------   CLUSTERED MESH IN Y-DIRECTION   -------------------
      DY1=YREGN(1)/DFLOAT(JREGN(1))/(RTOT-RTIN)
      IF(JREGN(1).LT.0) DY1=YREGN(2)/DFLOAT(JREGN(2))/(RTOT-RTIN)
      JL2=1
      DO 20 J=1,JL2
      AAJ(J)=DY1*DFLOAT(J-1)
   20 CONTINUE
      IF(NJREGN.LE.1) GO TO 123
      DO 122 N=1,NJREGN
      IF(N.GE.2) DY1=(AAJ(JL2)-AAJ(JL2-1))
      DYAV=YREGN(N)/(RTOT-RTIN)/DABS(DFLOAT(JREGN(N)))
      DYLAST=2.0*DYAV-DY1
      IF(DYLAST.LE.0.0) THEN
         WRITE(*,721) DYAV*(RTOT-RTIN),DY1*(RTOT-RTIN),N
  721    FORMAT(1X,'TOO MUCH COMPRESSION; DYAV & DY1=',2F10.6,'; N=',I2)
         ICODE=201
         RETURN
         END IF
      JSWT=0
      IF(JREGN(N).LT.0) THEN
                        JREGN(N)=-JREGN(N)
                        JSWT=1
                        END IF
      JL1=JL2+1
      JL2=JL1+JREGN(N)-1
      ADY=2.0*(DYAV-DY1)/DFLOAT(JREGN(N)-1)
      DO 22 J=JL1,JL2
      J2=JL2*JSWT+J*(1-JSWT)
      AAJ(J)=AAJ(JL1-1)+DY1*DFLOAT(J-JL1+1)
     1      +ADY*DFLOAT((J2-JL1)*(J2-JL1+1)-(J2-J-1)*(J2-J))/2.0
   22 CONTINUE
  122 CONTINUE
  123 CONTINUE
      NN=LJ-JL2
      IF(NN.LE.0) THEN
                  IF(AAJ(LJ).LT.0.99999) WRITE(*,722) AAJ(LJ)
  722             FORMAT(5X,'LAST POINT IN Y DIRECTION HAS A VALUE = ',
     1                       F12.8)
                  GO TO 129
                  END IF
      DNN=DFLOAT(NN)
      AL=(1.0-AAJ(JL2))
      DY=(AAJ(JL2)-AAJ(JL2-1))
      DYRAT=(AL/DNN)/DY
      IF(DYRAT.LT.1.0001) THEN
             WRITE(*,723) DY*(RTOT-RTIN),AL*(RTOT-RTIN),NN
  723        FORMAT(1X,'BETTER USE UNIFORM MESH IN Y',2X,2F10.6,2X,I5)
             ICODE=201
             RETURN
             END IF
      FUN1=DTAN(BTMIN/DNN)-DY*DTAN(BTMIN)/AL
      FUN2=DTAN(BTMAX/DNN)-DY*DTAN(BTMAX)/AL
      IF(FUN1.LE.0.0) THEN
             WRITE(*,724) DY*(RTOT-RTIN),AL*(RTOT-RTIN),NN
  724        FORMAT(1X,'TOO BIG INITIAL DY',2X,2F10.6,2X,I5)
             ICODE=202
             RETURN
             END IF
      IF(FUN2.GE.0.0) THEN
             WRITE(*,725) DY*(RTOT-RTIN),AL*(RTOT-RTIN),NN
  725        FORMAT(1X,'TOO SMALL INITIAL DY',2X,2F10.6,2X,I5)
             ICODE=203
             RETURN
             END IF
      DBETA=(BTMAX-BTMIN)/2.0
      BETA=BTMIN
   24 CONTINUE
      FUN=DTAN((BETA+DBETA)/DNN)-DY*DTAN((BETA+DBETA))/AL
      IF(FUN.GT.0.0) BETA=BETA+DBETA
      IF(FUN.LT.0.0) DBETA=DBETA/2.0
      IF(DABS(FUN).LE.1.0D-10) GO TO 26
      IF(DBETA.LE.1.0D-10) GO TO 26
      GO TO 24
   26 CONTINUE
      BETAY=BETA
      DO 28 J=JL2+1,LJ
      AAJ(J)=AAJ(JL2)+AL*DTAN(BETA*DFLOAT(J-JL2)/DNN)/DTAN(BETA)
   28 CONTINUE
  129 CONTINUE
C-----------------------------------------------------------------------
C----------------   UNIFORM MESH IN X-DIRECTION   ----------------------
      IF(NIREGN.EQ.0) THEN
                      DO 140 I=1,LI
                      AAI(I)=DFLOAT(I-1)/DFLOAT(LI-1)
  140                 CONTINUE
                      GO TO 149
                      END IF
C----------------   CLUSTER14D MESH IN X-DIRECTION   -------------------
      DX1=XREGN(1)/DFLOAT(IREGN(1))/ALENG
      IF(IREGN(1).LT.0) DX1=XREGN(2)/DFLOAT(IREGN(2))/ALENG
      IL2=1
      AAI(IL2)=0.0
      IF(NIREGN.LE.1) GO TO 143
      DO 142 N=1,NIREGN
      IF(N.GE.2) DX1=(AAI(IL2)-AAI(IL2-1))
      DXAV=XREGN(N)/ALENG/DABS(DFLOAT(IREGN(N)))
      DXLAST=2.0*DXAV-DX1
      IF(DXLAST.LE.0.0) THEN
         WRITE(*,741) DXAV*ALENG,DX1*ALENG,N
  741    FORMAT(1X,'TOO MUCH COMPRESSION; DXAV & DX1=',2F10.6,'; N=',I2)
         ICODE=201
         RETURN
         END IF
      ISWT=0
      IF(IREGN(N).LT.0) THEN
                        IREGN(N)=-IREGN(N)
                        ISWT=1
                        END IF
      IL1=IL2+1
      IL2=IL1+IREGN(N)-1
      ADX=2.0*(DXAV-DX1)/DFLOAT(IREGN(N)-1)
      DO 42 I=IL1,IL2
      I2=IL2*ISWT+I*(1-ISWT)
      AAI(I)=AAI(IL1-1)+DX1*DFLOAT(I-IL1+1)
     1      +ADX*DFLOAT((I2-IL1)*(I2-IL1+1)-(I2-I-1)*(I2-I))/2.0
   42 CONTINUE
  142 CONTINUE
  143 CONTINUE
      NN=LI-IL2
      IF(NN.LE.0) THEN
                  IF(AAI(LI).LT.0.99999) WRITE(*,742) AAI(LI)
  742             FORMAT(5X,'LAST POINT IN X DIRECTION HAS A VALUE = ',
     1                       F12.8)
                  GO TO 149
                  END IF
      DNN=DFLOAT(NN)
      AL=(1.0-AAI(IL2))
      DX=(AAI(IL2)-AAI(IL2-1))
      DXRAT=(AL/DNN)/DX
      IF(DXRAT.LT.1.0001) THEN
             WRITE(*,743) DX*ALENG,AL*ALENG,NN
  743        FORMAT(1X,'BETTER USE UNIFORM MESH IN X',2X,2F10.6,2X,I5)
             ICODE=211
             RETURN
             END IF
      FUN1=DTAN(BTMIN/DNN)-DX*DTAN(BTMIN)/AL
      FUN2=DTAN(BTMAX/DNN)-DX*DTAN(BTMAX)/AL
      IF(FUN1.LE.0.0) THEN
             WRITE(*,744) DX*ALENG,AL*ALENG,NN
  744        FORMAT(1X,'TOO BIG INITIAL DX',2X,2F10.6,2X,I5)
             ICODE=212
             RETURN
             END IF
      IF(FUN2.GE.0.0) THEN
             WRITE(*,745) DX*ALENG,AL*ALENG,NN
  745        FORMAT(1X,'TOO SMALL INITIAL DX',2X,2F10.6,2X,I5)
             ICODE=213
             RETURN
             END IF
      DBETA=(BTMAX-BTMIN)/2.0
      BETA=BTMIN
   44 CONTINUE
      FUN=DTAN((BETA+DBETA)/DNN)-DX*DTAN((BETA+DBETA))/AL
      IF(FUN.GT.0.0) BETA=BETA+DBETA
      IF(FUN.LT.0.0) DBETA=DBETA/2.0
      IF(DABS(FUN).LE.1.0D-10) GO TO 46
      IF(DBETA.LE.1.0D-10) GO TO 46
      GO TO 44
   46 CONTINUE
      BETAX=BETA
      DO 48 I=IL2+1,LI
      AAI(I)=AAI(IL2)+AL*DTAN(BETA*DFLOAT(I-IL2)/DNN)/DTAN(BETA)
   48 CONTINUE
  149 CONTINUE
      WRITE(*,70) BETAX,BETAY
   70 FORMAT(5X,'BETA-X='F10.6,4X,'BETA-Y=',F10.6)
C-----------------------------------------------------------------------
      DO 72 I=1,LI
      XDIS=AAI(I)*ALENG
      DO 72 J=1,LJ                                                      
      X(J,I)=XDIS/ALSTR
      Y(J,I)=(RTIN+(RTOT-RTIN)*AAJ(J))/ALSTR
   72 CONTINUE
C--------------------------Inserted Body Lengths-------------------------
      IF(NBODY.GT.0) THEN
         DO 80 N=1,NBODY
         IBODYM=0
         JBODYM=0
         IBODYP=0
         JBODYP=0
         XDIS1=XBDM(N)/ALSTR+1.0D-10
         XDIS2=XBDP(N)/ALSTR+1.0D-10
         YDIS1=YBDM(N)/ALSTR+1.0D-10
         YDIS2=YBDP(N)/ALSTR+1.0D-10
         DO 82 I=1,LI
         IF(X(1,I).LE.XDIS1) IBODYM=I
         IF(X(1,I).LE.XDIS2) IBODYP=I
   82    CONTINUE
         DO 84 J=1,LJ
         IF(Y(J,1).LE.YDIS1) JBODYM=J
         IF(Y(J,1).LE.YDIS2) JBODYP=J
   84    CONTINUE
         WRITE(16,86) N,XBDM(N),IBODYM,YBDM(N),JBODYM,
     1                  XBDP(N),IBODYP,YBDP(N),JBODYP
   86    FORMAT(//10X,'Inner Body-',I2,' is Inserted'/
     1   5X,'XMIN = ',F6.4,' (',I3,');   YMIN = ',F6.4,' (',I3,');'/
     2   5X,'XMAX = ',F6.4,' (',I3,');   YMAX = ',F6.4,' (',I3,')'/)
         IBDM(N)=IBODYM
         IBDP(N)=IBODYP
         JBDM(N)=JBODYM
         JBDP(N)=JBODYP
   80    CONTINUE
         END IF
C------------------------Inserted Injection Lengths---------------------
      IF(NFINJ.GT.0) THEN
         DO 90 N=1,NFINJ
         IFINJM=0
         JFINJM=0
         IFINJP=0
         JFINJP=0
         XDIS1=XFIM(N)/ALSTR+1.0D-10
         XDIS2=XFIP(N)/ALSTR+1.0D-10
         YDIS1=YFIM(N)/ALSTR+1.0D-10
         YDIS2=YFIP(N)/ALSTR+1.0D-10
         DO 92 I=1,LI
         IF(X(1,I).LE.XDIS1) IFINJM=I
         IF(X(1,I).LE.XDIS2) IFINJP=I
   92    CONTINUE
         DO 94 J=1,LJ
         IF(Y(J,1).LE.YDIS1) JFINJM=J
         IF(Y(J,1).LE.YDIS2) JFINJP=J
   94    CONTINUE
         WRITE(16,96) N,XFIM(N),IFINJM,YFIM(N),JFINJM,
     1                  XFIP(N),IFINJP,YFIP(N),JFINJP
   96    FORMAT(//10X,'Injection Line-',I2,' is Inserted'/
     1   5X,'XMIN = ',F6.4,' (',I3,');   YMIN = ',F6.4,' (',I3,');'/
     2   5X,'XMAX = ',F6.4,' (',I3,');   YMAX = ',F6.4,' (',I3,')'/)
         IFIM(N)=IFINJM
         IFIP(N)=IFINJP
         JFIM(N)=JFINJM
         JFIP(N)=JFINJP
   90    CONTINUE
         END IF
C---------------------Mesh Spacings for Calculations--------------------
      DO 74 I=2,LI-1
      XXC(I)=0.5*(X(1,I+1)-X(1,I-1))
      XXS(I)=(X(1,I)-X(1,I-1))
   74 CONTINUE
      XXC(1)=(X(1,2)-X(1,1))
      XXC(LI)=(X(1,LI)-X(1,LI-1))
      XXS(1)=XXC(1)
      XXS(LI)=XXC(LI)
      DO 76 J=2,LJ-1
      YYC(J)=0.5*(Y(J+1,1)-Y(J-1,1))
      YYS(J)=(Y(J,1)-Y(J-1,1))
   76 CONTINUE
      YYC(1)=(Y(2,1)-Y(1,1))
      YYC(LJ)=(Y(LJ,1)-Y(LJ-1,1))
      YYS(1)=YYC(1)
      YYS(LJ)=YYC(LJ)
      RETURN                                                            
      END                                                               
C***********************************************************************
      SUBROUTINE SKIPDF
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)                                   
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB09/ IBOT(LI),ITOP(LI),JLFT(LJ),JRGT(LJ),
     1       FBOT(8+LSP,LI),FTOP(8+LSP,LI),FLFT(8+LSP,LJ),FRGT(8+LSP,LJ)
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/FINJ/NFINJ,IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30)
C---------------Define Skip variable for flow Calculations--------------
C--------------- ISKIP=0 no skip
C--------------- ISKIP=1  skip epcilon calculations
C--------------- ISKIP=2  skip DT calculations
C--------------- ISKIP=3  Injection locations
C--------------- ISKIP=10 skip every calculations
C--------------- ISKIP=-9 introduce flame in the calculations
      DO 200 I=1,LI
      DO 200 J=1,LJ
      ISKIP(J,I)=0
  200 CONTINUE
      IF(NBODY.GE.1) THEN
          DO 202 N=1,NBODY
          IBODYM=IBDM(N)
          IBODYP=IBDP(N)
          JBODYM=JBDM(N)
          JBODYP=JBDP(N)
          DO 204 I=IBODYM,IBODYP
          IF(ISKIP(JBODYM,I).LT.10) THEN
                                    ISKIP(JBODYM,I)=11
                                    ELSE
                                    ISKIP(JBODYM,I)=10
                                    END IF
          IF(ISKIP(JBODYP,I).LT.10) THEN
                                    ISKIP(JBODYP,I)=11
                                    ELSE
                                    ISKIP(JBODYP,I)=10
                                    END IF
  204     CONTINUE
          DO 206 J=JBODYM,JBODYP
          IF(ISKIP(J,IBODYM).LT.10) THEN
                                    ISKIP(J,IBODYM)=11
                                    ELSE
                                    ISKIP(J,IBODYM)=10
                                    END IF
          IF(ISKIP(J,IBODYP).LT.10) THEN
                                    ISKIP(J,IBODYP)=11
                                    ELSE
                                    ISKIP(J,IBODYP)=10
                                    END IF
  206     CONTINUE
          DO 208 I=IBODYM+1,IBODYP-1
          DO 208 J=JBODYM+1,JBODYP-1
          ISKIP(J,I)=10
  208     CONTINUE
  202     CONTINUE
          END IF
      IF(NFINJ.GT.0) THEN
          DO 210 N=1,NFINJ
          IFINJM=IFIM(N)
          IFINJP=IFIP(N)
          JFINJM=JFIM(N)
          JFINJP=JFIP(N)
          DO 212 I=IFINJM,IFINJP
          DO J=JFINJM,JFINJP
          ISKIP(J,I)=3
          ENDDO
  212     CONTINUE
  210     CONTINUE
          END IF
      DO 222 I=2,LI-1
      ISKIP(1,I)=ISKIP(1,I)+1
      ISKIP(LJ,I)=ISKIP(LJ,I)+1
      IF(IBOT(I).EQ.1) THEN
                       ISKIP(2,I)=ISKIP(1,I)
                       ISKIP(1,I)=ISKIP(1,I)+1
                       END IF
      IF(ITOP(I).EQ.1) THEN
                       ISKIP(LJ-1,I)=ISKIP(LJ,I)
                       ISKIP(LJ,I)=ISKIP(LJ,I)+1
                       END IF
  222 CONTINUE
      DO 224 J=2,LJ-1
      ISKIP(J,1)=ISKIP(J,1)+1
      ISKIP(J,LI)=ISKIP(J,LI)+1
      IF(JLFT(J).EQ.1) THEN
                       ISKIP(J,2)=ISKIP(J,1)
                       ISKIP(J,1)=ISKIP(J,1)+1
                       END IF
      IF(JRGT(J).EQ.1) THEN
                       ISKIP(J,LI-1)=ISKIP(J,LI)
                       ISKIP(J,LI)=ISKIP(J,LI)+1
                       END IF
  224 CONTINUE
      ISKIP(1,1)=1
      ISKIP(1,LI)=1
      ISKIP(LJ,1)=1
      ISKIP(LJ,LI)=1
      RETURN                                                            
      END                                                               
C***********************************************************************
      SUBROUTINE LUMATX(KBOUND)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,
     1 LB=LJ-2,LL=(LI-2)*LB,LL1=LB*LL-(LB+1)*LB/2)
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/LUMAT/ IBAND(LB),PSL(LL1),PSD(LL)
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
      DO 608 J=2,LJ-1
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
  608 CONTINUE
      END IF
      IF(KBOUND.EQ.0) THEN
C------------------------  DIRICHLET CONDITION  ------------------------
      DO 609 I=2,LI-1
      DO 609 J=2,LJ-1
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
C***********************************************************************
      SUBROUTINE INFLOW(ICODE,INSP,
     1   RTIN,RTOT,ALENG,NSEG,ISIDE,ITYPE,BCVAL,
     2   FO2IN)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)                                   
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),
     1  W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDFSP(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/CB08/ TPOL1,TPOL2,POLSP(14,LSP),CISP(12,LSP)
      COMMON/CB09/ IBOT(LI),ITOP(LI),JLFT(LJ),JRGT(LJ),
     1       FBOT(8+LSP,LI),FTOP(8+LSP,LI),FLFT(8+LSP,LJ),FRGT(8+LSP,LJ)
      COMMON/FINJ/NFINJ,IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30)
      COMMON/HFORM/ HFO(LSP)
      COMMON/KEPS/ AKREF,EPREF,C1KE,C2KE,CMKE,CMKEQ,EE,AVON,PCON
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
      DIMENSION ISIDE(NSEG),ITYPE(NSEG),BCVAL(8+LSP,NSEG)
      SMALL=1.0D-10
      ICODE=0
      QSTR=(RSTR*VSTR*VSTR*VSTR)
C0000000000000000000000000000000000000000000000000000000000000000000000
C               Initialize the variables
C0000000000000000000000000000000000000000000000000000000000000000000000
C--------------------------REFERENCE VALUES-----------------------------
      HFFU=HFO(INSP)
      HFO2=HFO(02)
      HFSU1=HFO(06)
      HFN2=HFO(LSP)
      ININJ=06
C
      FN2IN=1.0-FO2IN
      TKD=TREF
      TKD2=TKD*TKD
      TKD3=TKD*TKD2
      TKD4=TKD*TKD3
      TKD5=TKD*TKD4
      IPOLY=1
      IF(TKD.GT.TPOL1) IPOLY=8
      RBAR=FO2IN/WM(2)+FN2IN/WM(LSP)
      DENS=PREF/(GASC*TKD*RBAR/WMSTR)
      HO2=(POLSP(IPOLY,2)*TKD+POLSP(IPOLY+1,2)*TKD2/2.0
     1   +POLSP(IPOLY+2,2)*TKD3/3.0+POLSP(IPOLY+3,2)*TKD4/4.0
     2   +POLSP(IPOLY+4,2)*TKD5/5.0+POLSP(IPOLY+5,2))*GASC/WM(2)/WMSTR
      HN2=(POLSP(IPOLY,LSP)*TKD+POLSP(IPOLY+1,LSP)*TKD2/2.0
     1+POLSP(IPOLY+2,LSP)*TKD3/3.0+POLSP(IPOLY+3,LSP)*TKD4/4.0
     2+POLSP(IPOLY+4,LSP)*TKD5/5.+POLSP(IPOLY+5,LSP))*GASC/WM(LSP)/WMSTR
      HTIN=HO2*FO2IN+HN2*FN2IN-HSTR*(HFO2*FO2IN+HFN2*FN2IN)
      DO 120 J=1,LJ
      DO 120 I=1,LI
      RHO(J,I)=RREF/RSTR
      U(J,I)=0.0
      V(J,I)=0.0
      W(J,I)=0.0
      P(J,I)=0.0
      HT(J,I)=HTIN/HSTR
      TK(J,I)=1.0
      AK(J,I)=AKREF/AKSTR
      EPS(J,I)=EPREF/EPSTR
      STMF(J,I)=0.0
      STND(J,I)=0.0
      STDF(J,I)=1.0
      VMU(J,I)=1.0
      VTC(J,I)=1.0
      DO 121 ISP=1,LSP
      VDFSP(J,I,ISP)=1.0
      FSP(J,I,ISP)=0.0
  121 CONTINUE
      FSP(J,I,2)=FO2IN
      FSP(J,I,LSP)=FN2IN
  120 CONTINUE
      IB=0
      JB=0
      IF(ISIDE(1).NE.1) THEN
                        WRITE(16,710) ISIDE(1)
                        ICODE=100
                        RETURN
                        END IF
  710 FORMAT(4X,'INPUT ERROR  ISIDE(1)=',I3)
      DO 10 N=1,NSEG
      IF((N.GE.2).AND.(ISIDE(N).NE.ISIDE(N-1))) THEN
                                                IB=0
                                                JB=0
                                                END IF
C11111111111111111111111111111111111111111111111111111111111111111111111
C----------   Here Input Data will be transfered on to I=1&LI  --------
C11111111111111111111111111111111111111111111111111111111111111111111111
      IF(ISIDE(N).LE.2) THEN
                  IA=IB+1
                  IF(IA.GE.LI) THEN
                               WRITE(16,712) N,ISIDE(N)
                               ICODE=101
                               RETURN
                               END IF
  712 FORMAT(4X,'INPUT ERROR I> LI ',2I3)
                  DO 12 I=1,LI
                  IF((X(1,I)*ALSTR+SMALL).GE.BCVAL(1,N)) GO TO 14
   12             CONTINUE
   14             CONTINUE
                  IB=I
                  IF(IB.GT.LI) IB=LI
                  DO 20 I=IA,IB
                  IF(ISIDE(N).EQ.1) THEN
                         IBOT(I)=ITYPE(N)
                         FBOT(1,I)=BCVAL(2,N)
                         END IF
                  IF(ISIDE(N).EQ.2) THEN
                         ITOP(I)=ITYPE(N)
                         FTOP(1,I)=BCVAL(2,N)
                         END IF
   20             CONTINUE
                  IF(ITYPE(N).EQ.0) THEN
                      TKD=BCVAL(6,N)
                      TKD2=TKD*TKD
                      TKD3=TKD*TKD2
                      TKD4=TKD*TKD3
                      TKD5=TKD*TKD4
                      IPOLY=1
                      IF(TKD.GT.TPOL1) IPOLY=8
                  FN2=1.0
                  RBAR=0.0
                  BCVAL(6,N)=0.5*(BCVAL(3,N)**2+BCVAL(4,N)**2
     1                      +BCVAL(5,N)**2)
                  DO 31 ISP=1,LSP-1
                  FN2=FN2-BCVAL(8+ISP,N)
                  RBAR=RBAR+BCVAL(8+ISP,N)/WM(ISP)
                  BCVAL(6,N)=BCVAL(6,N)+((POLSP(IPOLY,ISP)*TKD
     1                +POLSP(IPOLY+1,ISP)*TKD2/2.0
     2                +POLSP(IPOLY+2,ISP)*TKD3/3.0
     3                +POLSP(IPOLY+3,ISP)*TKD4/4.0
     4                +POLSP(IPOLY+4,ISP)*TKD5/5.0
     5                +POLSP(IPOLY+5,ISP))*GASC/WM(ISP)/WMSTR
     6                -HSTR*HFO(ISP))*BCVAL(8+ISP,N)
   31             CONTINUE
                  RBAR=RBAR+FN2/WM(LSP)
                  DENST=PREF/(GASC*TKD*RBAR/WMSTR)
                  BCVAL(6,N)=BCVAL(6,N)+((POLSP(IPOLY,LSP)*TKD
     1                +POLSP(IPOLY+1,LSP)*TKD2/2.0
     2                +POLSP(IPOLY+2,LSP)*TKD3/3.0
     3                +POLSP(IPOLY+3,LSP)*TKD4/4.0
     4                +POLSP(IPOLY+4,LSP)*TKD5/5.0
     5                +POLSP(IPOLY+5,LSP))*GASC/WM(LSP)/WMSTR
     6                -HSTR*HFO(LSP))*FN2
                      DO 22 I=IA,IB
                      IF(ISIDE(N).EQ.1) THEN
                            FBOT(2,I)=DENST/RSTR
                            FBOT(3,I)=BCVAL(3,N)/VSTR
                            FBOT(4,I)=BCVAL(4,N)/VSTR
                            FBOT(5,I)=BCVAL(5,N)/VSTR
                            FBOT(6,I)=BCVAL(6,N)/HSTR
                            FBOT(7,I)=BCVAL(7,N)/AKSTR
                            FBOT(8,I)=BCVAL(8,N)/EPSTR
                            JJ=1
                            END IF
                      IF(ISIDE(N).EQ.2) THEN
                            FTOP(2,I)=DENST/RSTR
                            FTOP(3,I)=BCVAL(3,N)/VSTR
                            FTOP(4,I)=BCVAL(4,N)/VSTR
                            FTOP(5,I)=BCVAL(5,N)/VSTR
                            FTOP(6,I)=BCVAL(6,N)/HSTR
                            FTOP(7,I)=BCVAL(7,N)/AKSTR
                            FTOP(8,I)=BCVAL(8,N)/EPSTR
                            JJ=LJ
                            END IF
                      RHO(JJ,I)=DENST/RSTR
                      U(JJ,I)=BCVAL(3,N)/VSTR
                      V(JJ,I)=BCVAL(4,N)/VSTR
                      W(JJ,I)=BCVAL(5,N)/VSTR
                      HT(JJ,I)=BCVAL(6,N)/HSTR
                      AK(JJ,I)=BCVAL(7,N)/AKSTR
                      EPS(JJ,I)=BCVAL(8,N)/EPSTR
                      FTOT=0.0
                      DO 25 ISP=1,LSP-1
                      IF(ISIDE(N).EQ.1) FBOT(8+ISP,I)=BCVAL(8+ISP,N)
                      IF(ISIDE(N).EQ.2) FTOP(8+ISP,I)=BCVAL(8+ISP,N)
                      FSP(JJ,I,ISP)=BCVAL(8+ISP,N)
                      FTOT=FTOT+BCVAL(8+ISP,N)
   25                 CONTINUE
                      FSP(JJ,I,LSP)=1.0-FTOT
                      IF(ISIDE(N).EQ.1) FBOT(8+LSP,I)=1.0-FTOT
                      IF(ISIDE(N).EQ.2) FTOP(8+LSP,I)=1.0-FTOT
                      TK(JJ,I)=TKD/TSTR
   22                 CONTINUE
                      END IF
                  IF(ITYPE(N).EQ.1) THEN
                      DO 26 I=IA,IB
                      IF(ISIDE(N).EQ.1) THEN
                                        FBOT(2,I)=BCVAL(3,N)
                                        FBOT(3,I)=BCVAL(4,N)/TSTR
                                        FBOT(4,I)=BCVAL(5,N)/QSTR
                                        END IF
                      IF(ISIDE(N).EQ.2) THEN
                                        FTOP(2,I)=BCVAL(3,N)
                                        FTOP(3,I)=BCVAL(4,N)/TSTR
                                        FTOP(4,I)=BCVAL(5,N)/QSTR
                                        END IF
   26                 CONTINUE
                      END IF
           END IF
      IF(ISIDE(N).GE.3) THEN
C22222222222222222222222222222222222222222222222222222222222222222222222
C----------   Here Input Data will be transferred on to J=1&LJ  --------
C22222222222222222222222222222222222222222222222222222222222222222222222
                  JA=JB+1
                  IF(JA.GE.LJ) THEN
                               WRITE(16,712) N,ISIDE(N)
                               ICODE=102
                               RETURN
                               END IF
                  DO 16 J=1,LJ
                  IF((Y(J,1)*ALSTR+SMALL).GE.BCVAL(1,N)) GO TO 18
   16             CONTINUE
   18             CONTINUE
                  JB=J
                  IF(JB.GT.LJ) JB=LJ
                  DO 30 J=JA,JB
                  IF(ISIDE(N).EQ.3) THEN
                         JLFT(J)=ITYPE(N)
                         FLFT(1,J)=BCVAL(2,N)
                         END IF
                  IF(ISIDE(N).EQ.4) THEN
                         JRGT(J)=ITYPE(N)
                         FRGT(1,J)=BCVAL(2,N)
                         END IF
   30             CONTINUE
                  IF(ITYPE(N).EQ.0) THEN
                      TKD=BCVAL(6,N)
                      TKD2=TKD*TKD
                      TKD3=TKD*TKD2
                      TKD4=TKD*TKD3
                      TKD5=TKD*TKD4
                      IPOLY=1
                      IF(TKD.GT.TPOL1) IPOLY=8
                  FN2=1.0
                  RBAR=0.0
                  BCVAL(6,N)=0.5*(BCVAL(3,N)**2+BCVAL(4,N)**2
     1                      +BCVAL(5,N)**2)
                  DO 33 ISP=1,LSP-1
                  FN2=FN2-BCVAL(8+ISP,N)
                  RBAR=RBAR+BCVAL(8+ISP,N)/WM(ISP)
                  BCVAL(6,N)=BCVAL(6,N)+((POLSP(IPOLY,ISP)*TKD
     1                +POLSP(IPOLY+1,ISP)*TKD2/2.0
     2                +POLSP(IPOLY+2,ISP)*TKD3/3.0
     3                +POLSP(IPOLY+3,ISP)*TKD4/4.0
     4                +POLSP(IPOLY+4,ISP)*TKD5/5.0
     5                +POLSP(IPOLY+5,ISP))*GASC/WM(ISP)/WMSTR
     6                -HSTR*HFO(ISP))*BCVAL(8+ISP,N)
   33             CONTINUE
                  RBAR=RBAR+FN2/WM(LSP)
                  DENST=PREF/(GASC*TKD*RBAR/WMSTR)
                  BCVAL(6,N)=BCVAL(6,N)+((POLSP(IPOLY,LSP)*TKD
     1                +POLSP(IPOLY+1,LSP)*TKD2/2.0
     2                +POLSP(IPOLY+2,LSP)*TKD3/3.0
     3                +POLSP(IPOLY+3,LSP)*TKD4/4.0
     4                +POLSP(IPOLY+4,LSP)*TKD5/5.0
     5                +POLSP(IPOLY+5,LSP))*GASC/WM(LSP)/WMSTR
     6                -HSTR*HFO(LSP))*FN2
C---------------- FLAT-HAT PARABOLIC VELOCITY PROFILE ------------------
                      R1BR0=(1.0-DABS(BCVAL(2,N)))
                      R1BR0S=R1BR0*R1BR0
                      VA=BCVAL(3,N)
                      VB=VA
                      IF(BCVAL(2,N).NE.0.0) THEN
                        VA=2.0*BCVAL(3,N)/(1.0-R1BR0S*R1BR0S)
                        VB=VA*(1.0-R1BR0S)
                        END IF
                      DO 32 J=JA,JB
                      JBB=JB+1
                      IF(JBB.GT.LJ) JBB=LJ
                      YBLYR=0.5*(Y(JB,1)+Y(JBB,1))-Y(JA,1)
                      YBLYR=Y(JB,1)-Y(JA,1)
                      RBR0=(Y(J,1)-Y(JA,1))/YBLYR
                      IF(BCVAL(2,N).LT.0.0) THEN
                         JAA=JA-1
                         IF(JAA.LT.1) JAA=1
                         YBLYR=Y(JB,1)-0.5*(Y(JA,1)+Y(JAA,1))
                         RBR0=(Y(JB,1)-Y(J,1))/YBLYR
                         END IF
                      PBVEL=VB
                      IF(RBR0.GT.R1BR0) PBVEL=VA*(1.0-RBR0*RBR0)
                      HDYNA=0.5*(PBVEL**2+BCVAL(4,N)**2+BCVAL(4,N)**2)
                      IF(ISIDE(N).EQ.3) THEN
                            FLFT(2,J)=DENST/RSTR
                            FLFT(3,J)=PBVEL/VSTR
                            FLFT(4,J)=BCVAL(4,N)/VSTR
                            FLFT(5,J)=BCVAL(5,N)/VSTR
                            FLFT(6,J)=(BCVAL(6,N)+HDYNA)/HSTR
                            FLFT(7,J)=BCVAL(7,N)/AKSTR
                            FLFT(8,J)=BCVAL(8,N)/EPSTR
                            II=1
                            END IF
                      IF(ISIDE(N).EQ.4) THEN
                            FRGT(2,J)=DENST/RSTR
                            FRGT(3,J)=PBVEL/VSTR
                            FRGT(4,J)=BCVAL(4,N)/VSTR
                            FRGT(5,J)=BCVAL(5,N)/VSTR
                            FRGT(6,J)=(BCVAL(6,N)+HDYNA)/HSTR
                            FRGT(7,J)=BCVAL(7,N)/AKSTR
                            FRGT(8,J)=BCVAL(8,N)/EPSTR
                            II=LI
                            END IF
                       RHO(J,II)=DENST/RSTR
                       U(J,II)=PBVEL/VSTR
                       V(J,II)=BCVAL(4,N)/VSTR
                       W(J,II)=BCVAL(5,N)/VSTR
                       HT(J,II)=(BCVAL(6,N)+HDYNA)/HSTR
                       AK(J,II)=BCVAL(7,N)/AKSTR
                       EPS(J,II)=BCVAL(8,N)/EPSTR
                       FTOT=0.0
                       DO 35 ISP=1,LSP-1
                       IF(ISIDE(N).EQ.3) FLFT(8+ISP,J)=BCVAL(8+ISP,N)
                       IF(ISIDE(N).EQ.4) FRGT(8+ISP,J)=BCVAL(8+ISP,N)
                       FSP(J,II,ISP)=BCVAL(8+ISP,N)
                       FTOT=FTOT+BCVAL(8+ISP,N)
   35                  CONTINUE
                       FSP(J,II,LSP)=1.0-FTOT
                       IF(ISIDE(N).EQ.3) FLFT(8+LSP,J)=1.0-FTOT
                       IF(ISIDE(N).EQ.4) FRGT(8+LSP,J)=1.0-FTOT
                       TK(J,II)=TKD/TSTR
   32                 CONTINUE
                      END IF
                  IF(ITYPE(N).EQ.1) THEN
                      DO 36 J=JA,JB
                      IF(ISIDE(N).EQ.3) THEN
                                        FLFT(2,J)=BCVAL(3,N)
                                        FLFT(3,J)=BCVAL(4,N)/TSTR
                                        FLFT(4,J)=BCVAL(5,N)/QSTR
                                        END IF
                      IF(ISIDE(N).EQ.4) THEN
                                        FRGT(2,J)=BCVAL(3,N)
                                        FRGT(3,J)=BCVAL(4,N)/TSTR
                                        FRGT(4,J)=BCVAL(5,N)/QSTR
                                        END IF
   36                 CONTINUE
                      END IF
           END IF
   10 CONTINUE
C22222222222222222222222222222222222222222222222222222222222222222222222
C---------------    Here Initial Data is Generated   -------------------
C22222222222222222222222222222222222222222222222222222222222222222222222
      DO 130 J=1,LJ                                                   
      IA=LI
      IB=LI
      IF((JLFT(J).EQ.0).AND.(JRGT(J).EQ.0)) IA=LI/2
      IF((JLFT(J).NE.0).AND.(JRGT(J).EQ.0)) IA=0
      DO 132 I=2,IA
      U(J,I)=U(J,1)
      V(J,I)=V(J,1)
      W(J,I)=W(J,1)
      P(J,I)=0.0
      IF(JLFT(J).EQ.0) THEN
                       RHO(J,I)=RHO(J,1)
                       HT(J,I)=HT(J,1)
                       AK(J,I)=AK(J,1)
                       EPS(J,I)=EPS(J,1)
                       TK(J,I)=TK(J,1)
                       VMU(J,I)=VMU(J,1)
                       DO 131 ISP=1,LSP
                       FSP(J,I,ISP)=FSP(J,1,ISP)
  131                  CONTINUE
                       END IF
  132 CONTINUE
      DO 134 I=IA+1,IB-1
      RHO(J,I)=RHO(J,LI)
      U(J,I)=U(J,LI)
      V(J,I)=V(J,LI)
      W(J,I)=W(J,LI)
      P(J,I)=0.0
      HT(J,I)=HT(J,LI)
      AK(J,I)=AK(J,LI)
      EPS(J,I)=EPS(J,LI)
      TK(J,I)=TK(J,LI)
      VMU(J,I)=VMU(J,LI)
      DO 135 ISP=1,LSP
      FSP(J,I,ISP)=FSP(J,LI,ISP)
  135 CONTINUE
  134 CONTINUE
  130 CONTINUE
C--------------------------INFLOW FOR INJECTIONS------------------------
      IF(NFINJ.GE.1) THEN
            HFINJ=HFO(06)
            DO 200 N=1,NFINJ
            TKD=FFI(N,4)
            TKD2=TKD*TKD
            TKD3=TKD*TKD2
            TKD4=TKD*TKD3
            TKD5=TKD*TKD4
            IPOLY=1
            IF(TKD.GT.TPOL1) IPOLY=8
            FN2=1.0-(FFI(N,7)+FFI(N,8))
            RBAR=FFI(N,7)/WM(INSP)+FFI(N,8)/WM(2)+FN2/WM(LSP)
            DENS=PREF/(GASC*TKD*RBAR/WMSTR)
            HFU=(POLSP(IPOLY,INSP)*TKD+POLSP(IPOLY+1,INSP)*TKD2/2.0
     1       +POLSP(IPOLY+2,INSP)*TKD3/3.0+POLSP(IPOLY+3,INSP)*TKD4/4.0
     2       +POLSP(IPOLY+4,INSP)*TKD5/5.0
     3       +POLSP(IPOLY+5,INSP))*GASC/WM(INSP)/WMSTR
            HO2=(POLSP(IPOLY,2)*TKD+POLSP(IPOLY+1,2)*TKD2/2.0
     1       +POLSP(IPOLY+2,2)*TKD3/3.0+POLSP(IPOLY+3,2)*TKD4/4.0
     2       +POLSP(IPOLY+4,2)*TKD5/5.0
     3       +POLSP(IPOLY+5,2))*GASC/WM(2)/WMSTR
            HINJ=(POLSP(IPOLY,ININJ)*TKD+POLSP(IPOLY+1,ININJ)*TKD2/2.0
     1       +POLSP(IPOLY+2,ININJ)*TKD3/3.0
     2       +POLSP(IPOLY+3,ININJ)*TKD4/4.0
     3       +POLSP(IPOLY+4,ININJ)*TKD5/5.0
     4       +POLSP(IPOLY+5,ININJ))*GASC/WM(ININJ)/WMSTR
            HN2=(POLSP(IPOLY,LSP)*TKD+POLSP(IPOLY+1,LSP)*TKD2/2.0
     1       +POLSP(IPOLY+2,LSP)*TKD3/3.0+POLSP(IPOLY+3,LSP)*TKD4/4.0
     2       +POLSP(IPOLY+4,LSP)*TKD5/5.0
     3       +POLSP(IPOLY+5,LSP))*GASC/WM(LSP)/WMSTR
            ENTH=HFU*FFI(N,7)+HO2*FFI(N,8)+HINJ*FFI(N,9)+HN2*FN2
     1       -HSTR*(HFFU*FFI(N,7)+HFO2*FFI(N,8)+HFINJ*FFI(N,9)+HFN2*FN2)
            HDYNA=0.5*(FFI(N,1)**2+FFI(N,2)**2+FFI(N,3)**2)
            FFI(N,1)=FFI(N,1)/VSTR
            FFI(N,2)=FFI(N,2)/VSTR
            FFI(N,3)=FFI(N,3)/VSTR
            FFI(N,4)=FFI(N,4)/TSTR
            FFI(N,5)=FFI(N,5)/AKSTR
            FFI(N,6)=FFI(N,6)/EPSTR
            FFI(N,10)=FN2
            FFI(N,11)=DENS/RSTR
            FFI(N,12)=0.0
            FFI(N,13)=(ENTH+HDYNA)/HSTR
            IFINJM=IFIM(N)
            IFINJP=IFIP(N)
            JFINJM=JFIM(N)
            JFINJP=JFIP(N)
            DO 202 I=IFINJM,IFINJP
            DO 202 J=JFINJM,JFINJP
            U(J,I)=FFI(N,1)
            V(J,I)=FFI(N,2)
            W(J,I)=FFI(N,3)
            TK(J,I)=FFI(N,4)
            AK(J,I)=FFI(N,5)
            EPS(J,I)=FFI(N,6)
            DO 204 ISP=1,LSP
            FSP(J,I,ISP)=0.0
  204       CONTINUE
            FSP(J,I,INSP)=FFI(N,7)
            FSP(J,I,2)=FFI(N,8)
            FSP(J,I,ININJ)=FFI(N,9)
            FSP(J,I,LSP)=FFI(N,10)
            RHO(J,I)=FFI(N,11)
            P(J,I)=FFI(N,12)
            HT(J,I)=FFI(N,13)
  202       CONTINUE
C--------------------------for boundary conditions
            IF(FFI(N,2).GT.0.0) THEN
            DO 210 I=IFINJM,IFINJP
            DO 211 J=JFINJP+1,LJ-1
            IF(ISKIP(J,I).GT.0) GO TO 212
            U(J,I)=U(JFINJP,I)
            V(J,I)=V(JFINJP,I)
            W(J,I)=W(JFINJP,I)
            RHO(J,I)=RHO(JFINJP,I)
            HT(J,I)=HT(JFINJP,I)
            TK(J,I)=TK(JFINJP,I)
            AK(J,I)=AK(JFINJP,I)
            EPS(J,I)=EPS(JFINJP,I)
            DO 211 ISP=1,LSP
            FSP(J,I,ISP)=FSP(JFINJP,I,ISP)
  211       CONTINUE
  212       CONTINUE
  210       CONTINUE
            END IF
            IF(FFI(N,2).LT.0.0) THEN
            DO 220 I=IFINJM,IFINJP
            DO 221 J=JFINJM-1,2,-1
            IF(ISKIP(J,I).GT.0) GO TO 222
            U(J,I)=U(JFINJP,I)
            V(J,I)=V(JFINJP,I)
            W(J,I)=W(JFINJP,I)
            RHO(J,I)=RHO(JFINJP,I)
            HT(J,I)=HT(JFINJP,I)
            TK(J,I)=TK(JFINJP,I)
            AK(J,I)=AK(JFINJP,I)
            EPS(J,I)=EPS(JFINJP,I)
            DO 221 ISP=1,LSP
            FSP(J,I,ISP)=FSP(JFINJP,I,ISP)
  221       CONTINUE
  222       CONTINUE
  220       CONTINUE
            END IF
            IF(FFI(N,1).GT.0.0) THEN
            DO 230 J=JFINJM,JFINJP
            DO 231 I=IFINJP+1,LI-1
            IF(ISKIP(J,I).GT.0) GO TO 232
            U(J,I)=U(J,IFINJP)
            V(J,I)=V(J,IFINJP)
            W(J,I)=W(J,IFINJP)
            RHO(J,I)=RHO(J,IFINJP)
            HT(J,I)=HT(J,IFINJP)
            TK(J,I)=TK(J,IFINJP)
            AK(J,I)=AK(J,IFINJP)
            EPS(J,I)=EPS(J,IFINJP)
            DO 231 ISP=1,LSP
            FSP(J,I,ISP)=FSP(J,IFINJP,ISP)
  231       CONTINUE
  232       CONTINUE
  230       CONTINUE
            END IF
            IF(FFI(N,1).LT.0.0) THEN
            DO 240 J=JFINJM,JFINJP
            DO 241 I=IFINJM-1,2,-1
            IF(ISKIP(J,I).GT.0) GO TO 242
            U(J,I)=U(J,IFINJP)
            V(J,I)=V(J,IFINJP)
            W(J,I)=W(J,IFINJP)
            RHO(J,I)=RHO(J,IFINJP)
            HT(J,I)=HT(J,IFINJP)
            TK(J,I)=TK(J,IFINJP)
            AK(J,I)=AK(J,IFINJP)
            EPS(J,I)=EPS(J,IFINJP)
            DO 241 ISP=1,LSP
            FSP(J,I,ISP)=FSP(J,IFINJP,ISP)
  241       CONTINUE
  242       CONTINUE
  240       CONTINUE
            END IF
  200    CONTINUE
         END IF
      RETURN                                                            
      END 
C***********************************************************************
      SUBROUTINE RESET                 
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LIJ=LI+LJ,LSP=52,ID1=7+LE,
     1   ID2=4+LSP,ID3=2*LE+2*LIJ+1,ID5=(9+LSP)*LE,IDNT=LE,IDN=20,
     2   IDV=(3+LSP)*LE,IDB1=2*LIJ,IDB2=2*(8+LSP)*LIJ,IDF=2*4*LIJ,
     3   LB=LJ-2,LL=(LI-2)*LB,IDLU=LB*LL-(LB+1)*LB/2+LL,LPD=55*LE)
      COMMON/CB01/ IDAT(ID1)                      
      COMMON/CB02/ DAT2(ID2)        
      COMMON/CB03/ DAT3(ID3)         
      COMMON/CB04/ DAT5(ID5)    
      COMMON/CB05/ DATNT(IDNT)
      COMMON/CB06/ DATN(IDN)
      COMMON/CB07/ DATV(IDV)
      COMMON/CB08/ DATC(2+26*LSP)
      COMMON/CB09/ IDATB(IDB1),DATB(IDB2)
      COMMON/CB10/ DATF(IDF)
      COMMON/LUMAT/ IDATLU(LB),DATLU(IDLU)
      COMMON/DUMMY/ DATD(LPD)
      DO 101 I=1,ID1
      IDAT(I)=0
  101 CONTINUE
      DO 102 I=1,ID2
      DAT2(I)=0.0
  102 CONTINUE
      DO 103 I=1,ID3
      DAT3(I)=0.0
  103 CONTINUE
      DO 105 I=1,ID5
      DAT5(I)=0.0
  105 CONTINUE
      DO 108 I=1,IDN
      DATN(I)=0.0
  108 CONTINUE
      DO 109 I=1,IDV
      DATV(I)=0.0
  109 CONTINUE
      DO 112 I=1,IDB1
      IDATB(I)=0
  112 CONTINUE
      DO 113 I=1,IDB2
      DATB(I)=0.0
  113 CONTINUE
      DO 114 I=1,IDF
      DATF(I)=0.0
  114 CONTINUE
      DO 115 I=1,LB
      IDATLU(I)=0
  115 CONTINUE
      DO 116 I=1,IDLU
      DATLU(I)=0.0
  116 CONTINUE
      DO 117 I=1,IDNT
      DATNT(I)=0.0
  117 CONTINUE
      DO 120 I=1,LPD
      DATD(I)=0.0
  120 CONTINUE
      RETURN                                                              
      END                                                               
C***********************************************************************
      SUBROUTINE PLOTS(NOPR,NOINJ,IPANM,
     1                 X1ANM,Y1ANM,X2ANM,Y2ANM,NFSURF,KORNT)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LSP=52,
     1   IPMAX=500,JPMAX=200,LNPR=200,LNMX=2000)
      CHARACTER *2 BFL2A(10),BFL2B(100)
      CHARACTER *5 BFL1,BFL2*4,BFL3*4,BBFILE*13
      CHARACTER*14 FRMTSTR
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),W(LJ,LI),
     1             P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1 HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/TRACK/ IPICK(101),JPICK(101),XPR(LNPR),YPR(LNPR),
     1              PATHX(LNMX,LNPR),PATHY(LNMX,LNPR),
     2              PUVEL(LNMX,LNPR),PVVEL(LNMX,LNPR),
     3              GUVEL(LNMX,LNPR),GVVEL(LNMX,LNPR)
      COMMON/FLOW/ KSYM,LIP,LJP,IOFF,JOFF,NCLR(LNPR),
     1  XS(LI),YS(LJ+LJ-1),RASMN,RASMX,RAS2MN,RAS2MX
      COMMON/APIXL/ IPAGE,NCOLOR,NBACK,NAA(IPMAX,JPMAX)
      COMMON/FRAME/ IPIXL,JPIXL,XMIN,YMIN,XMAX,YMAX,XPIXL,YPIXL
      COMMON/PEN/NPEN,NAPEN,NBPEN
      DATA BFL1,BFL3 /'movie','.bmp'/
      DATA BFL2A / '-0','-1','-2','-3','-4','-5','-6','-7','-8','-9'/
      DATA BFL2B / '00','01','02','03','04','05','06','07','08','09',
     1             '10','11','12','13','14','15','16','17','18','19',
     2             '20','21','22','23','24','25','26','27','28','29',
     3             '30','31','32','33','34','35','36','37','38','39',
     4             '40','41','42','43','44','45','46','47','48','49',
     5             '50','51','52','53','54','55','56','57','58','59',
     6             '60','61','62','63','64','65','66','67','68','69',
     7             '70','71','72','73','74','75','76','77','78','79',
     8             '80','81','82','83','84','85','86','87','88','89',
     9             '90','91','92','93','94','95','96','97','98','99'/
      CHARACTER *1 BB(IPMAX,JPMAX)
      DIMENSION XBDM(10),YBDM(10),XBDP(10),YBDP(10)
      DIMENSION FS(LI,LJ+LJ-1),XFL(LI),YFL(LI)
      DIMENSION FS2(LI,LJ)
      DIMENSION XSTG(LI),YSTG(LJ)
      IPIXL=IPMAX
      JPIXL=JPMAX
      NCLOT=255
      NPANM=IPANM
      IF(NPANM.GE.1) THEN
         II=NPANM/10
         IDOT=NPANM-II*10
         NPANM=II
         END IF
      NAPEN=1
      NBPEN=0
      ISVEC=4
      JSVEC=4
      ARP=0.0
      ARL=0.7
      ARW=0.15
      ARMAX=1.8
      VOFF=0.0
      AMAG=0.7
C     XSHIFT=0.5*(X(1,1)+X(1,2))
C     YSHIFT=0.5*(Y(1,1)+Y(2,1))
      XSHIFT=0.0
      YSHIFT=0.0
      DO 14 I=1,LI-1
      XSTG(I+1)=0.5*(X(1,I)+X(1,I+1))-XSHIFT
   14 CONTINUE
      XSTG(1)=XSTG(2)
      DO 15 J=1,LJ-1
      YSTG(J+1)=0.5*(Y(J,1)+Y(J+1,1))-YSHIFT
   15 CONTINUE
      YSTG(1)=YSTG(2)
      IF(IPAGE.GT.1) GO TO 100
      IPAGE=IPAGE+2
      XBMN=X1ANM*1000.0
      XBMX=X2ANM*1000.0
      YBMN=Y1ANM*1000.0
      YBMX=Y2ANM*1000.0
      KCNTL=0
      NBACK=0
      IF(KSYM.EQ.1) YBMN=0.0
      AMUL=ALSTR*1000.0
      IOFF=0
      JOFF=0
              DO 16 I=1,LI
              IF(XSTG(I)*AMUL.LE.XBMN) IOFF=I-1
              IF(XSTG(I)*AMUL.LE.XBMX) ILI=I
   16         CONTINUE
              DO 17 J=1,LJ
              IF(YSTG(J)*AMUL.LE.YBMN) JOFF=J-1
              IF(YSTG(J)*AMUL.LE.YBMX) JLJ=J
   17         CONTINUE
              IF(ILI.LE.1) ILI=LI
              IF(JLJ.LE.1) JLJ=LJ
              IF(IOFF.GE.ILI) IOFF=0
              IF(JOFF.GE.JLJ) JOFF=0
              LIP=ILI-IOFF
              LJP=JLJ-JOFF
      WRITE(16,18) X(1,IOFF+1)*AMUL,X(1,IOFF+LIP)*AMUL,
     1            Y(JOFF+1,1)*AMUL,Y(JOFF+LJP,1)*AMUL
   18 FORMAT(1X,79('*')/1X,79('*')/1X,'**',75X,'**'/
     1 1X,'**',28X,'---PLOTTING AREA---',28X,'**'/1X,'**',75X,'**'/
     2 1X,'**',20X,'XMIN= ',F6.1,' mm',5X,'XMAX= ',F6.1,' mm',20X,'**'/
     3 1X,'**',20X,'YMIN= ',F6.1,' mm',5X,'YMAX= ',F6.1,' mm',20X,'**'/
     4 1X,'**',75X,'**'/1X,79('*')/1X,79('*')//)
      LJ2=LJP+LJP-1
      LJJ=LJP
      IF(KSYM.EQ.1) LJJ=LJ2
      DO 21 I=1,LIP
      XS(I)=XSTG(I+IOFF)
   21 CONTINUE
      DO 22 J=1,LJP
      YS(J)=YSTG(J+JOFF)
   22 CONTINUE
      IF(KSYM.EQ.1) THEN
                    DO 24 J=1,LJP
                    JA=LJP-J+1
                    YS(J)=-YSTG(JA)
   24               CONTINUE
                    DO 26 J=LJP+1,LJ2
                    JA=J-LJP+1
                    YS(J)=YSTG(JA)
   26               CONTINUE
                    END IF
      XDMIN=1.0D+10
      YDMIN=1.0D+10
      XDMAX=-1.0D+10
      YDMAX=-1.0D+10
      DO 30 I=1,LIP
      IF(XS(I).GT.XDMAX) XDMAX=XS(I)
      IF(XS(I).LT.XDMIN) XDMIN=XS(I)
   30 CONTINUE
      DO 31 J=1,LJJ
      IF(YS(J).GT.YDMAX) YDMAX=YS(J)
      IF(YS(J).LT.YDMIN) YDMIN=YS(J)
   31 CONTINUE
      XMINT=XDMIN-0.05*(XDMAX-XDMIN)/0.9
      XMAXT=XDMAX+0.05*(XDMAX-XDMIN)/0.9
      YMINT=YDMIN-0.05*(YDMAX-YDMIN)/0.9
      YMAXT=YDMAX+0.05*(YDMAX-YDMIN)/0.9
      XPIXLT=FLOAT(IPIXL)/DABS(XMAXT-XMINT)
      YPIXLT=FLOAT(JPIXL)/DABS(YMAXT-YMINT)
      A1=DABS(YMAXT-YMINT)*XPIXLT
      XPIXL=XPIXLT
      YPIXL=XPIXL
      IF(A1.GT.DFLOAT(JPIXL)) THEN 
                             YPIXL=YPIXLT
                             XPIXL=YPIXL
                             END IF
      IF(KCNTL.EQ.1) THEN
                     XPIXL=XPIXLT
                     YPIXL=YPIXLT
                     END IF
      XMIN=XMINT+0.5*(1.0-XPIXLT/XPIXL)*(XMAXT-XMINT)
      XMAX=XMINT+0.5*(1.0+XPIXLT/XPIXL)*(XMAXT-XMINT)
      YMIN=YMINT+0.5*(1.0-YPIXLT/YPIXL)*(YMAXT-YMINT)
      YMAX=YMINT+0.5*(1.0+YPIXLT/YPIXL)*(YMAXT-YMINT)
      NCLR(1)=255
      IF(NPANM.EQ.1) THEN
        DO 32 N=2,NOPR
        NCLR(N)=1+INT(254.0*FLOAT(N-1)/FLOAT(NOPR-1))
   32   CONTINUE
        END IF
      IF(NPANM.EQ.2) THEN
        DO 33 N=2,LNMX
        NCLR(N)=1+INT(254.0*FLOAT(N-1)/FLOAT(LNMX-1))
   33   CONTINUE
        END IF
      IF(NPANM.EQ.0) THEN
        DO 34 N=2,NOPR
        NCLR(N)=255
   34   CONTINUE
        END IF
      RASMX=0.0
      RASMN=1.0D+10
      RAS2MX=0.0
      RAS2MN=1.0D+10
      FVMAX=0.0
      DO 36 I=1,LIP-1
      DO 36 J=1,LJP-1
      IA=I+IOFF
      JA=J+JOFF
      TKSTG=0.25*(TK(JA,IA)+TK(JA+1,IA)+TK(JA,IA+1)+TK(JA+1,IA+1))
      IF(TKSTG.GT.RASMX) RASMX=TKSTG
      IF(TKSTG.LT.RASMN) RASMN=TKSTG
C---------------SOOT VOLUME FRACTION------------------------------------
C     FS2STG=0.25*(RHO(JA,IA)*STMF(JA,IA)+RHO(JA+1,IA)*STMF(JA+1,IA)
C    1 +RHO(JA,IA+1)*STMF(JA,IA+1)+RHO(JA+1,IA+1)*STMF(JA+1,IA+1))
      FS2STG=0.25*(FSP(JA,IA,05)+FSP(JA+1,IA,05)
     1 +FSP(JA,IA+1,05)+FSP(JA+1,IA+1,05))
      IF(FS2STG.GT.RAS2MX) RAS2MX=FS2STG
      IF(FS2STG.LT.RAS2MN) RAS2MN=FS2STG
      VV=0.5*DSQRT((U(JA,IA+1)+U(JA+1,IA+1)-2.0*VOFF)**2
     1       +(V(JA+1,IA)+V(JA+1,IA+1))**2
     2   +0.25*(W(JA,IA)+W(JA+1,IA)+W(JA,IA+1)+W(JA+1,IA+1))**2)
      IF(VV.GT.FVMAX) FVMAX=VV  
   36 CONTINUE
  100 CONTINUE
      RASCON=(RAS2MX-RAS2MN)/(RASMX-RASMN)
      LJ2=LJP+LJP-1
      LJJ=LJP
      IF(KSYM.EQ.1) LJJ=LJ2
C--------------------------LOCATE THE FLAME-----------------------------
C-------------------- FOR FINITE-RATE CHEMISTRY MODELS -----------------
      IF(NFSURF.GE.1) THEN
            XFLP=0.0
            YFLP=1.0
            DO 1326 I=1,LI-1
            RFLAME=0.0
            TMAX=0.0
            JJ=1
            DO 1328 J=1,LJ-1
            TKSTG=0.25*(TK(J,I)+TK(J+1,I)+TK(J,I+1)+TK(J+1,I+1))
            IF(TKSTG.GE.TMAX) THEN
                              TMAX=TKSTG
                              JJ=J
                              END IF
 1328       CONTINUE
            IF(JJ.GT.LJ) JJ=LJP
            XFLN=XSTG(I)
            IF(JJ.LE.1) THEN
                        YFLN=0.0
                        ELSE
                        Y2=YSTG(JJ)-YSTG(JJ-1)
                        Y3=YSTG(JJ+1)-YSTG(JJ-1)
                        TKSTGC=0.25*(TK(JJ,I)+TK(JJ+1,I)
     1                        +TK(JJ,I+1)+TK(JJ+1,I+1))
                        TKSTGP=0.25*(TK(JJ+1,I)+TK(JJ+2,I)
     1                        +TK(JJ+1,I+1)+TK(JJ+2,I+1))
                        TKSTGM=0.25*(TK(JJ-1,I)+TK(JJ,I)
     1                        +TK(JJ-1,I+1)+TK(JJ,I+1))
                        A1=((TKSTGP-TKSTGM)*Y2
     1                    -(TKSTGC-TKSTGM)*Y3)/(Y2*Y3*(Y3-Y2))
                        A2=((TKSTGP-TKSTGM)*Y2*Y2
     1                    -(TKSTGC-TKSTGM)*Y3*Y3)/(Y2*Y3*(Y2-Y3))
                        YFLN=YSTG(JJ-1)-0.5*A2/A1
                        END IF
            XFL(I)=XFLN
            YFL(I)=YFLN
            IF(YFLN.LE.0.0.AND.YFLP.LE.0.0) XFL(I)=9.9D+20
            IF(YFLN.GT.0.0.AND.YFLP.LE.0.0) THEN
                        XFL(I-1)=XFLP
                        XFL(I)=XFLN
                        END IF
            YFLP=YFLN
            XFLP=XFLN
 1326       CONTINUE
            END IF
      DO 590 I=1,IPIXL
      DO 590 J=1,JPIXL
      NAA(I,J)=0
  590 CONTINUE
      IFPLOT=2
C---------------------------FLOW VARIABLE PLOTS-------------------------
        NCOLOR=NCLOT
        NPEN=0
        CALL LINE(XS(1),YS(1),XS(LIP),YS(1))
        CALL LINE(XS(1),YS(LJJ),XS(LIP),YS(LJJ))
        CALL LINE(XS(1),YS(1),XS(1),YS(LJJ))
        CALL LINE(XS(LIP),YS(1),XS(LIP),YS(LJJ))
      IF(IFPLOT.EQ.1) THEN
C-------------------------   VELOCITY VECTORS  -------------------------
      FMAX=(1.0-ARP)*FVMAX*AMAG
      FAAY=ARMAX*(1.0-ARP)*(XS(LIP)-XS(1))/DFLOAT(LIP-1)*DFLOAT(ISVEC)
      DO 216 I=1,LIP
      IF(I.GE.2.AND.MOD(I-1,ISVEC).NE.0) GO TO 217
      IA=I+IOFF
      DO 218 J=1,LJP
      IF(J.GE.2.AND.MOD(J-1,JSVEC).NE.0) GO TO 218
      JA=J+JOFF
      UCOMP=0.5*(U(JA,IA+1)+U(JA+1,IA+1))-VOFF
      VCOMP=0.5*(V(JA+1,IA)+V(JA+1,IA+1))
      FLOC=DSQRT((UCOMP**2+VCOMP**2))
      NCOLOR=255
      NPEN=0
      IF(FLOC.GT.FMAX) THEN
                       UCOMP=UCOMP*FMAX/FLOC
                       VCOMP=VCOMP*FMAX/FLOC
                       AAB=(FLOC/FMAX-1.0)*AMAG/(1.0-AMAG)
                       NCOLOR=1+INT(AAB*253.0)
                       END IF
      FAA1=UCOMP*FAAY/FMAX
      FAA2=VCOMP*FAAY/FMAX
      IF(FAA1.GT.FAA2.AND.FAA1.GT.FAAY) THEN
                     FAA1=FAAY
                     FAA2=VCOMP*FAAY/FMAX
                     END IF
      IF(FAA2.GE.FAA1.AND.FAA2.GT.FAAY) THEN
                     FAA1=UCOMP*FAAY/FMAX
                     FAA2=FAAY
                     END IF
      FAL=ARL*DSQRT(FAA1*FAA1+FAA2*FAA2)
      FAW=ARW*FAL/ARL
      XX1=XS(I)-ARP*FAA1
      YY1=YS(LJJ-LJP+J)-ARP*FAA2
      XX2=XS(I)+(1.0-ARP)*FAA1
      YY2=YS(LJJ-LJP+J)+(1.0-ARP)*FAA2
      CALL ARROW(XX1,YY1,XX2,YY2,FAL,FAW) 
      IF(KSYM.EQ.1) THEN 
                    YY1=-YY1
                    YY2=-YY2
                    CALL ARROW(XX1,YY1,XX2,YY2,FAL,FAW) 
                    END IF
  218 CONTINUE
  217 CONTINUE
  216 CONTINUE
      END IF
      IF(IFPLOT.EQ.2) THEN
C-------------------  RASTER & CONOUR PLOTS ----------------------------
      RASLMT=(RASMX-RASMN)
      DO 600 I=1,LIP
      DO 600 J=1,LJP
      IA=I+IOFF
      JA=J+JOFF
      FS(I,J)=0.25*(TK(JA,IA)+TK(JA+1,IA)+TK(JA,IA+1)+TK(JA+1,IA+1))
C     FS2(I,J)=0.25*(RHO(JA,IA)*STMF(JA,IA)+RHO(JA+1,IA)*STMF(JA+1,IA)
C    1 +RHO(JA,IA+1)*STMF(JA,IA+1)+RHO(JA+1,IA+1)*STMF(JA+1,IA+1))
      FS2(I,J)=0.25*(FSP(JA,IA,05)+FSP(JA+1,IA,05)
     1        +FSP(JA,IA+1,05)+FSP(JA+1,IA+1,05))
  600 CONTINUE
      DO 601 I=LIP,2,-1
      DO 601 J=LJP,2,-1
      FS(I,J)=FS(I-1,J-1)
      FS2(I,J)=FS2(I-1,J-1)
  601 CONTINUE
      IF(KSYM.EQ.1) THEN
                    DO 602 I=1,LIP
                    DO 602 J=LJ2,LJP,-1
                    JA=J-LJP+1
                    FS(I,J)=FS(I,JA)
  602               CONTINUE
C------------------------ADD THESE LINES FOR FUEL CONCENTRATION PLOT----
                    DO 604 I=1,LIP
                    DO 604 J=1,LJP-1
                    JA=LJP-J+1
                    FSNEW=1.0*(FS2(I,JA)-RAS2MN)/RASCON
                    IF(FSNEW.GT.RASLMT) FSNEW=RASLMT
                    FS(I,J)=RASMN+FSNEW
  604               CONTINUE
C-------------------added
C----------use these lines for making temperature on both the sides-----
C                   DO 606 I=1,LIP
C                   DO 606 J=1,LJP-1
C                   JA=LJ2-J+1
C                   FS(I,J)=FS(I,JA)
C 606               CONTINUE
C-------------------added
                    END IF
C-------------------  RASTER PLOTS ----------------------------
      DO 610 I=1,LIP
      DO 610 J=1,LJJ
      IF(FS(I,J).GT.RASMX) FS(I,J)=RASMX
      IF(FS(I,J).LT.RASMN) FS(I,J)=RASMN
  610 CONTINUE
      IBACK=999999
      IBACKM=IBACK-10
      BACK=DFLOAT(IBACKM)/(RASMX-RASMN)
      DO 620 I=1,LIP
      DO 620 J=1,LJJ
      FS(I,J)=0.0+BACK*(FS(I,J)-RASMN)
  620 CONTINUE
C-----------------   INITIALIZE 'NAA' WITH A BIG NUMBER   -----
      DO 630 I=1,IPIXL
      DO 630 J=1,JPIXL
      NAA(I,J)=IBACK
  630 CONTINUE
      DXP=(XMAX-XMIN)/DFLOAT(IPIXL-1)
      DYP=(YMAX-YMIN)/DFLOAT(JPIXL-1)
C-----------------   IN  I&J-DIRECTIONS   ------------------------
      IA=1
      DO 632 I=1,IPIXL
      XREF=XMIN+DFLOAT(I-1)*DXP
      DO 633 II=IA,LIP
      IF(XS(II).GE.XREF) GO TO 634
  633 CONTINUE
      GO TO 632
  634 CONTINUE
      IF(II.GE.2) THEN
                  IA=II-1
                  XM=(XREF-XS(IA))/(XS(II)-XS(IA))
                  XP=(XS(II)-XREF)/(XS(II)-XS(IA))
                  ELSE
                  GO TO 632
                  END IF
      JA=1
      DO 642 J=1,JPIXL
      YREF=YMIN+DFLOAT(J-1)*DYP
      DO 643 JJ=JA,LJJ
      IF(YS(JJ).GE.YREF) GO TO 644
  643 CONTINUE
      GO TO 642
  644 CONTINUE
      IF(JJ.GE.2) THEN
                  JA=JJ-1
                  YM=(YREF-YS(JA))/(YS(JJ)-YS(JA))
                  YP=(YS(JJ)-YREF)/(YS(JJ)-YS(JA))
                  ELSE
                  GO TO 642
                  END IF
      BETA=YM*(XM*FS(II,JJ)+XP*FS(IA,JJ))
     1    +YP*(XM*FS(II,JA)+XP*FS(IA,JA))
      NAA(I,J)=INT(BETA)
  642 CONTINUE
  632 CONTINUE
C---------------------- FIND MAXIMA & MINIMA  -------------------------
      IBACKW=IBACK-2
      IBMAX=0
      IBMIN=IBACK+10
      DO 650 I=1,IPIXL
      DO 650 J=1,JPIXL
      NIJ=NAA(I,J)
      IF(NIJ.GE.IBACKW) GO TO 651
      IF(NIJ.GT.IBMAX) IBMAX=NIJ
      IF(NIJ.LT.IBMIN) IBMIN=NIJ
  651 CONTINUE
  650 CONTINUE
      IF(IBMAX.LT.IBACKM) IBMAX=IBACKM
      IF(IBMIN.GT.0) IBMIN=0
      IPLOW=1
      IPHIG=253
      FACT=DFLOAT(IPHIG-IPLOW)/DFLOAT(IBMAX-IBMIN)
      DO 652 I=1,IPIXL
      DO 652 J=1,JPIXL
      IF(NAA(I,J).GE.IBACKW) THEN
                            NAA(I,J)=NBACK
                            GO TO 653
                            END IF
      NAA(I,J)=IPLOW+INT(FACT*DFLOAT(NAA(I,J)-IBMIN))
  653 CONTINUE
  652 CONTINUE
      END IF
C---------------------PLOT ON THE SECOND HALF---------------------------
C     IF((KSYM.EQ.1).AND.(IPANM.GE.1)) THEN
C                   DO 300 I=1,IPIXL
C                   DO 300 J=JPIXL/2+1,JPIXL
C                   NAA(I,J)=0
C 300               CONTINUE
C                   END IF
C-------------------------------BODY PLOTS------------------------------
      XSHIFT=0.5*(X(1,1)+X(1,2))
      YSHIFT=0.5*(Y(1,1)+Y(2,1))
      IF(NBODY.GE.1) THEN
          DO 1320 N=1,NBODY
          IBODYM=IBDM(N)
          IBODYP=IBDP(N)
          JBODYM=JBDM(N)
          JBODYP=JBDP(N)
          JJ=JBODYP+1
          IF(JJ.GT.LJ) JJ=LJ
          II=IBODYP+1
          IF(II.GT.LI) II=LI
          XBDM(N)=X(JBODYM,IBODYM)
          YBDM(N)=Y(JBODYM,IBODYM)
          XBDP(N)=X(JBODYP,IBODYP)
          YBDP(N)=Y(JBODYP,IBODYP)
 1320     CONTINUE
          CALL DRAWBY(KSYM)
          END IF
C---------------------------PARTICLE PLOTS------------------------------
        CALL LINE(XS(1),YS(1),XS(LIP),YS(1))
        CALL LINE(XS(1),YS(LJJ),XS(LIP),YS(LJJ))
        CALL LINE(XS(1),YS(1),XS(1),YS(LJJ))
        CALL LINE(XS(LIP),YS(1),XS(LIP),YS(LJJ))
      IF(IPANM.GE.1) THEN
        DO 322 I=1,NOINJ
        DO 324 N=1,NOPR
        IF(PATHX(I,N).LT.-999.0) GO TO 324
        X1=PATHX(I,N)
        Y1=PATHY(I,N)
        NCOLOR=NCLR(N)
        IF((X1.GT.X(1,IOFF+LIP)).OR.(X1.LT.X(1,IOFF+1))) GO TO 324
        IF((Y1.GT.Y(JOFF+LJP,1)).OR.(Y1.LT.Y(JOFF+1,1))) GO TO 324
        IDD=IDOT
        CALL DOTPR(X1,Y1,IDD)
  324   CONTINUE
  322   CONTINUE
        END IF
C---------------------------FLAME PLOT----------------------------------
      IF(NFSURF.GE.1) THEN
          NCOLOR=255
          IDOTP=0
          JDOTP=0
          IDSYM=KSYM
          DO 670 I=1,LIP
          X1=XFL(I+IOFF)
          IF(X1.GE.1.0D+20) GO TO 670
          Y1=YFL(I+IOFF)
          IF(Y1.LE.YS(LJJ)) THEN
                 CALL DOT(X1,Y1,NFSURF,IDOTP,JDOTP,IDSYM,KSYM)
                 IF(IDSYM.EQ.1) THEN
C                 CALL DOT(X1,-Y1,NFSURF,IDOTP,JDOTP,IDSYM,KSYM)
                   IDSYM=0
                   END IF
                 END IF
  670     CONTINUE
          END IF
           DO 400 I=1,IPIXL
           DO 400 J=1,JPIXL
           IS=NAA(I,J)
           IF(IS.GT.255) IS=255
           IF(IS.LT.0) IS=0
           BB(I,J)=CHAR(IS)
  400      CONTINUE
           IF(IPAGE.GE.3) IPAGE=IPAGE+1
           IF(IPAGE.EQ.2) THEN
                          BFL2=BFL2A(1)//BFL2B(1)
                          ELSE
                          IPP=IPAGE-3
                          IPPA=(IPP-1)/100
                          IPPB=IPP-IPPA*100
                          BFL2=BFL2A(IPPA+1)//BFL2B(IPPB)
                          END IF
C          IF(RASMX.LE.1000.0/TSTR) THEN
C                        IPAGE=IPAGE-2
C                        IF(IPAGE.NE.0) IPAGE=1
C                        END IF
           BBFILE=BFL1//BFL2//BFL3
           MPP=82
           OPEN(MPP,FILE=BBFILE,STATUS='UNKNOWN')
           ITMP=IPIXL*JPIXL
           IF(KORNT.EQ.1) THEN
                CALL header(MPP,IPIXL,JPIXL)
                WRITE(FRMTSTR,'(''('',i8.8,''A,$)'')') ITMP
                WRITE(MPP,FMT=FRMTSTR)
     1            ((BB(I,J),I=1,IPIXL),J=JPIXL,1,-1) 
                END IF
           IF(KORNT.EQ.2) THEN
                CALL header(MPP,JPIXL,IPIXL)
                WRITE(FRMTSTR,'(''('',i8.8,''A,$)'')') ITMP
                IF(KSYM.EQ.1) THEN
                   WRITE(MPP,FMT=FRMTSTR)
     1                  ((BB(I,J),J=JPIXL,1,-1),I=1,IPIXL) 
                   ELSE
                   WRITE(MPP,FMT=FRMTSTR)
     1                  ((BB(I,J),J=1,JPIXL),I=1,IPIXL) 
                   END IF
                END IF
           IF(KORNT.EQ.3) THEN
                CALL header(MPP,IPIXL,JPIXL)
                WRITE(FRMTSTR,'(''('',i8.8,''A,$)'')') ITMP
                WRITE(MPP,FMT=FRMTSTR)
     1            ((BB(I,J),I=IPIXL,1,-1),J=1,JPIXL) 
                END IF
           IF(KORNT.EQ.4) THEN
                CALL header(MPP,JPIXL,IPIXL)
                WRITE(FRMTSTR,'(''('',i8.8,''A,$)'')') ITMP
                WRITE(MPP,FMT=FRMTSTR)
     1            ((BB(I,J),J=1,JPIXL),I=1,IPIXL) 
                END IF
           CLOSE(MPP)
      RETURN
      END
C***********************************************************************
       subroutine header(mpp,ipix,jpix)
C***********************************************************************
       character*1 cpl(4,256)
       character*4 byt4
       character*1 byta,bytb,bytc,bytd
       character*14 frmtstr
       character*54 headmsw
C-----------------rainbow color palette---------------------------------
C-----red
       do 110 i=1,142
       irgb=float(128*(29-i))/28.0
       if(irgb.GT.255) irgb=255
       if(irgb.LT.0) irgb=0
       cpl(3,i)=CHAR(irgb)
 110   continue
       do 111 i=143,256
       irgb=float(255*(i-142))/57.0
       if(irgb.GT.255) irgb=255
       if(irgb.LT.0) irgb=0
       cpl(3,i)=CHAR(irgb)
 111   continue
C-----green
       do 114 i=1,199
       irgb=float(255*(i-29))/56.0
       if(irgb.GT.255) irgb=255
       if(irgb.LT.0) irgb=0
       cpl(2,i)=CHAR(irgb)
 114   continue
       do 115 i=200,256
       irgb=256.0-float(255*(i-199))/57.0
       if(irgb.GT.255) irgb=255
       if(irgb.LT.0) irgb=0
       cpl(2,i)=CHAR(irgb)
 115   continue
C------blue
       do 118 i=1,142
       irgb=256.0-float(255*(i-85))/57.0
       if(irgb.GT.255) irgb=255
       if(irgb.LT.0) irgb=0
       cpl(1,i)=CHAR(irgb)
 118   continue
       do 119 i=143,256
       irgb=0
       cpl(1,i)=CHAR(irgb)
 119   continue
C------last black
       irgb=255
       cpl(1,1)=CHAR(irgb)
       cpl(2,1)=CHAR(irgb)
       cpl(3,1)=CHAR(irgb)
C------first white
       irgb=0
       cpl(1,256)=CHAR(irgb)
       cpl(2,256)=CHAR(irgb)
       cpl(3,256)=CHAR(irgb)
C------extra
       do 120 i=1,256
       irgb=0
       cpl(4,i)=CHAR(irgb)
 120   continue
C--------------------------------end palette---------------
* header 1 (file header ; 1--14 byte)
         headmsw( 1: 2) = 'BM'       ! declaring this is BMP file
         itmp=1024+54+ipix*jpix      ! total file size = header + data
         call num2bit4(itmp,byt4)
         headmsw( 3: 6) = byt4(1:4)
         itmp = 0                    ! may be 0
         call num2bit4(itmp,byt4)
         headmsw( 7: 8) = byt4(1:2)
         headmsw( 9:10) = byt4(1:2)
         itmp = 1078              ! must be 54 : total length of header
         call num2bit4(itmp,byt4)
         headmsw(11:14) = byt4(1:4)
* header 2 (bit-map header ; 13--54 byte)
         itmp = 40                   ! 40 : length of bit-map header
         call num2bit4(itmp,byt4)
         headmsw(15:18) = byt4(1:4)
         itmp = ipix                 ! width
         call num2bit4(itmp,byt4)
         headmsw(19:22) = byt4(1:4)
         itmp = jpix                 ! height
         call num2bit4(itmp,byt4)
         headmsw(23:26) = byt4(1:4)
         itmp = 1                    ! must be 1
         call num2bit4(itmp,byt4)
         headmsw(27:28) = byt4(1:2)
         itmp =  8                   ! 8 : color depth in bit.
         call num2bit4(itmp,byt4)
         headmsw(29:30) = byt4(1:2)
         itmp = 0                    ! 0 : no compression
         call num2bit4(itmp,byt4)
         headmsw(31:34) = byt4(1:4)
         headmsw(35:38) = byt4(1:4)
         headmsw(39:42) = byt4(1:4)
         headmsw(43:46) = byt4(1:4)
         itmp = 256                  ! 256 : num. of colors used
         call num2bit4(itmp,byt4)
         headmsw(47:50) = byt4(1:4)
         itmp = 0                    ! 0 : num. of important color
         call num2bit4(itmp,byt4)
         headmsw(51:54) = byt4(1:4)
* writing header part
         write(mpp,'(a54,$)') headmsw(1:54)
* writing color palette data
         itmp = 256 * 4
         write(frmtstr,'(''('',i8.8,''A,$)'')') itmp
         write(mpp,fmt=frmtstr)
     &      ((cpl(k,i),k=1,4),i=1,256) 
	     return
		 end
C***********************************************************************
       subroutine num2bit4(inum,byt4)
C***********************************************************************
C --------------------------------------
C convert number to 8-bit characters
C --------------------------------------
       character*4 byt4
       character*1 byta,bytb,bytc,bytd
       itmp1 = inum
       itmp2 = itmp1 / 256**3
       bytd = char(itmp2)
       itmp1 =-itmp2 * 256**3 +itmp1
       itmp2 = itmp1 / 256**2
       bytc = char(itmp2)
       itmp1 =-itmp2 * 256**2 +itmp1
       itmp2 = itmp1 / 256
       bytb = char(itmp2)
       itmp1 =-itmp2 * 256    +itmp1
       byta = char(itmp1)
	   byt4=byta//bytb//bytc//bytd
       return
       end
C***********************************************************************
      SUBROUTINE DRAWBY(KSYM)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (LI=361,LJ=131,IPMAX=500,JPMAX=200)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/APIXL/ IPAGE,NCOLOR,NBACK,NAA(IPMAX,JPMAX)
      COMMON/FRAME/ IPIXL,JPIXL,XMIN,YMIN,XMAX,YMAX,XPIXL,YPIXL
      NCOLOR=255
      DO 200 N=1,NBODY
      IBODYM=IBDM(N)
      IBODYP=IBDP(N)
      JBODYM=JBDM(N)
      JBODYP=JBDP(N)
      DY=(YMAX-YMIN)/FLOAT(JPIXL)
      XA=X(JBODYM,IBODYM)
      IF(XA.LT.XMIN) XA=XMIN
      XB=X(JBODYP,IBODYP)
      IF(XB.GT.XMAX) XB=XMAX
      IF(XB.LT.XMIN) XB=XMIN
      YA=Y(JBODYM,IBODYM)
      IF(YA.LT.YMIN) YA=YMIN
  100 CONTINUE
      CALL LINE(XA,YA,XB,YA)
      IF(KSYM.EQ.1) CALL LINE(XA,-YA,XB,-YA)
      YA=YA+DY
      IF(YA.LE.Y(JBODYP,IBODYP).AND.YA.LE.YMAX) GO TO 100
  200 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE LINE(X1,Y1,X2,Y2)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IPMAX=500,JPMAX=200)
      COMMON/APIXL/ IPAGE,NCOLOR,NBACK,NAA(IPMAX,JPMAX)
      COMMON/FRAME/ IPIXL,JPIXL,XMIN,YMIN,XMAX,YMAX,XPIXL,YPIXL
      COMMON/PEN/NPEN,NAPEN,NBPEN
C
      EPS=1.0D-08
      XBIG=(X1-XMIN)*XPIXL
      XEND=(X2-XMIN)*XPIXL
      YBIG=(Y1-YMIN)*YPIXL
      YEND=(Y2-YMIN)*YPIXL
C---------------- CONVERT X-Y VALUES IN TO PIXELS  ---------------------
      IBIG=NINT(XBIG)
      IEND=NINT(XEND)
      JBIG=NINT(YBIG)
      JEND=NINT(YEND)
      IDIF=IEND-IBIG
      IF(IDIF.LT.0) IDIF=-IDIF
      JDIF=JEND-JBIG
      IF(JDIF.LT.0) JDIF=-JDIF
      IF(IDIF+JDIF.LE.0) THEN
                         KVEC=1
                         DX=1.0
                         DY=1.0
                         GO TO 11
                         END IF
      IF(IDIF.GE.JDIF) THEN
                       KVEC=IDIF
                       DX=(XEND-XBIG)/DABS(XEND-XBIG)
                       DY=(YEND-YBIG)/DABS(XEND-XBIG)
                       END IF
      IF(JDIF.GT.IDIF) THEN
                       KVEC=JDIF
                       DX=(XEND-XBIG)/DABS(YEND-YBIG)
                       DY=(YEND-YBIG)/DABS(YEND-YBIG)
                       END IF
C---------------------------   DRAW LINE  ------------------------------
   11 CONTINUE
      DO 12 K=1,KVEC
      NNCLR=NCOLOR
      IF(NPEN.GE.1) THEN
             NCPEN=NAPEN
             NAPEN=NBPEN
             NBPEN=NCPEN
             IF(NAPEN.EQ.1) NNCLR=0
             END IF
      XNEW=XBIG+DFLOAT(K-1)*DX
      YNEW=YBIG+DFLOAT(K-1)*DY
      I=NINT(XNEW)
      J=NINT(YNEW)
      NAA(I,J)=NNCLR
   12 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE ARROW(X1,Y1,X2,Y2,X3,X4)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IPMAX=500,JPMAX=200)
      COMMON/APIXL/ IPAGE,NCOLOR,NBACK,NAA(IPMAX,JPMAX)
      COMMON/FRAME/ IPIXL,JPIXL,XMIN,YMIN,XMAX,YMAX,XPIXL,YPIXL
      COMMON/PEN/NPEN,NAPEN,NBPEN
C
      EPS=1.0D-08
      AW1=DABS(X2-X1)
      IF((AW1+ABS(Y2-Y1)).LE.EPS) RETURN
      IF(AW1.LT.EPS) AW1=DABS(Y2-Y1)
      IF(X3.LT.EPS) X3=0.2*AW1
      IF(X4.LT.EPS) X4=0.2*AW1
      XBIG=(X1-XMIN)*XPIXL
      XEND=(X2-XMIN)*XPIXL
      YBIG=(Y1-YMIN)*YPIXL
      YEND=(Y2-YMIN)*YPIXL
      HEADL=X3*XPIXL
      HEADW=X4*XPIXL
C----------------------  FIX ARROW DIMENSIONS  -------------------------
      PI=3.1415926
      EPS=1.0D-08
      ALPHA=ATAN(0.5*HEADW/HEADL)
      SIDE=SQRT(HEADL*HEADL+HEADW*HEADW/4.0)
      XDIF=XEND-XBIG
      YDIF=YEND-YBIG
      TETA=DATAN(ABS(YDIF/(XDIF+EPS)))
      IF(DABS(XDIF).LT.EPS) TETA=PI/2.0
      S1=1.0
      S2=1.0
      IF(XDIF.LT.0.0) S1=-1.0
      IF(YDIF.LT.0.0) S2=-1.0
      S3=S1*S2
      BETA=TETA-S3*ALPHA
      XLFT=XEND-S1*SIDE*DCOS(BETA)
      YLFT=YEND-S2*SIDE*DSIN(BETA)
      BETA=TETA+S3*ALPHA
      XRGT=XEND-S1*SIDE*DCOS(BETA)
      YRGT=YEND-S2*SIDE*DSIN(BETA)
      TETAD=TETA*180.0/PI
      ALPHAD=ALPHA*180.0/PI
      BETAD=BETA*180.0/PI
C---------------- CONVERT X-Y VALUES IN TO PIXELS  ---------------------
      IBIG=NINT(XBIG)
      IEND=NINT(XEND)
      JBIG=NINT(YBIG)
      JEND=NINT(YEND)
      IDIF=IEND-IBIG
      IF(IDIF.LT.0) IDIF=-IDIF
      JDIF=JEND-JBIG
      IF(JDIF.LT.0) JDIF=-JDIF
      IF(IDIF+JDIF.LE.0) THEN
                         KVEC=1
                         DX=1.0
                         DY=1.0
                         GO TO 11
                         END IF
      IF(IDIF.GE.JDIF) THEN
                       KVEC=IDIF
                       DX=(XEND-XBIG)/DABS(XEND-XBIG)
                       DY=(YEND-YBIG)/DABS(XEND-XBIG)
                       END IF
      IF(JDIF.GT.IDIF) THEN
                       KVEC=JDIF
                       DX=(XEND-XBIG)/DABS(YEND-YBIG)
                       DY=(YEND-YBIG)/DABS(YEND-YBIG)
                       END IF
C---------------------------   DRAW LINE  ------------------------------
   11 CONTINUE
      DO 12 K=1,KVEC
      XNEW=XBIG+DFLOAT(K-1)*DX
      YNEW=YBIG+DFLOAT(K-1)*DY
      I=NINT(XNEW)
      J=NINT(YNEW)
      NAA(I,J)=NCOLOR
   12 CONTINUE
C-----------------------   DRAW ARROW HEAD (LEFT)-----------------------
      IBIG=NINT(XLFT)
      JBIG=NINT(YLFT)
      IDIF=IEND-IBIG
      IF(IDIF.LT.0) IDIF=-IDIF
      JDIF=JEND-JBIG
      IF(JDIF.LT.0) JDIF=-JDIF
      IF(IDIF+JDIF.LE.0) THEN
                         KVEC=1
                         DX=1.0
                         DY=1.0
                         GO TO 21
                         END IF
      IF(IDIF.GE.JDIF) THEN
                       KVEC=IDIF
                       DX=(XEND-XLFT)/DABS(XEND-XLFT)
                       DY=(YEND-YLFT)/DABS(XEND-XLFT)
                       END IF
      IF(JDIF.GT.IDIF) THEN
                       KVEC=JDIF
                       DX=(XEND-XLFT)/DABS(YEND-YLFT)
                       DY=(YEND-YLFT)/DABS(YEND-YLFT)
                       END IF
C--------------------------   DRAW LINE  -------------------------------
   21 CONTINUE
      DO 22 K=1,KVEC
      XNEW=XLFT+FLOAT(K-1)*DX
      YNEW=YLFT+FLOAT(K-1)*DY
      I=NINT(XNEW)
      J=NINT(YNEW)
      NAA(I,J)=NCOLOR
   22 CONTINUE
C-----------------------   DRAW ARROW HEAD (RIGHT)----------------------
      IBIG=NINT(XRGT)
      JBIG=NINT(YRGT)
      IDIF=IEND-IBIG
      IF(IDIF.LT.0) IDIF=-IDIF
      JDIF=JEND-JBIG
      IF(JDIF.LT.0) JDIF=-JDIF
      IF(IDIF+JDIF.LE.0) THEN
                         KVEC=1
                         DX=1.0
                         DY=1.0
                         GO TO 31
                         END IF
      IF(IDIF.GE.JDIF) THEN
                       KVEC=IDIF
                       DX=(XEND-XRGT)/DABS(XEND-XRGT)
                       DY=(YEND-YRGT)/DABS(XEND-XRGT)
                       END IF
      IF(JDIF.GT.IDIF) THEN
                       KVEC=JDIF
                       DX=(XEND-XRGT)/DABS(YEND-YRGT)
                       DY=(YEND-YRGT)/DABS(YEND-YRGT)
                       END IF
C---------------------------   DRAW LINE  ------------------------------
   31 CONTINUE
      DO 32 K=1,KVEC
      XNEW=XRGT+DFLOAT(K-1)*DX
      YNEW=YRGT+DFLOAT(K-1)*DY
      I=NINT(XNEW)
      J=NINT(YNEW)
      NAA(I,J)=NCOLOR
   32 CONTINUE
C-----------------------   DRAW ARROW HEAD (BACK) ----------------------
      IBIG=NINT(XRGT)
      IEND=NINT(XLFT)
      JBIG=NINT(YRGT)
      JEND=NINT(YLFT)
      IDIF=IEND-IBIG
      IF(IDIF.LT.0) IDIF=-IDIF
      JDIF=JEND-JBIG
      IF(JDIF.LT.0) JDIF=-JDIF
      IF(IDIF+JDIF.LE.0) THEN
                         KVEC=1
                         DX=1.0
                         DY=1.0
                         GO TO 41
                         END IF
      IF(IDIF.GE.JDIF) THEN
                       KVEC=IDIF
                       DX=(XLFT-XRGT)/DABS(XLFT-XRGT)
                       DY=(YLFT-YRGT)/DABS(XLFT-XRGT)
                       END IF
      IF(JDIF.GT.IDIF) THEN
                       KVEC=JDIF
                       DX=(XLFT-XRGT)/DABS(YLFT-YRGT)
                       DY=(YLFT-YRGT)/DABS(YLFT-YRGT)
                       END IF
C---------------------------   DRAW LINE  ------------------------------
   41 CONTINUE
      DO 42 K=1,KVEC
      XNEW=XRGT+DFLOAT(K-1)*DX
      YNEW=YRGT+DFLOAT(K-1)*DY
      I=NINT(XNEW)
      J=NINT(YNEW)
      NAA(I,J)=NCOLOR
   42 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE DOT(X1,Y1,IDOT,IDOTP,JDOTP,IDSYM,KSYM)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IPMAX=500,JPMAX=200)
      COMMON/APIXL/ IPAGE,NCOLOR,NBACK,NAA(IPMAX,JPMAX)
      COMMON/FRAME/ IPIXL,JPIXL,XMIN,YMIN,XMAX,YMAX,XPIXL,YPIXL
      COMMON/PEN/NPEN,NAPEN,NBPEN
C
      XBIG=(X1-XMIN)*XPIXL
      YBIG=(Y1-YMIN)*YPIXL
C---------------------------   DRAW DOT  ------------------------------
      II=NINT(XBIG)
      JJ=NINT(YBIG)
      IC=II-IDOTP
      JC=JJ-JDOTP
      IF(IDOTP.EQ.0) GO TO 10
      IF(IDSYM.EQ.1) GO TO 11
      IDSYM=0
      IF(IC.LT.0) IC=-IC
      IF(JC.LT.0) JC=-JC
      IDOT1=IDOT*2
      IF(IC.LT.IDOT1.AND.JC.LT.IDOT1) RETURN 
      IDSYM=KSYM
   10 CONTINUE
      IDOTP=II
      JDOTP=JJ
   11 CONTINUE
      DO 12 KVV=1,IDOT
      DO 12 KHH=1,IDOT
      KV=KVV-1
      KH=KHH-1
      IR=NINT(DSQRT(DFLOAT(KV*KV+KH*KH)))+1
      IF(IR.GT.IDOT) GO TO 14
      I=II+KV
      J=JJ+KH
      NAA(I,J)=NCOLOR
      I=II+KV
      J=JJ-KH
      NAA(I,J)=NCOLOR
      I=II-KV
      J=JJ+KH
      NAA(I,J)=NCOLOR
      I=II-KV
      J=JJ-KH
      NAA(I,J)=NCOLOR
   14 CONTINUE
   12 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE DOTPR(X1,Y1,IDOT)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IPMAX=500,JPMAX=200)
      COMMON/APIXL/ IPAGE,NCOLOR,NBACK,NAA(IPMAX,JPMAX)
      COMMON/FRAME/ IPIXL,JPIXL,XMIN,YMIN,XMAX,YMAX,XPIXL,YPIXL
C
      XBIG=(X1-XMIN)*XPIXL
      YBIG=(Y1-YMIN)*YPIXL
C---------------------------   DRAW DOT  ------------------------------
   11 CONTINUE
      II=NINT(XBIG)
      JJ=NINT(YBIG)
      IDOTO=IDOT
      IDOTI=0
      IF(IDOT.LT.0) THEN
            IDOTO=-IDOT
            IDOTI=IDOTO-1
            END IF
      DO 12 KVV=1,IDOTO
      DO 12 KHH=1,IDOTO
      KV=KVV-1
      KH=KHH-1
      IR=NINT(SQRT(FLOAT(KV*KV+KH*KH)))+1
      IF(IR.GT.IDOTO) GO TO 14
      IF(IR.LE.IDOTI) GO TO 14
      I=II+KV
      J=JJ+KH
      NAA(I,J)=NCOLOR
      I=II+KV
      J=JJ-KH
      NAA(I,J)=NCOLOR
      I=II-KV
      J=JJ+KH
      NAA(I,J)=NCOLOR
      I=II-KV
      J=JJ-KH
      NAA(I,J)=NCOLOR
   14 CONTINUE
   12 CONTINUE
      RETURN
      END
C***********************************************************************
      BLOCK DATA
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (LRX=544)
      COMMON/REAC/ LREV,NRLIND,ITBEF(LRX),ALF(LRX),AKF(LRX),EXAF(LRX),
     1             ALLOW(LRX),AKLOW(LRX),EALOW(LRX),TROE(LRX,4)
      COMMON/CB08/ TPOL1,TPOL2,P001(14),P002(14),P003(14),P004(14),
     1  P005(14),P006(14),P007(14),P008(14),P009(14),P010(14),P011(14),
     2  P012(14),P013(14),P014(14),P015(14),P016(14),P017(14),P018(14),
     4  P019(14),P020(14),P021(14),P022(14),P023(14),P024(14),P025(14),
     5  P026(14),P027(14),P028(14),P029(14),P030(14),P031(14),P032(14),
     6  P033(14),P034(14),P035(14),P036(14),P037(14),P038(14),P039(14),
     7  P040(14),P041(14),P042(14),P043(14),P044(14),P045(14),P046(14),
     8  P047(14),P048(14),P049(14),P050(14),P051(14),P052(14),
     *  C001(12),C002(12),C003(12),C004(12),
     *  C005(12),C006(12),C007(12),C008(12),C009(12),C010(12),C011(12),
     *  C012(12),C013(12),C014(12),C015(12),C016(12),C017(12),C018(12),
     *  C019(12),C020(12),C021(12),C022(12),C023(12),C024(12),C025(12),
     *  C026(12),C027(12),C028(12),C029(12),C030(12),C031(12),C032(12),
     *  C033(12),C034(12),C035(12),C036(12),C037(12),C038(12),C039(12),
     *  C040(12),C041(12),C042(12),C043(12),C044(12),C045(12),C046(12),
     *  C047(12),C048(12),C049(12),C050(12),C051(12),C052(12)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    Curve Fitts for the Thermochemical Properties                   CC
CC      Cp/R = A1 + A2*T + A3*T^2 + A4*T^3 +A5*T^4                    CC
CC      H/RT = A1 + A2*T/2 + A3*T^2/3 + A4*T^3/4 +A5*T^4/5 + A6/T     CC
CC       S/R = A1*Ln(T) + A2*T + A3*T^2/2 + A4*T^3/3 +A5*T^4/4 + A7   CC
CC   R=8.3144 J/g-mole/K ==>> Cp - J/g-mole/K ; H - J/g-mole          CC
CC      H - is the total enthalpy i.e. it includes heat of formation  CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C----------------   TEMPERATURE LIMITS FOR THE TWO POLYNOMIALS   -------
C
      DATA TPOL1,TPOL2
     1    /  1.0D+03,  5.0D+03 /
C
C-------------------------  H2               -------------------------
C
      DATA P001
     1 / 2.34433112D+00, 7.98052075D-03,-1.94781510D-05, 2.01572094D-08,
     2  -7.37611761D-12,-9.17935173D+02, 6.83010238D-01, 3.33727920D+00,
     3  -4.94024731D-05, 4.99456778D-07,-1.79566394D-10, 2.00255376D-14,
     4  -9.50158922D+02,-3.20502331D+00/
C
C-------------------------  O2               -------------------------
C
      DATA P002
     1 / 3.78245636D+00,-2.99673416D-03, 9.84730201D-06,-9.68129509D-09,
     2   3.24372837D-12,-1.06394356D+03, 3.65767573D+00, 3.28253784D+00,
     3   1.48308754D-03,-7.57966669D-07, 2.09470555D-10,-2.16717794D-14,
     4  -1.08845772D+03, 5.45323129D+00/
C
C-------------------------  H                -------------------------
C
      DATA P003
     1 / 2.50000000D+00, 7.05332819D-13,-1.99591964D-15, 2.30081632D-18,
     2  -9.27732332D-22, 2.54736599D+04,-4.46682853D-01, 2.50000001D+00,
     3  -2.30842973D-11, 1.61561948D-14,-4.73515235D-18, 4.98197357D-22,
     4   2.54736599D+04,-4.46682914D-01/
C
C-------------------------  O                -------------------------
C
      DATA P004
     1 / 3.16826710D+00,-3.27931884D-03, 6.64306396D-06,-6.12806624D-09,
     2   2.11265971D-12, 2.91222592D+04, 2.05193346D+00, 2.56942078D+00,
     3  -8.59741137D-05, 4.19484589D-08,-1.00177799D-11, 1.22833691D-15,
     4   2.92175791D+04, 4.78433864D+00/
C
C-------------------------  OH               -------------------------
C
      DATA P005
     1 / 4.12530561D+00,-3.22544939D-03, 6.52764691D-06,-5.79853643D-09,
     2   2.06237379D-12, 3.38153812D+03,-6.90432960D-01, 2.86472886D+00,
     3   1.05650448D-03,-2.59082758D-07, 3.05218674D-11,-1.33195876D-15,
     4   3.71885774D+03, 5.70164073D+00/
C
C-------------------------  H2O              -------------------------
C
      DATA P006
     1 / 4.19864056D+00,-2.03643410D-03, 6.52040211D-06,-5.48797062D-09,
     2   1.77197817D-12,-3.02937267D+04,-8.49032208D-01, 3.03399249D+00,
     3   2.17691804D-03,-1.64072518D-07,-9.70419870D-11, 1.68200992D-14,
     4  -3.00042971D+04, 4.96677010D+00/
C
C-------------------------  HO2              -------------------------
C
      DATA P007
     1 / 4.30179801D+00,-4.74912051D-03, 2.11582891D-05,-2.42763894D-08,
     2   9.29225124D-12, 2.94808040D+02, 3.71666245D+00, 4.01721090D+00,
     3   2.23982013D-03,-6.33658150D-07, 1.14246370D-10,-1.07908535D-14,
     4   1.11856713D+02, 3.78510215D+00/
C
C-------------------------  H2O2             -------------------------
C
      DATA P008
     1 / 4.27611269D+00,-5.42822417D-04, 1.67335701D-05,-2.15770813D-08,
     2   8.62454363D-12,-1.77025821D+04, 3.43505074D+00, 4.16500285D+00,
     3   4.90831694D-03,-1.90139225D-06, 3.71185986D-10,-2.87908305D-14,
     4  -1.78617877D+04, 2.91615662D+00/
C
C-------------------------  CO               -------------------------
C
      DATA P009
     1 / 3.57953347D+00,-6.10353680D-04, 1.01681433D-06, 9.07005884D-10,
     2  -9.04424499D-13,-1.43440860D+04, 3.50840928D+00, 2.71518561D+00,
     3   2.06252743D-03,-9.98825771D-07, 2.30053008D-10,-2.03647716D-14,
     4  -1.41518724D+04, 7.81868772D+00/
C
C-------------------------  CO2              -------------------------
C
      DATA P010
     1 / 2.35677352D+00, 8.98459677D-03,-7.12356269D-06, 2.45919022D-09,
     2  -1.43699548D-13,-4.83719697D+04, 9.90105222D+00, 3.85746029D+00,
     3   4.41437026D-03,-2.21481404D-06, 5.23490188D-10,-4.72084164D-14,
     4  -4.87591660D+04, 2.27163806D+00/
C
C-------------------------  HCO              -------------------------
C
      DATA P011
     1 / 4.22118584D+00,-3.24392532D-03, 1.37799446D-05,-1.33144093D-08,
     2   4.33768865D-12, 3.83956496D+03, 3.39437243D+00, 2.77217438D+00,
     3   4.95695526D-03,-2.48445613D-06, 5.89161778D-10,-5.33508711D-14,
     4   4.01191815D+03, 9.79834492D+00/
C
C-------------------------  CH2O             -------------------------
C
      DATA P012
     1 / 4.79372315D+00,-9.90833369D-03, 3.73220008D-05,-3.79285261D-08,
     2   1.31772652D-11,-1.43089567D+04, 6.02812900D-01, 1.76069008D+00,
     3   9.20000082D-03,-4.42258813D-06, 1.00641212D-09,-8.83855640D-14,
     4  -1.39958323D+04, 1.36563230D+01/
C
C-------------------------  CH4              -------------------------
C
      DATA P013
     1 / 5.14987613D+00,-1.36709788D-02, 4.91800599D-05,-4.84743026D-08,
     2   1.66693956D-11,-1.02466476D+04,-4.64130376D+00, 7.48514950D-02,
     3   1.33909467D-02,-5.73285809D-06, 1.22292535D-09,-1.01815230D-13,
     4  -9.46834459D+03, 1.84373180D+01/
C
C-------------------------  CH3              -------------------------
C
      DATA P014
     1 / 3.67359040D+00, 2.01095175D-03, 5.73021856D-06,-6.87117425D-09,
     2   2.54385734D-12, 1.64449988D+04, 1.60456433D+00, 2.28571772D+00,
     3   7.23990037D-03,-2.98714348D-06, 5.95684644D-10,-4.67154394D-14,
     4   1.67755843D+04, 8.48007179D+00/
C
C-------------------------  T-CH2            -------------------------
C
      DATA P015
     1 / 3.76267867D+00, 9.68872143D-04, 2.79489841D-06,-3.85091153D-09,
     2   1.68741719D-12, 4.60040401D+04, 1.56253185D+00, 2.87410113D+00,
     3   3.65639292D-03,-1.40894597D-06, 2.60179549D-10,-1.87727567D-14,
     4   4.62636040D+04, 6.17119324D+00/
C
C-------------------------  S-CH2            -------------------------
C
      DATA P016
     1 / 4.19860411D+00,-2.36661419D-03, 8.23296220D-06,-6.68815981D-09,
     2   1.94314737D-12, 5.04968163D+04,-7.69118967D-01, 2.29203842D+00,
     3   4.65588637D-03,-2.01191947D-06, 4.17906000D-10,-3.39716365D-14,
     4   5.09259997D+04, 8.62650169D+00/
C
C-------------------------  C2H4             -------------------------
C
      DATA P017
     1 / 3.95920148D+00,-7.57052247D-03, 5.70990292D-05,-6.91588753D-08,
     2   2.69884373D-11, 5.08977593D+03, 4.09733096D+00, 2.03611116D+00,
     3   1.46454151D-02,-6.71077915D-06, 1.47222923D-09,-1.25706061D-13,
     4   4.93988614D+03, 1.03053693D+01/
C
C-------------------------  CH3O             -------------------------
C
      DATA P018
     1 / 3.71180502D+00,-2.80463306D-03, 3.76550971D-05,-4.73072089D-08,
     2   1.86588420D-11, 1.30772484D+03, 6.57240864D+00, 4.75779238D+00,
     3   7.44142474D-03,-2.69705176D-06, 4.38090504D-10,-2.63537098D-14,
     4   3.90139164D+02,-1.96680028D+00/
C
C-------------------------  C2H5             -------------------------
C
      DATA P019
     1 / 4.30646568D+00,-4.18658892D-03, 4.97142807D-05,-5.99126606D-08,
     2   2.30509004D-11, 1.28416265D+04, 4.70720924D+00, 1.95465642D+00,
     3   1.73972722D-02,-7.98206668D-06, 1.75217689D-09,-1.49641576D-13,
     4   1.28575200D+04, 1.34624343D+01/
C
C-------------------------  C2H6             -------------------------
C
      DATA P020
     1 / 4.29142492D+00,-5.50154270D-03, 5.99438288D-05,-7.08466285D-08,
     2   2.68685771D-11,-1.15222055D+04, 2.66682316D+00, 1.07188150D+00,
     3   2.16852677D-02,-1.00256067D-05, 2.21412001D-09,-1.90002890D-13,
     4  -1.14263932D+04, 1.51156107D+01/
C
C-------------------------  CH               -------------------------
C
      DATA P021
     1 / 3.48981665D+00, 3.23835541D-04,-1.68899065D-06, 3.16217327D-09,
     2  -1.40609067D-12, 7.07972934D+04, 2.08401108D+00, 2.87846473D+00,
     3   9.70913681D-04, 1.44445655D-07,-1.30687849D-10, 1.76079383D-14,
     4   7.10124364D+04, 5.48497999D+00/
C
C-------------------------  C2H2             -------------------------
C
      DATA P022
     1 / 8.08681094D-01, 2.33615629D-02,-3.55171815D-05, 2.80152437D-08,
     2  -8.50072974D-12, 2.64289807D+04, 1.39397051D+01, 4.14756964D+00,
     3   5.96166664D-03,-2.37294852D-06, 4.67412171D-10,-3.61235213D-14,
     4   2.59359992D+04,-1.23028121D+00/
C
C-------------------------  C2H3             -------------------------
C
      DATA P023
     1 / 3.21246645D+00, 1.51479162D-03, 2.59209412D-05,-3.57657847D-08,
     2   1.47150873D-11, 3.48598468D+04, 8.51054025D+00, 3.01672400D+00,
     3   1.03302292D-02,-4.68082349D-06, 1.01763288D-09,-8.62607041D-14,
     4   3.46128739D+04, 7.78732378D+00/
C
C-------------------------  CH2CHO           -------------------------
C
      DATA P024
     1 / 1.01340010D+00, 2.26814670D-02,-1.57339440D-05, 4.04915030D-09,
     2   2.95990120D-13, 3.80428530D+02, 1.93565520D+01, 5.16620060D+00,
     3   1.08478260D-02,-4.46583680D-06, 8.06285480D-10,-4.84101930D-14,
     4  -7.31993470D+02,-1.96333610D+00/
C
C-------------------------  C2H4O            -------------------------
C
      DATA P025
     1 / 4.72945950D+00,-3.19328580D-03, 4.75349210D-05,-5.74586110D-08,
     2   2.19311120D-11,-2.15728780D+04, 4.10301590D+00, 5.40411080D+00,
     3   1.17230590D-02,-4.22631370D-06, 6.83724510D-10,-4.09848630D-14,
     4  -2.25931220D+04,-3.48079170D+00/
C
C-------------------------  CH2CO            -------------------------
C
      DATA P026
     1 / 2.13583630D+00, 1.81188721D-02,-1.73947474D-05, 9.34397568D-09,
     2  -2.01457615D-12,-7.04291804D+03, 1.22156480D+01, 4.51129732D+00,
     3   9.00359745D-03,-4.16939635D-06, 9.23345882D-10,-7.94838201D-14,
     4  -7.55105311D+03, 6.32247205D-01/
C
C-------------------------  HCCO             -------------------------
C
      DATA P027
     1 / 2.25172140D+00, 1.76550210D-02,-2.37291010D-05, 1.72757590D-08,
     2  -5.06648110D-12, 2.00594490D+04, 1.24904170D+01, 5.62820580D+00,
     3   4.08534010D-03,-1.59345470D-06, 2.86260520D-10,-1.94078320D-14,
     4   1.93272150D+04,-3.93025950D+00/
C
C-------------------------  C2H              -------------------------
C
      DATA P028
     1 / 2.88965733D+00, 1.34099611D-02,-2.84769501D-05, 2.94791045D-08,
     2  -1.09331511D-11, 6.68393932D+04, 6.22296438D+00, 3.16780652D+00,
     3   4.75221902D-03,-1.83787077D-06, 3.04190252D-10,-1.77232770D-14,
     4   6.71210650D+04, 6.63589475D+00/
C
C-------------------------  CH2OH            -------------------------
C
      DATA P029
     1 / 4.47832317D+00,-1.35069687D-03, 2.78483707D-05,-3.64867397D-08,
     2   1.47906775D-11,-3.52476728D+03, 3.30911984D+00, 5.09312037D+00,
     3   5.94758550D-03,-2.06496524D-06, 3.23006703D-10,-1.88125052D-14,
     4  -4.05813228D+03,-1.84690613D+00/
C
C-------------------------  CH3OH            -------------------------
C
      DATA P030
     1 / 5.71539582D+00,-1.52309129D-02, 6.52441155D-05,-7.10806889D-08,
     2   2.61352698D-11,-2.56427656D+04,-1.50409823D+00, 1.78970791D+00,
     3   1.40938292D-02,-6.36500835D-06, 1.38171085D-09,-1.17060220D-13,
     4  -2.53748747D+04, 1.45023623D+01/
C
C-------------------------  C2H5OH           -------------------------
C
      DATA P031
     1 / 5.76535800D-01, 2.89451200D-02,-1.61002000D-05, 3.59164100D-09,
     2   0.00000000D+00,-2.96359500D+04, 2.27081300D+01, 4.34717120D+00,
     3   1.86288000D-02,-6.77946700D-06, 8.16592600D-10, 0.00000000D+00,
     4  -3.06615743D+04, 3.24247304D+00/
C
C-------------------------  CH3CHO           -------------------------
C
      DATA P032
     1 / 4.72945950D+00,-3.19328580D-03, 4.75349210D-05,-5.74586110D-08,
     2   2.19311120D-11,-2.15728780D+04, 4.10301590D+00, 5.40411080D+00,
     3   1.17230590D-02,-4.22631370D-06, 6.83724510D-10,-4.09848630D-14,
     4  -2.25931220D+04,-3.48079170D+00/
C
C-------------------------  CH3CHOH          -------------------------
C
      DATA P033
     1 / 1.83974631D+00, 1.87789371D-02,-4.60544253D-06,-2.13116990D-09,
     2   9.43772653D-13,-6.29595195D+03, 2.01446141D+01, 7.26570301D+00,
     3   1.09588926D-02,-3.63662803D-06, 5.53659830D-10,-3.17012322D-14,
     4  -8.64371441D+03,-1.06822851D+01/
C
C-------------------------  CH2CH2OH         -------------------------
C
      DATA P034
     1 / 1.17714711D+00, 2.48115685D-02,-1.50299503D-05, 4.79006785D-09,
     2  -6.40994211D-13,-4.95369043D+03, 2.20081586D+01, 7.52244726D+00,
     3   1.10492715D-02,-3.72576465D-06, 5.72827397D-10,-3.30061759D-14,
     4  -7.29337464D+03,-1.24960750D+01/
C
C-------------------------  CH3CO            -------------------------
C
      DATA P035
     1 / 4.16342570D+00,-2.32616100D-04, 3.42678200D-05,-4.41052270D-08,
     2   1.72756120D-11,-2.65745290D+03, 7.34682800D+00, 5.94477310D+00,
     3   7.86672050D-03,-2.88658820D-06, 4.72708750D-10,-2.85998610D-14,
     4  -3.78730750D+03,-5.01367510D+00/
C
C-------------------------  CH3CH2O          -------------------------
C
      DATA P036
     1 /-2.71296378D-01, 2.98839812D-02,-1.97090548D-05, 6.37339893D-09,
     2  -7.77965054D-13,-3.16397196D+03, 2.47706003D+01, 8.31182392D+00,
     3   1.03426319D-02,-3.39186089D-06, 5.12212617D-10,-2.91601713D-14,
     4  -6.13097954D+03,-2.13985581D+01/
C
C-------------------------  C3H4             -------------------------
C
      DATA P037
     1 / 2.61304450D+00, 1.21225750D-02, 1.85398800D-05,-3.45251490D-08,
     2   1.53350790D-11, 2.15415670D+04, 1.02261390D+01, 6.31687220D+00,
     3   1.11337280D-02,-3.96293780D-06, 6.35642380D-10,-3.78755400D-14,
     4   2.01174950D+04,-1.09957660D+01/
C
C-------------------------  C3H3             -------------------------
C
      DATA P038
     1 / 1.35110927D+00, 3.27411223D-02,-4.73827135D-05, 3.76309808D-08,
     2  -1.18540923D-11, 4.01057783D+04, 1.52058924D+01, 7.14221880D+00,
     3   7.61902005D-03,-2.67459950D-06, 4.24914801D-10,-2.51475415D-14,
     4   3.89087427D+04,-1.25848436D+01/
C
C-------------------------  C3H5             -------------------------
C
      DATA P039
     1 / 1.36318350D+00, 1.98138210D-02, 1.24970600D-05,-3.33555550D-08,
     2   1.58465710D-11, 1.92456290D+04, 1.71732140D+01, 6.50078770D+00,
     3   1.43247310D-02,-5.67816320D-06, 1.10808010D-09,-9.03638870D-14,
     4   1.74824490D+04,-1.12430500D+01/
C
C-------------------------  C3H6             -------------------------
C
      DATA P040
     1 / 1.49330700D+00, 2.09251800D-02, 4.48679400D-06,-1.66891200D-08,
     2   7.15814600D-12, 1.07482600D+03, 1.61453400D+01, 6.73225700D+00,
     3   1.49083400D-02,-4.94989900D-06, 7.21202200D-10,-3.76620400D-14,
     4  -9.23570300D+02,-1.33133500D+01/
C
C-------------------------  C3H8             -------------------------
C
      DATA P041
     1 / 9.28510930D-01, 2.64605660D-02, 6.03324460D-06,-2.19149530D-08,
     2   9.49615440D-12,-1.40579070D+04, 1.92255380D+01, 7.52441520D+00,
     3   1.88982820D-02,-6.29210410D-06, 9.21614570D-10,-4.86844780D-14,
     4  -1.65643940D+04,-1.78383750D+01/
C
C-------------------------  I-C3H7           -------------------------
C
      DATA P042
     1 / 1.44491990D+00, 2.09991120D-02, 7.70362220D-06,-1.84762530D-08,
     2   7.12829620D-12, 9.42237240D+03, 2.01163170D+01, 6.51927410D+00,
     3   1.72201040D-02,-5.73642170D-06, 8.41307320D-10,-4.45659130D-14,
     4   7.32271930D+03,-9.08302150D+00/
C
C-------------------------  N-C3H7           -------------------------
C
      DATA P043
     1 / 1.04911730D+00, 2.60089730D-02, 2.35425160D-06,-1.95951320D-08,
     2   9.37202070D-12, 1.03123460D+04, 2.11360340D+01, 7.70974790D+00,
     3   1.60314850D-02,-5.27202380D-06, 7.58883520D-10,-3.88627190D-14,
     4   7.97622360D+03,-1.55152970D+01/
C
C-------------------------  C4H6             -------------------------
C
      DATA P044
     1 / 1.11078700D+01,-6.30279400D-03, 5.36192000D-05,-5.91451900D-08,
     2   2.12386300D-11, 9.68687700D+03,-2.99286900D+01, 9.84386200D+00,
     3   1.54451700D-02,-5.71720000D-06, 1.01451600D-09,-6.86559300D-14,
     4   9.07722800D+03,-2.80034300D+01/
C
C-------------------------  CHO              -------------------------
C
      DATA P045
     1 / 4.22118584D+00,-3.24392532D-03, 1.37799446D-05,-1.33144093D-08,
     2   4.33768865D-12, 3.83956496D+03, 3.39437243D+00, 2.77217438D+00,
     3   4.95695526D-03,-2.48445613D-06, 5.89161778D-10,-5.33508711D-14,
     4   4.01191815D+03, 9.79834492D+00/
C
C-------------------------  C5H8             -------------------------
C
      DATA P046
     1 / 2.68981400D+00, 2.09545500D-03, 1.13036870D-04,-1.54080700D-07,
     2   6.27636580D-11, 2.31396630D+03, 1.52940560D+01, 7.72447920D+00,
     3   2.83223160D-02,-1.15452360D-05, 2.15408150D-09,-1.50541780D-13,
     4  -7.82615730D+02,-1.97696980D+01/
C
C-------------------------  C7H16            -------------------------
C
      DATA P047
     1 /-1.49861609D+00, 8.56446777D-02,-5.10778379D-05, 1.45322121D-08,
     2  -1.45719946D-12,-2.56089558D+04, 3.65725432D+01, 5.15327664D+00,
     3   6.58448144D-02,-2.99765849D-05, 5.12173194D-09, 0.00000000D+00,
     4  -2.73334877D+04, 2.64465624D+00/
C
C-------------------------  C4H8             -------------------------
C
      DATA P048
     1 / 1.18113800D+00, 3.08533800D-02, 5.08652400D-06,-2.46548800D-08,
     2   1.11101920D-11,-1.79040000D+03, 2.10624700D+01, 2.05358400D+00,
     3   3.43505000D-02,-1.58831960D-05, 3.30896600D-09,-2.53610400D-13,
     4  -2.13972300D+03, 1.55432010D+01/
C
C-------------------------  C5H10            -------------------------
C
      DATA P049
     1 /-1.36099451D+00, 5.87684877D-02,-3.93332993D-05, 1.36933176D-08,
     2  -1.93663167D-12,-4.42142034D+03, 3.36502744D+01, 4.62768794D+00,
     3   3.99163471D-02,-1.76584355D-05, 2.94527974D-09, 0.00000000D+00,
     4  -5.90930391D+03, 3.39515508D+00/
C
C-------------------------  AR               -------------------------
C
      DATA P050
     1 / 2.50000000D+00, 0.00000000D+00, 0.00000000D+00, 0.00000000D+00,
     2   0.00000000D+00,-7.45375000D+02, 4.36600000D+00, 2.50000000D+00,
     3   0.00000000D+00, 0.00000000D+00, 0.00000000D+00, 0.00000000D+00,
     4  -7.45375000D+02, 4.36600000D+00/
C
C-------------------------  HE               -------------------------
C
      DATA P051
     1 / 2.50000000D+00, 0.00000000D+00, 0.00000000D+00, 0.00000000D+00,
     2   0.00000000D+00,-7.45375000D+02, 9.28723974D-01, 2.50000000D+00,
     3   0.00000000D+00, 0.00000000D+00, 0.00000000D+00, 0.00000000D+00,
     4  -7.45375000D+02, 9.28723974D-01/
C
C-------------------------  N2               -------------------------
C
      DATA P052
     1 / 3.29867700D+00, 1.40824040D-03,-3.96322200D-06, 5.64151500D-09,
     2  -2.44485400D-12,-1.02089990D+03, 3.95037200D+00, 2.92664000D+00,
     3   1.48797680D-03,-5.68476000D-07, 1.00970380D-10,-6.75335100D-15,
     4  -9.22797700D+02, 5.98052800D+00/
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC    Curve Fitts for the Collision Integrals                         CC
CC     Omega(2,2)* = A1 + A2*T + A3*T^2 + A4*T^3 +A5*T^4+A6*T^5       CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C--H2             
      DATA C001/ 1.42331D+00,-4.51194D-03, 1.39762D-05,-2.21917D-08,
     1           1.71345D-11,-5.11475D-15, 8.82237D-01,-2.55571D-04,
     2           1.15843D-07,-3.16256D-11, 4.54502D-15,-2.61452D-19/
C--O2             
      DATA C002/ 1.83482D+00,-5.18351D-03, 1.28542D-05,-1.71840D-08,
     1           1.16977D-11,-3.18626D-15, 9.26460D-01,-8.61554D-05,
     2          -1.98416D-08, 1.52914D-11,-2.97272D-15, 1.93292D-19/
C--H              
      DATA C003/ 2.21020D+00,-6.97744D-03, 1.71600D-05,-2.27054D-08,
     1           1.53676D-11,-4.18063D-15, 1.05302D+00,-2.75259D-04,
     2           1.23285D-07,-3.45219D-11, 5.05909D-15,-2.93579D-19/
C--O              
      DATA C004/ 1.59521D+00,-4.21102D-03, 1.10383D-05,-1.57625D-08,
     1           1.14263D-11,-3.28263D-15, 1.01432D+00,-3.15588D-04,
     2           1.36807D-07,-3.40444D-11, 4.36962D-15,-2.24609D-19/
C--OH             
      DATA C005/ 1.59521D+00,-4.21102D-03, 1.10383D-05,-1.57625D-08,
     1           1.14263D-11,-3.28263D-15, 1.01432D+00,-3.15588D-04,
     2           1.36807D-07,-3.40444D-11, 4.36962D-15,-2.24609D-19/
C--H2O            
      DATA C006/ 4.10756D+00,-1.02196D-02, 1.72825D-05,-1.69771D-08,
     1           9.11498D-12,-2.07347D-15, 1.98442D+00,-1.22760D-03,
     2           6.16069D-07,-1.63179D-10, 2.16928D-14,-1.14101D-18/
C--HO2            
      DATA C007/ 1.83482D+00,-5.18351D-03, 1.28542D-05,-1.71840D-08,
     1           1.16977D-11,-3.18626D-15, 9.26460D-01,-8.61554D-05,
     2          -1.98416D-08, 1.52914D-11,-2.97272D-15, 1.93292D-19/
C--H2O2           
      DATA C008/ 1.83482D+00,-5.18351D-03, 1.28542D-05,-1.71840D-08,
     1           1.16977D-11,-3.18626D-15, 9.26460D-01,-8.61554D-05,
     2          -1.98416D-08, 1.52914D-11,-2.97272D-15, 1.93292D-19/
C--CO             
      DATA C009/ 1.73729D+00,-4.69277D-03, 1.16287D-05,-1.55874D-08,
     1           1.06386D-11,-2.90277D-15, 9.46466D-01,-1.39189D-04,
     2           1.34202D-08, 5.71067D-12,-1.67109D-15, 1.25787D-19/
C--CO2            
      DATA C010/ 3.15911D+00,-1.13014D-02, 2.71561D-05,-3.52182D-08,
     1           2.34783D-11,-6.31013D-15, 1.24871D+00,-4.17811D-04,
     2           1.60342D-07,-2.98842D-11, 2.21619D-15,-2.45645D-20/
C--HCO            
      DATA C011/ 4.08126D+00,-1.14450D-02, 2.15265D-05,-2.30986D-08,
     1           1.32885D-11,-3.17936D-15, 1.81437D+00,-1.04544D-03,
     2           5.20251D-07,-1.36883D-10, 1.80263D-14,-9.35709D-19/
C--CH2O           
      DATA C012/ 4.08126D+00,-1.14450D-02, 2.15265D-05,-2.30986D-08,
     1           1.32885D-11,-3.17936D-15, 1.81437D+00,-1.04544D-03,
     2           5.20251D-07,-1.36883D-10, 1.80263D-14,-9.35709D-19/
C--CH4            
      DATA C013/ 2.17416D+00,-6.80497D-03, 1.67459D-05,-2.21759D-08,
     1           1.50217D-11,-4.09057D-15, 1.03463D+00,-2.46098D-04,
     2           1.03015D-07,-2.79940D-11, 4.07597D-15,-2.37526D-19/
C--CH3            
      DATA C014/ 2.20054D+00,-6.93349D-03, 1.70616D-05,-2.25906D-08,
     1           1.53010D-11,-4.16577D-15, 1.04755D+00,-2.66510D-04,
     2           1.17229D-07,-3.25773D-11, 4.76662D-15,-2.76904D-19/
C--T-CH2          
      DATA C015/ 2.20054D+00,-6.93349D-03, 1.70616D-05,-2.25906D-08,
     1           1.53010D-11,-4.16577D-15, 1.04755D+00,-2.66510D-04,
     2           1.17229D-07,-3.25773D-11, 4.76662D-15,-2.76904D-19/
C--S-CH2          
      DATA C016/ 2.20054D+00,-6.93349D-03, 1.70616D-05,-2.25906D-08,
     1           1.53010D-11,-4.16577D-15, 1.04755D+00,-2.66510D-04,
     2           1.17229D-07,-3.25773D-11, 4.76662D-15,-2.76904D-19/
C--C2H4           
      DATA C017/ 3.11195D+00,-1.11072D-02, 2.67356D-05,-3.47073D-08,
     1           2.31483D-11,-6.22223D-15, 1.25413D+00,-4.41052D-04,
     2           1.79778D-07,-3.67497D-11, 3.29715D-15,-8.70183D-20/
C--CH3O           
      DATA C018/ 3.97418D+00,-1.25472D-02, 2.61990D-05,-3.06465D-08,
     1           1.89290D-11,-4.80152D-15, 1.66188D+00,-9.18599D-04,
     2           4.74783D-07,-1.31173D-10, 1.81965D-14,-9.93878D-19/
C--C2H5           
      DATA C019/ 3.18243D+00,-1.13684D-02, 2.72351D-05,-3.52387D-08,
     1           2.34497D-11,-6.29350D-15, 1.24781D+00,-4.08862D-04,
     2           1.52707D-07,-2.72650D-11, 1.82713D-15,-3.91177D-21/
C--C2H6           
      DATA C020/ 3.18243D+00,-1.13684D-02, 2.72351D-05,-3.52387D-08,
     1           2.34497D-11,-6.29350D-15, 1.24781D+00,-4.08862D-04,
     2           1.52707D-07,-2.72650D-11, 1.82713D-15,-3.91177D-21/
C--CH             
      DATA C021/ 1.59521D+00,-4.21102D-03, 1.10383D-05,-1.57625D-08,
     1           1.14263D-11,-3.28263D-15, 1.01432D+00,-3.15588D-04,
     2           1.36807D-07,-3.40444D-11, 4.36962D-15,-2.24609D-19/
C--C2H2           
      DATA C022/ 3.31662D+00,-1.18547D-02, 2.81472D-05,-3.61868D-08,
     1           2.39704D-11,-6.41106D-15, 1.26713D+00,-4.15356D-04,
     2           1.55219D-07,-2.90795D-11, 2.41477D-15,-6.16491D-20/
C--C2H3           
      DATA C023/ 3.31662D+00,-1.18547D-02, 2.81472D-05,-3.61868D-08,
     1           2.39704D-11,-6.41106D-15, 1.26713D+00,-4.15356D-04,
     2           1.55219D-07,-2.90795D-11, 2.41477D-15,-6.16491D-20/
C--CH2CHO         
      DATA C024/ 4.01903D+00,-1.24157D-02, 2.54360D-05,-2.93102D-08,
     1           1.78817D-11,-4.48819D-15, 1.70159D+00,-9.56072D-04,
     2           4.90776D-07,-1.34084D-10, 1.83536D-14,-9.88286D-19/
C--C2H4O          
      DATA C025/ 4.01903D+00,-1.24157D-02, 2.54360D-05,-2.93102D-08,
     1           1.78817D-11,-4.48819D-15, 1.70159D+00,-9.56072D-04,
     2           4.90776D-07,-1.34084D-10, 1.83536D-14,-9.88286D-19/
C--CH2CO          
      DATA C026/ 4.01903D+00,-1.24157D-02, 2.54360D-05,-2.93102D-08,
     1           1.78817D-11,-4.48819D-15, 1.70159D+00,-9.56072D-04,
     2           4.90776D-07,-1.34084D-10, 1.83536D-14,-9.88286D-19/
C--HCCO           
      DATA C027/ 2.27092D+00,-7.31260D-03, 1.80643D-05,-2.39862D-08,
     1           1.62837D-11,-4.44060D-15, 1.07970D+00,-3.17783D-04,
     2           1.52625D-07,-4.39208D-11, 6.47102D-15,-3.74110D-19/
C--C2H            
      DATA C028/ 3.31662D+00,-1.18547D-02, 2.81472D-05,-3.61868D-08,
     1           2.39704D-11,-6.41106D-15, 1.26713D+00,-4.15356D-04,
     2           1.55219D-07,-2.90795D-11, 2.41477D-15,-6.16491D-20/
C--CH2OH          
      DATA C029/ 3.97418D+00,-1.25472D-02, 2.61990D-05,-3.06465D-08,
     1           1.89290D-11,-4.80152D-15, 1.66188D+00,-9.18599D-04,
     2           4.74783D-07,-1.31173D-10, 1.81965D-14,-9.93878D-19/
C--CH3OH          
      DATA C030/ 4.07329D+00,-1.17392D-02, 2.26149D-05,-2.47481D-08,
     1           1.44576D-11,-3.49837D-15, 1.78383D+00,-1.01920D-03,
     2           5.09705D-07,-1.34869D-10, 1.78507D-14,-9.30097D-19/
C--C2H5OH         
      DATA C031/ 4.05895D+00,-1.18778D-02, 2.31763D-05,-2.56093D-08,
     1           1.50657D-11,-3.66277D-15, 1.76423D+00,-1.00445D-03,
     2           5.05256D-07,-1.34593D-10, 1.79330D-14,-9.40105D-19/
C--CH3CHO         
      DATA C032/ 4.01903D+00,-1.24157D-02, 2.54360D-05,-2.93102D-08,
     1           1.78817D-11,-4.48819D-15, 1.70159D+00,-9.56072D-04,
     2           4.90776D-07,-1.34084D-10, 1.83536D-14,-9.88286D-19/
C--CH3CHOH        
      DATA C033/ 3.82127D+00,-1.28354D-02, 2.82945D-05,-3.44746D-08,
     1           2.19658D-11,-5.70984D-15, 1.51519D+00,-7.37262D-04,
     2           3.68879D-07,-1.00820D-10, 1.40518D-14,-7.79520D-19/
C--CH2CH2OH       
      DATA C034/ 3.82127D+00,-1.28354D-02, 2.82945D-05,-3.44746D-08,
     1           2.19658D-11,-5.70984D-15, 1.51519D+00,-7.37262D-04,
     2           3.68879D-07,-1.00820D-10, 1.40518D-14,-7.79520D-19/
C--CH3CO          
      DATA C035/ 4.01903D+00,-1.24157D-02, 2.54360D-05,-2.93102D-08,
     1           1.78817D-11,-4.48819D-15, 1.70159D+00,-9.56072D-04,
     2           4.90776D-07,-1.34084D-10, 1.83536D-14,-9.88286D-19/
C--CH3CH2O        
      DATA C036/ 4.05895D+00,-1.18778D-02, 2.31763D-05,-2.56093D-08,
     1           1.50657D-11,-3.66277D-15, 1.76423D+00,-1.00445D-03,
     2           5.05256D-07,-1.34593D-10, 1.79330D-14,-9.40105D-19/
C--C3H4           
      DATA C037/ 3.65330D+00,-1.26077D-02, 2.85201D-05,-3.53121D-08,
     1           2.27052D-11,-5.92784D-15, 1.39999D+00,-5.76229D-04,
     2           2.62344D-07,-6.61699D-11, 8.67119D-15,-4.62150D-19/
C--C3H3           
      DATA C038/ 3.65330D+00,-1.26077D-02, 2.85201D-05,-3.53121D-08,
     1           2.27052D-11,-5.92784D-15, 1.39999D+00,-5.76229D-04,
     2           2.62344D-07,-6.61699D-11, 8.67119D-15,-4.62150D-19/
C--C3H5           
      DATA C039/ 3.61115D+00,-1.25468D-02, 2.85919D-05,-3.55925D-08,
     1           2.29754D-11,-6.01540D-15, 1.37551D+00,-5.42683D-04,
     2           2.39677D-07,-5.85540D-11, 7.44637D-15,-3.87484D-19/
C--C3H6           
      DATA C040/ 3.56976D+00,-1.24763D-02, 2.86234D-05,-3.58081D-08,
     1           2.31976D-11,-6.08921D-15, 1.35403D+00,-5.14170D-04,
     2           2.20481D-07,-5.20547D-11, 6.38858D-15,-3.22220D-19/
C--C3H8           
      DATA C041/ 3.55052D+00,-1.24694D-02, 2.87590D-05,-3.61377D-08,
     1           2.35002D-11,-6.18916D-15, 1.34149D+00,-4.97038D-04,
     2           2.08937D-07,-4.81664D-11, 5.75880D-15,-2.83483D-19/
C--I-C3H7         
      DATA C042/ 3.55052D+00,-1.24694D-02, 2.87590D-05,-3.61377D-08,
     1           2.35002D-11,-6.18916D-15, 1.34149D+00,-4.97038D-04,
     2           2.08937D-07,-4.81664D-11, 5.75880D-15,-2.83483D-19/
C--N-C3H7         
      DATA C043/ 3.55052D+00,-1.24694D-02, 2.87590D-05,-3.61377D-08,
     1           2.35002D-11,-6.18916D-15, 1.34149D+00,-4.97038D-04,
     2           2.08937D-07,-4.81664D-11, 5.75880D-15,-2.83483D-19/
C--C4H6           
      DATA C044/ 3.79676D+00,-1.28003D-02, 2.83182D-05,-3.45861D-08,
     1           2.20739D-11,-5.74544D-15, 1.49759D+00,-7.12790D-04,
     2           3.53030D-07,-9.57971D-11, 1.32928D-14,-7.35913D-19/
C--CHO            
      DATA C045/ 4.08126D+00,-1.14450D-02, 2.15265D-05,-2.30986D-08,
     1           1.32885D-11,-3.17936D-15, 1.81437D+00,-1.04544D-03,
     2           5.20251D-07,-1.36883D-10, 1.80263D-14,-9.35709D-19/
C--C5H8           
      DATA C046/ 3.96091D+00,-1.26867D-02, 2.68312D-05,-3.17260D-08,
     1           1.97775D-11,-5.05769D-15, 1.64166D+00,-8.97483D-04,
     2           4.64414D-07,-1.28770D-10, 1.79519D-14,-9.86111D-19/
C--C7H16          
      DATA C047/ 4.04598D+00,-1.20356D-02, 2.38278D-05,-2.66527D-08,
     1           1.58418D-11,-3.88517D-15, 1.74443D+00,-9.89285D-04,
     2           5.00718D-07,-1.34388D-10, 1.80458D-14,-9.53172D-19/
C--C4H8           
      DATA C048/ 3.79188D+00,-1.28228D-02, 2.84462D-05,-3.48203D-08,
     1           2.22640D-11,-5.80369D-15, 1.49199D+00,-7.05257D-04,
     2           3.48164D-07,-9.42461D-11, 1.30565D-14,-7.22244D-19/
C--C5H10          
      DATA C049/ 3.89078D+00,-1.27136D-02, 2.73518D-05,-3.27233D-08,
     1           2.05651D-11,-5.28956D-15, 1.58588D+00,-8.31550D-04,
     2           4.27515D-07,-1.18681D-10, 1.66486D-14,-9.23132D-19/
C--AR             
      DATA C050/ 2.11828D+00,-6.49925D-03, 1.59012D-05,-2.09257D-08,
     1           1.40746D-11,-3.80369D-15, 1.01108D+00,-2.08453D-04,
     2           7.62543D-08,-1.92118D-11, 2.73295D-15,-1.60018D-19/
C--HE             
      DATA C051/ 9.26058D-01,-1.34821D-03, 2.96811D-06,-3.86963D-09,
     1           2.64318D-12,-7.29700D-16, 6.54685D-01,-5.49840D-05,
     2          -2.02748D-08, 1.38008D-11,-2.72860D-15, 1.83695D-19/
C--N2             
      DATA C052/ 1.73041D+00,-4.65219D-03, 1.15074D-05,-1.53948D-08,
     1           1.04823D-11,-2.85205D-15, 9.49062D-01,-1.45201D-04,
     2           1.74761D-08, 4.43706D-12,-1.48140D-15, 1.14964D-19/
C-----------------------------------------------------------------------
C-- LREV = 1 for one-way reactions; LREV = 2 for reversible reactions 
C         AKF(IR) -k  ==> Pressure-Dependent Fall-Off Reaction    
C            TROE(IR,1) = 0.0 ==> Lindemann Form (RFACT=1.0)
C            TROE(IR,1) = +a  ==> Troe Form
C            TROE(IR,1) = -a  ==> SRI Form
C            TROE(IR,2) = 0.0  ==> Linear Form
C            AKLOW(IR) = 0.0  ==> TSENG Form (K0=Kinf*polynomial)
C-----------------------REACTION RATES (Maurice)------------------------
      DATA LREV / 2 /
C-----------------------------------------------------------------------     
      DATA AKF(   1),ALF(   1),EXAF(   1)
     1                         / 3.520D+16,-7.000D-01, 1.7070D+04/
      DATA AKF(   3),ALF(   3),EXAF(   3)
     1                         / 5.060D+04, 2.670D+00, 6.2906D+03/
      DATA AKF(   5),ALF(   5),EXAF(   5)
     1                         / 1.170D+09, 1.300D+00, 3.6352D+03/
      DATA AKF(   7),ALF(   7),EXAF(   7)
     1                         / 7.600D+00, 3.840D+00, 1.2780D+04/
      DATA AKF(   9),ALF(   9),EXAF(   9)
     1                         / 1.300D+18,-1.000D+00, 0.0000D+00/
      DATA AKF(  11),ALF(  11),EXAF(  11)
     1                         / 4.000D+22,-2.000D+00, 0.0000D+00/
      DATA AKF(  13),ALF(  13),EXAF(  13)
     1                         / 6.170D+15,-5.000D-01, 0.0000D+00/
      DATA AKF(  15),ALF(  15),EXAF(  15)
     1                         / 4.710D+18,-1.000D+00, 0.0000D+00/
      DATA AKF(  17),ALF(  17),EXAF(  17)
     1                         / 8.000D+15, 0.000D+00, 0.0000D+00/
      DATA AKF(  19),ALF(  19),EXAF(  19)
     1                         /-4.650D+12, 4.400D-01, 0.0000D+00/
           DATA AKLOW(  19),ALLOW(  19),EALOW(  19)                             
     1     /5.75D+19,   -1.4,  0.0 /                                            
           DATA TROE(  19,1),TROE(  19,2),TROE(  19,3),TROE(  19,4)             
     1     / 0.5,   1.0D-30,      1.0D+30, 0.0  /                                    
      DATA AKF(  21),ALF(  21),EXAF(  21)
     1                         / 7.080D+13, 0.000D+00, 2.9500D+02/
      DATA AKF(  23),ALF(  23),EXAF(  23)
     1                         / 1.660D+13, 0.000D+00, 8.2300D+02/
      DATA AKF(  25),ALF(  25),EXAF(  25)
     1                         / 3.100D+13, 0.000D+00, 1.7208D+03/
      DATA AKF(  27),ALF(  27),EXAF(  27)
     1                         / 2.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF(  29),ALF(  29),EXAF(  29)
     1                         / 2.890D+13, 0.000D+00,-4.9710D+02/
      DATA AKF(  31),ALF(  31),EXAF(  31)
     1                         /-7.400D+13,-3.700D-01, 0.0000D+00/
           DATA AKLOW(  31),ALLOW(  31),EALOW(  31)                             
     1     /  2.300D+18,  -0.900, -1701.72 /                                    
           DATA TROE(  31,1),TROE(  31,2),TROE(  31,3),TROE(  31,4)             
     1     /   0.735,      94.0,     1756.0,     5182.0 /                       
      DATA AKF(  33),ALF(  33),EXAF(  33)
     1                         / 3.020D+12, 0.000D+00, 1.3862D+03/
      DATA AKF(  35),ALF(  35),EXAF(  35)
     1                         / 4.790D+13, 0.000D+00, 7.9588D+03/
      DATA AKF(  37),ALF(  37),EXAF(  37)
     1                         / 1.000D+13, 0.000D+00, 3.5850D+03/
      DATA AKF(  39),ALF(  39),EXAF(  39)
     1                         / 7.080D+12, 0.000D+00, 1.4340D+03/
      DATA AKF(  41),ALF(  41),EXAF(  41)
     1                         / 9.630D+06, 2.000D+00, 3.9914D+03/
      DATA AKF(  43),ALF(  43),EXAF(  43)
     1                         / 4.400D+06, 1.500D+00,-7.4090D+02/
      DATA AKF(  45),ALF(  45),EXAF(  45)
     1                         / 6.000D+13, 0.000D+00, 2.2945D+04/
      DATA AKF(  47),ALF(  47),EXAF(  47)
     1                         / 1.000D+12, 0.000D+00, 4.7700D+04/
      DATA AKF(  49),ALF(  49),EXAF(  49)
     1                         / 1.860D+17,-1.000D+00, 1.7000D+04/
      DATA AKF(  51),ALF(  51),EXAF(  51)
     1                         / 5.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF(  53),ALF(  53),EXAF(  53)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF(  55),ALF(  55),EXAF(  55)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF(  57),ALF(  57),EXAF(  57)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF(  59),ALF(  59),EXAF(  59)
     1                         / 7.580D+12, 0.000D+00, 4.1000D+02/
      DATA AKF(  61),ALF(  61),EXAF(  61)
     1                         / 5.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF(  63),ALF(  63),EXAF(  63)
     1                         /-1.090D+12, 4.800D-01,-2.6000D+02/
           DATA AKLOW(  63),ALLOW(  63),EALOW(  63)                             
     1     /  1.350D+24,  -2.570,   425.00/                                     
           DATA TROE(  63,1),TROE(  63,2),TROE(  63,3),TROE(  63,4)             
     1     /   .7824, 271.00, 2755.00, 6570.00 /                                
      DATA AKF(  65),ALF(  65),EXAF(  65)
     1                         / 5.740D+07, 1.900D+00, 2.7486D+03/
      DATA AKF(  67),ALF(  67),EXAF(  67)
     1                         / 3.500D+13, 0.000D+00, 3.5133D+03/
      DATA AKF(  69),ALF(  69),EXAF(  69)
     1                         / 3.900D+10, 8.900D-01, 4.0630D+02/
      DATA AKF(  71),ALF(  71),EXAF(  71)
     1                         / 6.000D+13, 0.000D+00, 4.0674D+04/
      DATA AKF(  73),ALF(  73),EXAF(  73)
     1                         / 4.110D+04, 2.500D+00, 1.0210D+04/
      DATA AKF(  75),ALF(  75),EXAF(  75)
     1                         / 1.300D+04, 3.000D+00, 8.0377D+03/
      DATA AKF(  77),ALF(  77),EXAF(  77)
     1                         / 1.600D+07, 1.830D+00, 2.7820D+03/
      DATA AKF(  79),ALF(  79),EXAF(  79)
     1                         / 1.900D+09, 1.440D+00, 8.6759D+03/
      DATA AKF(  81),ALF(  81),EXAF(  81)
     1                         / 3.980D+13, 0.000D+00, 5.6891D+04/
      DATA AKF(  83),ALF(  83),EXAF(  83)
     1                         / 9.030D+12, 0.000D+00, 2.4641D+04/
      DATA AKF(  85),ALF(  85),EXAF(  85)
     1                         / 1.800D+14, 0.000D+00, 1.5105D+04/
      DATA AKF(  87),ALF(  87),EXAF(  87)
     1                         / 1.550D+14, 0.000D+00, 1.3480D+04/
      DATA AKF(  89),ALF(  89),EXAF(  89)
     1                         / 4.000D+13, 0.000D+00, 2.5023D+03/
      DATA AKF(  91),ALF(  91),EXAF(  91)
     1                         / 8.430D+13, 0.000D+00, 0.0000D+00/
      DATA AKF(  93),ALF(  93),EXAF(  93)
     1                         / 4.220D+13, 0.000D+00, 0.0000D+00/
      DATA AKF(  95),ALF(  95),EXAF(  95)
     1                         / 5.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF(  97),ALF(  97),EXAF(  97)
     1                         / 3.300D+11, 0.000D+00, 8.9412D+03/
      DATA AKF(  99),ALF(  99),EXAF(  99)
     1                         / 1.100D+13, 0.000D+00, 2.7820D+04/
      DATA AKF( 101),ALF( 101),EXAF( 101)
     1                         / 1.000D+14, 0.000D+00, 3.2003D+04/
      DATA AKF( 103),ALF( 103),EXAF( 103)
     1                         / 3.160D+13, 0.000D+00, 1.4699D+04/
      DATA AKF( 105),ALF( 105),EXAF( 105)
     1                         /-1.270D+16,-6.300D-01, 3.8300D+02/
           DATA AKLOW( 105),ALLOW( 105),EALOW( 105)                             
     1     /  2.470D+33,  -4.760,  2440.00/                                     
           DATA TROE( 105,1),TROE( 105,2),TROE( 105,3),TROE( 105,4)             
     1     /  0.7830,  74.00, 2941.00, 6964.00 /                                
      DATA AKF( 107),ALF( 107),EXAF( 107)
     1                         /-1.810D+13, 0.000D+00, 0.0000D+00/
           DATA AKLOW( 107),ALLOW( 107),EALOW( 107)                             
     1     /  1.270D+41,  -7.000,  2762.91 /                                    
           DATA TROE( 107,1),TROE( 107,2),TROE( 107,3),TROE( 107,4)             
     1     /    0.62,   73.00,  1180.00,  0.0 /                                
      DATA AKF( 109),ALF( 109),EXAF( 109)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 111),ALF( 111),EXAF( 111)
     1                         / 3.130D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 113),ALF( 113),EXAF( 113)
     1                         / 3.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 115),ALF( 115),EXAF( 115)
     1                         / 6.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 117),ALF( 117),EXAF( 117)
     1                         / 6.020D+12, 0.000D+00,-1.7877D+03/
      DATA AKF( 119),ALF( 119),EXAF( 119)
     1                         / 2.500D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 121),ALF( 121),EXAF( 121)
     1                         / 1.130D+07, 2.000D+00, 2.9995D+03/
      DATA AKF( 123),ALF( 123),EXAF( 123)
     1                         / 8.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 125),ALF( 125),EXAF( 125)
     1                         / 4.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 127),ALF( 127),EXAF( 127)
     1                         / 2.630D+12, 0.000D+00, 1.4914D+03/
      DATA AKF( 129),ALF( 129),EXAF( 129)
     1                         / 6.580D+12, 0.000D+00, 1.4914D+03/
      DATA AKF( 131),ALF( 131),EXAF( 131)
     1                         / 1.000D+14, 0.000D+00, 0.0000D+00/
      DATA AKF( 133),ALF( 133),EXAF( 133)
     1                         / 4.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 135),ALF( 135),EXAF( 135)
     1                         / 1.770D+11, 7.600D-01,-4.7800D+02/
      DATA AKF( 137),ALF( 137),EXAF( 137)
     1                         / 1.170D+15,-7.500D-01, 0.0000D+00/
      DATA AKF( 139),ALF( 139),EXAF( 139)
     1                         / 4.800D+01, 3.220D+00,-3.2265D+03/
      DATA AKF( 141),ALF( 141),EXAF( 141)
     1                         / 2.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 143),ALF( 143),EXAF( 143)
     1                         / 1.600D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 145),ALF( 145),EXAF( 145)
     1                         / 5.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 147),ALF( 147),EXAF( 147)
     1                         / 1.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 149),ALF( 149),EXAF( 149)
     1                         / 4.280D-13, 7.600D+00,-3.5372D+03/
      DATA AKF( 151),ALF( 151),EXAF( 151)
     1                         / 7.780D+13, 0.000D+00, 1.3513D+04/
      DATA AKF( 153),ALF( 153),EXAF( 153)
     1                         / 5.400D+02, 3.500D+00, 5.2103D+03/
      DATA AKF( 155),ALF( 155),EXAF( 155)
     1                         / 1.400D+00, 4.300D+00, 2.7724D+03/
      DATA AKF( 157),ALF( 157),EXAF( 157)
     1                         / 2.200D+07, 1.900D+00, 1.1233D+03/
      DATA AKF( 159),ALF( 159),EXAF( 159)
     1                         / 5.500D-01, 4.000D+00, 8.2935D+03/
      DATA AKF( 161),ALF( 161),EXAF( 161)
     1                         /-8.850D+20,-1.230D+00, 1.0222D+05/
           DATA AKLOW( 161),ALLOW( 161),EALOW( 161)                             
     1     /  4.900D+42,  -6.430, 107170.17 /                                    
           DATA TROE( 161,1),TROE( 161,2),TROE( 161,3),TROE( 161,4)             
     1     /    0.84,  125.00,  2219.00,  6882.00 /                             
      DATA AKF( 163),ALF( 163),EXAF( 163)
     1                         / 1.320D+13, 0.000D+00, 2.0470D+04/
      DATA AKF( 165),ALF( 165),EXAF( 165)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 167),ALF( 167),EXAF( 167)
     1                         / 3.060D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 169),ALF( 169),EXAF( 169)
     1                         / 4.240D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 171),ALF( 171),EXAF( 171)
     1                         / 2.000D+12, 0.000D+00, 4.9952D+03/
      DATA AKF( 173),ALF( 173),EXAF( 173)
     1                         /-1.110D+10, 1.037D+00, 3.6769D+04/
           DATA AKLOW( 173),ALLOW( 173),EALOW( 173)                             
     1     /  3.990D+33,  -4.990, 40000.00 /                                    
           DATA TROE( 173,1),TROE( 173,2),TROE( 173,3),TROE( 173,4)             
     1     /   0.168, 1203.00,     1.0D-30,  0.00 /                                
      DATA AKF( 175),ALF( 175),EXAF( 175)
     1                         / 4.490D+07, 2.120D+00, 1.3360D+04/
      DATA AKF( 177),ALF( 177),EXAF( 177)
     1                         / 5.530D+05, 2.310D+00, 2.9636D+03/
      DATA AKF( 179),ALF( 179),EXAF( 179)
     1                         / 2.250D+06, 2.080D+00, 0.0000D+00/
      DATA AKF( 181),ALF( 181),EXAF( 181)
     1                         / 1.210D+06, 2.080D+00, 0.0000D+00/
      DATA AKF( 183),ALF( 183),EXAF( 183)
     1                         / 5.010D+14, 0.000D+00, 6.4700D+04/
      DATA AKF( 185),ALF( 185),EXAF( 185)
     1                         / 4.220D+13, 0.000D+00, 5.7623D+04/
      DATA AKF( 187),ALF( 187),EXAF( 187)
     1                         / 2.230D+12, 0.000D+00, 1.7189D+04/
      DATA AKF( 189),ALF( 189),EXAF( 189)
     1                         / 4.000D+12, 0.000D+00, 1.7008D+04/
      DATA AKF( 191),ALF( 191),EXAF( 191)
     1                         / 2.600D+17, 0.000D+00, 9.6568D+04/
      DATA AKF( 193),ALF( 193),EXAF( 193)
     1                         / 3.500D+16, 0.000D+00, 7.1532D+04/
      DATA AKF( 195),ALF( 195),EXAF( 195)
     1                         / 4.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 197),ALF( 197),EXAF( 197)
     1                         /-6.380D+09, 1.000D+00, 3.7627D+04/
           DATA AKLOW( 197),ALLOW( 197),EALOW( 197)                             
     1     /  1.510D+14,   0.100, 32686.42 /                                    
           DATA TROE( 197,1),TROE( 197,2),TROE( 197,3),TROE( 197,4)             
     1     /   0.3, 1.0D+30,   1.0D-30, 0.0 /                                          
      DATA AKF( 199),ALF( 199),EXAF( 199)
     1                         / 1.700D+29,-5.312D+00, 6.5031D+03/
      DATA AKF( 201),ALF( 201),EXAF( 201)
     1                         / 7.000D+14,-6.110D-01, 5.2624D+03/
      DATA AKF( 203),ALF( 203),EXAF( 203)
     1                         / 5.190D+15,-1.260D+00, 3.3126D+03/
      DATA AKF( 205),ALF( 205),EXAF( 205)
     1                         / 4.000D+14, 0.000D+00, 1.0660D+04/
      DATA AKF( 207),ALF( 207),EXAF( 207)
     1                         / 1.600D+14, 0.000D+00, 9.8948D+03/
      DATA AKF( 209),ALF( 209),EXAF( 209)
     1                         / 4.600D+15,-5.400D-01, 4.4933D+04/
      DATA AKF( 211),ALF( 211),EXAF( 211)
     1                         / 1.900D+07, 1.700D+00, 9.9900D+02/
      DATA AKF( 213),ALF( 213),EXAF( 213)
     1                         / 3.370D+07, 2.000D+00, 1.4001D+04/
      DATA AKF( 215),ALF( 215),EXAF( 215)
     1                         / 1.500D+09, 1.430D+00, 2.6888D+03/
      DATA AKF( 217),ALF( 217),EXAF( 217)
     1                         / 2.000D+13, 0.000D+00, 2.2944D+03/
      DATA AKF( 219),ALF( 219),EXAF( 219)
     1                         / 1.000D+13, 0.000D+00, 2.0004D+03/
      DATA AKF( 221),ALF( 221),EXAF( 221)
     1                         / 9.000D+10, 0.000D+00, 0.0000D+00/
      DATA AKF( 223),ALF( 223),EXAF( 223)
     1                         / 1.500D+14, 0.000D+00, 0.0000D+00/
      DATA AKF( 225),ALF( 225),EXAF( 225)
     1                         / 2.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 227),ALF( 227),EXAF( 227)
     1                         / 9.640D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 229),ALF( 229),EXAF( 229)
     1                         / 2.880D+07, 1.700D+00, 1.0014D+03/
      DATA AKF( 231),ALF( 231),EXAF( 231)
     1                         / 1.400D+07, 1.700D+00, 1.0014D+03/
      DATA AKF( 233),ALF( 233),EXAF( 233)
     1                         / 2.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 235),ALF( 235),EXAF( 235)
     1                         / 1.020D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 237),ALF( 237),EXAF( 237)
     1                         / 6.020D+11, 0.000D+00, 0.0000D+00/
      DATA AKF( 239),ALF( 239),EXAF( 239)
     1                         / 4.500D+15, 0.000D+00, 2.5096D+04/
      DATA AKF( 241),ALF( 241),EXAF( 241)
     1                         / 2.410D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 243),ALF( 243),EXAF( 243)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 245),ALF( 245),EXAF( 245)
     1                         / 2.500D+17,-9.300D-01, 5.1268D+03/
      DATA AKF( 247),ALF( 247),EXAF( 247)
     1                         / 2.400D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 249),ALF( 249),EXAF( 249)
     1                         / 5.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 251),ALF( 251),EXAF( 251)
     1                         / 5.000D+13, 0.000D+00, 2.5120D+04/
      DATA AKF( 253),ALF( 253),EXAF( 253)
     1                         / 1.000D+14, 0.000D+00, 1.9120D+04/
      DATA AKF( 255),ALF( 255),EXAF( 255)
     1                         / 1.020D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 257),ALF( 257),EXAF( 257)
     1                         / 1.440D+06, 2.000D+00,-8.3890D+02/
      DATA AKF( 259),ALF( 259),EXAF( 259)
     1                         / 4.400D+06, 2.000D+00, 1.5057D+03/
      DATA AKF( 261),ALF( 261),EXAF( 261)
     1                         / 1.353D+03, 3.200D+00, 3.4907D+03/
      DATA AKF( 263),ALF( 263),EXAF( 263)
     1                         / 6.830D+00, 0.000D+00, 3.4000D+00/
      DATA AKF( 265),ALF( 265),EXAF( 265)
     1                         / 1.000D+13, 0.000D+00, 4.6845D+03/
      DATA AKF( 267),ALF( 267),EXAF( 267)
     1                         / 6.200D+12, 0.000D+00, 1.9383D+04/
      DATA AKF( 269),ALF( 269),EXAF( 269)
     1                         / 2.000D+13, 0.000D+00, 4.4933D+04/
      DATA AKF( 271),ALF( 271),EXAF( 271)
     1                         /-1.900D+16, 0.000D+00, 9.1730D+04/
           DATA AKLOW( 271),ALLOW( 271),EALOW( 271)                             
     1     /  2.95D+44,    -7.35,  95460./                                       
           DATA TROE( 271,1),TROE( 271,2),TROE( 271,3),TROE( 271,4)             
     1     /     0.414, 279.0,    5459.0, 0.0 /                                    
      DATA AKF( 273),ALF( 273),EXAF( 273)
     1                         / 1.047D+37,-7.189D+00, 4.4340D+04/
      DATA AKF( 275),ALF( 275),EXAF( 275)
     1                         / 5.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 277),ALF( 277),EXAF( 277)
     1                         / 2.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 279),ALF( 279),EXAF( 279)
     1                         / 1.000D+14, 0.000D+00, 0.0000D+00/
      DATA AKF( 281),ALF( 281),EXAF( 281)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 283),ALF( 283),EXAF( 283)
     1                         / 3.000D+10, 0.000D+00, 0.0000D+00/
      DATA AKF( 285),ALF( 285),EXAF( 285)
     1                         / 4.900D+14,-5.000D-01, 0.0000D+00/
      DATA AKF( 287),ALF( 287),EXAF( 287)
     1                         / 7.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 289),ALF( 289),EXAF( 289)
     1                         / 3.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 291),ALF( 291),EXAF( 291)
     1                         / 1.170D+43,-9.800D+00, 4.3800D+04/
      DATA AKF( 293),ALF( 293),EXAF( 293)
     1                         / 7.000D+15, 0.000D+00, 8.1700D+04/
      DATA AKF( 295),ALF( 295),EXAF( 295)
     1                         /-3.000D+12, 0.000D+00, 1.6700D+04/
           DATA AKLOW( 295),ALLOW( 295),EALOW( 295)                             
     1     /  1.2D+15,         0.0,    1.25D+04/                                
           DATA TROE( 295,1),TROE( 295,2),TROE( 295,3),TROE( 295,4)             
     1     / 0.0, 0.0, 0.0, 0.0/                                                
      DATA AKF( 297),ALF( 297),EXAF( 297)
     1                         / 3.370D+12, 0.000D+00,-6.2000D+02/
      DATA AKF( 299),ALF( 299),EXAF( 299)
     1                         / 3.370D+11, 0.000D+00,-6.2000D+02/
      DATA AKF( 301),ALF( 301),EXAF( 301)
     1                         / 1.770D+18,-1.900D+00, 2.9800D+03/
      DATA AKF( 303),ALF( 303),EXAF( 303)
     1                         / 3.720D+13,-2.000D-01, 3.5600D+03/
      DATA AKF( 305),ALF( 305),EXAF( 305)
     1                         / 4.660D+13,-3.000D-01, 2.9900D+03/
      DATA AKF( 307),ALF( 307),EXAF( 307)
     1                         / 1.850D+12, 4.000D-01, 5.3600D+03/
      DATA AKF( 309),ALF( 309),EXAF( 309)
     1                         / 3.900D-07, 5.800D+00, 2.2000D+03/
      DATA AKF( 311),ALF( 311),EXAF( 311)
     1                         / 2.450D+01, 3.100D+00, 5.7300D+03/
      DATA AKF( 313),ALF( 313),EXAF( 313)
     1                         / 3.600D+19,-2.200D+00, 1.4000D+04/
      DATA AKF( 315),ALF( 315),EXAF( 315)
     1                         / 2.320D+11, 4.000D-01, 1.4900D+04/
      DATA AKF( 317),ALF( 317),EXAF( 317)
     1                         / 1.000D+14, 0.000D+00, 4.2200D+04/
      DATA AKF( 319),ALF( 319),EXAF( 319)
     1                         /-5.000D+15, 0.000D+00, 8.2000D+04/
           DATA AKLOW( 319),ALLOW( 319),EALOW( 319)                             
     1     /3.0D+16,   0.0,   58000.0/                                             
           DATA TROE( 319,1),TROE( 319,2),TROE( 319,3),TROE( 319,4)             
     1     /       0.5,   1.0D-30,      1.0D+30, 0.0 /                              
      DATA AKF( 321),ALF( 321),EXAF( 321)
     1                         /-8.000D+13, 0.000D+00, 6.5000D+04/
           DATA AKLOW( 321),ALLOW( 321),EALOW( 321)                             
     1     /1.0D+17,   0.0,  54000.0/                                            
           DATA TROE( 321,1),TROE( 321,2),TROE( 321,3),TROE( 321,4)             
     1     /       0.5,   1.0D-30,      1.0D+30, 0.0 /                              
      DATA AKF( 323),ALF( 323),EXAF( 323)
     1                         / 1.810D+11, 4.000D-01, 7.1700D+02/
      DATA AKF( 325),ALF( 325),EXAF( 325)
     1                         / 3.090D+10, 5.000D-01,-3.8000D+02/
      DATA AKF( 327),ALF( 327),EXAF( 327)
     1                         / 1.050D+10, 8.000D-01, 7.1700D+02/
      DATA AKF( 329),ALF( 329),EXAF( 329)
     1                         / 1.900D+07, 1.800D+00, 5.1000D+03/
      DATA AKF( 331),ALF( 331),EXAF( 331)
     1                         / 2.580D+07, 1.600D+00, 2.8300D+03/
      DATA AKF( 333),ALF( 333),EXAF( 333)
     1                         / 1.500D+07, 1.600D+00, 3.0400D+03/
      DATA AKF( 335),ALF( 335),EXAF( 335)
     1                         / 9.410D+07, 1.700D+00, 5.4600D+03/
      DATA AKF( 337),ALF( 337),EXAF( 337)
     1                         / 1.880D+07, 1.900D+00, 1.8200D+03/
      DATA AKF( 339),ALF( 339),EXAF( 339)
     1                         / 1.580D+07, 2.000D+00, 4.4500D+03/
      DATA AKF( 341),ALF( 341),EXAF( 341)
     1                         / 2.190D+02, 3.200D+00, 9.6200D+03/
      DATA AKF( 343),ALF( 343),EXAF( 343)
     1                         / 7.280D+02, 3.000D+00, 7.9500D+03/
      DATA AKF( 345),ALF( 345),EXAF( 345)
     1                         / 1.450D+02, 3.000D+00, 7.6500D+03/
      DATA AKF( 347),ALF( 347),EXAF( 347)
     1                         / 8.200D+03, 2.500D+00, 1.0800D+04/
      DATA AKF( 349),ALF( 349),EXAF( 349)
     1                         / 2.430D+04, 2.500D+00, 1.5800D+04/
      DATA AKF( 351),ALF( 351),EXAF( 351)
     1                         / 3.800D+12, 0.000D+00, 2.4000D+04/
      DATA AKF( 353),ALF( 353),EXAF( 353)
     1                         / 2.410D+11, 0.000D+00,-2.3800D+03/
      DATA AKF( 355),ALF( 355),EXAF( 355)
     1                         / 4.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 357),ALF( 357),EXAF( 357)
     1                         / 5.600D+34,-5.900D+00, 2.5300D+04/
      DATA AKF( 359),ALF( 359),EXAF( 359)
     1                         / 5.350D+37,-7.000D+00, 2.3800D+04/
      DATA AKF( 361),ALF( 361),EXAF( 361)
     1                         / 4.000D+10, 0.000D+00, 1.1000D+03/
      DATA AKF( 363),ALF( 363),EXAF( 363)
     1                         / 4.680D+02, 3.200D+00, 5.3800D+03/
      DATA AKF( 365),ALF( 365),EXAF( 365)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 367),ALF( 367),EXAF( 367)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 369),ALF( 369),EXAF( 369)
     1                         / 1.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 371),ALF( 371),EXAF( 371)
     1                         / 4.820D+13, 0.000D+00, 5.0200D+03/
      DATA AKF( 373),ALF( 373),EXAF( 373)
     1                         / 1.000D+14, 0.000D+00, 0.0000D+00/
      DATA AKF( 375),ALF( 375),EXAF( 375)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 377),ALF( 377),EXAF( 377)
     1                         / 3.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 379),ALF( 379),EXAF( 379)
     1                         / 4.000D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 381),ALF( 381),EXAF( 381)
     1                         / 5.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 383),ALF( 383),EXAF( 383)
     1                         / 1.000D+14, 0.000D+00, 2.5000D+04/
      DATA AKF( 385),ALF( 385),EXAF( 385)
     1                         / 2.000D+07, 1.800D+00, 1.0000D+03/
      DATA AKF( 387),ALF( 387),EXAF( 387)
     1                         / 2.560D+09, 1.100D+00, 1.3644D+04/
      DATA AKF( 389),ALF( 389),EXAF( 389)
     1                         / 7.300D+12, 0.000D+00, 2.2500D+03/
      DATA AKF( 391),ALF( 391),EXAF( 391)
     1                         /-3.000D+13, 0.000D+00, 0.0000D+00/
           DATA AKLOW( 391),ALLOW( 391),EALOW( 391)                             
     1     /  9.000D+15,   1.000,     0.00 /                                    
           DATA TROE( 391,1),TROE( 391,2),TROE( 391,3),TROE( 391,4)             
     1     /     0.5, 1.0D+30,      1.0D-30, 0.0/                                   
      DATA AKF( 393),ALF( 393),EXAF( 393)
     1                         / 2.500D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 395),ALF( 395),EXAF( 395)
     1                         / 5.300D+06, 2.000D+00, 2.0000D+03/
      DATA AKF( 397),ALF( 397),EXAF( 397)
     1                         / 3.000D+10, 0.000D+00, 2.8680D+03/
      DATA AKF( 399),ALF( 399),EXAF( 399)
     1                         /-4.000D+13, 0.000D+00, 0.0000D+00/
           DATA AKLOW( 399),ALLOW( 399),EALOW( 399)                             
     1     /  3.000D+24,  -2.000,     0.00 /                                    
           DATA TROE( 399,1),TROE( 399,2),TROE( 399,3),TROE( 399,4)             
     1     /     0.8, 1.0D+30,      1.0D-30,  0.00/                                   
      DATA AKF( 401),ALF( 401),EXAF( 401)
     1                         / 1.800D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 403),ALF( 403),EXAF( 403)
     1                         / 4.990D+15,-1.400D+00, 2.2428D+04/
      DATA AKF( 405),ALF( 405),EXAF( 405)
     1                         / 3.000D+12,-3.200D-01,-1.3090D+02/
      DATA AKF( 407),ALF( 407),EXAF( 407)
     1                         /-6.000D+08, 0.000D+00, 0.0000D+00/
           DATA AKLOW( 407),ALLOW( 407),EALOW( 407)                             
     1     /  2.000D+09,   1.000,     0.00 /                                    
           DATA TROE( 407,1),TROE( 407,2),TROE( 407,3),TROE( 407,4)             
     1     /  0.5, 1.0D+30,      1.0D-30, 0.0 /                                   
      DATA AKF( 409),ALF( 409),EXAF( 409)
     1                         / 6.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 411),ALF( 411),EXAF( 411)
     1                         / 2.500D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 413),ALF( 413),EXAF( 413)
     1                         / 8.000D+11, 0.000D+00, 0.0000D+00/
      DATA AKF( 415),ALF( 415),EXAF( 415)
     1                         / 4.000D+14, 0.000D+00, 4.1826D+04/
      DATA AKF( 417),ALF( 417),EXAF( 417)
     1                         / 3.500D+07, 1.650D+00,-9.7270D+02/
      DATA AKF( 419),ALF( 419),EXAF( 419)
     1                         / 3.100D+06, 2.000D+00,-2.9820D+02/
      DATA AKF( 421),ALF( 421),EXAF( 421)
     1                         / 1.200D+08, 1.650D+00, 3.2740D+02/
      DATA AKF( 423),ALF( 423),EXAF( 423)
     1                         / 1.700D+05, 2.500D+00, 2.4928D+03/
      DATA AKF( 425),ALF( 425),EXAF( 425)
     1                         /-2.000D+14, 0.000D+00, 0.0000D+00/
           DATA AKLOW( 425),ALLOW( 425),EALOW( 425)                             
     1     /  1.330D+60, -12.000,  5967.97 /                                    
           DATA TROE( 425,1),TROE( 425,2),TROE( 425,3),TROE( 425,4)             
     1     /    0.02,    1097.0,      1097.0,      6860.0 /                             
      DATA AKF( 427),ALF( 427),EXAF( 427)
     1                         / 2.660D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 429),ALF( 429),EXAF( 429)
     1                         / 3.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 431),ALF( 431),EXAF( 431)
     1                         /-2.500D+13, 0.000D+00, 0.0000D+00/
           DATA AKLOW( 431),ALLOW( 431),EALOW( 431)                             
     1     /  4.270D+58, -11.940,  9770.55 /                                    
           DATA TROE( 431,1),TROE( 431,2),TROE( 431,3),TROE( 431,4)             
     1     /   0.175,    1341.0,     6.0D+04, 1.014D+04 /                             
      DATA AKF( 433),ALF( 433),EXAF( 433)
     1                         / 1.600D+22,-2.390D+00, 1.1185D+04/
      DATA AKF( 435),ALF( 435),EXAF( 435)
     1                         / 1.500D+24,-2.830D+00, 1.8619D+04/
      DATA AKF( 437),ALF( 437),EXAF( 437)
     1                         /-1.100D+17, 0.000D+00, 8.4393D+04/
           DATA AKLOW( 437),ALLOW( 437),EALOW( 437)                             
     1     /  7.830D+18,   0.000, 64978.49 /                                    
           DATA TROE( 437,1),TROE( 437,2),TROE( 437,3),TROE( 437,4)             
     1     /    0.76, 1946.00,    38.00,  0.0/                                
      DATA AKF( 439),ALF( 439),EXAF( 439)
     1                         / 4.000D+13, 0.000D+00, 4.7500D+04/
      DATA AKF( 441),ALF( 441),EXAF( 441)
     1                         / 4.000D+13, 0.000D+00, 5.0932D+04/
      DATA AKF( 443),ALF( 443),EXAF( 443)
     1                         / 1.300D+06, 2.400D+00, 4.4710D+03/
      DATA AKF( 445),ALF( 445),EXAF( 445)
     1                         / 1.330D+06, 2.540D+00, 6.7614D+03/
      DATA AKF( 447),ALF( 447),EXAF( 447)
     1                         / 4.760D+04, 2.710D+00, 2.1073D+03/
      DATA AKF( 449),ALF( 449),EXAF( 449)
     1                         / 1.900D+05, 2.680D+00, 3.7184D+03/
      DATA AKF( 451),ALF( 451),EXAF( 451)
     1                         / 1.400D+03, 2.660D+00, 5.2720D+02/
      DATA AKF( 453),ALF( 453),EXAF( 453)
     1                         / 2.700D+04, 2.390D+00, 3.9310D+02/
      DATA AKF( 455),ALF( 455),EXAF( 455)
     1                         / 9.640D+03, 2.600D+00, 1.3917D+04/
      DATA AKF( 457),ALF( 457),EXAF( 457)
     1                         / 4.760D+04, 2.550D+00, 1.6491D+04/
      DATA AKF( 459),ALF( 459),EXAF( 459)
     1                         / 8.400D-03, 4.200D+00, 8.6759D+03/
      DATA AKF( 461),ALF( 461),EXAF( 461)
     1                         /-1.330D+13, 0.000D+00, 1.5607D+03/
           DATA AKLOW( 461),ALLOW( 461),EALOW( 461)                             
     1     /  8.700D+42,  -7.500,  4732.31 /                                    
           DATA TROE( 461,1),TROE( 461,2),TROE( 461,3),TROE( 461,4)             
     1     /       1.0,    1000.0,    645.4,     6844. /                          
      DATA AKF( 463),ALF( 463),EXAF( 463)
     1                         / 1.300D+11, 0.000D+00, 0.0000D+00/
      DATA AKF( 465),ALF( 465),EXAF( 465)
     1                         /-1.230D+13,-1.000D-01, 3.0210D+04/
           DATA AKLOW( 465),ALLOW( 465),EALOW( 465)                             
     1     /  5.490D+49, -10.000, 35779.16 /                                    
           DATA TROE( 465,1),TROE( 465,2),TROE( 465,3),TROE( 465,4)             
     1     /   -1.17,     251.0,    1.0D-30,      1185. /                           
      DATA AKF( 467),ALF( 467),EXAF( 467)
     1                         /-1.330D+13, 0.000D+00, 3.2600D+03/
           DATA AKLOW( 467),ALLOW( 467),EALOW( 467)                             
     1     /  6.260D+38,  -6.660,  7000.48 /                                    
           DATA TROE( 467,1),TROE( 467,2),TROE( 467,3),TROE( 467,4)             
     1     /       1.0,    1000.0,     1310.0, 4.81D+04 /                          
      DATA AKF( 469),ALF( 469),EXAF( 469)
     1                         / 9.000D+10, 0.000D+00, 0.0000D+00/
      DATA AKF( 471),ALF( 471),EXAF( 471)
     1                         / 1.580D+16, 0.000D+00, 1.1000D+05/
      DATA AKF( 473),ALF( 473),EXAF( 473)
     1                         / 1.800D+13, 0.000D+00, 8.5127D+04/
      DATA AKF( 475),ALF( 475),EXAF( 475)
     1                         / 1.260D+13, 0.000D+00, 0.0000D+00/
      DATA AKF( 477),ALF( 477),EXAF( 477)
     1                         / 5.000D+11, 0.000D+00, 0.0000D+00/
      DATA AKF( 479),ALF( 479),EXAF( 479)
     1                         / 6.300D+10, 7.000D-01, 6.0014D+03/
      DATA AKF( 481),ALF( 481),EXAF( 481)
     1                         / 5.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 483),ALF( 483),EXAF( 483)
     1                         / 7.000D+13, 0.000D+00, 1.8413D+04/
      DATA AKF( 485),ALF( 485),EXAF( 485)
     1                         / 5.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 487),ALF( 487),EXAF( 487)
     1                         / 1.000D+16, 0.000D+00, 7.2933D+04/
      DATA AKF( 489),ALF( 489),EXAF( 489)
     1                         / 3.160D+12, 0.000D+00, 5.7075D+04/
      DATA AKF( 491),ALF( 491),EXAF( 491)
     1                         / 3.160D+12, 0.000D+00, 5.7075D+04/
      DATA AKF( 493),ALF( 493),EXAF( 493)
     1                         / 3.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 495),ALF( 495),EXAF( 495)
     1                         / 3.000D+12, 0.000D+00, 0.0000D+00/
      DATA AKF( 497),ALF( 497),EXAF( 497)
     1                         / 1.000D+14, 0.000D+00, 0.0000D+00/
      DATA AKF( 499),ALF( 499),EXAF( 499)
     1                         / 1.000D+14, 0.000D+00, 0.0000D+00/
      DATA AKF( 501),ALF( 501),EXAF( 501)
     1                         / 1.000D+18, 0.000D+00, 8.5325D+04/
      DATA AKF( 503),ALF( 503),EXAF( 503)
     1                         / 2.600D+06, 2.400D+00, 4.4694D+03/
      DATA AKF( 505),ALF( 505),EXAF( 505)
     1                         / 2.080D+06, 2.400D+00, 4.4694D+03/
      DATA AKF( 507),ALF( 507),EXAF( 507)
     1                         / 1.320D+06, 2.540D+00, 6.7638D+03/
      DATA AKF( 509),ALF( 509),EXAF( 509)
     1                         / 1.300D+06, 2.400D+00, 4.4694D+03/
      DATA AKF( 511),ALF( 511),EXAF( 511)
     1                         / 5.460D+06, 2.000D+00,-1.3145D+03/
      DATA AKF( 513),ALF( 513),EXAF( 513)
     1                         / 4.380D+06, 2.000D+00,-1.3145D+03/
      DATA AKF( 515),ALF( 515),EXAF( 515)
     1                         / 4.750D+06, 2.000D+00,-5.9750D+02/
      DATA AKF( 517),ALF( 517),EXAF( 517)
     1                         / 2.180D+07, 1.800D+00, 9.7990D+02/
      DATA AKF( 519),ALF( 519),EXAF( 519)
     1                         / 1.000D+16, 0.000D+00, 7.2897D+04/
      DATA AKF( 521),ALF( 521),EXAF( 521)
     1                         / 3.160D+12, 0.000D+00, 5.7122D+04/
      DATA AKF( 523),ALF( 523),EXAF( 523)
     1                         / 7.080D+07, 1.900D+00, 1.6730D+02/
      DATA AKF( 525),ALF( 525),EXAF( 525)
     1                         / 1.300D+06, 2.400D+00, 4.4694D+03/
      DATA AKF( 527),ALF( 527),EXAF( 527)
     1                         / 7.230D+12, 0.000D+00, 1.2906D+03/
      DATA AKF( 529),ALF( 529),EXAF( 529)
     1                         / 7.230D+12, 0.000D+00, 1.2906D+03/
      DATA AKF( 531),ALF( 531),EXAF( 531)
     1                         / 6.600D+05, 2.540D+00, 6.7638D+03/
      DATA AKF( 533),ALF( 533),EXAF( 533)
     1                         / 1.150D+05, 2.500D+00, 2.4856D+03/
      DATA AKF( 535),ALF( 535),EXAF( 535)
     1                         / 1.000D+16, 0.000D+00, 7.2897D+04/
      DATA AKF( 537),ALF( 537),EXAF( 537)
     1                         / 7.230D+12, 0.000D+00, 1.2906D+03/
      DATA AKF( 539),ALF( 539),EXAF( 539)
     1                         / 7.230D+12, 0.000D+00, 1.2906D+03/
      DATA AKF( 541),ALF( 541),EXAF( 541)
     1                         / 6.600D+05, 2.540D+00, 6.7638D+03/
      DATA AKF( 543),ALF( 543),EXAF( 543)
     1                         / 2.080D+06, 2.000D+00,-2.8680D+02/
C-----------------------------------------------------------------------     
      END 

