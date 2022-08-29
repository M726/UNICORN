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
      END PROGRAM UNICORN

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
      INCLUDE 'lib\ptrelx.f'
      INCLUDE 'lib\readff.f'
      INCLUDE 'lib\reset.f'
      INCLUDE 'lib\skipdf.f'
      INCLUDE 'lib\solveh.f'
      INCLUDE 'lib\solvesp.f'
      INCLUDE 'lib\spgsolv.f'
      INCLUDE 'lib\stdata.f'
      INCLUDE 'lib\swlsolv.f'


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