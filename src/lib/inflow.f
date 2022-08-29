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
	  
C------------Initialize every grid point
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
C-----------Initialize all species
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