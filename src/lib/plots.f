      
      INCLUDE 'lib\arrow.f'
      INCLUDE 'lib\dot.f'
      INCLUDE 'lib\dotpr.f'
      INCLUDE 'lib\drawby.f'
      INCLUDE 'lib\header.f'
      INCLUDE 'lib\line.f'

      SUBROUTINE PLOTS(NOPR,NOINJ,IPANM,
     1                 X1ANM,Y1ANM,X2ANM,Y2ANM,NFSURF,KORNT)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LSP=52,
     1   IPMAX=500,JPMAX=200,LNPR=200,LNMX=2000)
      CHARACTER *2 BFL2A(10)
      CHARACTER *5 BFL2*4
      CHARACTER*14 FRMTSTR
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),U(LJ,LI),V(LJ,LI),W(LJ,LI),
     1             P(LJ,LI),HT(LJ,LI),TK(LJ,LI),AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1 HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
      COMMON/BODY/IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10),NBODY
      COMMON/TRACK/ IPICK(101),JPICK(101),XPR(LNPR),YPR(LNPR),
     1              PATHX(LNMX,LNPR),PATHY(LNMX,LNPR),
     2              PUVEL(LNMX,LNPR),PVVEL(LNMX,LNPR),
     3              GUVEL(LNMX,LNPR),GVVEL(LNMX,LNPR)
      COMMON/FLOW/ RASMN,RASMX,RAS2MN,RAS2MX,XS(LI),YS(LJ+LJ-1),
     1             KSYM,LIP,LJP,IOFF,JOFF,NCLR(LNPR)
      COMMON/APIXL/ IPAGE,NCOLOR,NBACK,NAA(IPMAX,JPMAX)
      COMMON/FRAME/ IPIXL,JPIXL,XMIN,YMIN,XMAX,YMAX,XPIXL,YPIXL
      COMMON/PEN/NPEN,NAPEN,NBPEN
      DATA BFL2A / '-0','-1','-2','-3','-4','-5','-6','-7','-8','-9'/
      CHARACTER *1 BB(IPMAX,JPMAX)
      CHARACTER *25 FILENAME
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
      DO 37 J=1,LJP-1
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
   37 CONTINUE
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
      DO 591 J=1,JPIXL
      NAA(I,J)=0
  591 CONTINUE
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
      DO 599 J=1,LJP
      IA=I+IOFF
      JA=J+JOFF
      FS(I,J)=0.25*(TK(JA,IA)+TK(JA+1,IA)+TK(JA,IA+1)+TK(JA+1,IA+1))
C     FS2(I,J)=0.25*(RHO(JA,IA)*STMF(JA,IA)+RHO(JA+1,IA)*STMF(JA+1,IA)
C    1 +RHO(JA,IA+1)*STMF(JA,IA+1)+RHO(JA+1,IA+1)*STMF(JA+1,IA+1))
      FS2(I,J)=0.25*(FSP(JA,IA,05)+FSP(JA+1,IA,05)
     1        +FSP(JA,IA+1,05)+FSP(JA+1,IA+1,05))
  599 CONTINUE
  600 CONTINUE
      DO 601 I=LIP,2,-1
      DO 603 J=LJP,2,-1
      FS(I,J)=FS(I-1,J-1)
      FS2(I,J)=FS2(I-1,J-1)
  603 CONTINUE
  601 CONTINUE
      IF(KSYM.EQ.1) THEN
                    DO 607 I=1,LIP
                    DO 602 J=LJ2,LJP,-1
                    JA=J-LJP+1
                    FS(I,J)=FS(I,JA)
  602               CONTINUE
  607               CONTINUE
C------------------------ADD THESE LINES FOR FUEL CONCENTRATION PLOT----
                    DO 604 I=1,LIP
                    DO 605 J=1,LJP-1
                    JA=LJP-J+1
                    FSNEW=1.0*(FS2(I,JA)-RAS2MN)/RASCON
                    IF(FSNEW.GT.RASLMT) FSNEW=RASLMT
                    FS(I,J)=RASMN+FSNEW
  605               CONTINUE
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
      DO 611 I=1,LIP
      DO 610 J=1,LJJ
      IF(FS(I,J).GT.RASMX) FS(I,J)=RASMX
      IF(FS(I,J).LT.RASMN) FS(I,J)=RASMN
  610 CONTINUE
  611 CONTINUE
      IBACK=999999
      IBACKM=IBACK-10
      BACK=DFLOAT(IBACKM)/(RASMX-RASMN)
      DO 621 I=1,LIP
      DO 620 J=1,LJJ
      FS(I,J)=0.0+BACK*(FS(I,J)-RASMN)
  620 CONTINUE
  621 CONTINUE
C-----------------   INITIALIZE 'NAA' WITH A BIG NUMBER   -----
      DO 631 I=1,IPIXL
      DO 630 J=1,JPIXL
      NAA(I,J)=IBACK
  630 CONTINUE
  631 CONTINUE
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
      DO 649 I=1,IPIXL
      DO 650 J=1,JPIXL
      NIJ=NAA(I,J)
      IF(NIJ.GE.IBACKW) GO TO 651
      IF(NIJ.GT.IBMAX) IBMAX=NIJ
      IF(NIJ.LT.IBMIN) IBMIN=NIJ
  651 CONTINUE
  650 CONTINUE
  649 CONTINUE
      IF(IBMAX.LT.IBACKM) IBMAX=IBACKM
      IF(IBMIN.GT.0) IBMIN=0
      IPLOW=1
      IPHIG=253
      FACT=DFLOAT(IPHIG-IPLOW)/DFLOAT(IBMAX-IBMIN)
      DO 652 I=1,IPIXL
      DO 653 J=1,JPIXL
      IF(NAA(I,J).GE.IBACKW) THEN
                            NAA(I,J)=NBACK
                            GO TO 654
                            END IF
      NAA(I,J)=IPLOW+INT(FACT*DFLOAT(NAA(I,J)-IBMIN))
  654 CONTINUE
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
           DO 401 J=1,JPIXL
           IS=NAA(I,J)
           IF(IS.GT.255) IS=255
           IF(IS.LT.0) IS=0
           BB(I,J)=CHAR(IS)
  401      CONTINUE
  400      CONTINUE
           IF(IPAGE.GE.3) IPAGE=IPAGE+1
           IF(IPAGE.NE.2) IPP=IPAGE-3

           WRITE(FILENAME, 10) IPP
   10      FORMAT('movie-',(I6.6),'.bmp')
           MPP=82
           OPEN(MPP,FILE="images\"//FILENAME,STATUS='UNKNOWN')
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