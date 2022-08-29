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