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