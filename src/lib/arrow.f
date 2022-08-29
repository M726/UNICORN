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