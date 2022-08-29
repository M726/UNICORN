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