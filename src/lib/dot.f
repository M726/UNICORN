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