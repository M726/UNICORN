      SUBROUTINE DOTPR(X1,Y1,IDOT)
            
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IPMAX=500,JPMAX=200)
      COMMON/APIXL/ IPAGE,NCOLOR,NBACK,NAA(IPMAX,JPMAX)
      COMMON/FRAME/ IPIXL,JPIXL,XMIN,YMIN,XMAX,YMAX,XPIXL,YPIXL
      
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