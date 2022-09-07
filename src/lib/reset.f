      SUBROUTINE RESET                 
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LIJ=LI+LJ,LSP=52,ID1=7+LE,
     1   ID2=4+LSP,ID3=1+2*LE+2*LI+2*LJ,ID5=(9+LSP)*LE,IDNT=LE,IDN=20,
     2   IDV=(3+LSP)*LE,IDB1=2*LIJ,IDB2=2*(8+LSP)*LIJ,IDF=2*4*LIJ,
     3   LB=LJ-2,LL=(LI-2)*LB,LL1=LB*LL-(LB+1)*LB/2,LPD=55*LE)
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
      COMMON/LUMAT/ PSL(LL1),PSD(LL),IDATLU(LB)
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
      PSL(I)=0.0
  116 CONTINUE
      DO 117 I=1,IDNT
      DATNT(I)=0.0
  117 CONTINUE
      DO 120 I=1,LPD
      DATD(I)=0.0
  120 CONTINUE


C-----New - Maintains correct /lumat/ structure and still clears everything.
      DO 118 I=1, LL
      PSD(I)=0.0
  118 CONTINUE
      RETURN                                                              
      END     