      SUBROUTINE PRNTDK(MP,INFUEL)                      
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)                                   
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LE,LSP-1),FSN2(LE),
     6  U(LJ,LI),V(LJ,LI),W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),
     7  AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/HRLS/ QDOT(LJ,LI)
      COMMON/BODY/IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10),NBODY
      COMMON/SOOT/ RSOOT,STMF(LE),STND(LE),STDF(LE)
C---add 200000 (for swirl), 400000 (for soot), 100000 (for bodies)
      IFLAME=200000+LSP*100+INFUEL
      IFLAME=IFLAME+400000
      IF(NBODY.GE.1) IFLAME=IFLAME+100000
      INERT=0
      WRITE(MP,100) IFLAME
      WRITE(MP,102) LI,LJ,ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,INERT,ITR,
     1              DT
      WRITE(MP,104) ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      IF(NBODY.GE.1) 
     1 WRITE(MP,103) NBODY,(IBDM(N),JBDM(N),IBDP(N),JBDP(N),N=1,NBODY)
      WRITE(MP,104) X,Y,RHO,FSP,U,V,W,P,QDOT,TK
      WRITE(MP,104) STMF,STND
      IF(IFLOW.EQ.2) WRITE(MP,104) AK,EPS
  100 FORMAT(I6,2('-'),'FLAME DATA FROM UNICOND-HEPT-SD(',
     1    ' 52 SPECIES & 544 REACTIONS OF SANDIEGO MECH) DATA',2('-'))
  102 FORMAT(10I6,F12.10)
  103 FORMAT(21I6)
  104 FORMAT(8(1PE14.7,1X))
      RETURN
      END