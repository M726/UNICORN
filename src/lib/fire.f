      SUBROUTINE FIRE(INSP,IGKEEP)
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
      IFIRE=0
      DO 12 I=2,LI-1
      DO 12 J=2,LJ-1
      IF(FSP(J,I,INSP).GE.0.6.AND.FSP(J,I,INSP).LE.0.95) THEN
             IFIRE=IFIRE+1
             ISKIP(J,I)=-9
             END IF
   12 CONTINUE
      IF(IFIRE.EQ.0) WRITE(*,812) 
  812 FORMAT(5X,'FIRE IS NOT INTRODUCED INTO THE FLOW')
      IF(IFIRE.GE.1) WRITE(*,814) IFIRE,ITR,ITR+IGKEEP
  814 FORMAT(5X,50('&'),
     1  /5X,'&&',46X,'&&',
     2  /5X,'&&',5X,'FIRE IS INITIATED AT ',I4,' GRID POINTS',4X,'&&',
     3  /5X,'&&',8X,'BETWEEN ITERATIONS ',I3,' AND ',I3,8X,'&&',
     4  /5X,'&&',46X,'&&'/5X,50('&')/)
      RETURN
      END