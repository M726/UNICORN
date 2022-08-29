      SUBROUTINE PRNTFF(INSP,ITRANS)
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)                                   
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LE)                  
      COMMON/CB02/ RREF,VREF,TREF,PREF,WM(LSP)
      COMMON/CB03/ DT,X(LJ,LI),Y(LJ,LI),XXC(LI),XXS(LI),YYC(LJ),YYS(LJ)           
      COMMON/CB04/ RHO(LJ,LI),FSP(LJ,LI,LSP),
     *  U(LJ,LI),V(LJ,LI),W(LJ,LI),P(LJ,LI),HT(LJ,LI),TK(LJ,LI),
     *  AK(LJ,LI),EPS(LJ,LI)
      COMMON/CB06/ ALSTR,VSTR,RSTR,PSTR,TSTR,TMSTR,AMSTR,ACSTR,DSTR,  
     1  HSTR,CPSTR,AKSTR,EPSTR,WMSTR,GASC,BETA1,BETA2,BETA3,BETA4,REI         
      COMMON/CB07/ VMU(LJ,LI),VTC(LJ,LI),VDFSP(LJ,LI,LSP),TMU(LJ,LI)
      COMMON/CB08/ TPOL1,TPOL2,POLSP(14,LSP),CISP(12,LSP)
      COMMON/HRLS/ QDOT(LJ,LI)
      COMMON/SOOT/ RSOOT,STMF(LJ,LI),STND(LJ,LI),STDF(LJ,LI)
      DIMENSION IP(16),A(16)
      CHARACTER CAPT1(LSP+16)*18                                             
      DATA CAPT1 
     1 /'     PRESSURE     ','      DENSITY     ','   TEMPERATURE    ',    
     2  '   U VELOCITY     ','    V VELOCITY    ','  SWIRL VELOCITY  ',
     3  '  KINETIC ENERGY  ','  KE DISSIPATION  ',' N2 MASS FRACTION ',
     *  'SOOT MASS FRACTION','SOOT NUMBE DENSITY',
     *  'LAMINAR VISCOSITY ','TURBULT VISCOSITY ','HEAT CONDUCTIVITY ',
     *  'CH4-DIFFUSION CO. ','O2-DIFFUSION CO.  ','   HRR (J/CC/S)   ',
     *  'H2                ','O2                ','H                 ',
     *  'O                 ','OH                ','H2O               ',
     *  'HO2               ','H2O2              ','CO                ',
     *  'CO2               ','HCO               ','CH2O              ',
     *  'CH4               ','CH3               ','T-CH2             ',
     *  'S-CH2             ','C2H4              ','CH3O              ',
     *  'C2H5              ','C2H6              ','CH                ',
     *  'C2H2              ','C2H3              ','CH2CHO            ',
     *  'C2H4O             ','CH2CO             ','HCCO              ',
     *  'C2H               ','CH2OH             ','CH3OH             ',
     *  'C2H5OH            ','CH3CHO            ','CH3CHOH           ',
     *  'CH2CH2OH          ','CH3CO             ','CH3CH2O           ',
     *  'C3H4              ','C3H3              ','C3H5              ',
     *  'C3H6              ','C3H8              ','I-C3H7            ',
     *  'N-C3H7            ','C4H6              ','CHO               ',
     *  'C5H8              ','C7H16             ','C4H8              ',
     *  'C5H10             ','AR                ','HE                '/
      KA=LI
      KB=LJ
      IF(ITRANS.EQ.1) THEN
                      KA=LJ
                      KB=LI
                      END IF
      IVAL=8
      IP(1)=1                                                           
      DO 10 I=2,IVAL-1                                                      
      IP(I)=INT(DFLOAT(KA*(I-1))/DFLOAT(IVAL-1))                             
      IF(IP(I).EQ.IP(I-1)) IP(I)=IP(I)+1
      IF(IP(I).GT.LI) IP(I)=KA
   10 CONTINUE
      IP(IVAL)=KA
      WRITE(16,810)                                                      
  810 FORMAT(/20X,'GEOMETRICAL DATA IN millimeters'/)                          
      WRITE(16,812)                                                      
  812 FORMAT(2X,'NO.',4X,'X-Dis.',5X,'DX',7X,'R-outer ',     
     1       5X,'NO.',3X,'Y-Dis.',5X,'DY',7X,'Velocity')     
      AA=1000.0*ALSTR
      IJ=LI
      IF(LJ.GT.LI) IJ=LJ
      DO 12 I=1,IJ
      DXX=0.0
      DYY=0.0
      IF((I.LE.LI).AND.(I.LE.LJ)) THEN
          IF(I.GE.2) DXX=(X(1,I)-X(1,I-1))*AA
          IF(I.GE.2) DYY=(Y(I,1)-Y(I-1,1))*AA
          WRITE(16,814) I,X(1,I)*AA,DXX,
     1    Y(LJ,I)*AA,I,Y(I,1)*AA,DYY,U(I,1)*VSTR
          END IF
      IF((I.LE.LI).AND.(I.GT.LJ)) THEN
          DXX=(X(1,I)-X(1,I-1))*AA
          WRITE(16,814) I,X(1,I)*AA,DXX,Y(LJ,I)*AA
          END IF
      IF((I.GT.LI).AND.(I.LE.LJ)) THEN
          DYY=(Y(I,1)-Y(I-1,1))*AA
          WRITE(16,815) I,Y(I,1)*AA,DYY,U(I,1)*VSTR
          END IF
  814 FORMAT(2X,I3,F10.3,1X,F10.7,F10.5,5X,I3,F10.5,1X,F10.7,1X,F10.6)                    
  815 FORMAT(42X,I3,F10.5,2X,F10.7,1X,F10.6)                    
   12 CONTINUE                                                          
      DO 20 II=1,LSP+16                                                     
      IF(II.LE.17) THEN
                   WRITE(16,818) ITR,CAPT1(II)
                   ELSE
                   WRITE(16,819) ITR,II-17,CAPT1(II)
                   END IF
  818 FORMAT(/1X,14('@'),'  ITR=',I7,7X,A18,10X,14('@')/)                            
  819 FORMAT(/1X,14('@'),'  ITR=',I7,3X,I3,'-',A18,10X,14('@')/)                            
      IF(ITRANS.EQ.0) WRITE(16,830) (IP(I),I=1,IVAL)                                       
      IF(ITRANS.EQ.1) WRITE(16,831) (IP(I),I=1,IVAL)                                       
  830 FORMAT(1X,'(J)   ',16('I=',I3,4X))                         
  831 FORMAT(1X,'(I)   ',16('J=',I3,4X))                         
      DO 22 K=1,KB                                                      
      DO 24 IA=1,IVAL
      IB1=IP(IA)                                                        
      J=K*(1-ITRANS)+IB1*ITRANS
      I=K*ITRANS+IB1*(1-ITRANS)
      IF(II.EQ.1) A(IA)=P(J,I)*PSTR                       
      IF(II.EQ.2) A(IA)=RHO(J,I)*RSTR                                 
      IF(II.EQ.3) A(IA)=TK(J,I)*TSTR                       
      IF(II.EQ.4) A(IA)=U(J,I)*VSTR                       
      IF(II.EQ.5) A(IA)=V(J,I)*VSTR                          
      IF(II.EQ.6) A(IA)=W(J,I)*VSTR                          
      IF(II.EQ.7) A(IA)=AK(J,I)*AKSTR                          
      IF(II.EQ.8) A(IA)=EPS(J,I)*EPSTR                          
      IF(II.EQ.9) A(IA)=FSP(J,I,LSP)
      IF(II.EQ.10) A(IA)=STMF(J,I)                     
      IF(II.EQ.11) A(IA)=STND(J,I)                       
      IF(II.EQ.12) A(IA)=VMU(J,I)*AMSTR                       
      IF(II.EQ.13) A(IA)=TMU(J,I)*AMSTR                        
      IF(II.EQ.14) A(IA)=VTC(J,I)*ACSTR                       
      IF(II.EQ.15) A(IA)=VDFSP(J,I,INSP)*DSTR                   
      IF(II.EQ.16) A(IA)=VDFSP(J,I,2)*DSTR                   
      IF(II.EQ.17) A(IA)=QDOT(J,I)                   
      IF(II.GE.18) A(IA)=FSP(J,I,II-17)                     
   24 CONTINUE                                                          
      IF(II.EQ.1) WRITE(16,824) K,(A(IA),IA=1,IVAL)                                
      IF(II.EQ.2) WRITE(16,821) K,(A(IA),IA=1,IVAL)                                
      IF(II.EQ.3) WRITE(16,822) K,(A(IA),IA=1,IVAL)                                
      IF((II.GE.4).AND.(II.LE.6)) WRITE(16,823) K,(A(IA),IA=1,IVAL)                                
      IF((II.GE.7).AND.(II.LE.8)) WRITE(16,826) K,(A(I),IA=1,IVAL)                                
      IF(II.EQ.9) WRITE(16,825) K,(A(IA),IA=1,IVAL)                                
      IF(II.GE.10) WRITE(16,826) K,(A(IA),IA=1,IVAL)                                
  821 FORMAT(1X,I3,16F9.5)                                      
  822 FORMAT(1X,I3,16F9.2)                                      
  823 FORMAT(1X,I3,16F9.5)                                      
  824 FORMAT(1X,I3,1P16D9.2)                                      
  825 FORMAT(1X,I3,16F9.5)                                      
  826 FORMAT(1X,I3,1P16D9.2)                                      
   22 CONTINUE                                                          
   20 CONTINUE                                                          
      RETURN                                                            
      END