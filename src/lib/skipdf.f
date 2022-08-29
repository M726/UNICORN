      SUBROUTINE SKIPDF
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)                                        
      PARAMETER(LI=711,LJ=131,LE=LI*LJ,LSP=52)                                   
      COMMON/CB01/ ISYM,IFLOW,ISWRL,ITHRM,ICHEM,IGRAV,ITR,ISKIP(LJ,LI)                  
      COMMON/CB09/ IBOT(LI),ITOP(LI),JLFT(LJ),JRGT(LJ),
     1       FBOT(8+LSP,LI),FTOP(8+LSP,LI),FLFT(8+LSP,LJ),FRGT(8+LSP,LJ)
      COMMON/BODY/NBODY,IBDM(10),JBDM(10),IBDP(10),JBDP(10),TBD(10)
      COMMON/FINJ/NFINJ,IFIM(10),JFIM(10),IFIP(10),JFIP(10),FFI(10,30)
C---------------Define Skip variable for flow Calculations--------------
C--------------- ISKIP=0 no skip
C--------------- ISKIP=1  skip epcilon calculations
C--------------- ISKIP=2  skip DT calculations
C--------------- ISKIP=3  Injection locations
C--------------- ISKIP=10 skip every calculations
C--------------- ISKIP=-9 introduce flame in the calculations
      DO 200 I=1,LI
      DO 200 J=1,LJ
      ISKIP(J,I)=0
  200 CONTINUE
      IF(NBODY.GE.1) THEN
          DO 202 N=1,NBODY
          IBODYM=IBDM(N)
          IBODYP=IBDP(N)
          JBODYM=JBDM(N)
          JBODYP=JBDP(N)
          DO 204 I=IBODYM,IBODYP
          IF(ISKIP(JBODYM,I).LT.10) THEN
                                    ISKIP(JBODYM,I)=11
                                    ELSE
                                    ISKIP(JBODYM,I)=10
                                    END IF
          IF(ISKIP(JBODYP,I).LT.10) THEN
                                    ISKIP(JBODYP,I)=11
                                    ELSE
                                    ISKIP(JBODYP,I)=10
                                    END IF
  204     CONTINUE
          DO 206 J=JBODYM,JBODYP
          IF(ISKIP(J,IBODYM).LT.10) THEN
                                    ISKIP(J,IBODYM)=11
                                    ELSE
                                    ISKIP(J,IBODYM)=10
                                    END IF
          IF(ISKIP(J,IBODYP).LT.10) THEN
                                    ISKIP(J,IBODYP)=11
                                    ELSE
                                    ISKIP(J,IBODYP)=10
                                    END IF
  206     CONTINUE
          DO 208 I=IBODYM+1,IBODYP-1
          DO 208 J=JBODYM+1,JBODYP-1
          ISKIP(J,I)=10
  208     CONTINUE
  202     CONTINUE
          END IF
      IF(NFINJ.GT.0) THEN
          DO 210 N=1,NFINJ
          IFINJM=IFIM(N)
          IFINJP=IFIP(N)
          JFINJM=JFIM(N)
          JFINJP=JFIP(N)
          DO 212 I=IFINJM,IFINJP
          DO J=JFINJM,JFINJP
          ISKIP(J,I)=3
          ENDDO
  212     CONTINUE
  210     CONTINUE
          END IF
      DO 222 I=2,LI-1
      ISKIP(1,I)=ISKIP(1,I)+1
      ISKIP(LJ,I)=ISKIP(LJ,I)+1
      IF(IBOT(I).EQ.1) THEN
                       ISKIP(2,I)=ISKIP(1,I)
                       ISKIP(1,I)=ISKIP(1,I)+1
                       END IF
      IF(ITOP(I).EQ.1) THEN
                       ISKIP(LJ-1,I)=ISKIP(LJ,I)
                       ISKIP(LJ,I)=ISKIP(LJ,I)+1
                       END IF
  222 CONTINUE
      DO 224 J=2,LJ-1
      ISKIP(J,1)=ISKIP(J,1)+1
      ISKIP(J,LI)=ISKIP(J,LI)+1
      IF(JLFT(J).EQ.1) THEN
                       ISKIP(J,2)=ISKIP(J,1)
                       ISKIP(J,1)=ISKIP(J,1)+1
                       END IF
      IF(JRGT(J).EQ.1) THEN
                       ISKIP(J,LI-1)=ISKIP(J,LI)
                       ISKIP(J,LI)=ISKIP(J,LI)+1
                       END IF
  224 CONTINUE
      ISKIP(1,1)=1
      ISKIP(1,LI)=1
      ISKIP(LJ,1)=1
      ISKIP(LJ,LI)=1
      RETURN                                                            
      END                                                               
      