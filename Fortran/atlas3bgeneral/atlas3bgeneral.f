C ATLAS OF THREE BODY RESONANCES
C restricted case
C TABARE GALLARDO 2016
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KGAUS,KG2
C MAXIMUM NUMBER OF RESONANCES
      PARAMETER (MAXRES=100000)
      DIMENSION SE(10),SM(10),SN(10)
      DIMENSION NPLA1(MAXRES),NPLA2(MAXRES),KAST(MAXRES),KPLA1(MAXRES)
      DIMENSION KPLA2(MAXRES),AAST(MAXRES),INDICE(MAXRES)
      DIMENSION NKAST(MAXRES),NKPLA1(MAXRES),NKPLA2(MAXRES)
      DIMENSION NNPLA1(MAXRES),NNPLA2(MAXRES),SEMI(MAXRES)
      DIMENSION II(MAXRES)
      COMMON SE,SM,SN,STARM
      TWOPI = 8.0D0*DATAN(1.0D0)
      CERO  = 0.0D0
      UNO   = 1.0D0
      PI=TWOPI/2.D0
      G2R=PI/180.D0
      DOST=2.D0/3.D0
      UNT=1.D0/3.D0
      KGAUS=0.01720209895D0
      KG2=KGAUS**2

C THE UNIT OF STRENGHT IS THE VALUE FOR 490 VERITAS
C IN THE RESONANCE  2Veritas - 5Jupiter + 2Saturn

      OPEN(1,FILE="massivebodies.inp")
      read(1,1004) APESOS
C MASS OF THE central body (star or planet) IN SOLAR MASSES
      READ(1,*) STARM
      read(1,1004) APESOS
C NUMBER OF PLANETS
      READ(1,*) NTOT
C SEMIMAJOR AXES AND MASSES
      read(1,1004) APESOS
      DO I=1,NTOT
        READ(1,*) SE(I),SM(I)
      ENDDO
      CLOSE(1)

      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)' Given a planetary system this program calculates all '
      WRITE(*,*)' THREE BODY RESONANCES of a massless particle in a '
      WRITE(*,*)' given range of semimajor axes with all the planets '
      WRITE(*,*)' taken by pairs. Gallardo (2014, Icarus 231, 273).'
      WRITE(*,*)'          www.fisica.edu.uy/~gallardo/atlas'
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'  '
      WRITE(*,*)'Put the data for the central star and planets in the '
      WRITE(*,*)'input file              massivebodies.inp  '
      WRITE(*,*)'  '
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)' '

C max p
      WRITE(*,*)'max p = |k0| + |k1| + |k2| (less than 50)? '
      READ(*,*)SUPERIOR
      ISUP=IDINT(SUPERIOR)

C max q
      WRITE(*,*)'max order q = |k0 + k1 + k2| (less than 30)? '
      READ(*,*)QSUPERIOR

C PLANETS
      WRITE(*,*)'initial and final planet? (ex: 2,7) '
      READ(*,*) IPLA1,IPLA2
C LIMITS IN SEMIMAJOR AXIS
      WRITE(*,*)'interval in au: a_min,a_max? (ex: 3.1,3.3)'
      READ(*,*)A2MIN,A2MAX
C ORBITAL ELEMENTS OF THE test PARTICLE
      WRITE(*,*)'e,i,argper for the test particle? (ex: 0.3,35,127)'
      READ(*,*)EXC,YNCLINA,ARG
C TYPICAL OF MAIN BELT ASTEROIDS
C      exc=0.15
C      ynclina=6.0
C      arg=60.0
      DO I=1,MAXRES
        II(I)=0
      ENDDO
C PLANETARY MEAN MOTIONS
      DO I=1,NTOT
        SN(I)=DSQRT(KG2*(STARM+SM(I))/SE(I)**3)
      ENDDO
C LOOKING FOR ALL POSSIBLE RESONANCES
      I=0
      DO 100 IP1=IPLA1,IPLA2
      DO 101 IP2=IP1+1,IPLA2
      DO 102 IAST=1,ISUP
        CAST=IAST*1.D0
        DO 103 IPL1=1,2*ISUP+1
          CPL1=(IPL1-ISUP-1)*1.D0
          IF(CPL1.EQ.0.D0) THEN
            GOTO 103
          END IF
          DO 104 IPL2=1,2*ISUP+1
            CPL2=(IPL2-ISUP-1)*1.D0
          IF(CPL2.EQ.0.D0) THEN
            GOTO 104
          END IF
C THREE CONSTANTS CANNOT BE ALL POSITIVE
          IF(CPL1.GT.0.D0.AND.CPL2.GT.0.D0) THEN
            GOTO 104
          END IF
          HINDICE=CAST+DABS(CPL1)+DABS(CPL2)
C ORDER OF THE RESONANCE
          HORDER=DABS(CAST+CPL1+CPL2)
          
          AMM=-(CPL1*SN(IP1)+CPL2*SN(IP2))/CAST
C MEAN MOTIONS MUST BE POSITIVE
          IF(AMM.LE.0.D0) THEN
            GOTO 104
          ENDIF
        SEM=(KG2*STARM/AMM**2)**UNT
            IF (SEM.LE.A2MAX.AND.SEM.GE.A2MIN)  THEN
              IF(HINDICE.LE.SUPERIOR.AND.HORDER.LE.QSUPERIOR) THEN
C VALID RESONANCE
                I=I+1
                IF(I.GT.MAXRES) THEN
                  WRITE(*,*)"ERROR: TOO MUCH RESONANCES"
                  READ(*,*) NOTHING
                ENDIF
                KAST(I)=INT(CAST)
                KPLA1(I)=INT(CPL1)
                KPLA2(I)=INT(CPL2)
                NPLA1(I)=IP1
                NPLA2(I)=IP2
                AAST(I)=SEM
              ENDIF
            ENDIF
  104     CONTINUE
  103   CONTINUE
  102 CONTINUE
  101 CONTINUE
  100 CONTINUE
C NUMBER OF POSSIBLE RESONANCES, MUST BE LESS THAN 100000
      IDAT=I
C ELIMINATE REPEATED
      DO I=1,IDAT
        DO J=I+1,IDAT
          IF(NPLA1(J).EQ.NPLA1(I).AND.NPLA2(J).EQ.NPLA2(I)) THEN
            DIF=DABS(AAST(J)-AAST(I))
            IF (DIF.LT.0.000001D0.AND.ABS(KAST(J)).GT.ABS(KAST(I))) THEN
C 1 MEANS ELIMINATED
              II(J)=1
            ENDIF
            IF (DIF.LT.0.000001D0.AND.ABS(KAST(J)).LT.ABS(KAST(I))) THEN
              II(I)=1
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      NCONT=0
      DO I=1,IDAT
        IF (II(I).NE.1) THEN
C VALID RESONANCES
          NCONT=NCONT+1
          NKAST(NCONT)=KAST(I)
          NKPLA1(NCONT)=KPLA1(I)
          NKPLA2(NCONT)=KPLA2(I)
          NNPLA1(NCONT)=NPLA1(I)
          NNPLA2(NCONT)=NPLA2(I)
          SEMI(NCONT)=AAST(I)
          INDICE(NCONT)=ABS(KAST(I))+ABS(KPLA1(I))+ABS(KPLA2(I))
        ENDIF
      ENDDO
      WRITE(*,*)'CALCULATING ',NCONT,' RESONANCES'
C =================================================================
C STARTING CALCULATIONS FOR ALL VALID RESONANCES
        OPEN(1,FILE='atlas3bres.dat',STATUS='UNKNOWN',ACCESS='APPEND')
      WRITE(1,49) "the unit of strenght is the value for 490 Veritas ",
     *"in the resonance  2Ver - 5Ju + 2Sa"
        WRITE(1,49) "        n      a (au)  k0    k1 P1   k2 P2    p ",
     *"   q    e       i     w      strength "
        CLOSE(1)


      DO I=1,NCONT
        A0=SEMI(I)
        IP1=NNPLA1(I)
        IP2=NNPLA2(I)
        CTE0=DFLOAT(NKAST(I))
        CTE1=DFLOAT(NKPLA1(I))
        CTE2=DFLOAT(NKPLA2(I))
        ARGPER=ARG*G2R
        YNC=YNCLINA*G2R
        E0=EXC
        IORDEN=NKAST(I)+NKPLA1(I)+NKPLA2(I)
        CALL FUERZA3(I,A0,IP1,IP2,CTE0,CTE1,CTE2,ARGPER,YNC,E0,FME,FUE)
        OPEN(1,FILE='atlas3bres.dat',STATUS='UNKNOWN',ACCESS='APPEND')
        WRITE(1,46) I,SEMI(I),NKAST(I),NKPLA1(I),NNPLA1(I),NKPLA2(I),
     *NNPLA2(I),INDICE(I),abs(IORDEN),EXC,YNCLINA,ARG,FUE
        WRITE(*,48) I,NCONT,SEMI(I),NKAST(I),NKPLA1(I),NNPLA1(I),
     *NKPLA2(I),NNPLA2(I),FUE
        CLOSE(1)
      ENDDO
C 46   FORMAT(I10,F12.5,I4,I8,I2,I8,I2,2I8,2F6.1,F6.3,F25.20,E20.10)
 46   FORMAT(I10,F12.6,I4,I6,I2,I6,I2,2I6,F8.3,2F6.1,E14.4)
 48   FORMAT(2I6,F12.6,I4,I8,I2,I8,I2,F25.20)
 49   FORMAT(A49,A39)
1004  FORMAT(A76)
      STOP
      END
C ==================================================================
C CALCULATION OF THE STRENGTHS ACCORDING TO TABARE GALLARDO 2013
      SUBROUTINE FUERZA3(IRESO,A0,IP1,IP2,CTE0,CTE1,CTE2,ARGPER,J0,E0,
     *FME,FUE)
      IMPLICIT REAL*8 (A-H,J-Z)
      DIMENSION SE(10),SM(10),SN(10)
      DIMENSION RP(400)
      COMMON SE,SM,SN,STARM
      CERO  = 0.0D0
      UNO   = 1.0D0
      KGAUS=0.01720209895D0
      KG2=KGAUS**2
      TWOPI = 8.0D0*DATAN(1.0D0)
      PI=TWOPI/2.D0
      G2R=PI/180.D0
      ERROR=1.D-12
      IPOR=50
      POR=DFLOAT(IPOR)
C ORDER
      SUC=CTE0+CTE1+CTE2
      A1=SE(IP1)
      A2=SE(IP2)
C NUMBER OF EVALUATIONS OF RHO(SIGMA) BETWEEN O AND 360
C FOR MORE PRECISTION USE MORE POINTS
      ISIMAX=36
C STEPS OF THE NUMERICAL INTEGRATION
C FOR MORE PRECITION USE MORE POINTS
      IPASOS=180
C SIGMA = CRITICAL ANGLE
      DO 900 ISI=1,ISIMAX
C TETA IS SIGMA
        TETA=DFLOAT(ISI-1)/DFLOAT(ISIMAX)*360.D0
C TO RADIANS
        TETAR=TETA*G2R
C 0 IS THE PARTICLE
C 1 IS PLANET P1
C 2 IS PLANET P2
C INTEGRAL IN LAMBDA0, INTEGRAL IN LAMBDA1
        RTOT=CERO
        DELTAAM0=TWOPI/DFLOAT(IPASOS)
        AM0=-DELTAAM0
        DO 901 IMEDIA0=1,IPASOS
C PARTICLE IN A FIXED POSITION
          AM0=AM0+DELTAAM0
C ASSUME LONG NODE=0 THEN LONGPER=ARGPER
          L0=0.D0
          P0=ARGPER
          LAMBDA0=AM0+P0
 222      IF (LAMBDA0.GT.TWOPI) THEN
            LAMBDA0=LAMBDA0-TWOPI
            GOTO 222
          END IF
C SOLVING KEPLER FOR PARTICLE (ECCENTRIC)
          CALL SOLKEP(E0,AM0,AEX0,IN)
C TRUE ANOMALY
          AVE0=DACOS((DCOS(AEX0)-E0)/(UNO-E0*DCOS(AEX0)))
          IF(AM0.GT.PI) AVE0=TWOPI-AVE0
C HELIOCENTRIC DISTANCE
          R0=A0*(UNO-E0*DCOS(AEX0))
C RADIUS VECTOR FOR particle
          X0=R0*(DCOS(L0)*DCOS(P0-L0+AVE0)-DSIN(L0)*DSIN(P0-L0+AVE0)
     **DCOS(J0))
          Y0=R0*(DSIN(L0)*DCOS(P0-L0+AVE0)+DCOS(L0)*DSIN(P0-L0+AVE0)
     **DCOS(J0))
          Z0=R0*DSIN(P0-L0+AVE0)*DSIN(J0)
C PLANET 1
          DELTAAM1=TWOPI/DFLOAT(IPASOS)
          AM1=-DELTAAM1
          DO 902 IMEDIA1=1,IPASOS
C LAMBDA PLANET 1
            AM1=AM1+DELTAAM1
C ASSUMING LONG PERI=0  LONG NODE=0
            P1=0.D0
            L1=0.D0
            LAMBDA1=AM1+P1
C PLANET 1 IN CIRCULAR ORBIT
            AVE1=AM1
C HELIOCENTRIC DISTANCE
            R1=A1
C RADIUS VECTOR FOR PLANET 1
            X1=R1*DCOS(P1+AVE1)
            Y1=R1*DSIN(P1+AVE1)
            Z1=0.D0
C DISTANCE ASTEROID-P1
            DELTA01=DSQRT((X1-X0)**2+(Y1-Y0)**2+(Z1-Z0)**2)
C SCALAR PRODUCT R1*R0
            R1R0=X1*X0+Y1*Y0
C  DERIV R01 / DX1  and  DERIV R01 / DY1
            DR01DX1=(X0-X1)/DELTA01**3 - X0/R1**3 + 3.D0*X1*R1R0/R1**5
            DR01DY1=(Y0-Y1)/DELTA01**3 - Y0/R1**3 + 3.D0*Y1*R1R0/R1**5
C  DERIV R01 / DX0
            DR01DX0 = (X1-X0)/DELTA01**3 - X1/R1**3
C  DERIV R01 / DY0
            DR01DY0 = (Y1-Y0)/DELTA01**3 - Y1/R1**3
C  DERIV R01 / DZ0
            DR01DZ0 = (Z1-Z0)/DELTA01**3
C PLANET 2
C LAMBDA2 IS FIXED
            LAMBDA20=(TETAR + SUC*P0 -CTE0*LAMBDA0-CTE1*LAMBDA1)/CTE2
 502        IF (LAMBDA20.GT.TWOPI) THEN
              LAMBDA20=LAMBDA20-TWOPI
              GOTO 502
            END IF
 503        IF (LAMBDA20.LT.CERO) THEN
              LAMBDA20=LAMBDA20+TWOPI
              GOTO 503
            END IF
C K2 POSSIBLE VALUES F0R LAMBDA2
            DO 903 IL=1,INT(ABS(CTE2))
              LAMBDA2=LAMBDA20 + TWOPI*DFLOAT(IL-1)/DABS(CTE2)
C ASSUMING LONG PERI=0  LONG NODE=0
              L2=0.D0
              P2=0.D0
              AM2=LAMBDA2-P2
 202          IF (AM2.GT.TWOPI) THEN
                AM2=AM2-TWOPI
                GOTO 202
              END IF
 203          IF (AM2.LT.CERO) THEN
                AM2=AM2+TWOPI
                GOTO 203
              END IF
C PLANET IN CIRCULAR ORBIT
              AVE2=AM2
C HELIOCENTRIC DISTANCE
              R2=A2
C RADIUS VECTOR FOR PLANET  2
              X2=R2*DCOS(P2+AVE2)
              Y2=R2*DSIN(P2+AVE2)
              Z2=0.D0
C DISTANCE ASTEROID-P2
              DELTA02=DSQRT((X2-X0)**2+(Y2-Y0)**2+(Z2-Z0)**2)
C SCALAR PRODUCT R2*R0
              R2R0=X2*X0+Y2*Y0
C  DERIV R02 / DX2  Y  DERIV R02 / DY2
      DR02DX2=(X0-X2)/DELTA02**3 - X0/R2**3 + 3.D0*X2*R2R0/R2**5
      DR02DY2=(Y0-Y2)/DELTA02**3 - Y0/R2**3 + 3.D0*Y2*R2R0/R2**5
C DISTANCE P1-P2
              DELTA12=DSQRT((X2-X1)**2+(Y2-Y1)**2)
C  DERIV R21 / DX2  Y  DERIV R21 / DY2
              DR21DX2=(X1-X2)/DELTA12**3 - X1/R1**3
              DR21DY2=(Y1-Y2)/DELTA12**3 - Y1/R1**3
C  DERIV R12 / DX1  Y  DERIV R12 / DY1
              DR12DX1=(X2-X1)/DELTA12**3 - X2/R2**3
              DR12DY1=(Y2-Y1)/DELTA12**3 - Y2/R2**3
C  DERIV R02 / DX0
              DR02DX0 = (X2-X0)/DELTA02**3 - X2/R2**3
C  DERIV R02 / DY0
              DR02DY0 = (Y2-Y0)/DELTA02**3 - Y2/R2**3
C  DERIV R02 / DZ0
              DR02DZ0 = (Z2-Z0)/DELTA02**3
      COSO = 2.D0*(DR02DX0*DR01DX0 + DR02DY0*DR01DY0 + DR02DZ0*DR01DZ0)
              DIRECTA = COSO
              YNDIRECTA = DR02DX2*DR21DX2 + DR02DY2*DR21DY2
      YNDIRECTA = YNDIRECTA + DR01DX1*DR12DX1 + DR01DY1*DR12DY1
              RTOT=RTOT + YNDIRECTA + DIRECTA
 903        CONTINUE
 902      CONTINUE
 901    CONTINUE
C END OF LOOP FOR LAMBDAS
C VALUE OF THE INTEGRAL
      RTOT=SM(IP1)*SM(IP2)*RTOT/DABS(CTE2)/DFLOAT(IPASOS)**2
C DELTA T SQR OVER 2
      DTSQ2=0.5D0*A0*A1*A2*(2.D0*PI/KGAUS)**2/STARM
C RHO(SIGMA):
      RP(ISI)=RTOT*DTSQ2
 900  CONTINUE
C END OF LOOP FOR SIGMA
C ====================================================================
C CALCULATION OF MEAN VALUE <rho>
      VMEDIO=0.D0
      DO I2=1,ISIMAX
        VMEDIO=VMEDIO+RP(I2)
      ENDDO
      VMEDIO=VMEDIO/DFLOAT(ISIMAX)
C FIND rho_max AND rho_min
      VRMAX=-999999.D0
      VRMIN=999999.D0
      DO I2=1,ISIMAX
        IF(RP(I2).LT.VRMIN) THEN
          VRMIN=RP(I2)
        ENDIF
        IF(RP(I2).GT.VRMAX) THEN
          VRMAX=RP(I2)
        ENDIF
      ENDDO
      FME=VMEDIO
C      FUE=VMEDIO-VRMIN
C STRENGTH
C THE UNIT OF STRENGHT IS THE VALUE FOR 490 VERITAS = 0.5358E-04
C RESONANCE  2Veritas - 5Jupiter + 2Saturn
      FUE=(VRMAX-VRMIN)/2.d0/0.5358D-04
      RETURN
 40   FORMAT(I10,F6.1,E20.10)
      END
C ========================================================================
C SOLVING KEPLER EQUATION (Adrian Brunini)
	SUBROUTINE SOLKEP(EX,M,E,NITER)
C SOLUCION ITERATIVA DE LA ECUACION DE KEPLER
C ENTRA:EX   EXCENTRICIDAD            (<1)
C       M    ANOMALIA MEDIA           (RADIANES)
C SALE: E    ANOMALIA EXCENTRICA      (RADIANES)
C
	 IMPLICIT REAL*8 (A-H,O-Z)
	 REAL*8 M,MK
	 TOLE=1.D-11
	 DOSPI=8.D0*DATAN(1.D0)
	 M=DMOD(M,DOSPI)
	 E=M
	 NITER=0
 100     E0=E
	 SE=DSIN(E0)
	 CE=DCOS(E0)
	 ES=EX*SE
	 EC=1.D0-EX*CE
	 MK=E0-ES
	 U=(MK-M)/EC
	 XPRI=E0-U
	 XSEG=E0-U/(1.D0-U*ES)
	 E=(XPRI+XSEG)/2.D0
	 DEX=DABS(E-E0)
	 NITER=NITER+1
	 IF(NITER.GT.20)GOTO 200
	 IF(DEX.GT.TOLE)GOTO 100
	 RETURN
200      CONTINUE
	 NDIC=0
	 E0=-DOSPI
	 DE0=DOSPI/10.D0
400      DE=DE0/(10.D0**NDIC)
	 SE=DSIN(E0)
	 CE=DCOS(E0)
	 ES=EX*SE
	 EM0=E0-ES-M
	 NITER=0
300      E1=E0+DE
	 NITER=NITER+1
	 IF(NITER.GT.100)THEN
	 WRITE(*,*)'ERROR IN KEPLER'
	 RETURN
	 ENDIF
	 SE=DSIN(E1)
	 CE=DCOS(E1)
	 ES=EX*SE
	 EM1=E1-ES-M
	 IF(EM1*EM0.GT.0.D0)THEN
	 E0=E1
	 EM0=EM1
	 GOTO 300
	 ELSE
	 NDIC=NDIC+1
	 IF(NDIC.EQ.3)THEN
	 E=E1
	 RETURN
	 ENDIF
	 GOTO 400
	 ENDIF
	 RETURN
	 END
C===================================================================



