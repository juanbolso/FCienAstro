C PROGRAM RESONALYZERv2.F  - TABARE GALLARDO, JANUARY 2020
C gallardo@fisica.edu.uy
C CALCULATES RESONANT DISTURBING FUNCTION R(sigma), THE EQUILIBRIUM POINTS,
C LIBRATION PERIODS, RESONANCE STRENGTH AND WIDTH OF THE STABLE RESONANCE IN AU
C PLANETS AND STAR GIVEN   BY PLASYS.INP
C REFERENCE: GALLARDO 2019, CMDA
      IMPLICIT REAL*8 (A-H,J-Z)
      DIMENSION RP(400),SE(10),SM(10),EX(10),RH(10),DERIV1(400),
     *DERIV2(400),TLIB(400),MD(400),IRH(400)
C THIS 50000 DIMENSION IS FOR STORING DATA PARTICLE, CAN BE MODIFIED
C MAX NUMBER OF POINTS IN THE INTEGRAL   = 1000*MAX(KK,KP)
      DIMENSION VX2(50000),VY2(50000),VZ2(50000),VR2(50000),VLA2(50000)
C FOR PLOTTING WITH DISLIN  LIBRARIES
      REAL*4 RPG(400),SMA(400)
      CHARACTER*9 AEQ(400),ACE(400)
      character*40 zcomment
      TWOPI = 8.0D0*DATAN(1.0D0)
      CERO  = 0.0D0
      UNO   = 1.0D0
      PI=TWOPI/2.D0
      G2R=PI/180.D0
      KGAUS=0.01720209895D0
      KG2=KGAUS**2

      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'               RESONANCE ANALYZER version 2 '
      WRITE(*,*)'  for eccentric planetary system given by plasys.inp'
      WRITE(*,*)'  This program calculates R(sigma), its strength'
      WRITE(*,*)'SR = <R>-Rmin, total stable width of resonance in au,'
      WRITE(*,*)'    equilibrium points and libration periods. '
      WRITE(*,*)'    Output files: rsigma.dat and resonance.dat '
      WRITE(*,*)'Tabare Gallardo                gallardo@fisica.edu.uy '
      WRITE(*,*)'Reference: Gallardo 2020, CMDA.'
      WRITE(*,*)'               version January 22, 2020'
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'              I assume the archive plasys.inp is ready'
      WRITE(*,*)'press ENTER'
      READ(*,*)
C INPUT
      OPEN(1,FILE='plasys.inp',STATUS='old')
      READ (1,*)  ZCOMMENT
C MASS CENTRAL STAR
      READ(1,*) MST
      READ (1,*)  ZCOMMENT
C READ PLANETS A E MASS
      DO  I=1,10
        READ (1,*,ERR=199) SE(I),EX(I),SM(I)
      ENDDO

      READ (1,*)  ZCOMMENT
  199 ITOTPLA=I-1
C THERE ARE ITOTPLA PLANETS
      CLOSE (1)
C HILL RADIUS
      DO I=1,ITOTPLA
        RH(I)=SE(I)*(SM(I)/(MST+SM(I))/3.D0)**(1.D0/3.D0)
      ENDDO
C NUMBER OF RHILL ADMITTED TO CALCULATE DISTURBING FUNCTION
C THIS ADMIT A FINE TUNNING ACCORDING TO THE ORBITAL INCLINATION
C FOR DIRECT ORBITS 3<RHTOL<4
C FOR RETROGRADE ORBITS 2<RHTOL<3
      RHTOL=3.D0
C PLANETS TO CONSIDER
      WRITE(*,*)'you have defined ',itotpla,' planets'
      WRITE(*,*)'now indicate which planet generates the resonance'
      READ(*,*) Inpla
      WRITE(*,*)'RESONANCE kp:k ?'
      WRITE(*,*)'(ex. 2,3 for plutinos; 3,2 for hildas, <50,50)'
      READ(*,*) KP,KK
      IF(KP.GE.KK) THEN
        WRITE(*,*)'sigma =',INT(KK),'*lambda',INT(-KP),'*lambda_pla +',
     *INT(kp-kk),'*varpi'
        ELSE
        WRITE(*,*)'sigma =',INT(KK),'*lambda',INT(-KP),'*lambda_pla',
     *INT(kp-kk),'*varpi'
      ENDIF
C MAXIMUM FOR CALCULATION OF THE INTEGRAL WITH ENOUGH PRECISION
      MAXFAC=KP
      IF(KK.GT.KP) THEN
        MAXFAC=KK
      ENDIF
      WRITE(*,*) '   '
 112  WRITE(*,*)'ASTEROID''S ORBITAL ELEMENTS:  e,i,argper,lonNode ?'
C angles IN DEGREES
      READ(*,*)EXC,YNC,ARGPER,nnpla
C lon node planet = 0
C NODE ASTEROID - NODE PLANET
      LONOD=nnpla
      LOPER=ARGPER+LONOD
      IF (LOPER.GT.360.0) THEN
            LOPER=LOPER-360.0d0
      ENDIF

C THE PLANET IS 1
      A1=SE(INPLA)
      E1=EX(INPLA)
      J1G=0.D0
      L1G=0.D0
      P1G=0.D0
      J1=J1G*G2R
      L1=L1G*G2R
      P1=P1G*G2R
C PARTICLE IS 2
      A2=A1*((DABS(KK/KP))**2*MST/(MST+SM(INPLA)))**(1.D0/3.D0)
      SALSEM=A2
      E2=EXC
      J2GRA=YNC
      J2=J2GRA*G2R
      L2G=LONOD
      AR2=ARGPER*G2R
C      P2G=LOPER
      L2=L2G*G2R
C      P2=P2G*G2R
C NUMBER OF EVALUATIONS OF R(sigma) BETWEEN O AND 360 DEGREES
      ISIMAX=360
C R(360)=R(0)
C STEPS OF THE NUMERICAL INTEGRATION FROM 0 TO 2PI*KP IN lambda PARTICLE
C CAN BE AS SMALL AS 100*INT(MAXFAC), AT YOUR OWN RISK
      IPASOS=1000*INT(MAXFAC)
C ==============================================================
C FIRST CALCULATE THE DATA FOR PARTICLE, THEY WILL BE USED ISIMAX TIMES
        DO I2=1,IPASOS
C ALL POSSIBLE CONFIGURATIONS OCCUR AFTER KP REVOLUTIONS OF THE PARTICLE
C LAMBDA PARTIC
          LA2=DFLOAT(I2)/DFLOAT(IPASOS)*TWOPI*DABS(KP)
C STORE VECTOR DATA LAMBDA PART
          VLA2(I2)=LA2
C MEAN ANOMALY PARTIC
          AM2=LA2-(L2+AR2)
 200      IF (AM2.GT.TWOPI) THEN
            AM2=AM2-TWOPI
            GOTO 200
          END IF
 201      IF (AM2.LT.CERO) THEN
            AM2=AM2+TWOPI
            GOTO 201
          END IF
C SOLVING KEPLER FOR PARTICLE
          CALL SOLKEP(E2,AM2,AEX2,INITER)
C TRUE ANOMALY
          AVE2=DACOS((DCOS(AEX2)-E2)/(UNO-E2*DCOS(AEX2)))
          IF(AM2.GT.PI) AVE2=TWOPI-AVE2
C HELIOCENTRIC DISTANCE
          R2=A2*(UNO-E2*DCOS(AEX2))
C VECTOR DATA RPART
          VR2(I2)=R2
C RADIUS VECTOR FOR PARTICLE
C EC 2.122 MURRAY AND DERMOT
      X2=R2*(DCOS(L2)*DCOS(AR2+AVE2)-DSIN(L2)*DSIN(AR2+AVE2)*DCOS(J2))
      Y2=R2*(DSIN(L2)*DCOS(AR2+AVE2)+DCOS(L2)*DSIN(AR2+AVE2)*DCOS(J2))
          Z2=R2*DSIN(AR2+AVE2)*DSIN(J2)
C VECTOR DATA XYZ
          VX2(I2)=X2
          VY2(I2)=Y2
          VZ2(I2)=Z2
        ENDDO
C I HAVE CALCULATED ALL THE POSITIONS FOR PARTICLE THAT I WILL USE LATER
C ================================================================
C NOW DEFINE SIGMA
      DO ISI=1,ISIMAX
C ACRIT IS SIGMA FROM 0 TO 359
C CLOSE ENCOUNTER INDICATOR
        ACE(ISI)='         '
C LIBRATION PERIOD
        TLIB(ISI)=0.D0
        ACRIT=DFLOAT(ISI-1)/DFLOAT(ISIMAX)*360.D0
        TETA=ACRIT+(KK-KP)*LOPER
 901    IF (TETA.GT.360.D0) THEN
          TETA=TETA-360.D0
          GOTO 901
        END IF
 902    IF (TETA.LT.0.D0) THEN
          TETA=TETA+360.D0
          GOTO 902
        END IF
C TO RADIANS
        TETAR=TETA*G2R
C DISTURBING FUNCTION
        RTOT=CERO
C NUMBER OF ENCOUNTERS INSIDE THE HILL SPHERE
        IRH(ISI)=0
C MINIMUM DISTANCE PARTICLE-PLANET
        MINDIS=9999999.D0
C GIVEN SIGMA CALCULATE THE INTEGRAL
C CALCULATION OF THE INTEGRAL IN LAMBDA PARTICLE
C (IN THE PAPER IS IN LAMBDA PLANET)
        DO ILAMBDA2=1,IPASOS
C ALL POSSIBLE CONFIGURATIONS OCCUR AFTER KP REVOLUTIONS OF THE PARTICLE
C HELIOCENTRIC DISTANCE
          R2=VR2(ILAMBDA2)
C MEAN LONG BETWEEN 0 AND TWOPI*KP
          LA2=VLA2(ILAMBDA2)
C RADIUS VECTOR FOR PARTICLE
          X2=VX2(ILAMBDA2)
          Y2=VY2(ILAMBDA2)
          Z2=VZ2(ILAMBDA2)
C GIVEN LAMBDA PARTICLE CALCULATE LA1  PLANET
          LA1=(KK*LA2-TETAR)/KP
 300      IF (LA1.GT.TWOPI) THEN
            LA1=LA1-TWOPI
            GOTO 300
          END IF
 301      IF (LA1.LT.CERO) THEN
            LA1=LA1+TWOPI
            GOTO 301
          END IF
C LA1=AM1 BECAUSE NODE PLANET = 0
C SOLVING KEPLER FOR PLANET
          CALL SOLKEP(E1,LA1,AEX1,INITER)
C TRUE ANOMALY
          AVE1=DACOS((DCOS(AEX1)-E1)/(UNO-E1*DCOS(AEX1)))
          IF(LA1.GT.PI) AVE1=TWOPI-AVE1
C HELIOCENTRIC DISTANCE
          R1=A1*(UNO-E1*DCOS(AEX1))
C RADIUS VECTOR FOR PLANET
          X1=R1*DCOS(AVE1)
          Y1=R1*DSIN(AVE1)
          Z1=0.D0
C MUTUAL DISTANCE
          DELTA=DSQRT((X1-X2)**2+(Y1-Y2)**2+(Z2)**2)
C REGISTER CLOSE ENCOUNTERS AND NUMBER
          IF(DELTA.LT.RH(INPLA)*RHTOL) THEN
            ACE(ISI)='CLOSE ENC'
            IRH(ISI)=IRH(ISI)+1
          ENDIF
C REGISTER MINIMUM DISTANCE PLANET-PARTICLE IN RHILL
          DELRH=DELTA/RH(INPLA)
          IF(DELRH.LT.MINDIS)THEN
            MINDIS=DELRH
          ENDIF
C CALCULATION DISTURBING FUNCTION
C SCALAR PRODUCT R1*R2
          R1R2=X1*X2+Y1*Y2
C DIRECT PART
          RDIR=UNO/DELTA
C INDIRECT PART (PROBABLY THIS IS A WASTE OF TIME)
          RIND=-R1R2/R1**3
C DISTURBING FUNCTION
          RPER=RDIR+RIND
C SUMATION FOR A CRUDE INTEGRAL
C I TRYED SOFISTICATED INTEGRATIONS METHODS THAT GIVE ONLY TRASH
          RTOT=RTOT+RPER
        ENDDO
C END OF CALCULATION OF THE INTEGRAL
        RTOT=RTOT/DFLOAT(IPASOS)
C NOW MULTIPLY BY MASS_PLANET*k**2*MASSTAR
        RP(ISI)=RTOT*SM(INPLA)*KG2*mst
        SMA(ISI)=ACRIT
        MD(ISI)=MINDIS
      ENDDO
C END OF CALCULATION FOR ALL SIGMAS
C CALCULATION OF MEAN VALUE <R>  FOR CALCULATING STRENGTH SR
      VMEDIO=0.D0
      DO I2=1,ISIMAX
        VMEDIO=VMEDIO+RP(I2)
      ENDDO
      VMEDIO=VMEDIO/DFLOAT(ISIMAX)
C FIND R_max AND R_min
      VRMAX=-999999.D0
      PRMAX=-999999.D0
      VRMIN=999999.D0
      DO I2=1,ISIMAX
        IF(RP(I2).LT.VRMIN) THEN
          VRMIN=RP(I2)
        ENDIF
C ABSOLUTE RMAX
        IF(RP(I2).GT.PRMAX) THEN
             PRMAX=RP(I2)
        ENDIF
C FOR CALCULATING RMAX WE DISCARD CLOSE ENCOUNTERS TO LESS THAN RHTOL*RHILL
        IF(RP(I2).GT.VRMAX.AND.MD(I2).GT.RHTOL) THEN
          VRMAX=RP(I2)
        ENDIF
      ENDDO
C READY, RMAX AND RMIN CALCULATED DISCARDING CLOSE ENCOUNTERS
C CALCULATION OF DELTA R
         DELR=VRMAX-VRMIN
C JUST IN CASE...
         IF(DELR.LT.0.D0)   THEN
           DELR=0.D0
         ENDIF
C IF DELTA R IS VERY SMALL WE CANNOT DETECT ANY RESONANCE
C 1.D-11 IS A VERY REASONABLE LIMIT FOR RELATIVE VARIATION
C IT IS REALLY A SMALL VARIATION
C      DETEC=DELR/VMEDIO
      DETEC2=DABS((VMEDIO-VRMIN)/VMEDIO)
C      DETEC=DELR/VMEDIO
      IF(DETEC2.LT.1.D-11) THEN
         delr=0.d0
C        WRITE(*,*)'        '
C        WRITE(*,*)'*** UNDETECTABLE RESONANCE ***'
C        WRITE(*,*)'        '
        GOTO 111
      ENDIF
C CALCULATION OF R - R_min FOR PLOTTING WITH DISLIN LIBRARIES
C REQUIRES SIMPLE PRECISION
      MAXDELR=PRMAX-VRMIN
      DO I=1,ISIMAX
        RPG(I)=(RP(I)-VRMIN)/MAXDELR
      ENDDO
C FIRST DERIVATIVES IN RADIANS FOR EACH SIGMA
      DO I=1,ISIMAX
        ISIMA1=I+1
        ISIME1=I-1
        IF(ISIMA1.GT.360) ISIMA1=1
        IF(ISIME1.LT.1) ISIME1=360
        DERIV1(I)=(RP(ISIMA1)-RP(ISIME1))/2.D0/G2R
      ENDDO
C SECOND DERIVATIVES
      DO I=1,ISIMAX
        ISIMA1=I+1
        ISIME1=I-1
        IF(ISIMA1.GT.360) ISIMA1=1
        IF(ISIME1.LT.1) ISIME1=360
        DERIV2(I)=(DERIV1(ISIMA1)-DERIV1(ISIME1))/2.D0/G2R
      ENDDO
C DEFINING EQUILIBRIUM POINTS, NO CLOSE ENCOUNTERS ALLOWED
      DO I=1,ISIMAX
        AEQ(I)='         '
        ISIME1=I-1
        IF(ISIME1.LT.1) ISIME1=360
C CONDITION: NO CLOSE ENCOUNTERS AT LESS THAN 0.5 RHILL
C I AM INTERESTED IN DETECTING THE LOCATION OF THE STABLE EQ POINTS
C EVEN IF THEY ARE INSIDE 3*RHILL, BUT NOT INSIDE 0.5RHILL
        IF(MD(I).GT.0.5D0) THEN
C IF THE FIRST DERIVATIVE CHANGES ITS SIGN
          VALOR=DERIV1(I)
          IF(DERIV1(ISIME1)*VALOR.LT.0.D0)THEN
C DECIDE IF I OR I-1 IS THE POINT
            IEXTREMO=I
            IF(DERIV2(I).GT.0.)  THEN
              IF(RP(ISIME1).LT.RP(I)) THEN
                IEXTREMO=ISIME1
              ENDIF
C IT IS A STABLE EQUILIBRIUM POINT AND LIBRATION PERIOD IN YEARS
              AEQ(IEXTREMO)='E. STABLE'
      TLIB(IEXTREMO)=A2*TWOPI/DABS(KK)/DSQRT(3.D0*DABS(DERIV2(IEXTREMO
     *)))/365.25D0
              ELSE
C UNSTABLE
              IF(RP(ISIME1).GT.RP(I))THEN
                IEXTREMO=ISIME1
              ENDIF
            AEQ(IEXTREMO)=' UNSTABLE'
          ENDIF
        ENDIF
      ENDIF
      ENDDO
C ANCHO2 IS RESONANCE FULL stable WIDTH
 111     ANCHO2=8.D0/3.D0/(KG2*MST)*(DELR)*A2**3
      ANCHO2=2.D0*DSQRT(ANCHO2)
C FILE OUTPUT
      OPEN(2,FILE="resonance.dat",STATUS="UNKNOWN",ACCESS="APPEND")
      WRITE(2,*)'Star mass =  ',mst
      WRITE(2,*)'Planet mass=  ',sm(INPLA)
      WRITE(2,*)'Planet a =  ',se(INPLA)
      WRITE(2,*)'Planet e =  ',ex(INPLA)
      WRITE(2,*)'Resonance =  ',INT(KP),':',INT(KK)
      WRITE(2,*)'a = ',SALSEM, ' au'
      WRITE(2,*)'e = ',EXC
      WRITE(2,*)'i = ',YNC
      WRITE(2,*)'arg perihelion = ',ARGPER
      WRITE(2,*)'long node = ',nnpla
      WRITE(2,*)'<R> = ',VMEDIO, '  [M_sun,au,days]'
      WRITE(2,*)'Strength (<R>-Rmin) = ',(VMEDIO-VRMIN)
      WRITE(2,*)'Rmax-Rmin = ',DELR
      WRITE(2,*)'Close enc. dist. = ',rhtol, ' RHill'
      WRITE(2,*)'stable width = ',ANCHO2, ' au'
C SCREEN
      WRITE(*,*)'        '
      WRITE(*,*)'Star mass =  ',mst
      WRITE(*,*)'Planet mass=  ',sm(INPLA)
      WRITE(*,*)'Planet a =  ',se(INPLA)
      WRITE(*,*)'Planet e =  ',ex(INPLA)
      WRITE(*,*)'Resonance =  ',INT(KP),':',INT(KK)
      WRITE(*,*)'a = ',SALSEM, ' au'
      WRITE(*,*)'e = ',EXC
      WRITE(*,*)'i = ',YNC
      WRITE(*,*)'arg perihelion = ',ARGPER
      WRITE(*,*)'long node = ',nnpla
      WRITE(*,*)'<R> = ',VMEDIO, '  [M_sun,au,days]'
      WRITE(*,*)'Strength (<R>-Rmin) = ',(VMEDIO-VRMIN)
      WRITE(*,*)'Rmax-Rmin = ',DELR
      WRITE(*,*)'Close enc. dist. = ',rhtol, ' RHill'
      WRITE(*,*)'stable width = ',ANCHO2, ' au'
C WRITE THE FILE WITH R(SIGMA) AND OTHER DATA
      OPEN(1,FILE="rsigma.dat",STATUS="UNKNOWN")
      WRITE(1,*)' sigma   R-R_min '
      DO I=1,ISIMAX
C FOR STABLE EQUILIBRIUM POINTS WRITE LIBRATION PERIOD
        IF(TLIB(I).NE.0.D0) THEN
          WRITE(1,1001) SMA(I),RP(I)-VRMIN,ACE(I),AEQ(I),TLIB(I)
          WRITE(*,1004) 'sigma=',SMA(I),AEQ(I),'  P=',TLIB(I)," years"
          WRITE(2,1004) 'sigma=',SMA(I),AEQ(I),'  P=',TLIB(I)," years"
          ELSE
          WRITE(1,1002) SMA(I),RP(I)-VRMIN,ACE(I),AEQ(I)
        ENDIF
      ENDDO
      close(1)
      WRITE(*,*)'        '
      WRITE(2,*)'        '
      CLOSE(2)
C  DISLIN  START HERE
C IF YOU DONT WANT TO COMPILE WITH DISLIN LIBRARIES CUT HERE -----------
      call METAFL('XWIN')
      CALL DISINI
      CALL PAGHDR('Resonant Disturbing Function','T. Gallardo',4,0)
      CALL COMPLX
      CALL NAME('sigma','X')
      CALL NAME('R(sigma) - Rmin','Y')
      CALL TITLIN('Shape of the Resonant Disturbing Function',3)
      CALL GRAF(0.,360.,0.,60.,0.,1.,0.,0.2)
      CALL QPLOT(SMA,RPG,ISIMAX)
      CALL ENDGRF
      CALL DISFIN
C ----------------END DISLIN-------------------------------------------
      GOTO 112
 1001 FORMAT(F7.2,E16.8,2A10,E16.8)
 1002 FORMAT(F7.2,E16.8,2A10)
 1004 FORMAT(A7,F4.0,A10,A4,E16.8,A10)
      END

C ========================================================================
C SOLVING KEPLER
      SUBROUTINE SOLKEP(EX,M,E,NITER)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M,MK
      TOLE=1.D-12
      TWOPI=8.D0*DATAN(1.D0)
      M=DMOD(M,TWOPI)
      E=M
      NITER=0
 100  E0=E
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
C	 IF(NITER.GT.20)GOTO 200
      IF(DEX.GT.TOLE)GOTO 100
      RETURN
      END

