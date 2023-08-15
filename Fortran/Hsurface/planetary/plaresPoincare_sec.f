C PROGRAM plarespoincare  - TABARE GALLARDO, july 2020
C gallardo@fisica.edu.uy
C CALCULATES RESONANT DISTURBING FUNCTION R(sigma), THE EQUILIBRIUM POINTS,
C LIBRATION PERIODS, RESONANCE STRENGTH AND WIDTH OF THE RESONANCE IN AU
C PLANETS AND STAR GIVEN   BY PLARES.INP
      IMPLICIT REAL*8 (A-H,J-Z)
      INTEGER NEQ,NEQ_max,Ne,Nw,CONT,I_1,I_2,IAUX,IMAX,IMIN,ISIMAX,ISIMI
     *N,IND,I_EQ,AM_mode
      DIMENSION SE(2),SM(2),EX(2),YN(2),LN(2),LP(2)
      DIMENSION RP(1441),DERIV1(1441),DERIV2(1441),TLIB(1441),MD(1441)
C THIS 50000 DIMENSION IS FOR STORING DATA PARTICLE, CAN BE MODIFIED
C MAX NUMBER OF POINTS IN THE INTEGRAL   = 1000*MAX(K1,K2)
      DIMENSION VX2(50000),VY2(50000),VZ2(50000),VR2(50000),VLA2(50000)
      DIMENSION VXP2(50000),VYP2(50000),VZP2(50000)
C FOR PLOTTING WITH DISLIN  LIBRARIES
      REAL*4 SMA(1441)
      CHARACTER*250 RHTOL_str, PREC_str, AM_str
      CHARACTER*40 zcomment, name1, name2
      LOGICAL PRIMERO, IGUAL, EQUILIBRIO, CONFIRMADO, NOISE_mode
      
      TWOPI = 8.0D0*DATAN(1.0D0)
      CERO  = 0.0D0
      UNO   = 1.0D0
      PI=TWOPI/2.D0
      G2R=PI/180.D0
      KGAUS=0.01720209895D0
      KG2=KGAUS**2

      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'           *** plaresPoincare_secV4 ***               '
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'                                                      '
      WRITE(*,*)'This program allows the study of the secular evolution'
      WRITE(*,*)'of 2 eccentric coplanar planets locked in resonance.  '
      WRITE(*,*)'                                                      '
      WRITE(*,*)'The system info is excpected to be in plares.inp with '
      WRITE(*,*)'the following convention:                             '
      WRITE(*,*)'Planet 1 is the internal and planet 2 is the external.'
      WRITE(*,*)'In this file, lines 1 and 3 are skipped.              '
      WRITE(*,*)'The 2nd line should contain the central body mass and '
      WRITE(*,*)'the 4th and 5th lines the planets orbital elements:   '
      WRITE(*,*)' a(au) e i(ø) node(ø) perihelia(ø) mass(solar masses) '
      WRITE(*,*)'                                                      '
      WRITE(*,*)'PROCEDURE: calculate R(sigma), then its local minimums'
      WRITE(*,*)'in order to find the resonant equilibrium points.     '
      WRITE(*,*)'Also calculates minimum distances to find encounters. '
      WRITE(*,*)'This procedure is looped in 2 variables:              '
      WRITE(*,*)'  1) delta_varpi = varpi1 - varpi2 = varpi1           '
      WRITE(*,*)'  2) e1 or e2 depending of their variation range.     '
      WRITE(*,*)'Note: e1 and e2 are related through AM(=constant).    '
      WRITE(*,*)'                                                      '
      WRITE(*,*)'OUTPUT FILES:                                         '
      WRITE(*,*)'HSUP1_prec_RHtol=_XX.dat                              '
      WRITE(*,*)'HSUP2_prec_RHtol=_XX.dat                              '
      WRITE(*,*)'These files contain the data to plot the H surfaces.  '
      WRITE(*,*)'                                                      '
      WRITE(*,*)'Dont hesitate in contact us for any doubt/problem!    '
      WRITE(*,*)'Tabare Gallardo - gallardo@fisica.edu.uy'
      WRITE(*,*)'Juan Pons - juan.pons.93@gmail.com'
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)' P.S.: I assume the input file plares.inp is ready!   '
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'                                                      '
      
C INPUT
      OPEN(1,FILE='plares.inp',STATUS='old')
      READ (1,*)  ZCOMMENT
C MASS CENTRAL STAR
      READ(1,*) MST
      READ (1,*)  ZCOMMENT
C READ 2 PLANETS A E i lonod, loper, MASS
      DO  I=1,2
        READ (1,*) SE(I),EX(I),YN(I),LN(I),LP(I),SM(I)
      ENDDO
      CLOSE (1)

C WE ASSUME COPLANAR CASE (for now):
      YN(1) = 0.D0
      YN(2) = 0.D0
      LN(1) = 0.D0
      LN(2) = 0.D0
C AND EXTERIOR VARPI2 = 0 (-> deltavarpi = varpi1):
      LP(2) = 0.D0

  77  WRITE(*,*)'RESONANCE K2:K1 ? (K2 >= K1)'
      WRITE(*,*)'(ex. 3,2 for hildas or plutinos, <50,50)'
      READ(*,*) K2,K1
      IF(K2.lt.K1) THEN
        WRITE(*,*)'I SAID K2 >= K1'
        GOTO 77
      ENDIF
      IF(K2.GE.K1) THEN
        WRITE(*,*)'sigma =',INT(K1),'*lambda1',INT(-K2),'*lambda2 +',
     *INT(K2-K1),'*varpi1'
      ENDIF
C MAXIMUM FOR CALCULATION OF THE INTEGRAL WITH ENOUGH PRECISION
      MAXFAC=K2
      WRITE(*,*)'                                                      '
C NUMBER OF RHILL ADMITTED TO CALCULATE DISTURBING FUNCTION
C THIS ADMIT A FINE TUNNING ACCORDING TO THE ORBITAL INCLINATION
      WRITE(*,*)'CRITERIA FOR ENCOUNTER: min{|r-rp|}<RHtol*RHill'
      WRITE(*,*)'Usual values for direct orbits: 3<RHtol<4 '
      WRITE(*,*)'Usual values for retrograde orbits: 2<RHtol<3'
      WRITE(*,*)'ENTER RHtol DESIRED VALUE PLEASE:'
      READ(*,*) RHTOL
	  IF (RHTOL.LT.0) THEN
        RHTOL = 3.D0
        WRITE(*,*)'Invalid value. RHTOL set to 3.'
      ENDIF
      WRITE(RHTOL_str, '(F4.2)') RHTOL
      WRITE(*,*)'                                                      '
C REGISTER MIN DISTANCES TO WRITE IN A FILE
      WRITE(*,*)'WANT TO REGISTER MIN. DISTANCES? (0=NO, 1=YES)'
      READ(*,*) MIN_flag
      IF ((MIN_flag.NE.0).AND.(MIN_flag.NE.1)) THEN
        MIN_flag = 0
        WRITE(*,*)'Invalid value. By default min. distances are not calc
     *ulated.'
      ENDIF
      WRITE(*,*)'                                                      '
C The precision of calculation defines number of points in e,varpi y sigma:
      WRITE(*,*)'PRECISION OF CALCULATION? (0=LOW, 1=MIDDLE, 2=HIGH)'
      READ(*,*) PREC
      PREC = INT(PREC)
      IF ((PREC.NE.0).AND.(PREC.NE.1).AND.(PREC.NE.2)) THEN
        PREC = 0
        WRITE(*,*)'Invalid value. PRECISION set to 0.'
      ENDIF
      WRITE(*,*)'                                                      '
C WE SWEPT e AND varpi OF THE PARTICLE.
C DELTA_SMA defines de precision of sigma variable.
      IF (PREC.EQ.0) THEN
        PREC_str = 'low'
        DELTA_SMA = 3.D0
        eMIN = 0.01D0
        eMAX = 0.96D0
        deltaE = 0.05D0
        wMIN = 0.D0
        wMAX = 350.D0
        deltaW = 10.D0
      ENDIF
      IF (PREC.EQ.1) THEN
        PREC_str = 'middle'
        DELTA_SMA = 1.D0
        eMIN = 0.01D0
        eMAX = 0.99D0
        deltaE = 0.01D0
        wMIN = 0.D0
        wMAX = 357.D0
        deltaW = 3.D0
      ENDIF
      IF (PREC.EQ.2) THEN
        PREC_str = 'high'
        DELTA_SMA = 0.25D0
        eMIN = 0.01D0
        eMAX = 0.99D0
        deltaE = 0.005D0
        wMIN = 0.D0
        wMAX = 359.D0
        deltaW = 1.D0
      ENDIF

C From min,max,deltas, we obtain the number of indexes for main DO loops.
      Ne = NINT((eMAX - eMIN)/deltaE + 1)
      Nw = NINT((wMAX - wMIN)/deltaW + 1)

C EXTERIOR PLANET IS 2
      A2=SE(2)

C MEAN MOTION N2
      N2=DSQRT(KG2*(MST+SM(2))/A2**3)
      E2=EX(2)
      B2=A2*DSQRT(UNO-E2*E2)
      J2G=YN(2)
      J2=J2G*G2R
      L2G=LN(2)
      P2G=LP(2)
      L2=L2G*G2R
      P2=P2G*G2R
      ARGPER=P2G-L2G
      AR2=ARGPER*G2R
      
C THE INTERIOR PLANET IS 1
      N1=K2*N2/K1
      A1=(KG2*(MST+SM(1))/N1**2)**(UNO/3.D0)
      E1=EX(1)
      B1=A1*DSQRT(UNO-E1*E1)
      J1G=YN(1)
      L1G=LN(1)
      P1G=LP(1)
      J1=J1G*G2R
      L1=L1G*G2R
      P1=P1G*G2R
      AR1=P1-L1
      
C SECOND DERIVATIVE HII
      BETA1=SM(1)*MST/(SM(1)+MST)
      BETA2=SM(2)*MST/(SM(2)+MST)
      HII=3.D0*(K1/A1)**2/BETA1 + 3.D0*(K2/A2)**2/BETA2

C RHILL MUTUAL
      RHS=(A1+A2)/2.D0*((SM(1)+SM(2))/(3.D0*MST))**(1.D0/3.D0)

C Constants used for the angular momentum definition
      C1 = (SM(1)*MST/(SM(1)+MST))*(N1*A1**2)
      C2 = (SM(2)*MST/(SM(2)+MST))*(N2*A2**2)
      AM_min = 0
      AM_max = C1 + C2
        
C The precision of calculation defines number of points in e,varpi y sigma:
      WRITE(*,*)'WANT TO CALCULATE ANGULAR MOMENTUM FROM SYSTEM DATA IN
     *INPUT FILE (0) OR ENTER ANGULAR MOMENTUM VALUE MANUALLY (1)?'
      READ(*,*) AM_mode
      AM_mode = INT(AM_mode)
      IF ((AM_mode.NE.0).AND.(AM_mode.NE.1)) THEN
        AM_mode = 0
        WRITE(*,*)'Invalid value. Angular momentum will be calculated fr
     *om system data in input file.'
      ENDIF
      WRITE(*,*)'                                                      '
      
C ANGULAR MOMENTUM (CTE.):
      IF (AM_mode.EQ.1) THEN
        WRITE(*,*)'Minimum angular momentum is:',AM_min
        WRITE(*,*)'Maximum angular momentum is:',AM_max
 135    WRITE(*,*)'Please, enter the fraction of maximum angular momentu
     *m to be considered (a number between 0 and 1)'
        READ(*,*) AM
        IF((AM.LT.0).OR.(AM.GT.1)) GOTO 135
        AM = AM*AM_max
      WRITE(*,*)'Angular momentum = ',AM
      ELSE
        AM = C1*DSQRT(1-E1**2) + C2*DSQRT(1-E2**2)
        WRITE(*,*)'Angular momentum = ',AM
        WRITE(*,*)'which corresponds to the',(100*AM/AM_max),'% of the m
     *maximum angular momentum.'
      ENDIF
C Lets put in the output filename the normalized AM
      WRITE(AM_str, '(F6.4)') AM/AM_max
      WRITE(*,*)'                                                      '
      WRITE(*,*)'We have that: AM = C1*SQRT(1-e1^2) + C2*SQRT(1-e2^2)'
      WRITE(*,*)'where C1 =',C1,'and C2=',C2,'. So, we are in the...'
      WRITE(*,*)'                                                      '
      
C We inform the user min and max values of e1,e2.
C Therer are 4 different cases which determine the way of calculation
      IF ((AM.GT.C1).AND.(AM.GT.C2)) THEN
        E1min = 0
        E2min = 0
        E1max = DSQRT(1-((AM-C2)/C1)**2)
        E2max = DSQRT(1-((AM-C1)/C2)**2)
        WRITE(*,*)'...case 1: AM > C1,C2'
      ENDIF
      IF ((AM.GT.C1).AND.(AM.LT.C2)) THEN
        E1min = 0
        E2min = DSQRT(1-(AM/C2)**2)
        E1max = 1
        E2max = DSQRT(1-((AM-C1)/C2)**2)
        WRITE(*,*)'...case 2: AM > C1 ; AM < C2'
      ENDIF
      IF ((AM.LT.C1).AND.(AM.GT.C2)) THEN
        E1min = DSQRT(1-(AM/C1)**2)
        E2min = 0
        E1max = DSQRT(1-((AM-C2)/C1)**2)
        E2max = 1
        WRITE(*,*)'...case 3: AM < C1 ; AM > C2'
      ENDIF
      IF ((AM.LT.C1).AND.(AM.LT.C2)) THEN
        E1min = DSQRT(1-(AM/C1)**2)
        E2min = DSQRT(1-(AM/C2)**2)
        E1max = 1
        E2max = 1
        WRITE(*,*)'...case 4: AM < C1,C2'
      ENDIF
      
      WRITE(*,*)'                                                      '
      WRITE(*,*)'(MIN e1, MAX e1) = (', E1min,',',E1max,')'
      WRITE(*,*)'(MIN e2, MAX e2) = (', E2min,',',E2max,')'
      WRITE(*,*)'                                                      '
C From max and mins we obtain the range of variaton of e1 and e2
      deltaE1 = E1max - E1min
      deltaE2 = E2max - E2min
      
C Output files name definition (one per H surface):
      name1='HSUP1_'//TRIM(PREC_str)//'_RHtol='//TRIM(RHTOL_str)//'_AM='
     *//TRIM(AM_str)//'.dat'
      name2='HSUP2_'//TRIM(PREC_str)//'_RHtol='//TRIM(RHTOL_str)//'_AM='
     *//TRIM(AM_str)//'.dat'
      OPEN(1,FILE=name1,STATUS="UNKNOWN",ACCESS="APPEND")
      OPEN(2,FILE=name2,STATUS="UNKNOWN",ACCESS="APPEND")
      
      IF (MIN_flag.EQ.0) THEN
C       OPEN(3,FILE='min_distances.dat',STATUS="UNKNOWN",ACCESS="APPEND")
        WRITE(1,*)' e1   deltavarpi   sigma1   R '
        WRITE(1,*)'                              '
        WRITE(2,*)' e2   deltavarpi   sigma1   R '
        WRITE(2,*)'                              '
      ELSE
        WRITE(1,*)' e1   deltavarpi   sigma1   R   mindis '
        WRITE(1,*)'                                       '
        WRITE(2,*)' e2   deltavarpi   sigma1   R   mindis '
        WRITE(2,*)'                                       '
      ENDIF
      
C NUMBER OF EVALUATIONS OF R(sigma) BETWEEN O AND 360 DEGREES
      ISIMIN = 1
      ISIMAX=INT(360/DELTA_SMA)
C R(360)=R(0)
C STEPS OF THE NUMERICAL INTEGRATION FROM 0 TO 2PI*K2 IN lambda PARTICLE
C CAN BE AS SMALL AS 100*INT(MAXFAC), AT YOUR OWN RISK
      IPASOS=1000*INT(MAXFAC)
      
C DOUBLE DO LOOP TO SWEEP ECCENTRICITY AND LONGITUDE OF PERIASTRO DIFFERENCE:
C Lets assume varpi2 = 0 and just sweep varpi1 from 0 to 360.
C Lets sweep only e1. e2 is determined from angular momentum conservation.
C **********************************************************************
C *************ASTEROID'S EXCENTRICITY AND PERI-CENTER LOOP:************
C **********************************************************************
      DO ind_e = 1, Ne
      EXC = (ind_e-1)*deltaE + eMIN
      DO ind_w = 1, Nw
      loper = (ind_w-1)*deltaW + wMIN
C Flag to check if we are in the first case of the (e,w) grid.
      PRIMERO = (ind_e.EQ.1).AND.(ind_w.EQ.1)
C MODIFY VARIABLES AFFECTED BY LOOPS:
      P1G=loper
      P1=P1G*G2R
      AR1=P1-L1
C Choose which excentricty is more convenient to sweep
      IF (deltaE1.GT.deltaE2) THEN
        E1=EXC
        E2 = 1-((AM-C1*DSQRT(1-E1**2))/C2)**2
C Check if given e1 and AM, e2 exists:
        IF ((E2.LT.0.D0).OR.(AM.LT.C1*DSQRT(1-E1**2))) GOTO 321
        E2 = DSQRT(E2)
      ELSE
        E2=EXC
        E1 = 1-((AM-C2*DSQRT(1-E2**2))/C1)**2
C Check if given e1 and AM, e2 exists:
        IF ((E1.LT.0.D0).OR.(AM.LT.C2*DSQRT(1-E2**2))) GOTO 321
        E1 = DSQRT(E1)
      ENDIF
      
C TO AVOID NUMERICAL PROBLEMS:
C      IF ((E1.LT.0.01).OR.(E2.LT.0.01)) GO TO 312
      IF ((E1.GT.0.99).OR.(E2.GT.0.99)) GO TO 312
      
      B1=A1*DSQRT(UNO-E1*E1)
      B2=A2*DSQRT(UNO-E2*E2)
      
C CONSTANTS L1,M1,N1 FOR PLANET 2 ROY PAGE  93
      L12=DCOS(L2)*DCOS(AR2)-DSIN(L2)*DSIN(AR2)*DCOS(J2)
      M12=DSIN(L2)*DCOS(AR2)+DCOS(L2)*DSIN(AR2)*DCOS(J2)
      N12=DSIN(AR2)*DSIN(J2)

C CONSTANTS L2,M2,N2 FOR PLANET 2, ROY PAGE  93
      L22=-DCOS(L2)*DSIN(AR2)-DSIN(L2)*DCOS(AR2)*DCOS(J2)
      M22=-DSIN(L2)*DSIN(AR2)+DCOS(L2)*DCOS(AR2)*DCOS(J2)
      N22=DCOS(AR2)*DSIN(J2)

C CONSTANTS L1,M1,N1 FOR PLANET 1 ROY PAGE  93
      L11=DCOS(L1)*DCOS(AR1)-DSIN(L1)*DSIN(AR1)*DCOS(J1)
      M11=DSIN(L1)*DCOS(AR1)+DCOS(L1)*DSIN(AR1)*DCOS(J1)
      N11=DSIN(AR1)*DSIN(J1)
C CONSTANTS L2,M2,N2 FOR PLANET 1  ROY PAGE  93
      L21=-DCOS(L1)*DSIN(AR1)-DSIN(L1)*DCOS(AR1)*DCOS(J1)
      M21=-DSIN(L1)*DSIN(AR1)+DCOS(L1)*DCOS(AR1)*DCOS(J1)
      N21=DCOS(AR1)*DSIN(J1)

C ==============================================================
C FIRST CALCULATE THE DATA EXT PLANET, THEY WILL BE USED ISIMAX TIMES
        DO I2=1,IPASOS
C ALL POSSIBLE CONFIGURATIONS OCCUR AFTER K1 REVOLUTIONS OF THE EXT PLA
C LAMBDA PLA 2
          LA2=DFLOAT(I2)/DFLOAT(IPASOS)*TWOPI*DABS(K1)
C STORE VECTOR DATA LAMBDA PLA2
          VLA2(I2)=LA2
C SE PODRIA CALCULAR EN ANOM EXCENTRICA
C MEAN ANOMALY PLA2
          AM2=LA2-P2
          AM2=DMOD(AM2,TWOPI)
          
C SOLVING KEPLER FOR PLA2
          CALL SOLKEP(E2,AM2,AEX2,INITER)
          COSE=DCOS(AEX2)
          SINE=DSIN(AEX2)
C HELIOCENTRIC DISTANCE
          R2=A2*(UNO-E2*COSE)
C FORMULAS ROY PAGE 93
          X2=A2*L12*COSE+B2*L22*SINE-A2*E2*L12
          Y2=A2*M12*COSE+B2*M22*SINE-A2*E2*M12
          Z2=A2*N12*COSE+B2*N22*SINE-A2*E2*N12
C VELOCITIES
          XP2=(B2*L22*COSE-A2*L12*SINE)*N2*A2/R2
          YP2=(B2*M22*COSE-A2*M12*SINE)*N2*A2/R2
          ZP2=(B2*N22*COSE-A2*N12*SINE)*N2*A2/R2


C VECTOR DATA RPLA2
          VR2(I2)=R2

C VECTOR DATA XYZ
          VX2(I2)=X2
          VY2(I2)=Y2
          VZ2(I2)=Z2
C VECTOR DATA XPYPZP
          VXP2(I2)=XP2
          VYP2(I2)=YP2
          VZP2(I2)=ZP2

        ENDDO
C I HAVE CALCULATED ALL THE POSITIONS FOR PLA2 THAT I WILL USE LATER
C ================================================================
C NOW DEFINE SIGMA1
      DO ISI=ISIMIN,ISIMAX
C ACRIT IS SIGMA1 FROM 0 TO 359
C LIBRATION PERIOD
        TLIB(ISI)=0.D0
        ACRIT=DFLOAT(ISI-1)/DFLOAT(ISIMAX)*360.D0
C THETA IS THE COMBINATION OF LAMBDAS
        TETA=ACRIT+(K1-K2)*P1G
        TETA=DMOD(TETA,360.D0)
C TO RADIANS
        TETAR=TETA*G2R
C DISTURBING FUNCTION
        RTOT=CERO
C MINIMUM DISTANCE PARTICLE-PLANET
        MINDIS=9999999.D0
C GIVEN SIGMA1 CALCULATE THE INTEGRAL
C CALCULATION OF THE INTEGRAL IN LAMBDA PLA2
        DO ILAMBDA2=1,IPASOS
C ALL POSSIBLE CONFIGURATIONS OCCUR AFTER K2 REVOLUTIONS OF THE PLANET 1
C HELIOCENTRIC DISTANCE
          R2=VR2(ILAMBDA2)
C MEAN LONG BETWEEN 0 AND TWOPI*K2
          LA2=VLA2(ILAMBDA2)
C RADIUS VECTOR FOR PLA EXT 2
          X2=VX2(ILAMBDA2)
          Y2=VY2(ILAMBDA2)
          Z2=VZ2(ILAMBDA2)
C VELOCITY VECTOR FOR PLA EXT 2
          XP2=VXP2(ILAMBDA2)
          YP2=VYP2(ILAMBDA2)
          ZP2=VZP2(ILAMBDA2)

C GIVEN LAMBDA PLA EXTERIOR 2 CALCULATE LA1  PLANET  INT 1
          LA1=(K2*LA2+TETAR)/K1
C MEAN ANOMALY PLA 1
          AM1=LA1-P1
          AM1=DMOD(AM1,TWOPI)

C SOLVING KEPLER FOR PLANET
          CALL SOLKEP(E1,AM1,AEX1,INITER)
          COSE=DCOS(AEX1)
          SINE=DSIN(AEX1)
C HELIOCENTRIC DISTANCE
          R1=A1*(UNO-E1*COSE)
C FORMULAS ROY PAGE 93
          X1=A1*L11*COSE+B1*L21*SINE-A1*E1*L11
          Y1=A1*M11*COSE+B1*M21*SINE-A1*E1*M11
          Z1=A1*N11*COSE+B1*N21*SINE-A1*E1*N11
C BARICENTRIC VELOCITIES
          XP1=(B1*L21*COSE-A1*L11*SINE)*N1*A1/R1
          YP1=(B1*M21*COSE-A1*M11*SINE)*N1*A1/R1
          ZP1=(B1*N21*COSE-A1*N11*SINE)*N1*A1/R1

C MUTUAL DISTANCE
          DELTA=DSQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)

C REGISTER MINIMUM DISTANCE PLANET-PARTICLE IN RHILL (FOR ENCOUNTERS)
          DELRH=DELTA/RHS
          IF(DELRH.LT.MINDIS)THEN
            MINDIS=DELRH
          ENDIF
C CALCULATION DISTURBING FUNCTION
C SCALAR PRODUCT V1*V2
          V1V2=XP1*XP2+YP1*YP2+ZP1*ZP2

C INDIRECT PART
          RIND=-SM(1)*SM(2)*V1V2/MST
C DIRECT PART
          RDIR=KG2*SM(1)*SM(2)/DELTA
C DISTURBING FUNCTION
          RPER=RDIR+RIND
C SUMATION FOR A CRUDE INTEGRAL
          RTOT=RTOT+RPER
        ENDDO
C END OF CALCULATION OF THE INTEGRAL
        RTOT=RTOT/DFLOAT(IPASOS)
        RP(ISI)=RTOT
        SMA(ISI)=ACRIT
        MD(ISI)=MINDIS
      ENDDO
C *******END OF CALCULATION OF DISTURBING FUNCTION FOR ALL SIGMAS*******

C CALCULATION OF SOME PARAMETERS, STARTING WITH:
C MEAN VALUE <R>  FOR CALCULATING STRENGTH SR
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

C if delta R is very small the resonance is unrealistic
C we compare delta R with the central term which is of the order of kg2*mst*(m1+m2)/a1
         cterm=kg2*mst*(sm(1)+sm(2))/a1
         if(delr/cterm.lt.1.d-20)   then
           delr=0.d0
         endif
C if delta R is very small we cannot detect any resonance
C 1.d-11 is a very reasonable limit for relative variation
C it is really a small variation
      detec2=dabs((vmedio-vrmin)/vmedio)
      if(delr.eq.0.d0.or.detec2.lt.1.d-11) then
        delr=0.d0
        write(*,*)'Undetectable resonance at (e1,w1)=(',E1,',',P1G,')'
C        goto 111
      endif

C *******DEFINING EQUILIBRIUM POINTS (NO CLOSE ENCOUNTERS ALLOWED)******
C Previously, we start looking eq. points in normal mode.
C If NEQ is too big, we activate noisy mode to avoid false eq. points due to a noisy R(sigma) function.
      NOISE_mode = .false.
      GOTO 200
 201  NOISE_mode = .true.
      WRITE(*,*)'Noisy mode activated for (e1,w1)=(',E1,',',P1G,')'
 200  CONTINUE

C Some initialization before loop
      NEQ = 0
      CONT = 0
      IGUAL = .false.
      IMAX = ISIMAX
      IMIN = ISIMIN
 
C     LOOP:
      DO I=IMIN,IMAX
      
      IF ((NEQ.GT.10).AND.(.NOT.NOISE_mode)) GOTO 201
C     Assume beforehand there is no eq. in each point
      EQUILIBRIO = .false.
C     Define index + 1 and index + 2
      IF (I.LE.(IMAX-2)) THEN
        I_1 = I + 1
        I_2 = I + 2
      ELSE
        IF (I.EQ.IMAX-1) THEN
          I_1 = I + 1
          I_2 = IMIN
        ELSE
          I_1 = IMIN
          I_2 = IMIN + 1
        ENDIF
      ENDIF
C     Algorithm:
      IF ((RP(I).GT.RP(I_1)).AND.(.NOT.IGUAL)) THEN
        IF (RP(I_1).LT.RP(I_2)) THEN
C     Normal equilibrium point:
          EQUILIBRIO = .true.
          I_EQ = I_1
        ELSE
          IF (RP(I_1).EQ.RP(I_2)) IGUAL = .true.
        ENDIF
      ELSE
        IF (IGUAL.AND.(RP(I).EQ.RP(I_1))) THEN
          CONT = CONT + 1
        ELSE
          IF ((RP(I).LT.RP(I_1)).AND.(CONT.GT.0).AND.(I.GT.CONT+1)) THEN
C     Equilibrium point in a equal RP valley:
            EQUILIBRIO = .true.
            I_EQ = I - INT(AINT(DFLOAT(CONT/2)))
          ENDIF
          CONT = 0
          IGUAL = .false.
        ENDIF
      ENDIF

C     BORDER CASE:
      IF (IGUAL.AND.(I.EQ.IMAX)) THEN
        IAUX = IMIN
        DO WHILE (IGUAL)
          IAUX = IAUX + 1
          IGUAL = (RP(IAUX).EQ.RP(IAUX+1))
        ENDDO
        IF (RP(IAUX).LT.RP(IAUX+1)) THEN
          IF (IAUX.GT.CONT) THEN
            I_EQ = INT(AINT(DFLOAT( (IAUX - CONT)/2 )))
            IF (I_EQ.EQ.0) I_EQ = 1
          ELSE
            I_EQ = IMAX -  INT(AINT(DFLOAT( (CONT - IAUX)/2 )))
          ENDIF
          EQUILIBRIO = .true.
        ENDIF
      ENDIF

C     WRITE IN OUTPUT FILES THE 2 H-SURFACES:
      IF (EQUILIBRIO) THEN
        SMA_EQ = SMA(I_EQ)
        RP_EQ = RP(I_EQ)
        MD_EQ = MD(I_EQ)
C     Check for close encounters
        IF (MD_EQ.GT.RHTOL) THEN
        
C     Confirm equilibrium point (dismiss those due to noise in R function):
          CONFIRMADO = .True.
          IF (.NOT.NOISE_mode) GOTO 202
C     In the noisy mode, we check R value in 2 points at "ventana" degrees of
C     distance in sigma and compare them with the possible libration center
C     to decide if this is a real eq. point.
          ventana = 5
          I_vent = INT(ventana/DELTA_SMA)
          IF (CONT.GE.3) I_vent = I_vent + INT(CONT/2)
          Ieq_men_vent = I_EQ - I_vent
          IF (Ieq_men_vent.LT.IMIN) Ieq_men_vent = IMAX + Ieq_men_vent
          Ieq_mas_vent = I_EQ + I_vent
          IF (Ieq_mas_vent.GT.IMAX) Ieq_mas_vent = IMAX - Ieq_mas_vent
     *+ IMIN
          IF((RP(Ieq_men_vent).LT.RP_EQ).OR.(RP(Ieq_mas_vent).LT.RP_EQ))
     * THEN
            CONFIRMADO = .False.
          ENDIF

 202      CONTINUE
          IF (CONFIRMADO) THEN
C     Lets register the total hamiltonian:
          H_EQ = Ha - RP_EQ
            IF(MIN_flag.EQ.0) THEN
              WRITE(1,1004)E1,loper,SMA_EQ,RP_EQ
              WRITE(2,1004)E2,loper,SMA_EQ,RP_EQ
            ELSE
              WRITE(1,1005)E1,loper,SMA_EQ,RP_EQ,MD_EQ
              WRITE(2,1005)E2,loper,SMA_EQ,RP_EQ,MD_EQ
            ENDIF
          NEQ = NEQ + 1
          ENDIF
        ENDIF
      ENDIF
      

      ENDDO

C Show an approximately percentaje of the program completion:
      IF (PRIMERO) WRITE(*,*)'PROGRESS:'
C      porcentaje = 100*EXC + (loper - wMIN)/360
      porcentaje = (((DFLOAT(ind_e)-1)*Nw + ind_w)/(DFLOAT(Ne*Nw)))*100
      WRITE(*,1000)porcentaje,'% completed...'
      
C ancha1 is resonance full stable width  planet interior 1
C 111  ff= 2.d0*dsqrt(2.d0*delr/hii/kg2)
C      ancha1=k1*dsqrt(a1*(mst+sm(1)))/mst/sm(1)*ff
C      ancha1=2.d0*ancha1
C width planet exterior 2
C      ancha2=k2*dsqrt(a2*(mst+sm(2)))/mst/sm(2)*ff
C      ancha2=2.d0*ancha2

      IF (NEQ.EQ.0) THEN
        WRITE(*,*)'No equilibrium points found for (e1,dw)=(',E1,',',P1G
     *,')'
      ENDIF

C **********************************************************************
C ********END OF ASTEROID'S EXCENTRICITY AND PERI-CENTER LOOPS:*********
C **********************************************************************
      ENDDO
      GOTO 312
 321  IF (deltaE1.GT.deltaE2) THEN
        WRITE(*,*)'e2 does not exists for (e1,AM)=(',E1,',',AM,')'
      ELSE
        WRITE(*,*)'e1 does not exists for (e2,AM)=(',E2,',',AM,')'
      ENDIF
 312  CONTINUE
      ENDDO

      CLOSE(1)
      CLOSE(2)
 
C     FORMATS USED:
 1000 FORMAT(F7.3,A16)
 1001 FORMAT(F7.2,E16.8,2A10,E16.8)
 1002 FORMAT(F7.2,E16.8,2A10)
 1003 FORMAT(A7,F4.0,A10,A4,E16.8,A10)
 1004 FORMAT(F5.2,F6.0,F8.2,E20.12)
 1005 FORMAT(F5.2,F6.0,F8.2,E20.12,F14.8)
C 1110 WRITE(*,*)'More than 25 equilibrium points found! -> Program abort
C     *ed'
      GOTO 1111

C 1100 WRITE(*,*)'No equilibrium points found for (e1,w1)=(',E1,',',P1G,'
C     *)'
C      WRITE(*,*)'Ejecution aborted!'
C Useful endless loop to debug in Windows (if not, terminal shutdowns):
C 888  CONTINUE
C      GOTO 888
 1111 END
C     ********************** END OF MAIN PROGRAM ***********************

C     *** SUBROUTINES SECTION: ***
C ======================================================================
C SOLVING KEPLER
      SUBROUTINE SOLKEP(EX,M,E,NITER)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 M,MK
      TOLE=1.D-12
C      TWOPI=8.D0*DATAN(1.D0)
C      M=DMOD(M,TWOPI)
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
C     *** END OF SOURCE FILE ***
