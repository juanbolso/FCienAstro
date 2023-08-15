C     PROGRAM calculo de superficies 3D H. Auths.: Tabare Gallardo and Juan Pons

      IMPLICIT REAL*8 (A-H,J-Z)
      REAL INF
      INTEGER NEQ,NEQ_max,Ne,Nw,CONT,I_1,I_2,IAUX,IMAX,IMIN,ISIMAX,ISIMI
     *N,IND,I_EQ
      CHARACTER*250 RHTOL_str, PREC_str
      CHARACTER*50 APESO, nombre
      LOGICAL PRIMERO, EQUILIBRIO, ENCUENTRO, IGUAL, EXTRA
      DIMENSION SMA(1441),RP(1441),Htot(1441),ENC(1441),MINDIS(1441)
C     1441 = 360/DELTA_SMA where DELTA_SMA = 0.25 as better case.
      COMMON KPL,KPA,TETAR,A1,E1,J1,L1,P1,A2nom,E2,J2,L2,P2
      TWOPI = 8.0D0*DATAN(1.0D0)
      CERO  = 0.0D0
      UNO   = 1.0D0
      PI=TWOPI/2.D0
      G2R=PI/180.D0
      KGAUS=0.01720209895D0
      KG2=KGAUS**2
C     Boolean for adding extra info at the start of the output file:
      EXTRA = .false.

C     Terminal message when starting the program:
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'This program calculates the maximum H(a0, sigma) to'
      WRITE(*,*)'obtain the centers of resonant libration in the planar'
      WRITE(*,*)'case. It sweeps in eccentricity and perihelion'
      WRITE(*,*)'obtaining a 2D surface in the 3D space (e,varpi,sigma)'
      WRITE(*,*)'called H surface. This is a modification of Hasigma.f'
      WRITE(*,*)'program made by Juan Pons (juan.pons.93@gmail.com).'
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'H is calculated according to Gallardo 2020, CMDA 132 9'
      WRITE(*,*)'Author of the original program:'
      WRITE(*,*)'Tabar‚ Gallardo (gallardo@fisica.edu.uy).'
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'The PLANET is assumed with i = lonod = loper = 0ø but'
      WRITE(*,*)'no restriction for any eccentricity exists.'
      WRITE(*,*)'The input file must be planet.inp containing the'
      WRITE(*,*)'central body mass (in line 2), and this planet info:'
      WRITE(*,*)'   a_p   e_p   m_p     (in line 4)'
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'                                                      '
C     Read planet.inp input file:
      open(1,file="planet.inp",status="old")
      READ(1,*) APESO
      READ(1,*) STARM
      READ(1,*) APESO
      READ(1,*) A1,E1,SM1
      CLOSE(1)
      
      WRITE(*,*)'NOTE: This program searchs for ALL the equilibrium'
      WRITE(*,*)'points, calculating the hamiltonianïs local maximums.'
      WRITE(*,*)'                                                      '

      WRITE(*,*)'RESONANCE Kp,K ?'
      WRITE(*,*)'(ex. 2,3 for plutinos; 3,2 for hildas)'
      READ(*,*) KPL,KPA
C     KPL : K Planet
C     KPA : K Particle
C     In order to assure a valid pair of numbers:
      KPL = NINT(DABS(KPL))
      KPA = NINT(DABS(KPA))
C     Resonance order:
      ORD=KPL-KPA

C     THE PLANET IS 1
      J1=0.D0
      L1=0.D0
      P1=0.D0

C     COPLANAR CASE, SO:
      YNC = 0.D0
      lonod = 0.D0

C     planet mean motion
      mmp = DSQRT(kg2*(STARM+sm1)/a1**3)
      
C     PARTICLE IS 2
C     resonant particle mean motion
      mmr = (kpl/kpa)*mmp
C     A2 NOMINAL
      A2NOM = ((1+SM1/STARM)**(-1.D0/3.D0))*(kg2/mmr**2)**(1.D0/3.D0)
      WRITE(*,*)'Particle semi-major axis: ',A2NOM,'au'
C      A2NOM=A1/(DABS(NPQ/NP))**(2.D0/3.D0)
      A2 = A2NOM
      
C     Hamiltonian part that depends only on the semi-major axis:
      Ha = -KG2/2.D0/A2 - mmp*DSQRT(kg2*A2)*(kpl/kpa)

C     Hill radius to detect close encounters:
      RHS = A1*(SM1/(3.D0*STARM))**(1.D0/3.D0)
C     NUMBER OF RHILL ADMITTED TO CALCULATE DISTURBING FUNCTION
C     THIS ADMIT A FINE TUNNING ACCORDING TO THE ORBITAL INCLINATION
      WRITE(*,*)'                                                      '
      WRITE(*,*)'CRITERIA FOR ENCOUNTER: min{|r-rp|}<RHtol*RHill'
      WRITE(*,*)'Usual values for direct orbits: 3<RHtol<4 '
      WRITE(*,*)'Usual values for retrograde orbits: 2<RHtol<3'
      WRITE(*,*)'ENTER RHtol DESIRED VALUE PLEASE:'
      READ(*,*) RHTOL
      IF (RHTOL.LT.0) THEN
        RHTOL = DABS(RHTOL)
        WRITE(*,*)'Absolute value was considered.'
      ENDIF
      IF (RHTOL.GE.10.D0) THEN
        RHTOL = 9.99
        WRITE(*,*)'Set to max. value allowed: 9.99.'
      ENDIF
      WRITE(RHTOL_str, '(F4.2)') RHTOL   
      WRITE(*,*)'                                                      '
      
C     REGISTER MIN DISTANCES TO WRITE IN A FILE
      WRITE(*,*)'WANT TO REGISTER MIN. DISTANCES? (0=NO, 1=YES)'
      READ(*,*) MIN_flag
      IF ((MIN_flag.NE.0).AND.(MIN_flag.NE.1)) THEN
        MIN_flag = 0
        WRITE(*,*)'Invalid value. By default min. distances are not regi
     *stered.'
      ENDIF
      
C     The precision of calculation defines number of points in e,varpi y sigma:
      WRITE(*,*)'                                                      '
      WRITE(*,*)'PRECISION OF CALCULATION? (0=LOW, 1=HIGH)'
      READ(*,*) PREC

      PREC = INT(PREC)
      
      IF ((PREC.NE.0).AND.(PREC.NE.1)) THEN
        PREC = 0
        WRITE(*,*)'Invalid value. PRECISION set to 0.'
      ENDIF

C     WE SWEPT e AND varpi OF THE PARTICLE.
C     DELTA_SMA defines de precision of sigma variable.
      IF (PREC.EQ.0) THEN
        PREC_str = 'low'
        DELTA_SMA = 1.D0
        eMIN = 0.01D0
        eMAX = 0.96D0
        deltaE = 0.05D0
        wMIN = 0.D0
        wMAX = 355.D0
C        wMAX = 270.D0
        deltaW = 5.D0
C        deltaW = 90.D0
      ELSE
        PREC_str = 'high'
        DELTA_SMA = 0.25D0
        eMIN = 0.01D0
        eMAX = 0.99D0
        deltaE = 0.01D0
        wMIN = 0.D0
        wMAX = 359.D0
        deltaW = 1.D0
      ENDIF

C     From min,max,deltas, we obtain the number of indexes for main DO loops.
      Ne = NINT((eMAX - eMIN)/deltaE + 1)
      Nw = NINT((wMAX - wMIN)/deltaW + 1)
      
C     Max number of equilibrium points counter is initialized:
      NEQ_max = 1

C     NUMBER OF EVALUATIONS OF R(sigma) BETWEEN O AND 360
      ISIMIN = 1
      ISIMAX = NINT(360.D0/DELTA_SMA)

C     The output file name "nombre" has some info related to the particular run:
      nombre='HSUP_'//TRIM(PREC_str)//'_RHtol='//TRIM(RHTOL_str)//'.dat'
      OPEN(1,FILE=nombre,STATUS="UNKNOWN",ACCESS="APPEND")
      WRITE(1,*)'                                                      '
      WRITE(*,*)'                                                      '
      WRITE(*,*)'Writing in output file: '//nombre
      WRITE(*,*)'                                                      '
      IF (EXTRA) THEN
        WRITE(1,*)'Particle semi-major axis (au): ', A2
        WRITE(1,*)'Non-perturvative hamiltonian: ', Ha
        WRITE(1,*)'Hillïs radius (au): ', RHS
        WRITE(1,*)'                                                    '
      ENDIF
C     Write header:
      IF(MIN_flag.EQ.0) THEN
        WRITE(1,*)' e_ast   w_ast   sigma   H '
      ELSE
        WRITE(1,*)' e_ast   w_ast   sigma   H   mindis '
      ENDIF
      WRITE(1,*)'                                                      '
      
C     ******************************************************************
C     *********ASTEROID'S EXCENTRICITY AND PERI-CENTER LOOP:************
C     ******************************************************************
      DO ind_e = 1, INT(Ne)
      EXC = (ind_e-1)*deltaE + eMIN

      DO ind_w = 1, INT(Nw)
      loper = (ind_w-1)*deltaW + wMIN

C     Flag to check if we are in the first case of the (e,w) grid.
      PRIMERO = (ind_e.EQ.1).AND.(ind_w.EQ.1)
 
      E2 = EXC
      P2G = loper 
C     Pass the angles to radians:
      J2GRA=YNC
      J2=J2GRA*G2R
      L2G=LONOD
      L2=L2G*G2R
      P2=P2G*G2R

C     *************STARTS FUNCTION R(sigma) CALCULATION:****************
C     Before calculating R(sigma), lets initialize the ENCounter flags:
      ENC(:) = 0

C     Sigma equivalent definitions:
C      sigma = (p+q)*lambda_p - p*lambda - q*varpi
C      sigma = k1*lambda_p - k2*lambda + (k2-k1)*varpi

      DO ISI=ISIMIN,ISIMAX
        IND = ISI - ISIMIN + 1
C     ACRIT IS SIGMA
        ACRIT = DFLOAT(ISI-1)*DELTA_SMA       
        TETA = ACRIT - ORD*P2G

C     To left the angle in the (0,360) domain: (is this necessary?)
        IF (TETA.LT.0.D0) THEN
          TETA = DMOD(TETA,360.D0) + 360.D0
        ELSE
          TETA = DMOD(TETA,360.D0)
        ENDIF

C     TO RADIANS
        TETAR = TETA*G2R
        RTOT = CERO
C        DELTMIN=999999.0
        TINTEGRAL = 0.D0
C     INITIAL STEP 1 DEGREE
        PASO = 1.D0
        LA1G = -1.D0
C     LAMBDA PLANET FROM 0 TO 360*NP
        LIMIT = 360.D0*KPA
   88   LA1G = LA1G + PASO
   
C     CALCULATE DISTURBING FUNCTION A DISTANCE
        CALL RCALC(LA1G,RPER,DELTA)

C     REGISTER CLOSE ENCOUNTERS FOR THAT FIXED SIGMA:
        ENCUENTRO = (DELTA.LT.(RHTOL*RHS))
        IF (ENCUENTRO) ENC(IND) = 1
C     CALCULATE MIN DISTANCE NORMALIZED IN RHILL:
	IF (LA1G.EQ.0.D0) THEN
          MINDIS(IND) = DELTA/RHS
	ELSE
	  IF (DELTA.LT.MINDIS(IND)*RHS) MINDIS(IND) = DELTA/RHS
	ENDIF
C     WE GOT 360*NP, WE DO NOT COMPUTE THE INTEGRAL (IS THE SAME THAN 0)
        IF (LA1G.GE.LIMIT) THEN
          GOTO 99
        ENDIF
        TINTEGRAL=TINTEGRAL+RPER*PASO
        GOTO 88
C     CRUDE INTEGRAL
  99    RTOT = TINTEGRAL/(360.D0*KPA)
  
C     NOW MULTIPLY BY MASS AND CONSTANT OF GRAVITATION
        RP(IND)= RTOT*SM1*KG2
C         RPlog(IND) = DLOG10(DABS(RP(IND)))

C     To left the angle in the (0,360) domain...:
        IF (ACRIT.LT.0.D0) THEN
          SMA(IND) = DMOD(ACRIT,360.D0) + 360.D0
        ELSE
          SMA(IND) = DMOD(ACRIT,360.D0)
        ENDIF
      
C     ****************END OF R(sigma) CALCULATION***********************
      ENDDO

C     **********STARTS THE EQUILIBRIUM POINTS SEARCH:*******************
C     Some initialization first
      CONT = 0
      IGUAL = .false.
      IMAX = ISIMAX
      IMIN = ISIMIN
C     LOOP:
      DO I=IMIN,IMAX
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
          SMA_EQ = SMA(I_EQ)
          RP_EQ = RP(I_EQ)
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
            SMA_EQ = SMA(I_EQ)
            RP_EQ = RP(I_EQ)
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
          SMA_EQ = SMA(I_EQ)
          RP_EQ = RP(I_EQ)
        ENDIF
      ENDIF

C     WRITE IN OUTPUT FILES THE 2 H-SURFACES:
      IF (EQUILIBRIO) THEN
        IF (INT(ENC(I_EQ)).EQ.0) THEN
C     Lets register the total hamiltonian:
        H_EQ = Ha - RP_EQ
          IF(MIN_flag.EQ.0) THEN
            WRITE(1,1001)E2,P2G,SMA_EQ,H_EQ
          ELSE
            WRITE(1,1002)E2,P2G,SMA_EQ,H_EQ,MINDIS(I_EQ)
          ENDIF
        ENDIF
      NEQ = NEQ + 1
      ENDIF

      ENDDO

C     Show an approximately percentaje of the program completion:
      IF (PRIMERO) WRITE(*,*)'PROGRESS:'
C      porcentaje = 100*EXC + (loper - wMIN)/360
      porcentaje = (((DFLOAT(ind_e)-1)*Nw + ind_w)/(DFLOAT(Ne*Nw)))*100
      WRITE(*,1000)porcentaje,'% completed...'

C     ******************************************************************
C     ****END OF ASTEROID'S EXCENTRICITY AND PERI-CENTER LOOPS:*********
C     ******************************************************************
      ENDDO
      ENDDO

      CLOSE(1)
      
      WRITE(*,*)'End of the execution! :)'
      GOTO 1111
C      GOTO 888

C     FORMATS USED:
 1000 FORMAT(F7.3,A16)
 1001 FORMAT(F5.2,F6.0,F8.2,E20.12)
 1002 FORMAT(F5.2,F6.0,F8.2,E20.12,F14.8)
 
C 1100 WRITE(*,*)'No equilibrium points found -> Program aborted'
C     Useful endless loop to debug in Windows (if not, terminal shutdowns):
 888  CONTINUE
      GOTO 888
 1111 END
 
C     ******************************************************************
C     ******************************************************************
C     ******************END OF THE MAIN PROGRAM*************************
C     ******************************************************************
C     ******************************************************************

C     ***************** SUBROUTINES SECTION: ***************************

      SUBROUTINE RCALC(LA1G,RPER,DELTA)
      IMPLICIT REAL*8 (A-H,J-Z)
      COMMON KPL,KPA,TETAR,A1,E1,J1,L1,P1,A2nom,E2,J2,L2,P2
      TWOPI = 8.0D0*DATAN(1.0D0)
      CERO  = 0.0D0
      UNO   = 1.0D0
      PI=TWOPI/2.D0
      G2R=PI/180.D0
      ERROR=1.D-12

      a2=a2nom

      LA1=LA1G*g2r
C     MEAN ALY PLANET assuming loper=0
      AM1=LA1
 200  IF (AM1.GT.TWOPI) THEN
        AM1=AM1-TWOPI
        GOTO 200
      END IF
 201  IF (AM1.LT.CERO) THEN
        AM1=AM1+TWOPI
        GOTO 201
      END IF
C     PLANET IN excentric orbit
C     CRUDE METHOD FOR SOLVING KEPLER FOR planet (ECCENTRIC)
      AEX=AM1
 210  AEXN=AM1+E1*DSIN(AEX)
      IF(DABS(AEXN-AEX).GT.ERROR) THEN
        AEX=AEXN
        GOTO 210
      ENDIF
      AEX1=AEXN
C     TRUE ANOMALY
      AVE1=DACOS((DCOS(AEX1)-E1)/(UNO-E1*DCOS(AEX1)))
      IF(AM1.GT.PI) AVE1=TWOPI-AVE1
C     HELIOCENTRIC DISTANCE
      R1=A1*(UNO-E1*DCOS(AEX1))
C     RADIUS VECTOR FOR planet
      X1=R1*(DCOS(L1)*DCOS(P1-L1+AVE1)-DSIN(L1)*DSIN(P1-L1+AVE1)
     **DCOS(J1))
      Y1=R1*(DSIN(L1)*DCOS(P1-L1+AVE1)+DCOS(L1)*DSIN(P1-L1+AVE1)
     **DCOS(J1))
      Z1=R1*DSIN(P1-L1+AVE1)*DSIN(J1)

C     GIVEN LAMBDA PLANET CALCULATE LAMBDA PARTICLE
C      LA2 = (KPL*LA1-TETAR)/KPA
      LA2 = (KPL*LA1+TETAR)/KPA
      
C     MEAN ANOMALY PARTICLE
      AM2=LA2-P2
 202  IF (AM2.GT.TWOPI) THEN
        AM2=AM2-TWOPI
        GOTO 202
      END IF
 203  IF (AM2.LT.CERO) THEN
        AM2=AM2+TWOPI
        GOTO 203
      END IF
C     CRUDE METHOD FOR SOLVING KEPLER FOR PARTICLE (ECCENTRIC)
      AEX=AM2
 211  AEXN=AM2+E2*DSIN(AEX)
      IF(DABS(AEXN-AEX).GT.ERROR) THEN
        AEX=AEXN
        GOTO 211
      ENDIF
      AEX2=AEXN
C     TRUE ANOMALY
      AVE2=DACOS((DCOS(AEX2)-E2)/(UNO-E2*DCOS(AEX2)))
      IF(AM2.GT.PI) AVE2=TWOPI-AVE2
C     HELIOCENTRIC DISTANCE
      R2=A2*(UNO-E2*DCOS(AEX2))
C     RADIUS VECTOR FOR PARTICLE
      X2=R2*(DCOS(L2)*DCOS(P2-L2+AVE2)-DSIN(L2)*DSIN(P2-L2+AVE2)
     **DCOS(J2))
      Y2=R2*(DSIN(L2)*DCOS(P2-L2+AVE2)+DCOS(L2)*DSIN(P2-L2+AVE2)
     **DCOS(J2))
      Z2=R2*DSIN(P2-L2+AVE2)*DSIN(J2)
C     DISTANCE
      DELTA=DSQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
C     SCALAR PRODUCT R1*R2
      R1R2=X1*X2+Y1*Y2+Z1*Z2
C     DIRECT PART
      RDIR=UNO/DELTA
C     INDIRECT PART
      RIND=-R1R2/R1**3
C     DISTURBING FUNCTION
      RPER=RDIR+RIND
      RETURN
      END
C     **********END OF SUBROUTINE RCALC(LA1G,RPER,DELTA)****************

C     ************** END OF SUBROUTINES SECTION ************************

C     ******************************************************************
C     ********************** END OF FILE *******************************
C     ******************************************************************
