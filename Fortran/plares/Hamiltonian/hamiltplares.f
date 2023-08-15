c program hamiltplares.f  - gallardo, beauge, giuppone 2020
c gallardo@fisica.edu.uy
c similar to atlaspr.f but restricted to 2 planets
c and including calculation of Hamiltonian
C CALCULATES HAMILTONIAN (level curves)
C CALCULATES RESONANT DISTURBING FUNCTION R(sigma), THE EQUILIBRIUM POINTS,
C LIBRATION PERIODS, RESONANCE STRENGTH AND WIDTH OF THE RESONANCE IN AU
C PLANETS AND STAR GIVEN BY PLASYS.INP in JACOBI SYSTEM
      IMPLICIT REAL*8 (A-H,J-Z)
      DIMENSION SE(2),SM(2),EX(2),YN(2),LN(2),LP(2),RH(10)
      DIMENSION RP(400),DERIV1(400),DERIV2(400),TLIB(400),MD(400),
     *IRH(400),sma(400)
      dimension h(150,400)
C THIS 50000 DIMENSION IS FOR STORING DATA PARTICLE, CAN BE MODIFIED
C MAX NUMBER OF POINTS IN THE INTEGRAL   = 1000*MAX(K1,K2)
      DIMENSION VX2(50000),VY2(50000),VZ2(50000),VR2(50000),VLA2(50000)
C      DIMENSION VXP2(50000),VYP2(50000),VZP2(50000)
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
      WRITE(*,*)'               hamiltplares      '
      WRITE(*,*)'  calculates level curves of the Hamiltonian'
      WRITE(*,*)'        for 2 planets in resonance'
      WRITE(*,*)'Reference: Gallardo, Beauge, Giuppone 2020'
      WRITE(*,*)'               version January 26, 2022'
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'             I assume the archive plasysj.inp is ready'
C      WRITE(*,*)'press ENTER'
C      READ(*,*)
C INPUT
      OPEN(1,FILE='plasysj.inp',STATUS='old')
      READ (1,*)  ZCOMMENT
C MASS CENTRAL STAR
      READ(1,*) MST
      READ (1,*)  ZCOMMENT
C READ 2 PLANETS A E i lonod, loper, MASS
      DO  I=1,2
        READ (1,*) SE(I),EX(I),YN(I),LN(I),LP(I), SM(I)
      ENDDO
      READ (1,*)  ZCOMMENT
      CLOSE (1)
C NUMBER OF RHILL ADMITTED TO CALCULATE DISTURBING FUNCTION
C THIS ADMIT A FINE TUNNING ACCORDING TO THE ORBITAL INCLINATION
C FOR DIRECT ORBITS 3<RHTOL<4
C FOR RETROGRADE ORBITS 2<RHTOL<3
C PUT 0 FOR OBTAINING MAXIMUM TEORETHICAL WIDTH (SEPARATRICES)
      RHTOL=0.D0
      
      M3=MST+SM(1)+SM(2)
      
      
  10    WRITE(*,*)'RESONANCE K2:K1 ? (K2 >= K1)'
      WRITE(*,*)'(ex. 3,2 for hildas or plutinos, <50,50)'
      READ(*,*) K2,K1
      IF(K2.GE.K1) THEN
        WRITE(*,*)'sigma =',INT(K1),'*lambda1',INT(-K2),'*lambda2 +',
     *INT(K2-K1),'*varpi1'
C        ELSE
C        WRITE(*,*)'sigma =',INT(K1),'*lambda',INT(-K2),'*lambda_pla',
C     *INT(K2-K1),'*varpi'
      ELSE
        WRITE(*,*)'please... K2 cant be < K1'
        goto 10
      ENDIF
C MAXIMUM FOR CALCULATION OF THE INTEGRAL WITH ENOUGH PRECISION
      MAXFAC=K2
      WRITE(*,*) '   '

C FACTORS TO BE USED
      FA01=SM(1)/(MST+SM(1))

      A1=SE(1)
C MEAN MOTION OF M1 WITH GIVEN A1
      N1=DSQRT(KG2*(MST+SM(1))/A1**3)
C RESONANT N2
      N2=K1*N1/K2
C SEMIMAJOR AXIS A2
      A2=(KG2*(MST+SM(1)+SM(2))/N2**2)**(1.D0/3.D0)

C EXTERIOR PLANET IS 2
C HILL
      RH(2)=A2*(SM(2)/(MST+SM(2))/3.D0)**(1.D0/3.D0)

      E2=EX(2)
      B2=A2*DSQRT(UNO-E2*E2)
      J2G=YN(2)
      J2=J2G*G2R
      L2G=LN(2)
      P2G=LP(2)
      L2=L2G*G2R
      P2=P2G*G2R
      AR2=P2-L2
C CONSTANTS L1,M1,N1 FOR PLANET 2 ROY PAGE  93
      L12=DCOS(L2)*DCOS(AR2)-DSIN(L2)*DSIN(AR2)*DCOS(J2)
      M12=DSIN(L2)*DCOS(AR2)+DCOS(L2)*DSIN(AR2)*DCOS(J2)
      N12=DSIN(AR2)*DSIN(J2)

C CONSTANTS L2,M2,N2 FOR PLANET 2, ROY PAGE  93
      L22=-DCOS(L2)*DSIN(AR2)-DSIN(L2)*DCOS(AR2)*DCOS(J2)
      M22=-DSIN(L2)*DSIN(AR2)+DCOS(L2)*DCOS(AR2)*DCOS(J2)
      N22=DCOS(AR2)*DSIN(J2)



C THE INTERIOR PLANET IS 1
C HILL
      RH(1)=A1*(SM(1)/(MST+SM(1))/3.D0)**(1.D0/3.D0)
      E1=EX(1)
      B1=A1*DSQRT(UNO-E1*E1)
      J1G=YN(1)
      L1G=LN(1)
      P1G=LP(1)
      J1=J1G*G2R
      L1=L1G*G2R
      P1=P1G*G2R
      AR1=P1-L1
C CONSTANTS L1,M1,N1 FOR PLANET 1 ROY PAGE  93
      L11=DCOS(L1)*DCOS(AR1)-DSIN(L1)*DSIN(AR1)*DCOS(J1)
      M11=DSIN(L1)*DCOS(AR1)+DCOS(L1)*DSIN(AR1)*DCOS(J1)
      N11=DSIN(AR1)*DSIN(J1)
C CONSTANTS L2,M2,N2 FOR PLANET 1  ROY PAGE  93
      L21=-DCOS(L1)*DSIN(AR1)-DSIN(L1)*DCOS(AR1)*DCOS(J1)
      M21=-DSIN(L1)*DSIN(AR1)+DCOS(L1)*DCOS(AR1)*DCOS(J1)
      N21=DCOS(AR1)*DSIN(J1)

C SECOND DERIVATIVE HII
      BETA1=SM(1)*MST/(SM(1)+MST)
      BETA2=SM(2)*(MST+SM(1))/(SM(2)+SM(1)+MST)
      HII=3.D0*(K1/A1)**2/BETA1 +  3.D0*(K2/A2)**2/BETA2
      
      
C RHILL mutual
      rhmass=  (sm(1)+ sm(2))/mst
      rhs=(a2+a1)/2.d0*(rhmass/3.d0)**(1.d0/3.d0)

      
C NUMBER OF EVALUATIONS OF R(sigma) BETWEEN O AND 360 DEGREES
      ISIMAX=360
C R(360)=R(0)
C STEPS OF THE NUMERICAL INTEGRATION FROM 0 TO 2PI*K2 IN lambda PARTICLE
C CAN BE AS SMALL AS 100*INT(MAXFAC), AT YOUR OWN RISK
      IPASOS=1000*INT(MAXFAC)
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
C          XP2=(B2*L22*COSE-A2*L12*SINE)*N2*A2/R2
C          YP2=(B2*M22*COSE-A2*M12*SINE)*N2*A2/R2
C          ZP2=(B2*N22*COSE-A2*N12*SINE)*N2*A2/R2


C VECTOR DATA RPLA2
          VR2(I2)=R2

C VECTOR DATA XYZ
          VX2(I2)=X2
          VY2(I2)=Y2
          VZ2(I2)=Z2
C VECTOR DATA XPYPZP
C          VXP2(I2)=XP2
C          VYP2(I2)=YP2
C          VZP2(I2)=ZP2

          
          
        ENDDO
C I HAVE CALCULATED ALL THE POSITIONS FOR PLA2 THAT I WILL USE LATER
C ================================================================
C NOW DEFINE SIGMA1
      DO ISI=1,ISIMAX
C ACRIT IS SIGMA1 FROM 0 TO 359
C CLOSE ENCOUNTER INDICATOR
        ACE(ISI)='         '
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
C NUMBER OF ENCOUNTERS INSIDE THE HILL SPHERE
        IRH(ISI)=0
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
C          XP2=VXP2(ILAMBDA2)
C          YP2=VYP2(ILAMBDA2)
C          ZP2=VZP2(ILAMBDA2)


          
          
          
C GIVEN LAMBDA PLA EXTERIOR 2 CALCULATE LA1  PLANET  INT 1
          LA1=(K2*LA2+TETAR)/K1
C MEAN ANOMALY PLA 1
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
C VELOCITIES
C          XP1=(B1*L21*COSE-A1*L11*SINE)*N1*A1/R1
C          YP1=(B1*M21*COSE-A1*M11*SINE)*N1*A1/R1
C          ZP1=(B1*N21*COSE-A1*N11*SINE)*N1*A1/R1

C VECTOR R02
          X02=FA01*X1+X2
          Y02=FA01*Y1+Y2
          Z02=FA01*Z1+Z2
C VECTOR R12
          X12=X02-X1
          Y12=Y02-Y1
          Z12=Z02-Z1

C MUTUAL DISTANCE
      DELTA02=DSQRT(X02**2+Y02**2+Z02**2)
      DELTA12=DSQRT(X12**2+Y12**2+Z12**2)


      
C REGISTER CLOSE ENCOUNTERS AND NUMBER
          IF(DELTA12.LT.RHS*RHTOL) THEN
            ACE(ISI)='CLOSE ENC'
            IRH(ISI)=IRH(ISI)+1
          ENDIF
C REGISTER MINIMUM DISTANCE PLANET-PARTICLE IN RHILL
          DELRH=DELTA12/RHS
          IF(DELRH.LT.MINDIS)THEN
            MINDIS=DELRH
          ENDIF
C CALCULATION DISTURBING FUNCTION

          RPER=MST/DELTA02 + SM(1)/DELTA12 - (MST+SM(1))/R2
          
C SUMATION FOR A CRUDE INTEGRAL
          RTOT=RTOT+RPER
        ENDDO
          RTOT=RTOT*KG2*SM(2)
        
C END OF CALCULATION OF THE INTEGRAL
        RTOT=RTOT/DFLOAT(IPASOS)
        RP(ISI)=RTOT
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
C        IF(RP(I2).GT.VRMAX) THEN
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
C        IF(MD(I).GT.0.5D0) THEN
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
              
C LIBRATION PERIOD
      TLIB(IEXTREMO)=TWOPI/DSQRT(HII*DABS(DERIV2(IEXTREMO
     *)))/365.25D0
              ELSE
C UNSTABLE
              IF(RP(ISIME1).GT.RP(I))THEN
                IEXTREMO=ISIME1
              ENDIF
            AEQ(IEXTREMO)=' UNSTABLE'
          ENDIF
C        ENDIF
      ENDIF
      ENDDO
C ANCHA1 IS RESONANCE FULL stable WIDTH  PLANET INTERIOR 1
 111  ANCHA1=2.D0*K1/MST/SM(1)*DSQRT(2.D0*DELR*A1*(SM(1)+MST)/HII/KG2)
      ANCHA1=2.D0*ANCHA1
C WIDTH PLANET EXTERIOR 2
      ANCHA2=FA01*MST/SM(2)*dsqrt((SM(2)+SM(1)+MST)*a2/(SM(1)+MST)/a1)
      ancha2=K2/K1*ancha2*ancha1
      
C FILE OUTPUT
      OPEN(2,FILE="resonance.dat",STATUS="UNKNOWN",ACCESS="APPEND")
      WRITE(2,*)'Star mass =  ',mst
      WRITE(2,*)'Planet masses=  ',sm(1),sm(2)
      WRITE(2,*)'Planet a1, a2 =  ',a1,a2,' Jacobi system'
      WRITE(2,*)'Planet e1, e2 =  ',e1,e2
      WRITE(2,*)'Planet i1, i2 =  ',J1G,J2G
      WRITE(2,*)'Planet LN1, LN2 =  ',L1G,L2G
      WRITE(2,*)'Planet LP1, LP2 =  ',P1G,P2G

      
      WRITE(2,*)'Resonance =  ',INT(K2),':',INT(K1)
      WRITE(2,*)'<R> = ',VMEDIO, '  [M_sun,au,days]'
      WRITE(2,*)'Strength (<R>-Rmin) = ',(VMEDIO-VRMIN)
      WRITE(2,*)'Rmax-Rmin = ',DELR
      WRITE(2,*)'Close enc. dist. = ',rhtol, ' RHill'
      WRITE(2,*)'stable full width Da1, Da2= ',ANCHA1, ANCHA2
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
        WRITE(2,*)'sigma =',INT(K1),'*lambda1',INT(-K2),'*lambda2 +',
     *INT(K2-K1),'*varpi1'

      WRITE(2,*)'        '
      CLOSE(2)
      
C here is the interesting part********************************************
C CONSTANT ACTION I2
      ACI2=(MST+SM(1))*SM(2)*KGAUS*DSQRT(A2/M3)
      ACI2=ACI2+K2*MST*SM(1)*KGAUS/K1*DSQRT(A1/(MST+SM(1)))
C FACTORS IN HAMILTONIAN
      HA1=KG2*(MST*SM(1))**3.d0/2.D0/(MST+SM(1))
      HA2=KG2*((MST+SM(1))*SM(2))**3.d0/2.D0/M3



C HAMILTONIAN
      aini=a2-ancha2*0.5d0
      deltaa=ancha2*1.2d0
      suma=0.d0
      do ia=1,101
        a=aini+deltaa/100.d0*dfloat(ia-1)
C L2
        L2= (MST+SM(1))*SM(2)*KGAUS*DSQRT(A/M3)
C ACTION I1
        ACI1= (ACI2-L2)/K2
        do is=1,360
C HAMILTONIAN
          h(ia,is)=-HA1/(K1*ACI1)**2.d0-HA2/(ACI2-K2*ACI1)**2.d0
C multiply by gauss2 and substract R
          h(ia,is)=h(ia,is)*kg2-rp(is)
          suma=suma+h(ia,is)
        enddo
      enddo
C HAMILTONIAN
      OPEN(2,FILE="hamiltplares.dat",STATUS="UNKNOWN",ACCESS="APPEND")
      aini=a2-ancha2*0.5d0
      deltaa=ancha2*1.2d0
      do ia=1,101
        a=aini+deltaa/100.d0*dfloat(ia-1)
        do is=1,360
          h(ia,is)=h(ia,is)-SUMA/DFLOAT(360*101)
        
          write(2,202)a,is-1,h(ia,is)
        enddo
        write(2,*)'        '
      enddo
      CLOSE (2)
      
      
  202 format(f10.6,i6,E30.20)
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
      IF(DEX.GT.TOLE)GOTO 100
      RETURN
      END

