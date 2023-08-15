C PROGRAM SUPERATLASv2.F - TABARE GALLARDO, January 2020
C CALCULATES AN ATLAS OF RESONANCE'S STRENGTHS, WIDTHS, LIBRATION CENTERS
C AND PERIODS
C IS THE OLD ATLAS.F SUPERUPDATED
C CONSIDERING PLANETARY ECCENTRICITIES
C REFERENCE: GALLARDO 2019, CMDA
C gallardo@fisica.edu.uy
      IMPLICIT REAL*8 (A-H,O-Z)
C PLANETARY SEMIMAJOR AXES, ECCENTRicities, MASSes, RADIUS HILL
      DIMENSION XSE(10),XEC(10),XMA(10),RH(10)
C 30000 IS THE MAXIMUM NUMBER OF RESONANCES TO CALCULATE BY ONCE
C CAN BE CHANGED
      DIMENSION NPLA(30000),NPQ(30000),NP(30000)
      DIMENSION A2(30000),NCPLA(30000),NCPQ(30000),NCP(30000)
      COMMON XM0,XSE,XEC,XMA,RH
      CHARACTER*49 cabezal1
      CHARACTER*55 cabezal2
      CHARACTER*110 cabezal4
      CHARACTER*48 cabsigma
      character*40 zcomment
      cabezal1='pla  kp:k    a(au)    e     i     w  ln      <R> '
      cabezal2='     <R>-R_min  width(au)  sigma_0          periods(yr)'
      cabezal4=cabezal1//cabezal2
      cabsigma='   sigma = k*lambda - kp*lambda_p + (kp-k)*varpi'
C INPUT
      OPEN(1,FILE='plasys.inp',STATUS='old')
      READ (1,*)  ZCOMMENT
C MASS CENTRAL STAR
      READ(1,*) XM0
      READ (1,*)  ZCOMMENT
C READ PLANETS A E MASS
      DO  I=1,10
        READ (1,*,ERR=199) XSE(I),XEC(I),XMA(I)
      ENDDO
      READ (1,*)  ZCOMMENT
  199 NUMPLA=I-1
C THERE ARE NUMPLA PLANETS
      close (1)
C PLANETARY HILL RADIUS
      DO I=1,NUMPLA
        RH(I)=XSE(I)*(XMA(I)/(XM0+XMA(I))/3.D0)**(1.D0/3.D0)
      ENDDO
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'SUPERATLASv2 OF MEAN MOTIONS RESONANCES kp:k'
      WRITE(*,*)'This program calculates all resonance''s strengths SR,'
      WRITE(*,*)'widths in au, libration centers and periods'
      WRITE(*,*)'in a given semimajor axis interval with a set '
      WRITE(*,*)'of given planets described in plasys.inp'
      WRITE(*,*)'Output: superatlasv2.out'
      WRITE(*,*)'Gallardo 2020, CMDA'
      WRITE(*,*)'Tabare Gallardo                gallardo@fisica.edu.uy '
      WRITE(*,*)'------------------------------------------------------'
      WRITE(*,*)'         '
      WRITE(*,*)'Do you need to know which resonances dominate in some '
      WRITE(*,*)'interval of semimajor axis? --> Go ahead!    '
      WRITE(*,*)'         '
      WRITE(*,*)'--> For resonances  kp:k indicate:'
C MAX |P+Q|
      WRITE(*,*)'max kp ?  (less than 50)'
      READ(*,*)NLIMPQ
C MAX DEGREE
      WRITE(*,*)'max k ?  (less than 50)'
      READ(*,*)IPF
C MAX ORDER
      WRITE(*,*)'max order q=|kp-k| ?  (less than 50)'
      READ(*,*)IQF
C LIMITS IN SEMIMAJOR AXIS
      WRITE(*,*)'interval to be analyzed in au: a_min,a_max?'
      WRITE(*,*)'(example: 0.3,25.4)'
      READ(*,*)A2MIN,A2MAX
C PLANETS TO CONSIDER
      WRITE(*,*)'you have defined ',numpla,' planets'
      WRITE(*,*)'now indicate first and last planet to be considered'
      WRITE(*,*)'(example: 1,8)'
      READ(*,*)IFIRSTP,ILASTP
C READ ELEMENTS TEST PARTICLE
      write(*,*) '  '
      WRITE(*,*)'indicate orbital elements of test particle e,i,w,ln'
      WRITE(*,*)'example: 0.2,15,60,0'
      read (*,*) exc,ync,arg,ome
C e,i,argper and ome is longitude of the node asuming that the perihelion
C of the planet is at the origin
C MASS OF THE PARTICLE =0
      XM2=0.0D0
      ZM2=XM0+XM2
C COUNTER I
      I=0
C FOR EACH PLANET.....
      DO 200 JPLA=IFIRSTP,ILASTP
      XM1=XMA(JPLA)
      ZM1=XM0+XM1
      A1=XSE(JPLA)
C FOR EACH ORDER...
      DO IQ=0,IQF
        Q=IQ*1.D0
C FOR EACH DEGREE FOR INTERIOR RESONANCES....
        DO IP=1,IPF
          P=IP*1.D0
C LIMIT |P+Q|
          isum=abs(iq+ip)
          if(isum.le.nlimpq) then
C CALCULATE SEMIMAJOR AXIS semc
          semC=((P/(P+Q))**2*ZM2/ZM1)**(1.D0/3.D0)*A1
            IF (semC.LE.A2MAX.AND.semC.GE.A2MIN)  THEN
              i=i+1
              NPLA(I)=JPLA
              NPQ(I)=INT(ABS(P+Q))
              NP(I)=INT(ABS(P))
              A2(I)=semC
            endif
          endif
        ENDDO
      ENDDO
C NOW THE SAME FOR EXTERIOR RESONANCES
      DO IQ=0,IQF
        Q=IQ*1.D0
        DO IP=IQ,IPF-1
          P=-(IP+1)*1.D0
C LIMIT |P+Q|
          isum=abs(iq-ip-1)
          if(isum.le.nlimpq) then
          semC=((P/(P+Q))**2*ZM2/ZM1)**(1.D0/3.D0)*A1
            IF (semC.LE.A2MAX.AND.semC.GE.A2MIN)  THEN
              I=I+1
              NPLA(I)=JPLA
              NPQ(I)=INT(ABS(P+Q))
              NP(I)=INT(ABS(P))
              A2(I)=semC
            endif
          endif
        ENDDO
      ENDDO
 200  CONTINUE
C NUMBER OF POSSIBLE (redundant) RESONANCES, MUST BE LESS THAN 30000
      IDAT=I
      if(idat.ge.30000) then
        write(*,*)'error: too much resonances'
        stop
      end if
C ELIMINATE REPEATED, EXAMPLE: 4:2 IS THE SAME AS 2:1
      DO I=1,IDAT
        DO J=I+1,IDAT
          IF(NPLA(J).EQ.NPLA(I)) THEN
            IF (A2(J).EQ.A2(I).AND.NPQ(J).GE.NPQ(I)) THEN
              A2(J)=0.D0
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      NCONT=0
      DO I=1,IDAT
        IF (A2(I).NE.0.D0) THEN
C LIMIT |P+Q|
C          IF (NPQ(I).LE.NLIMPQ)  THEN
            NCONT=NCONT+1
            NCPQ(NCONT)=NPQ(I)
            NCP(NCONT)=NP(I)
            NCPLA(NCONT)=NPLA(I)
C            SEMI(NCONT)=A2(I)
C          ENDIF
        ENDIF
      ENDDO
C RESONANCES REPEATED WERE ELIMINATED
      WRITE(*,*)'I WILL CALCULATE ',NCONT,' RESONANCES'
C OK. WE NEED TO CALCULATE  NCONT RESONANCES
      OPEN(1,FILE='superatlasv2.out',STATUS='UNKNOWN',ACCESS='APPEND')
      WRITE(1,*)'I WILL CALCULATE ',NCONT,' RESONANCES',cabsigma
      write(1,*) cabezal4
      DO I=1,NCONT

C IF YOU WANT TO CALCULATE WIDTHS IN (a,e) YOU CAN USE A SIMPLE LOOP
C LIKE THIS ONE
C        exc=-0.02d0
C        do iex=1,50
C          exc=exc+0.02d0

C VPQ IS k_p AND VP IS k
        VPQ=DFLOAT(NCPQ(I))
        VP=DFLOAT(NCP(I))
        write(*,*) I
        CALL FUERZAW(VPQ,VP,NCPLA(I),EXC,YNC,ARG,OME)

C END OF THE LOOP FOR CALCULATING MAPS (a,e)
C        enddo
      ENDDO
      CLOSE(1)
      STOP
      END
C --------------------------------------------------------------------
C BASED ON RESONALYZER.F  - TABARE GALLARDO, July 26, 2019
      SUBROUTINE FUERZAW(KP,KK,INPLA,EXC,YNC,ARGPER,OME)
      IMPLICIT REAL*8 (A-H,J-Z)
      DIMENSION XSE(10),XEC(10),XMA(10),RH(10)
      DIMENSION RP(400),DERIV1(400),DERIV2(400),MD(400),IRH(400)
      DIMENSION TLIB(4),ISIGEQ(4)
C THIS 50000 DIMENSION IS FOR STORING DATA PARTICLE, CAN BE MODIFIED
C MAX NUMBER OF POINTS IN THE INTEGRAL   = 1000*MAX(KK,KP)
      DIMENSION VX2(50000),VY2(50000),VZ2(50000),VR2(50000),VLA2(50000)

      COMMON XM0,XSE,XEC,XMA,RH


      TWOPI = 8.0D0*DATAN(1.0D0)
      CERO  = 0.0D0
      UNO   = 1.0D0
      PI=TWOPI/2.D0
      G2R=PI/180.D0
      KGAUS=0.01720209895D0
      KG2=KGAUS**2
C MASS OF THE STAR IN SOLAR MASSES
      MST=XM0
C NUMBER OF RHILL ADMITTED TO CALCULATE DISTURBING FUNCTION
C THIS ADMIT A FINE TUNNING ACCORDING TO THE ORBITAL INCLINATION
C FOR DIRECT ORBITS 3<RHTOL<4
C FOR RETROGRADE ORBITS 2<RHTOL<3
      RHTOL=3.D0
C MAXIMUM FOR CALCULATION OF THE INTEGRAL WITH ENOUGH PRECISION
      MAXFAC=KP
      IF(KK.GT.KP) THEN
        MAXFAC=KK
      ENDIF
C INC AND ARG PER IN DEGREES
C sigma is invariant to changes in aries then I take aries at the node
C so longnode(particle)=0
C The planet is in ECCENTRIC ORBIT WITH i=0 orbit AND longper=0
C so I take \omega=\Omega=0 for the planet
C LON NODE ASTEROID
      LONOD=OME
      LOPER=ARGPER+LONOD
      IF (LOPER.GT.360.0) THEN
            LOPER=LOPER-360.0d0
      ENDIF

C THE PLANET IS 1
      A1=XSE(INPLA)
      E1=XEC(INPLA)
C INCLINATION, LONNODE, LONPER
      J1G=0.D0
      L1G=0.D0
      P1G=0.D0
C TO RADIANS...
      J1=J1G*G2R
      L1=L1G*G2R
      P1=P1G*G2R
C PARTICLE IS 2
      A2=A1*((DABS(KK/KP))**2*MST/(MST+XMA(INPLA)))**(1.D0/3.D0)
      SALSEM=A2
      E2=EXC
C INCLINATION
      J2GRA=YNC
      J2=J2GRA*G2R
C LONG NODE
      L2G=LONOD
      L2=L2G*G2R
C ARG PER
      AR2=ARGPER*G2R


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
C MEAN ANOMALY PARTICLE
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
C INDEX TO DEFINE EQUILIBRIUM POINTS
        ICEQ=0
C NOW DEFINE SIGMA
      DO ISI=1,ISIMAX
C ACRIT IS SIGMA FROM 0 TO 359
C sigma = kk*lambda - kp*lambda_pla + (kp-k)*lonper
C BORRAR ESTO        TLIB(ISI)=0.D0
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
C GIVEN SIGMA CALCULATE THE INTEGRAL  ++++++++++++++++++++++++++
C CALCULATION OF THE INTEGRAL IN LAMBDA PARTICLE
C (IN THE PAPER IS IN LAMBDA PLANET)
        DO ILAMBDA2=1,IPASOS
C ALL POSSIBLE CONFIGURATIONS OCCUR AFTER KP REVOLUTIONS OF THE PARTICLE
C HELIOCENTRIC DISTANCE particle
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
C TRUE ANOMALY planet
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
C END OF CALCULATION OF THE INTEGRAL ++++++++++++++++++++++++++++
        RTOT=RTOT/DFLOAT(IPASOS)
C NOW MULTIPLY BY MASS_PLANET*k**2*MASSTAR
        RP(ISI)=RTOT*XMA(INPLA)*KG2*mst
        MD(ISI)=MINDIS
      ENDDO
C END OF CALCULATION FOR ALL SIGMAS
C NOW WE HAVE THE COMPLETE FUNCTION R(SIGMA)
C CALCULATION OF MEAN VALUE <R>  FOR CALCULATING STRENGTH SR
      VMEDIO=0.D0
      DO I2=1,ISIMAX
        VMEDIO=VMEDIO+RP(I2)
      ENDDO
      VMEDIO=VMEDIO/DFLOAT(ISIMAX)
C FIND R_max AND R_min
      VRMAX=-999999.D0
      VRMIN=999999.D0
      DO I2=1,ISIMAX
        IF(RP(I2).LT.VRMIN) THEN
          VRMIN=RP(I2)
        ENDIF
C FOR CALCULATING RMAX WE DISCARD CLOSE ENCOUNTERS TO LESS THAN RHTOL*RHILL
        IF(RP(I2).GT.VRMAX.AND.MD(I2).GT.RHTOL) THEN
          VRMAX=RP(I2)
        ENDIF
      ENDDO
C READY, RMAX AND RMIN CALCULATED DISCARDING CLOSE ENCOUNTERS
C =====================================================
C SEARCHING FOR STABLE EQUILIBRIUM POINTS, MAX 4
      DO II=1,4
C LIBRATION PERIOD=0 AND CENTER AT 999 MEANS THAT NO STABLE EQUILIBRIUM POINT WAS DETECTED
        TLIB(II)=0.0D0
        ISIGEQ(II)=999
      ENDDO
C CALCULATION OF DELTA R
      DELR=VRMAX-VRMIN
C JUST IN CASE...
      IF(DELR.LT.0.D0)   THEN
        DELR=0.D0
      ENDIF
C IF DELTA R IS VERY SMALL WE CANNOT DETECT ANY RESONANCE
C 1.D-11 IS A VERY REASONABLE LIMIT FOR RELATIVE VARIATION
C IT IS REALLY A SMALL VARIATION
      DETEC2=DABS((VMEDIO-VRMIN)/VMEDIO)
C      DETEC=DELR/VMEDIO
      IF(DETEC2.LT.1.D-11) THEN
         delr=0.d0
C        WRITE(*,*)'        '
C        WRITE(*,*)'*** UNDETECTABLE RESONANCE ***'
C        WRITE(*,*)'        '
        GOTO 111
      ENDIF

C FOR EACH SIGMA CALCULATE FIRST DERIVATIVES IN RADIANS
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
C        AEQ(I)='         '
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
C SO I INCREASE THE COUNTER
              ICEQ=ICEQ+1
C CORRECTION
              ISIGEQ(ICEQ)=IEXTREMO - 1
      TLIB(ICEQ)=A2*TWOPI/DABS(KK)/DSQRT(3.D0*DABS(DERIV2(IEXTREMO
     *)))/365.25D0
              ELSE
C UNSTABLE
              IF(RP(ISIME1).GT.RP(I))THEN
                IEXTREMO=ISIME1
              ENDIF
          ENDIF
        ENDIF
      ENDIF
      ENDDO
C OUTPUT
      IARG=INT(ARGPER)
      IOME=INT(OME)
C ANCHO2 IS RESONANCE FULL stable WIDTH
 111     ANCHO2=8.D0/3.D0/(KG2*MST)*(DELR)*A2**3
      ANCHO2=2.D0*DSQRT(ANCHO2)
      WRITE(1,7)INPLA,INT(KP),INT(KK),SALSEM,EXC,YNC,IARG,IOME,VMEDIO,
     *VMEDIO-VRMIN,ANCHO2,ISIGEQ(1),ISIGEQ(2),ISIGEQ(3),ISIGEQ(4),
     *TLIB(1),TLIB(2),TLIB(3),TLIB(4)
      WRITE(*,*)INPLA,KP,KK,SALSEM,ANCHO2
    7 FORMAT(I2,I5,I3,F10.5,F6.3,F6.1,2I4,1p,3E12.4,4I4,4E11.3)
      RETURN
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

