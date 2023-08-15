C CORRECCION EN ECUACION DE KEPLER 25 DE JULIO 2009
C REVISAR EPHEM1, SOLKEP, MOVREL, KEPREL PARA ORBITAS MUY EXCENTRICAS
C EFECTO RELATIVISTA OPTATIVO
C PROGRAMA EVORB14.FOR 1 de setiembre de 2007
C VERSION DIC 2005
C JUNIO 2005 INCLUYO EFECTOS RELATIVISTAS DEBIDOS AL CUERPO CENTRAL
C VERSION 3 DE ABRIL DE 2005 (MURIO EL PAPA)
C CORRECCION DE CALCULO DE ENERGIA DE SISTEMA
C CORRECCION DE CALCULO DE DISTANCIA MINIMA DE ENCUENTROS
C PROGRAMA QUE CALCULA LA EVOLUCION ORBITAL DE N PARTICULAS SIN MASA
C EN EL CAMPO DE M PLANETAS.
C UTILIZA UN INTEGRADOR  SIMPLECTICO DE SEGUNDO ORDEN.
C CONSIDERA ENCUENTROS CERCANOS CON LOS PLANETAS pero NO entre planetas
C ULTIMA ACTUALIZACION: 10-08-2000, La Plata
C                                       Adrian Brunini
C DIMENSIONADO PARA MAXIMO DE 500 CUERPOS
C LOGICA DE ENCUENTROS ACTUALIZADA 5-MAYO-2001    
C CORRECCIONES MENORES 10 MARZO 2002      
c mejoras en la llevada de planetas a la epoca de particulas
c intento de que integre para atras en tiempo

	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 MPLA
	REAL*8 M0
	CHARACTER*50 FSAL0,FSAL1,FSAL2,FSAL3,FSAL4,
     *FSAL5,FSAL6,FSAL7,APESOS
     
      PARAMETER (NTT=500)
C      include 'numcuerp.inc'
     
C coordenadas y velocidades:
      DIMENSION XPLA(NTT,3),VPLA(NTT,3)
      DIMENSION XPLA0(NTT,3),VPLA0(NTT,3)
C elementos orbitales a,e,
      DIMENSION APLA(NTT),EPLA(NTT),TAST0(NTT)
C elementos orbitales previos a los encuentros:
      DIMENSION APREV(NTT),EPREV(NTT)
C masas, y distancias heliocentricas:
      DIMENSION MPLA(NTT),RPLA(NTT),RHILL(NTT)
      DIMENSION M0(NTT),RADIO(NTT)
      DIMENSION RPH(NTT)
C vectores de uso circunstancial:
      DIMENSION POS(3),VEL(3)
C vectores de indices:
      DIMENSION INDH(500),IPRIM(NTT),NOMBRE(NTT)
      DIMENSION IND(500),INDS(500)
      DIMENSION INDA(500)
      dimension indrp(NTT)

C  BLOCK DE COMMONS:
      COMMON /DDPI/DOSPI
	COMMON/ARCHI/FSAL3,FSAL4
	COMMON/CTE/GM,GM0,UAPYR
	COMMON/PMASS/MPLA
	COMMON/UNIR/NPLA,NBOD   
C SI LOS PLANETAS NO SON MUY EXCENTRICOS SE PODRIA ASUMIR HILL CTE
C	COMMON/PRE/RHILL
	COMMON/PRIMOR/IPRIM
	COMMON/ELPLA/APLA,EPLA
	COMMON/NPLAINI/NOMBRE
	COMMON/RADIOSOL/RADSOL,RCRIT         
      COMMON/CODIGO/RCRIT1,RCRIT2
	COMMON/RELATIV/CLUZ2,H,CMASS
	COMMON/INDRELATIV/IEFREL



C ***********************USUARIO: CAMBIAR A GUSTO*********************
c datos de entrada
      FSAL0='datos.dat'
c salida elem orbit
      FSAL1='orbitas.sal'
c archivo de objetos eliminados
      FSAL2='finales.sal'
c detalle de encuentros
      FSAL3='encuent.sal'
c particulas colisionadas
      FSAL4='colisio.sal'
c control energia 
      FSAL5='energia.sal'
c control cte Jacobi
      FSAL6='cjacobi.sal'     
c archivo para rearranque      
      FSAL7='condin.xxx'
C MAXIMO SEMIEJE PERMITIDO
C       AMAX=9000.D0
C DISTANCIA A PARTIR DE LA CUAL SE ELIMINA SI a>amax
C       RCRIT2=60.D0
c distancia a menos de la cual si Q<RADSOL se decide eliminar por sung
       RCRIT=0.1D0
C ************************FIN ZONA DE CAMBIOS**************************

C UNIDADES Y CONSTANTES:    
c unidad de tiempo = año juliano
c unidad de distancia = UA
c unidad de masa = masa sol
c yerba preferida = la herboristeria para nerviosos y hepaticos
c mejor yerba = salud y flora
c vuelvo a cambiar, ahora es La Selva

	PI=4.D0*DATAN(1.D0)
	DOSPI=2.D0*PI
      gau=0.01720209895d0
      anio=365.25d0
      GM=(gau*anio)**2
C	DGM=2.D0*GM
	GAR=PI/180.D0
C	TDOS=DLOG(2.D0)
	UAKM=1.D0/1.495978707D+08   
      YRSE=365.25D0*24.D0*3600.D0
       UAPYR=UAKM*YRSE
C     VELOC LUZ
      CLUZ2=(2.99792458D10*UAKM/1.D5*3600.D0*24.D0*ANIO)**2
C POR DEFECTO SUPONGO RPUNTO NEGATIVO INICIALMENTE
	do 100 l=1,NTT
        indrp(l)=-1
 100  continue  

C INSTANTE INICIAL                    
      CTIEMPO=0.d0
      TINIC=0.D0
      T=TINIC
 
      write(*,*)'O---------------------------------------------------O'       
      write(*,*)'|                  programa EVORB15                 |'
      write(*,*)'|                VERSION 15 (JUL 2009)              |'
      write(*,*)'|       calcula evolucion orbital de sistema        |'
      write(*,*)'|      para fines didacticos y experimentales       |'
      write(*,*)'|   (mejora en computo deorbitas muy excentricas)   |'
      write(*,*)'|   integrador de Brunini, adaptacion de Gallardo   |'
      write(*,*)'|              gallardo@fisica.edu.uy               |'
      write(*,*)'O---------------------------------------------------O'
       
C--------------------------------------------------------------------------
C LECTURA DE DATOS:   los READ()APESOS son lineas de explicacion en DATOS
	
	OPEN(20,FILE=FSAL0)
	READ(20,1004)APESOS
	READ(20,*) IEFREL         !INDICADOR PARA CORRECCION RELAT
	READ(20,1004)APESOS
	READ(20,*) IXYZ         !INDICADOR PARA GUARDAR XYZVEL
	READ(20,1004)APESOS
	READ(20,*) RCRIT1,RCRIT2   ! LIMITES PARA ELIMINACION DE PARTICULAS
	READ(20,1004)APESOS
	READ(20,*) CMASS,RSTAR    ! MASA CENTRAL EN MASAS SOLARES y radio
	READ(20,1004)APESOS
	READ(20,*) TMAX         ! tiempo maximo de simulacion en años
	READ(20,1004)APESOS     ! espaciado de las salidas
	READ(20,*) TSAL	
	READ(20,1004)APESOS
	READ(20,*) H            ! paso de integracion en años	        
c obligo que tsal sea un numero entero de pasos h	
	tsal=dnint(tsal/h)*h
	
c si tmax es negativo para que integre hacia atras
      if(tmax.lt.0.) H=-H
	
	
	READ(20,1004)APESOS
	READ(20,*) ISTAT        ! indice que indica si es o no re-arranque

c radio roche del sol EN UAs
      RADSOL=RSTAR*1.14D0*UAKM         
C MASA CENTRAL POR CTE GRAVITACIONAL
      GM0=gm*CMASS

C ISTAT: si no es 0, arranca desde el archivo que genero la corrida
C previa y no genera condiciones inciales.
      
	IF(ISTAT.NE.0)THEN
      
      write(*,*)'------- CONTINUACION DE INTEGRACION ANTERIOR --------'
      write(*,*)'- Lectura de condin.xxx'
      
	   OPEN(80,FILE=FSAL7,STATUS='OLD',FORM='UNFORMATTED')
	   READ(80)T,NBOD,NPLA
	    TINIC=T
	    DO I=1,NBOD
	    READ(80)XPLA(I,1),XPLA(I,2),XPLA(I,3)
	    READ(80)VPLA(I,1),VPLA(I,2),VPLA(I,3)
	    READ(80)IPRIM(I),RPH(I),MPLA(I),RPLA(I),NOMBRE(IPRIM(I))
	    ENDDO
	   CLOSE(80)
      
c si tmax es menor que el instante en que se grabo condin.xxx debo
c integrar para atras
         if(tmax.lt.t) H=-dabs(h)
      
	  GOTO 4800

C Es una corrida que comienza y por lo tanto debe abrir archivos nuevos
C
	  ELSE

	    OPEN(10,FILE=FSAL1,STATUS='unknown',access='append')                                            
	write(10,*)'    t       sobrevivientes'
	write(10,*)'  N         a       e       i     nodo  ',
     #'argper anomed'  
	    OPEN(30,FILE=FSAL2,STATUS='unknown',access='append')
      write(30,*)'     t elim      r elim     a pre   e pre ', 
     #'    i    nodo  argper    N'
	    OPEN(40,FILE=FSAL3,STATUS='unknown',access='append')  
	write(40,*)'       t enc  Dm/R  Vp(K/s)  energ   Pla  Par' 
	    OPEN(60,FILE=FSAL4,STATUS='unknown',access='append') 
	write(60,*)'   t colision    Dmin/Rad  Vcol(K/s) ',
     #' energia   Pla  Par'
      close (10)
      close (30)
      close (40)
      close (60)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(IXYZ.EQ.1) THEN

	OPEN(11,FILE='xyzv.sal',STATUS='unknown',access='append')
	write(11,*)'    t       sobrevivientes'
      write(11,*)'   N        x            y            z           v',
     %'x           vy           vz'
	close(11)
      ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++




	ENDIF

C ******************************************************************
C *******************************************************************

C corrida nueva: lectura de los elementos orbitales de los planetas

	NBOD=0
      write(*,*)'- Voy a leer elementos de planetas:'
	READ(20,1004)APESOS              !DATOS PLANETAS
      READ(20,*) NPLA                  !NUMERO PLANETAS      
      write(*,*)'  son',npla,'  planetas.'      
	READ(20,1004)APESOS            
       read(20,*) T0PLA                ! EPOCA PLANETAS
	READ(20,1004)APESOS               ! ELEMENTOS     
	
      DO 101 IU=1,NPLA    
        READ(20,*)L,A0,E0,XI0,XN0,W0,EME0,VAMA,XRADI       
	NBOD=NBOD+1
C EVITO INC=0 PUES DA ERROR
      IF(XI0.LT.0.00001) XI0=0.00001D0      
	XI0=XI0*GAR
	I=NBOD
	NOMBRE(I)=L
C EVITO MASA PLANETA=0 (NO DEBERIA SER CERO) PUES DA ERROR      
      IF(VAMA.LT.1.D-20) THEN
        XM0=1.D20  
        else
        xm0=1.d0/vama
      END IF      
      M0(I)=XM0
      RADIO(I)=XRADI    
      eme0=eme0*gar 
	XN0=XN0*GAR
	W0=W0*GAR
	Q0=A0*(1.D0-E0)   
c n calculado incluyendo masa planeta
	ENE=DSQRT(GM*(CMASS+vama)/A0**3)	
	TAST0(I) = T0PLA     
c JD del instante del pasaje por perihelio
	TPER0=TAST0(I)-EME0/ENE*anio 
C	TPER=(TPER0-TAST0(I))/ANIO
      MPLA(I)=GM/M0(I)
      alpm=mpla(i)
	TPAN=TPER0/ANIO
	T0AN=T0PLA/ANIO               
C calculo de las posiciones y velocidades:      
	CALL EPHEM1(GM0,Q0,E0,XI0,XN0,W0,TPAN,T0AN,POS,VEL,alpm)
      RPLA(I)=0.D0
	DO J=1,3
	RPLA(I)=RPLA(I)+POS(J)**2
	XPLA(I,J)=POS(J)
	VPLA(I,J)=VEL(J)
	ENDDO
	RPLA(I)=DSQRT(RPLA(I))
C      RHILL(I)=RPLA(I)*(MPLA(I)/(3.D0*GM0))**0.33333
c radio planeta en UAs
	RPH(I)=RADIO(I)*UAKM

101     CONTINUE

C-----------------------------------------------------------------
C lectura de los elementos orbitales de las PARTICULAS:
	NBOD=npla                   
              write(*,*)'- Voy a leer elementos particulas:'
	READ(20,1004)APESOS              !DATOS particulas
      READ(20,*) NPARTI                  !NUMERO PARTICULAS      
              write(*,*)'  son',nparti,'  particulas.'
	READ(20,1004)APESOS            
       read(20,*) TJUL                ! EPOCA PARTICULAS
	READ(20,1004)APESOS               ! ELEMENTOS	
      DO 106 IU=1,NPARTI          
        READ(20,*)L,A0,E0,XI0,XN0,W0,EME0
C EVITO INC=0 PUES DA ERROR
      IF(XI0.LT.0.00001) XI0=0.00001D0      
	NBOD=NBOD+1
	XI0=XI0*GAR
      eme0=eme0*gar
	XN0=XN0*GAR
	W0=W0*GAR
	Q0=A0*(1.D0-E0)
c mov medio en radianes por dia
	ENE=DSQRT(GM0/A0**3)
	I=NBOD
	NOMBRE(I)=L
	TAST0(I) = TJUL 
	TPER0=TAST0(I)-EME0/ENE*anio
C	TPER=(TPER0-TAST0(I))/ANIO
	TPAN=TPER0/ANIO
	TJAN=TJUL/ANIO    
C calculo de las posiciones y velocidades:      
	CALL EPHEM1(GM0,Q0,E0,XI0,XN0,W0,TPAN,TJAN,POS,VEL,0.d0)
      RPLA(I)=0.D0
	DO J=1,3
	  RPLA(I)=RPLA(I)+POS(J)**2
	  XPLA(I,J)=POS(J)
	  VPLA(I,J)=VEL(J)
	ENDDO
	RPLA(I)=DSQRT(RPLA(I))
	MPLA(I)=0.D0
	RPH(I)=0.D0
	RHILL(I)=0.D0
106   CONTINUE

C inicializacion del vector de indice primordial
	DO I=1,NBOD
	IPRIM(I)=I              
	ENDDO
		
C calcula las condiciones iniciales de los objetos todos al mismo 
C instante: lleva a los planetas al instante en que vienen dadas las
C condiciones inciales de los asteroides (para todos el mismo t).  
c el tiempo inicial pasa a ser la epoca de los asteroides.
                
        write(*,*)'- Llevo los planetas a la epoca de las particulas.'
                
	CALL LLEVA(XPLA,VPLA,RPLA,TAST0,T0PLA)
C graba en el archivo las condiciones iniciales:        
	OPEN(80,FILE=FSAL7,FORM='UNFORMATTED')
	WRITE(80)T,NBOD,NPLA
	DO I=1,NBOD
	WRITE(80)XPLA(I,1),XPLA(I,2),XPLA(I,3)
	WRITE(80)VPLA(I,1),VPLA(I,2),VPLA(I,3)
	WRITE(80)IPRIM(I),RPH(I),MPLA(I),RPLA(I),NOMBRE(IPRIM(I))
	ENDDO
	CLOSE(80)

c si es una integracion nueva guardo los elementos en T=0
	OPEN(10,FILE=FSAL1,STATUS='unknown',access='append')	
	WRITE(10,1003)T,Nbod
	DO KK=1,Nbod
	CALL PLANORB(XPLA,VPLA,KK,CLINA,dospi,gm0,mpla,on,w,amed,se,sa)
	WRITE(10,1001)NOMBRE(IPRIM(KK)),sa,se,CLINA,on,w,amed	
c si npla=1 calculo jacobi
      if(npla.eq.1) then
c identifico a, masa, planeta
        if(kk.eq.1) then
           aj=sa
           uj=mpla(1)/gm0
        endif
c identificacion,a,e,i particula
        if(kk.eq.2) then
          a0jacob=sa
          e0jacob=se
          x0jacob=clina           
        endif
      endif
	ENDDO
      close(10)      

c+++++++++++++++++++++++++++++++++++++++++++++
      IF(IXYZ.EQ.1) THEN
      open(11,file='xyzv.sal',status='unknown',access='append')
      write(11,1003) t,nbod
      do i=1,nbod
      write(11,999)nombre(iprim(i)),XPLA(I,1),XPLA(I,2),XPLA(I,3),
     %VPLA(I,1),VPLA(I,2),VPLA(I,3)
      enddo
      close(11)
      ENDIF
  999 format(i5,1P6e13.5)
c+++++++++++++++++++++++++++++++++++++++++++++

c calculo jacobi
      if(npla.eq.1) then
        OPEN(21,FILE=FSAL6,access='append',status='unknown')
        WRITE(21,*) 'paso'
        write(21,707)h
        WRITE(21,*) 'a planeta'
        write(21,707)aj
        WRITE(21,*) 'masa planeta'
        write(21,707)uj
        WRITE(21,*) 'a particula'
        write(21,707)a0jacob
        WRITE(21,*) 'e particula'
        write(21,707)e0jacob
        WRITE(21,*) 'i particula'
        write(21,707)x0jacob
C PLANETA
        X11=XPLA(1,1)
        X12=XPLA(1,2)
        X13=XPLA(1,3)
        V11=VPLA(1,1)
        V12=VPLA(1,2)
        V13=VPLA(1,3)
C PARTICULA        
        X21=XPLA(2,1)
        X22=XPLA(2,2)
        X23=XPLA(2,3)
        V21=VPLA(2,1)
        V22=VPLA(2,2)
        V23=VPLA(2,3)

        CALL JACOBO(X11,X12,X13,X21,X22,X23,
     *V11,V12,V13,V21,V22,V23,C,uj,aj,CMASS)
        WRITE(21,422) T,C,NOMBRE(IPRIM(2))
        close(21)
      end if

c si es una corrida que continua salta aqui
4800    CONTINUE
	HK = 0.5D0*H

        write(*,*)'- Comienza integracion numerica.'
        write(*,*)'- El tiempo esta en anios julianos de 365.25 dias.'
        write(*,*)'- Instante inicial: epoca de las particulas.'
      write(*,*)'-----------------------------------------------------'
        write(*,*)'Tiempo         Sobrevivientes'


 707  format(f20.10)

      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C %%%%%%%%%%%%%%%%%%%%% INTEGRACION NUMERICA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c
c SUBRUTINAS UTILIZADAS:
c
c
c  PERTURB: evalua la perturbacion planetaria. 
c  MOVIKEP: avanza un objeto en orbita heliocentrica, si su orbita es 
c           hiperbolica, avisa con ICOD=1 
c
c  VECTORES DE INDICES:
c
c  INDA:   indica cuales pares de cuerpos acretan
c  INDH:   indica cuales cuerpos se eliminan del total de objetos
c          incluidos los acretados.
c  inds   indica cuales continuan en encuentro y vivos luego de pasar por perturb
  
       NCL=0
       NCLS=0
  
C LAZO DE ESCRITURA:

200     TGRID=0.D0

C LAZO DE UN PASO DE INTEGRACION:
	
210     TGRID=TGRID+H
	
C ESTE CONTADOR ES MAS PRECISO QUE SUMAR MUCHAS VECES H
       CTIEMPO=CTIEMPO+1.d0
       T=TINIC+CTIEMPO*H

C	T=T+H
	NHYP=0
	NA=0
	
C GUARDA CANTIDADES PREVIAS AL PASO DE INTEGRACION, QUE USA EN CASO DE
C QUE SE PRODUZCA UN ENCUENTRO ENTRE UN PLANETA Y UNA PARTICULA:

	DO I=1,NBOD
	APREV(I)=APLA(I)
	EPREV(I)=EPLA(I)
	DO L=1,3
	XPLA0(I,L)=XPLA(I,L)
	VPLA0(I,L)=VPLA(I,L)
	ENDDO
	ENDDO

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%% PASO DE INTEGRACION NUMERICA %%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c avanza keplerianamente H/2 los planetas, en el caso de particulas
c no avanzara aquellas que esten muy proximas a algun planeta (INDS,NCLS)
c pues la orbita heliocentrica puede ser una locura (hiperbola, con eyeccion)
c no hay drama pues peturb detectara la situacion e intergrara con BS un paso H
	CALL MOVIKEP(NBOD,HK,XPLA,VPLA,RPLA,APLA,EPLA,NHYP,INDH,NCL,IND,
     #mpla,ncls,inds,1,t) 

	CALL PERTURB(XPLA,VPLA,XPLA0,VPLA0,RPLA,APLA,EPLA,
     1NBOD,NPLA,RPH,MPLA,INDA,NA,IND,NCL,H,t,indrp,INDS,NCLS,NHYP,INDH)

c avanza keplerianamente H/2 los planetas, en el caso de particulas
c no avanzara aquellas que perturb haya detectado en encuentro (IND,NCL)
c pues peturb ya las intergraro con BS
	CALL MOVIKEP(NBOD,HK,XPLA,VPLA,RPLA,APLA,EPLA,NHYP,INDH,NCL,IND,
     #mpla,ncls,inds,2,t)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%% FIN DE LA INTEGRACION NUMERICA %%%%%%%%%%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C ACRECION DE CUERPOS:
	
	DO M=1,NA,2       
	I=INDA(M)
	J=INDA(M+1)
	NHYP=NHYP+1                                ! el J lo tiene que   ! 
	INDH(NHYP)=J                               ! eliminar            !
	ENDDO       
                    
                    
C ************************************************************************
C ELIMINA LAS PARTICULAS EYECTADAS EN ORBITAS HYPERBOLICAS
C Y LAS COLISIONADAS:        

	DO J=1,NHYP
	 I=INDH(J)

	CALL PLANORB(XPLA0,VPLA0,I,CLINA,dospi,gm0,mpla,on,w,amed,se,sa)
	OPEN(30,FILE=FSAL2,STATUS='unknown',access='append')

	WRITE(30,1002)T,RPLA(I),APREV(I),EPREV(I),
C	WRITE(30,1002)T,RPLA(I),SA,SE,
     $CLINA,on,w,NOMBRE(IPRIM(I))
      close(30)
	 CALL ORDENA(XPLA,VPLA,RPLA,APLA,EPLA,RPH,MPLA,I,IPRIM)   
	
	  DO K=J+1,NHYP                       ! corre los indices hiperboli-!
	   L=INDH(K)                          ! cos de los cuerpos que      !
	   IF(L.GT.I)INDH(K)=INDH(K)-1        ! estan mas arriba que el I   !
	  ENDDO                               ! para tirar los correctos    !
	 
	ENDDO
    
C ************************************************************************
C *******************          ESCRITURA        **************************
C ************************************************************************
	
	IF(DABS(TGRID-dsign(TSAL,h)).LT.1.D-05)THEN    
	
	WRITE(*,1012)T,NBOD

c CALCULOS REFERENTES A LOS PLANETAS EXISTENTES:            
	OPEN(10,FILE=FSAL1,STATUS='unknown',access='append')
	
	WRITE(10,1003)T,Nbod
	DO KK=1,Nbod
	CALL PLANORB(XPLA,VPLA,KK,CLINA,dospi,gm0,mpla,on,w,amed,se,sa)
	WRITE(10,1001)NOMBRE(IPRIM(KK)),sa,se,CLINA,on,w,amed
	ENDDO
      close(10)

c+++++++++++++++++++++++++++++++++++++++++++++
      IF(IXYZ.EQ.1) THEN
      open(11,file='xyzv.sal',status='unknown',access='append')
      write(11,1003) t,nbod
      do i=1,nbod
      write(11,999)nombre(iprim(i)),XPLA(I,1),XPLA(I,2),XPLA(I,3),
     %VPLA(I,1),VPLA(I,2),VPLA(I,3)
      enddo
      close(11)
      ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++


c calculo jacobi  si npla=1    
      if(npla.eq.1) then
           aj=APLA(1)
           uj=mpla(1)/gm0
C PLANETA i=1
        X11=XPLA(1,1)
        X12=XPLA(1,2)
        X13=XPLA(1,3)
        V11=VPLA(1,1)
        V12=VPLA(1,2)
        V13=VPLA(1,3)
C PARTICULA i=2       
        X21=XPLA(2,1)
        X22=XPLA(2,2)
        X23=XPLA(2,3)
        V21=VPLA(2,1)
        V22=VPLA(2,2)
        V23=VPLA(2,3)

        CALL JACOBO(X11,X12,X13,X21,X22,X23,
     *V11,V12,V13,V21,V22,V23,C,uj,aj,CMASS)
      
        OPEN(21,FILE=FSAL6,access='append',status='unknown')
        WRITE(21,422) T,C,NOMBRE(IPRIM(2))
        close(21)
      end if


C Guarda las condiciones iniciales para arrancar en caso de parada:        
	 OPEN(80,FILE=FSAL7,STATUS='OLD',FORM='UNFORMATTED')
	 WRITE(80)T,NBOD,NPLA
	 DO I=1,NBOD
	 WRITE(80)XPLA(I,1),XPLA(I,2),XPLA(I,3)
	 WRITE(80)VPLA(I,1),VPLA(I,2),VPLA(I,3)
	 WRITE(80)IPRIM(I),RPH(I),MPLA(I),RPLA(I),NOMBRE(IPRIM(I))
	 ENDDO
	 CLOSE(80)


C CHEQUEO ENERGIA-------------------------------------------------
C CALCULO VEL BARICENTRO
C UNIDADES UA Y ANIO, MASA SOLAR=1
      TOTMAS=CMASS
      DO J=1,NPLA
        TOTMAS=TOTMAS+MPLA(J)/GM
      ENDDO
      VXB=0.D0
      VYB=0.D0
      VZB=0.D0
      DO 40 J=1,NPLA
        VXB=VXB+MPLA(J)/GM*VPLA(J,1)
        VYB=VYB+MPLA(J)/GM*VPLA(J,2)
        VZB=VZB+MPLA(J)/GM*VPLA(J,3)
  40  CONTINUE 
      VXB=VXB/TOTMAS
      VYB=VYB/TOTMAS
      VZB=VZB/TOTMAS 
C ENERGIA CINETICA      
      ENCINET=0.D0 
C CONTRIBUCION DE LOS PLANETAS
      DO 43 J=1,NPLA
        VRB2=(VPLA(J,1)-VXB)**2+(VPLA(J,2)-VYB)**2+(VPLA(J,3)-VZB)**2
        ENCINET=ENCINET+MPLA(J)/GM/2.D0*VRB2
  43  CONTINUE 
C CONTRIBUCION DEL SOL
      ENCINET=ENCINET+1.D0/2.D0*(VXB*VXB+VYB*VYB+VZB*VZB)*CMASS
C ENERGIA POT      
      SUMPOT=0.D0    
C TERMINOS ENTRE PLANETAS
      DO K=1,NPLA-1
        DO J=K+1,NPLA 
          SX=(XPLA(J,1)-XPLA(K,1))**2
          SY=(XPLA(J,2)-XPLA(K,2))**2
          SZ=(XPLA(J,3)-XPLA(K,3))**2
          SUMPOT=SUMPOT-GM*MPLA(K)/GM*MPLA(J)/GM/DSQRT(SX+SY+SZ)
        ENDDO
      ENDDO
C TERMINOS CON EL SOL
      DO J=1,NPLA
        SX=XPLA(J,1)**2
        SY=XPLA(J,2)**2
        SZ=XPLA(J,3)**2
        SUMPOT=SUMPOT-GM*CMASS*MPLA(J)/GM/DSQRT(SX+SY+SZ)
      ENDDO
C ENERGIA TOTAL
      ENTOTAL=ENCINET+SUMPOT
  
      OPEN(37,FILE=FSAL5,STATUS='UNKNOWN',ACCESS='APPEND')
      write(37,420) T,ENTOTAL,encinet,sumpot,NPLA
      CLOSE(37) 
      
 420  FORMAT(1PE13.6,E16.8,E16.8,E16.8,I4)
 
C AGREGADO    SROLAND
 422  FORMAT(1PE13.6,E16.8,I4)




	GOTO 200                                   ! vuelve al principio !   

	ENDIF

C SI no SE ALCANZO EL TIEMPO MAXIMO, sigue LA EJECUCION:        
                                  
                                  
                                  
	IF((tmax-t)*dsign(1.d0,h).gt.0.)GOTO 210

C FORMATOS DE ESCRITURA:

1001  FORMAT(I5,f12.6,f9.6,4F7.2)
1002  FORMAT(F14.3,2f10.4,f9.6,3f7.2,I5)
1003  FORMAT(1PE13.6,I5)
1004  FORMAT(A76)
1012  FORMAT(F14.3,1X,'|',I6,'|') 
C 1013  FORMAT(F7.3,1X,I5)
C 1020  FORMAT(F13.3,1X,E14.10)   


C      close (10)
      close (20)
      close (30)
      close (40)
C      close (50)
      close (60)
      close (70)
      close (80)

	STOP        
	END
C ========================================================================
C ======================== FIN DEL MAIN PROGRAM ==========================
C ========================================================================      

C ========================================================================      
C ========================================================================
C subrutina que calcula el movimiento kepleriano de los objetos
	SUBROUTINE MOVIKEP(NBOD,H,XPLA,VPLA,RPLA,APLA,EPLA,NHYP,INDH,
     1                     NCL,IND,mpla,ncls,inds,momo,t)
	IMPLICIT REAL*8(A-H,O-Z)
	real*8 mpla
      CHARACTER*50 FSAL3,FSAL4
      PARAMETER (NTT=500)
C      include 'numcuerp.inc'
	DIMENSION XPLA(NTT,3),VPLA(NTT,3)
	DIMENSION APLA(NTT),EPLA(NTT)
	DIMENSION RPLA(NTT),mpla(NTT),nombre(NTT)
	DIMENSION POS(3),VEL(3)
	DIMENSION INDH(500),IND(500),inds(500)
	DIMENSION IPRIM(NTT)
	COMMON/ARCHI/FSAL3,FSAL4              
	COMMON/CTE/GM,GM0,UAPYR
	COMMON/RADIOSOL/RADSOL,RCRIT  

      COMMON/CODIGO/RCRIT1,RCRIT2
      COMMON/PRIMOR/IPRIM
	COMMON/NPLAINI/NOMBRE
C	UAKM=1.D0/1.495978707D+08
C      YRSE=365.25D0*24.D0*3600.D0
c indice del sol = 0      
      isol=0
      
	DO 410 I=1,NBOD
	
c agregado 3 de marzo 2001
c se fija si la particula fue eliminada asi no la vuelve a eliminar
c en realidad eliminaria al siguiente injustamente!!!
c esto solo corre en el segundo llamado de movikep pues en el primero
c tenemos nhyp=0

        do ip=1,nhyp   
          if(i.eq.indh(ip)) goto 410
        enddo                                          
                                                     

C SE FIJA SI EL CUERPO I ESTABA EN UN ENCUENTRO, ASI NO LO 
C AVANZA EN LA ORBITA KEPLERIANA, PUES YA LO ESTA HACIENDO 
C LA INTEGRACION POR B&S.  
C Y ADEMAS PORQUE PUEDE DAR ORBITA HELIOCENTRICA HIPERBOLICA EYECTANDOLO  

      if(momo.eq.1) then
	DO J=1,NCLs
	IF(IPRIM(I).EQ.INDs(J))GOTO 410
	ENDDO
      end if


      if(momo.eq.2) then
	DO J=1,NCL
	IF(IPRIM(I).EQ.IND(J))GOTO 410
	ENDDO
	end if         
	

	
	R0=RPLA(I)
	V02=VPLA(I,1)*VPLA(I,1)+VPLA(I,2)*VPLA(I,2)+VPLA(I,3)*VPLA(I,3)
c	A0=GM*R0/(2.D0*GM-R0*V02)                   ! calcula el semieje      !
c mejor lo calculo considerando la masa

      ua0=2.d0/r0-v02/(gm0+mpla(i))

c           write(*,*) ua0


      a0=1.d0/ua0



c--------------------------------------------------------------------------     
	IF(A0.LT.5.D-06)THEN
	APLA(I)=A0                              ! si corresponde a una    !
	NHYP=NHYP+1                             ! orbita hiperbolica     !   
	INDH(NHYP)=I                            !  lo          !
	GOTO 410                                ! eliminara luego         !
	ENDIF
c--------------------------------------------------------------------------

c--------------------------------------------------------------------------     
	IF(RPLA(I).LT.RADSOL)THEN
	APLA(I)=A0                              ! si corresponde a una    !
	NHYP=NHYP+1                             ! orbita     !   
	INDH(NHYP)=I                            ! Sun-Grazer, lo          !

	DO L=1,3
	POS(L)=XPLA(I,L)
	VEL(L)=VPLA(I,L)
	ENDDO      
	ALPM=MPLA(I)
c	CALL MOVREL(GM0,H,POS,VEL,R0,E0,Q0,ICOD,alpm)   
        energ=0.d0
        vcol=dsqrt(v02)/UAPYR

	  OPEN(60,FILE=FSAL4,STATUS='unknown',access='append')
	  WRITE(60,1007)t,r0/radsol,vcol,energ,isol,NOMBRE(IPRIM(i))  
        close(60)

	GOTO 410                                ! eliminara luego         !
	ENDIF
c--------------------------------------------------------------------------



c--------------------------------------------------------------------------     
C	IF(A0.GT.amax.AND.RPLA(I).GT.RCRIT2)THEN
	IF(RPLA(I).GT.RCRIT2.OR.RPLA(I).LT.RCRIT1)THEN
	APLA(I)=A0                              ! si corresponde a una    !
	NHYP=NHYP+1                             ! orbita hiperbolica o    !   
	INDH(NHYP)=I                            ! Sun-Grazer, lo          !
	GOTO 410                                ! eliminara luego         !
	ENDIF
c--------------------------------------------------------------------------





	DO L=1,3
	POS(L)=XPLA(I,L)
	VEL(L)=VPLA(I,L)
	ENDDO      
	ALPM=MPLA(I)
	CALL MOVREL(GM0,H,POS,VEL,R0,E0,Q0,ICOD,alpm)   
	Q0=A0*(1.D0-E0)

c si el perihelio es proximo al sol y esta muy cerca del sol lo elimino       
      mato=0
      if(q0.lt.radsol.and.r0.lt.rcrit) then
        mato=1  
        energ=0.d0
        vcol=dsqrt(v02)/UAPYR
	  OPEN(60,FILE=FSAL4,STATUS='unknown',access='append')
	  WRITE(60,1007)t,q0/radsol,vcol,energ,isol,NOMBRE(IPRIM(i))  
        close(60)
      end if

1007    FORMAT(F14.3,F11.3,1P2E11.3,2I5)   


c si vino de movrel con codigo de eyeccion tambien mato
	IF(ICOD.GT.0.OR.mato.eq.1)THEN

c agregado 5 de marzo
      rpla(i)=R0

	EPLA(I)=E0
	NHYP=NHYP+1                                
	INDH(NHYP)=I                            
	GOTO 410                                
	ENDIF
c--------------------------------------------------------------------------
	
	RPLA(I)=R0
	APLA(I)=A0
	EPLA(I)=E0
       
	DO L=1,3
	XPLA(I,L)=POS(L)
	VPLA(I,L)=VEL(L)
	ENDDO
		       
410     CONTINUE

	RETURN
	END

C ========================================================================      
C ========================================================================
	SUBROUTINE FYG(X,V,H,R,A0,EX0,I,ICOD,mpla)
C CALCULA LAS COORDENADAS Y VELOCIDADES EN T=T0+H DADAS LAS
C COORDENADAS Y VELOCIDADES EN T0 USANDO LAS F Y G DE LAGRANGE
C EN FORMA CERRADA.
C ENTRA COMO ARGUMENTO: X(3)=     POSICION EN T0
C                       V(3)=     VELOCIDAD EN T0
C                       H   =     INTERVALO DE TIEMPO QUE SE AVANZA
C                       I   =     NUMERO DE ORDEN DEL CUERPO
C SALE COMO ARGUMENTO:  X(3)=     POSICION EN T0+H
C                       V(3)=     VELOCIDAD EN T0+H
C ENTRA POR COMMON:  COMMON/CTE/GM,2*GM,SQRT(GM) DONDE GM
C                    ES LA CTE DE GRAVITACION UNIVERSAL QUE DEFINE
C                    ASI LAS UNIDADES (DEBEN SER COMPAT. CON LAS DE X Y V
c-------------------------------------------------------------------------
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 M,N,mpla
      PARAMETER (NTT=500)
C      include 'numcuerp.inc'
	DIMENSION X(NTT,3),V(NTT,3)
	COMMON/CTE/GM,GM0,UAPYR
	dimension mpla(NTT)


	ICOD=0
	R0=R
	RGMA=DSQRT(A0*(GM0+mpla(i)))
	CE0=1.D0-(R0/A0)
	SE0=(X(I,1)*V(I,1)+X(I,2)*V(I,2)+X(I,3)*V(I,3))/RGMA
	EX02=SE0*SE0+CE0*CE0
	EX0=DSQRT(EX02)
      
c dejo limite de exc=0.99 pues por aca solo pasan los planetas
c que se supune tienen exc<0.99      
      
	IF(EX0.GT.0.99D0)THEN
	   ICOD=1
	   RETURN
	ENDIF

	R0E0=R0*EX0
	E0=DATAN2(SE0,CE0)
	N=RGMA/A0**2
	M=E0-SE0+N*H

	CALL SOLKEP(EX0,M,E)
	
	SE=DSIN(E)
	CE=DCOS(E)
	F=A0*((CE-EX0)*CE0+SE*SE0)/R0E0
	G=((CE0-EX02)*SE-(CE-EX0)*SE0)/(N*EX0)
	X1=F*X(I,1)+G*V(I,1)
	X2=F*X(I,2)+G*V(I,2)
	X3=F*X(I,3)+G*V(I,3)
	R=DSQRT(X1*X1+X2*X2+X3*X3)
	FP=-RGMA*(SE*CE0-CE*SE0)/(R*R0E0)
	GP=(1.D0+G*FP)/F
	V1=FP*X(I,1)+GP*V(I,1)
	V2=FP*X(I,2)+GP*V(I,2)
	V3=FP*X(I,3)+GP*V(I,3)
	X(I,1)=X1
	X(I,2)=X2
	X(I,3)=X3
	V(I,1)=V1
	V(I,2)=V2
	V(I,3)=V3
	RETURN
	END
	
C ========================================================================      
C ========================================================================        
	SUBROUTINE ORDENA(X,V,R,A,E,RPH,EME,I,IPRIM)
	IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (NTT=500)
C      include 'numcuerp.inc'
	DIMENSION X(NTT,3),V(NTT,3),R(NTT),EME(NTT)
	DIMENSION RPH(NTT),A(NTT),E(NTT)
	DIMENSION IPRIM(NTT)

	COMMON/UNIR/NPLA,NBOD
	IF(I.LE.NPLA)NPLA=NPLA-1
	DO J=I,NBOD-1
	
	R(J)=R(J+1)
	A(J)=A(J+1)
	E(J)=E(J+1)
	RPH(J)=RPH(J+1)
	EME(J)=EME(J+1)
	IPRIM(J)=IPRIM(J+1)

	DO L=1,3
	X(J,L)=X(J+1,L)
	V(J,L)=V(J+1,L)
	ENDDO

	ENDDO
	NBOD=NBOD-1

	RETURN
	END
C ========================================================================      
C ========================================================================        
	SUBROUTINE PERTURB(XPLA,VPLA,XPLA0,VPLA0,RPLA,APLA,EPLA,NBOD,
     1NPLA,RPH,MPLA,INDA,NA,IND,NCL,H,tiempo,indrp,INDS,NCLS,NHYP,INDH)

	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 MPLA
      CHARACTER*50 FSAL3,FSAL4
      PARAMETER (NTT=500)
C      include 'numcuerp.inc'
      dimension indrp(NTT)
	DIMENSION XPLA(NTT,3),VPLA(NTT,3),RPLA(NTT)
	DIMENSION XPLA0(NTT,3),VPLA0(NTT,3)
	DIMENSION NOMBRE(NTT),IND(500),INDH(500)
	DIMENSION RPH(NTT),MPLA(NTT),INDA(500),IPRIM(NTT)
	DIMENSION APLA(NTT),EPLA(NTT)
	DIMENSION Y(12)
	DIMENSION FINDX(10),FINDY(10),FINDZ(10)
	dimension INDS(500)
      dimension pertx(NTT),perty(NTT),pertz(NTT)
      dimension posx(NTT),posy(NTT),posz(NTT)
	DIMENSION FREL(NTT,3)
C	COMMON/PRE/RHILL
	COMMON/ARCHI/FSAL3,FSAL4
	COMMON/CTE/GM,GM0,UAPYR
	COMMON/NPLAINI/NOMBRE
	COMMON/PRIMOR/IPRIM
	COMMON/A/EME1
	COMMON/INDRELATIV/IEFREL
C	UAKM=1.D0/1.495978707D+08
C      YRSE=365.25D0*24.D0*3600.D0

C ********************************************************************
C EFECTOS RELATIVISTAS JUNIO 2005 TABARE
      IF(IEFREL.GE.1) THEN
      CALL EFREL(XPLA,VPLA,RPLA,APLA,EPLA,NBOD,IEFREL,FREL)
      ENDIF

c importanteeee!!!!!
c si todas las particulas son clones de un objeto puede haber centenares
c de encuentros en un instante, hay que cerar muchos indices!!!
c cerar hasta 300 es demasiado, tal vez alcance con 30 aunque
c dificilmente se use mas de 10
        do ii=1,500
          ind(ii)=0
          inds(ii)=0
          enddo                          
                                    
                                    
c atenti le meto valores prefijados de facrh0, facrh1 y facrh2      
          facrh0=3.d0
c limite para no hacer calculos heliocentricos
          facrh1=2.d0
c limite para distancia anterior
          facrh2=facrh0                     
                     
      NCLS=0       
	NA=0
	DT0=H/10.D0        
	TOL=1.D-10
	NCL=0
	 
	    do 100 l=1,nbod+1
	      pertx(l)=0.d0  
	      perty(l)=0.d0  
	      pertz(l)=0.d0  
 100      continue  
	        
	DO 2000 I=1,NPLA
C I: PLANETA
	  XGM=MPLA(I)  
	  EME1=XGM
c radio hill	
	  RH=RPLA(I)*(XGM/(3.D0*GM0))**(0.3333333333333D0)

c         rh=3.d0*(xgm/gm/0.001d0)**0.5/3.d0



C TERMINOS INDIRECTOS DE LA PERTURBACION PLANETARIA:
	  RP3=XGM/RPLA(I)**3
	  FINDX(I)=XPLA(I,1)*RP3    
	  FINDY(I)=XPLA(I,2)*RP3    
	  FINDZ(I)=XPLA(I,3)*RP3    

	  DO 1000 J=NPLA+1,NBOD

C J: PARTICULA

C CRITERIO MALHOTRA
C      if(rpla(j).lt.0.22d0) goto 1000

C PARA NO CALCULAR LOS QUE YA SE FUERON
          DO IP=1,NHYP
            IF(J.EQ.INDH(IP)) GOTO 1000
          ENDDO


c A partir de aqui evalua a ver si esta dentro de la esfera de Hill
c y si se alejan o se acercan mutuamente:

	    DX=XPLA(I,1)-XPLA(J,1)                         
	    DY=XPLA(I,2)-XPLA(J,2)   
	    DZ=XPLA(I,3)-XPLA(J,3)   
	    DRQ=DX*DX+DY*DY+DZ*DZ
	    DIJ=DSQRT(DRQ)                 ! distancia mutua actual ! 
	    DX0=XPLA0(I,1)-XPLA0(J,1)                         
	    DY0=XPLA0(I,2)-XPLA0(J,2)   
	    DZ0=XPLA0(I,3)-XPLA0(J,3)   
	    DRQ0=DX0*DX0+DY0*DY0+DZ0*DZ0
	    DIJ0=DSQRT(DRQ0)               ! distancia mutua previa ! 

	    DvX=vPLA(I,1)-vPLA(J,1)                         
	    DvY=vPLA(I,2)-vPLA(J,2)   
	    DvZ=vPLA(I,3)-vPLA(J,3)   

C CORRECCION 31/12/02
          rpunto=(dx*dvx+dy*dvy+dz*dvz)/dij*DSIGN(1.D0,H)

c criterio para no perder encuentros en caso de RH chico
          facrh=facrh0+dabs(rpunto*H/rh)

	    IF(DIJ.GT.RH*facrh.and.DIJ0.gt.rh*facrh2)GOTO 900   
c osea: si ahora (habiendo avanzado H/2) esta en RH*facrh o si en el final
c del paso H anterior estaba a menos de RH*facrh2 entra en la BS
c                  (pensadlo)

c             write(*,*) facrh


C----------------------------------------------------------------------------
C LA PARTICULA ESTARIA DENTRO DE facrh veces el RHill DEL PLANETA.
C EL ENCUENTRO LO CALCULA INTEGRANDO CON B&S CON PASO VARIABLE, PERO
C COMO EL MOVIMIENTO HELIOCENTRICO LO INTEGRA TAMBIEN HAY QUE 
C DESHACER EL AVANCE EN LAS ORBITAS HELIOCENTRICAS QUE YA HIZO
C LA SUBRUTINA MOVIKEP EN EL PROGRAMA PRINCIPAL.
C
	    DO L=1,3
	    Y(L)=XPLA0(I,L)         !  1  a  3  POSICION  PLANETA
	    Y(L+3)=VPLA0(I,L)       !  4  a  6  VELOCIDAD PLANETA
	    Y(L+6)=XPLA0(J,L)       !  7  a  9  POSICION  PARTICULA
	    Y(L+9)=VPLA0(J,L)       !  10 a 12  VELOCIDAD PARTICULA
	    ENDDO
	
C Ya estamos como antes de avanzar el paso H en las orbitas
C keplerianas.
	 
	    DT=DT0     
	    T=0.0D0   
	    RMERGE=RPH(I)  

301       CALL SEXTRA (12,DT,TOL,T,Y,DIJ)        
                               
          drx=y(7)-y(1)                         
          dry=y(8)-y(2)                         
          drz=y(9)-y(3)                         
          dvx=y(10)-y(4)                         
          dvy=y(11)-y(5)                         
          dvz=y(12)-y(6)                         

	    DRQ=DrX*DrX+DrY*DrY+DrZ*DrZ
	    DIJ=DSQRT(DRQ)

C PARA TESTEAR DETECCION DE ENCUENTROS 31/12/02 ----------
C       if(dij.lt.0.01d0) then
C                tival4=tiempo-h+t+dt
C	OPEN(40,FILE=FSAL3,STATUS='unknown',access='append')
C	WRITE(40,*) tival4,dij,i
C      close(40)
C        endif
C---------------------------------------------------------

          rpunto=(drx*dvx+dry*dvy+drz*dvz)/dij*DSIGN(1.D0,H)
c      si en el calculo anterior se estaba acercando......
          if(indrp(iprim(j)).eq.-1) then
c y si ahora se esta alejando...
            if(rpunto.ge.0.d0)then
              dvq=dvx*dvx+dvy*dvy+dvz*dvz
              dr=dij
C MODIFICACIONES 15 FEBRERO 2004********************************
c es hora de calcular detalles de encuentro pues acaba de pasar por
c el punto mas proximo.
c energia del sistema doble (solo tiene sentido si dr<Rhill)
	        E0=0.5D0*DVQ-XGM/DR
C SI DR>RHILL ASUMO Q=DR, SI NO ASUMO QUE HUBO UN MOVIMIENTO PLANETOCENTRICO
C Y CALCULO SU PERIASTRO
                Q=DR
              IF(DR.LT.RH) THEN
	        HX=dry*dvz-dvy*drz
	        HY=dvx*drz-drx*dvz
	        HZ=drx*dvy-dry*dvx
	        HQ=HX*HX+HY*HY+HZ*HZ
c distancia pericentrica 
	        QRA=XGM**2+2.D0*HQ*E0       ! conservacion de la
	        Q=-XGM+DSQRT(QRA)           ! energia y momento
	        Q=0.5D0*Q/E0                ! angular (Danby)
              ENDIF

c guardo si q menor que 2RHill (un poquito mas)                         
              if(q.le.rh*2.1d0) then
c calculo velocidad de pericentro
                vperq=(e0+xgm/q)*2.d0
                vper=dsqrt(vperq)
c paso velocidad a km/s
                vper=vper/UAPYR
c instante corresp a estos valores: CORREGIDO 31/12/02
                tival=tiempo-h+t+dt
c arreglo salidas
                if(vper.gt.99.99) vper=99.99
                if(e0.lt.-99.9) e0=-99.999
                if(e0.gt.99.9) e0=99.999
                                   
	OPEN(40,FILE=FSAL3,STATUS='unknown',access='append')
c guarda datos del encuentro:t,
c pericentro, vel planetoc, energia planetoc, planeta, particula 
	WRITE(40,1107)tival,q/rmerge,
     *vper,e0,Nombre(IPRIM(i)),NOMBRE(IPRIM(J))
      close(40)


C PREGUNTA SI HAY COLISION
C PARA EVITAR PROBLEMAS, SOLO GUARDA LA PRIMERA ACRECION         
C SI ES QUE EL MISMO PLANETA TUVO MAS DE UNA EN ESE MISMO PASO
C COSA QUE ES ALTAMENTE IMPROBABLE
c si esta a menos de RH/3  ASUME QUE EL CALCULO DE Q ESTABA BIEN Y
c decide si hay colision
C                if (q.lt.rmerge) then
C CORRECCION DE FEBRERO 2004
                if (DIJ.LT.RH/3.D0.AND.q.lt.rmerge) then
                          
	            DO NN=1,NA                                   
	            IF(I.EQ.INDA(NN).OR.J.EQ.INDA(NN))GOTO 1000  
	            ENDDO                                        
	            NA=NA+1                                        
	            INDA(NA)=I                                     
	            NA=NA+1                                       
	            INDA(NA)=J  

c calculo velocidad de colision
                  vcolq=(e0+xgm/rmerge)*2.d0
                  vcol=dsqrt(vcolq)
c paso velocidad a km/s
                  vcol=vcol/UAPYR


	OPEN(60,FILE=FSAL4,STATUS='unknown',access='append')
	WRITE(60,1007)tival,q/rmerge,vcol,e0,
     #NOMBRE(IPRIM(I)),NOMBRE(IPRIM(J))  
      close(60)

	            GOTO 1001
	          ENDIF             
              end if        
            end if 
          end if                       
c si ya tuvo el pasaje por el periastro 
c o sea rpunto>0
c ppppppppp           
          indrp(iprim(j))=1
          if(rpunto.lt.0.d0) indrp(iprim(j))=-1

C AJUSTA EL PASO DE INTEGRACION:
	    DT=H*(DIJ/RPLA(I))   
	    IF(dabs(T+DT).LT.dabs(H))GOTO 301   
	    DT=H-T
	    CALL SEXTRA (12,DT,TOL,T,Y,DIJ) 

C TERMINO DE INTEGRAR EL ENCUENTRO un paso H.     
c aqui deberia averiguar si por casualidad justo paso por el periastro    
c y calcular detalles del encuentro pero si no lo hago aca no
c importa pues lo hara al pricipio del proximo paso H.    
    
C AVERIGUO SI se encuentra muy proximo al planeta  pues en ese caso
c la orbita heliocentrica que calculara moviekep (en la primer llamada)
c podria ser espureamente hiperbolica y lo ejectaria.
C NCLS=NUM ENCUENTROS SOBREVIVIENTES peligrosos (QUE CONTINUAN)
c defino peligroso cuando dist<RH*facrh1 (es arbitrario)
          IF(DIJ.LT.RH*facrh1) THEN
            NCLS=NCLS+1
            INDS(NCLS)=IPRIM(J) 
          ENDIF
C ESTOS SON LOS ENCUENTROS QUE MOVIEKEP DEBE EVITAR CALCULAR
C EN LA PRIMERA LLAMADA DEL LAZO PARA EVITAR OBTENER POSIBLES
C ORBITAS HIPERBOLICAS
c nota: no puedo poner exactamente el limite rh*facrh pues el
c planeta no estara exactamente donde lo deja BS (pues hay
c que sumar otras perturbaciones  entonces cuando se evalue
c la DIJ tal vez sea > rh*facrh y no lo avanzara ni perturb ni movikep
c produciendo un error feo


C AL PlANETA LO DEJA DONDE ESTABA ANTES DEL ENCUENTRO:        

1001	    posx(J)=Y(7)
	    posy(J)=Y(8)
	    posz(J)=Y(9)
	    pertx(j)=pertx(j)+Y(10)
	    perty(j)=perty(j)+Y(11)
	    pertz(j)=pertz(j)+Y(12)
	    RPLA(J)=DSQRT(Y(7)**2+Y(8)**2+Y(9)**2)  
c comentarios viejos de Adrian:
C AQUI GUARDA LOS INDICES DE CUERPOS QUE TUVIERON UN ENCUENTRO:
C ASI NO LOS VUELVE A AVANZAR EN ORBITA KEPLERIANA.
C DE ESTE MODO LA LOGICA ES:
C CUANDO DETECTA EL ENCUENTRO POR PRIMERA VEZ, VUELVE ATRAS
C (COORDENADAS PREVIAS), INTEGRA
C CON B&S UN PASO H, Y LO IDENTIFICA. EL PROXIMO AVANCE NO LO HACE ASI
C QUE CUANDO VUELVE A ESTA SUBRUTINA, SI EL ENCUENTRO CONTINUA,
C VUELVE A LAS COORDENADAS PREVIAS, PERO COMO LA ORBITA KEPLERIANA NO AVANZO, 
C LAS COORDENADAS PREVIAS SON LAS DEL FINAL DE ESTA INTEGRACION.
C ASI HASTA QUE SALE DE LA ESFERA DE HILL. SOLO LO HACE CON LA
C PARTICULA PUES AL PLANETA LO INTEGRA INDEPENDIENTEMENTE
C EN FORMA SIMPLECTICA SIEMPRE.

	    NCL=NCL+1
	    IND(NCL)=IPRIM(J)
	    GOTO 1000

C NO HAY ENCUENTRO PROXIMO ENTRE EL PLANETA I Y LA PARTICULA J
C ASI QUE APLICA LA PERTURBACION A LA VELOCIDAD DE ACUERDO AL
C METODO SIMPLECTICO QUE ESTAMOS UTILIZANDO:

900       CONTINUE
	    RRP3=XGM/DIJ**3
	    PPLAX=DX*RRP3-FINDX(I)
	    PPLAY=DY*RRP3-FINDY(I)
	    PPLAZ=DZ*RRP3-FINDZ(I)

C CALCULO DE LAS PERTURBACIONES A LA VELOCIDAD:

	    pertx(J)=pertx(J)+H*PPLAX
	    perty(J)=perty(J)+H*PPLAY 
	    pertz(J)=pertz(J)+H*PPLAZ 

1000    CONTINUE
2000  CONTINUE

c suma todas las perturbaciones a las particulas
	DO 200 J=NPLA+1,NBOD 

	  VPLA(J,1)=VPLA(J,1)+pertx(j)
	  VPLA(J,2)=VPLA(J,2)+perty(j)
	  VPLA(J,3)=VPLA(J,3)+pertz(j)       

        do lclo=1,ncl
          if(ind(lclo).eq.IPRIM(j)) then
c si la part tuvo encuentro, pert incluye la velocidad
	      VPLA(J,1)=pertx(j)
	      VPLA(J,2)=perty(j)
	      VPLA(J,3)=pertz(j)       
c si la particula tuvo un encuentro queda donde la dejo BS
            xpla(j,1)=posx(j)
            xpla(j,2)=posy(j)
            xpla(j,3)=posz(j)
          end if
        enddo
 200  continue








C FREL, 15 DE JUNIO 2005, TABARE
      IF(IEFREL.GE.1) THEN
	DO  J=1,NBOD
	  VPLA(J,1)=VPLA(J,1)+FREL(J,1)*H
	  VPLA(J,2)=VPLA(J,2)+FREL(J,2)*H
	  VPLA(J,3)=VPLA(J,3)+FREL(J,3)*H
        ENDDO
      ENDIF
 
 
C AHORA PERTURBA MUTUAMENTE A LOS PLANETAS:
	DO I=1,NPLA
	DO J=I+1,NPLA
	DX=XPLA(I,1)-XPLA(J,1)   
	DY=XPLA(I,2)-XPLA(J,2)   
	DZ=XPLA(I,3)-XPLA(J,3)   
	DRQ=DX*DX+DY*DY+DZ*DZ
	DIJ3=DSQRT(DRQ)**3 
	DX=DX/DIJ3
	DY=DY/DIJ3
	DZ=DZ/DIJ3
	PIJX=MPLA(I)*DX-FINDX(I)
	PIJY=MPLA(I)*DY-FINDY(I)
	PIJZ=MPLA(I)*DZ-FINDZ(I)
	VPLA(J,1)=VPLA(J,1)+H*PIJX
	VPLA(J,2)=VPLA(J,2)+H*PIJY 
	VPLA(J,3)=VPLA(J,3)+H*PIJZ
	PJIX=-MPLA(J)*DX-FINDX(J)
	PJIY=-MPLA(J)*DY-FINDY(J)
	PJIZ=-MPLA(J)*DZ-FINDZ(J)
	VPLA(I,1)=VPLA(I,1)+H*PJIX
	VPLA(I,2)=VPLA(I,2)+H*PJIY 
	VPLA(I,3)=VPLA(I,3)+H*PJIZ
	ENDDO
	ENDDO

1007    FORMAT(F14.3,F11.3,1P2E11.3,2I5)   
c1017    FORMAT(2F13.3,F11.3,1P2E11.3,2I5)   
1107    FORMAT(F15.4,F7.1,F6.2,F9.4,2I5)

	RETURN
	END
C ========================================================================      
C ========================================================================        
C SOLUCION DE LA ECUACION DE KEPLER ELIPTICA:
	SUBROUTINE SOLKEP(EX,M,E)
C
C SOLUCION ITERATIVA DE LA ECUACION DE KEPLER
C ENTRA:EX   EXCENTRICIDAD            (<1)
C       M    ANOMALIA MEDIA           (RADIANES)
C SALE: E    ANOMALIA EXCENTRICA      (RADIANES)
C
	 IMPLICIT REAL*8 (A-H,O-Z)
	 REAL*8 M,MK
      COMMON /DDPI/DOSPI
	 TOLE=1.D-14
C	 DPI=8.D0*DATAN(1.D0)
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
	 IF(NITER.GT.20) THEN
C	 WRITE(*,*)'MUCHAS ITERACIONES LA ECUACION DE KEPLER'
           GOTO 200
         ENDIF
	 IF(DEX.GT.TOLE)GOTO 100
	 RETURN


C SI EL NUMERO DE ITERACIONES ES > 20 PRUEBA CON BISECCION:

C METODO DICOTOMICO:
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
C	 WRITE(*,*)'ERROR EN LA SOLUCION DE LA ECUACION DE KEPLER'
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

C ========================================================================      
C ========================================================================        
C CALCULO DE LA INCLINACION ORBITAL y otras cosas:
	SUBROUTINE PLANORB(XPLA,VPLA,I,INCLI,dospi,gm0,mpla,on,w,amed,
     *se,sa)    
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 INCLI, mpla
      PARAMETER (NTT=500)
C      include 'numcuerp.inc'
	DIMENSION XPLA(NTT,3),VPLA(NTT,3),mpla(NTT)

C calculo metido por tabare      
        
           
      xx=xpla(i,1)     
      yy=xpla(i,2)
      zz=xpla(i,3) 

      r=dsqrt(xx*xx+yy*yy+zz*zz)

      vx=vpla(i,1)   
      vy=vpla(i,2) 
      vz=vpla(i,3) 

      v2=vx*vx+vy*vy+vz*vz

      hx=yy*vz-zz*vy
      hy=zz*vx-xx*vz
      hz=xx*vy-yy*vx
      hh=dsqrt(hx*hx+hy*hy+hz*hz)
       
       
      
      
      on=datan2(hx,-hy)
      yn=dacos(hz/hh) 
      incli=yn

      IF(ON.LT.0D0) THEN
        ON=ON+dosPI
      END IF

      CMU=gm0+mpla(i)

      ex=(vy*hz-vz*hy)/CMU-xx/r
      ey=(vz*hx-vx*hz)/CMU-yy/r
      ez=(vx*hy-vy*hx)/CMU-zz/r

      ee=dsqrt(ex*ex+ey*ey+ez*ez)
      
      se=ee
      
      sew=ez/dsin(yn)
      cow=ex*dcos(on)+ey*dsin(on)
      w=datan2(sew,cow)
      IF(W.LT.0D0) THEN
        W=W+dosPI
      END IF

      a=1.d0/(2.d0/r-v2/CMU)   
      
      sa=a
      
      AE=0D0
      IF(A.GT.0.D0) THEN
         see=(xx*vx+yy*vy+zz*vz)/dsqrt(a*CMU)
         coe=1.d0-r/a
         ae=datan2(see,coe)
           IF(AE.LT.0D0) THEN
              AE=AE+dosPI
           END IF
        
         amed=ae-ee*dsin(ae)  
         amed=dmod(amed,dospi)
           IF(Amed.LT.0D0) THEN
              Amed=Amed+dosPI
           END IF        
      END IF
                       
c si es hiperbola no calculo anomalia media y le meto=0
      if(a.lt.0.d0) then
         amed=0.d0
      end if                 
                       
                       
                       
      incli=incli/dospi*36.d1      
      on=on/dospi*36.d1      
      w=w/dospi*36.d1      
      amed=amed/dospi*36.d1      
            

	RETURN
	END
 
C ========================================================================      
C ========================================================================        
C SOLO SE LLAMA DESDE MOVREL CUANDO LA EXC ES MUY ALTA
      SUBROUTINE ELEM(GM0,T,X,V,Q,E,FI,FN,PER,T0,ALPM)
CC
CC
CC PURPOSE    :  THIS SUBROUTINE CALCULATES THE (OSCULATING)
CC               ORBITAL ELEMENTS OF A CELESTIAL BODY , WHOSE
CC               CARTESIAN POSITION- AND VELOCITY COMPONENTS ARE
CC               GIVEN IN THE ARRAYS X(K) , V(K) , K=1,2,3
CC
CC PARAMETERS :
CC         IN :  GM    : GRAVITY-CONSTANT (GAUSSIAN CONSTANT R*8
CC                       PLANETARY APPLICATIONS)
CC                       IN  (UA**1/3)  / DAYS * SOLARMASS**1/2
CC               T     : EPOCH TO WHICH X,V (AND ELEMENTS)   R*8
CC                       REFER (MJD)
CC               X,V   : POSITION AND VELOCITY OF CELESTIAL R*8(3)
CC                       BODY AT TIME T
CC        OUT :  Q      : PERICENTRIC DISTANCE                R*8
CC               E      : NUMERICAL EXCENTRICITY              R*8
CC               FI     : INCLINATION OF ORBITAL PLANE WITH   R*8
CC                        RESPECT TO FUNDAMENTAL PLANE OF
CC                        COORDINATE SYSTEM
CC               FN     : LONGITUDE OF ASCENDING NODE         R*8
CC               PER    : ARGUMENT OF PERIGEE/HELION          R*8
CC               T0     : TIME OF PERIGEE/HELION PASSING      R*8
CC
CC                                           D.M.C Y A.B. 14/2/90
C*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(3),V(3),XX(3),H(3)
      COMMON /CODIGO/ RCRIT1,RCRIT2
      
C CONSIDERO MASA CENTRAL Y MASA CUERPO
      GGM=GM0+ALPM      
      
C-- MOMENTO ANGULAR  HX,-HY,HZ
      H(1)= X(2)*V(3)-X(3)*V(2)
      H(2)=-X(3)*V(1)+X(1)*V(3)
      H(3)= X(1)*V(2)-X(2)*V(1)
C-- LONG. OF NODE ,INCLINATION AND ARG OF LATITUDE
      FN=DATAN2(H(1),H(2))
      CK=DCOS(FN)
      SK=DSIN(FN)
      FI=DATAN2(DSQRT(H(1)**2+H(2)**2),H(3))
      CI=DCOS(FI)
      SI=DSIN(FI)
C
      XX(1)=    CK*X(1)   +SK*X(2)
      XX(2)=-CI*SK*X(1)+CI*CK*X(2)+SI*X(3)
      XX(3)= SI*SK*X(1)-SI*CK*X(2)+CI*X(3)
      U=DATAN2(XX(2),XX(1))
C-- P/GM, R, V Y COS(ANGULO ENTRE R Y V)
      P=(H(1)**2+H(2)**2+H(3)**2)/GGM
      R= DSQRT(X(1)**2+X(2)**2+X(3)**2)
      VE=DSQRT(V(1)**2+V(2)**2+V(3)**2)
      RRP=X(1)*V(1)+X(2)*V(2)+X(3)*V(3)
      CT=RRP/(R*VE)
C-- ENERGIA POR UNIDAD DE MASA,Y EXCENTRICIDAD
      ALFA= VE * VE / 2.D0 - GGM / R
      E= DSQRT( 1.D0 + 2.D0*P*ALFA/GGM)

c attenti: voy a permitir exc mas proximas de 1
c	IF(DABS(E-1.D0).LT.1.D-02)GOTO 200



	IF(DABS(E-1.D0).LT.1.D-06)GOTO 200





	IF(E.LT.1.D0)GOTO 100
	IF(E.GT.1.D0)GOTO 300


C-- CASO ELIPTICO
100   A= P / (1.D0 - E*E )
      IF(R.GT.RCRIT2.OR.R.LT.RCRIT1)ICOD=1
      CEX= (A-R) / (A*E)
      SEX= RRP/(E*DSQRT(GGM*A))



      EX= DATAN2(SEX,CEX)




      V1= 2.D0 * DATAN( DSQRT((1.D0+E)/(1.D0-E)) * DTAN(EX/2.D0) )
      T0= T- (EX -E*DSIN(EX)) * DSQRT(A**3/GGM)
      PER= U-V1
      Q= A*(1.D0-E)
      RETURN
      
C-- CASO PARABOLICO
200   E=1.D0
      Q= P / 2.D0
      ICOD=1

	CC=  (P - R) / R
	SS=DSQRT(1.D0-CC**2)
      V1 = DATAN(SS/CC)
	IF (CT.LT.0) V1=-V1
      EX=DTAN(V1/2.D0)
      T0= T - (EX +(EX**3)/ 3.D0) * DSQRT(2.D0 * Q**3 / GGM)
      PER=U-V1
      RETURN
      
C-- CASO HIPERBOLICO
300   A= P /(E**2 -1.D0)
      
      CHEX= (A+R)/(A*E)
      SHEX= DSQRT(CHEX**2-1.D0)
      EX=DLOG(CHEX+SHEX)
      IF (CT.LT.0) THEN
      EX=-EX
      SHEX=-SHEX
      ENDIF
      Z=DSQRT(A**3/GGM)
      ZE=E*SHEX
      T0= T - (E * SHEX - EX) * DSQRT(A**3/GGM)
      Q= A *(E-1.D0)
      V1=2.D0*DATAN(DTANH(EX/2.D0)/DSQRT(-(1.D0-E)/(1.D0+E)))
      PER=U-V1
      RETURN
      
      END

C ========================================================================      
C ========================================================================        
	SUBROUTINE EPHEM1(GM0,Q,E,I,ANODE,PERI,T0,T,X,XP,alpm)
CC
CC NAME       :  EPHEM1
CC
CC      CALL EPHEM1(GM,Q,E,I,ANODE,PERI,T0,T,X,XP)
CC
CC PURPOSE    :  COMPUTE POSITION X AND VELOCITY XP OF A CELESTIAL
CC               BODY WHOSE OSCULATING ELEMENTS ARE GIVEN
CC               (PARABOLAS, HYPERBOLAS, ELYPSES)
CC
CC PARAMETERS :
CC         IN :  GM     : GRAVITY - CONSTANT                  R*8
CC               Q      : PERICENTRIC DISTANCE                R*8
CC               E      : NUMERICAL EXCENTRICITY              R*8
CC               I      : INCLINATION (RADIAN)                R*8
CC               ANODE   : RIGHT ASCENSION OF ASCENDING NODE   R*8
CC                        (RADIAN)
CC               PERI   : ARGUMENT OF PERICENTRE (RADIAN)     R*8
CC               T0     : PERICENTRE-PASSING-TIME             R*8
CC               T      : EPOCH OF EPHEMERIS COMPUTATION      R*8
CC               X      : POSITION OF SATELLITE               R*8
CC               XP     : VELOCITY OF SATELLITE               R*8
CC               X(K), XP(K),K=1,2,3: POSITION AND VELOCITY   R*8
CC
CC
CC                                                D.M.C  14/02/90
CC
C*
c yo (tabare) le agrego la masa alpm en el calculo de n
c los tiempos son en ANIOS julianos
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 I,M,N
	DIMENSION X(3),XP(3)
c	COMMON/SIGP/SIG,NSIG,NCOM
	DPI=8.D0*DATAN(1.D0)


c atenti: permito que la exc sea mas proxima a 1
c	IF(DABS(E-1.D0).LT.1.D-02)GOTO 200
	
	IF(DABS(E-1.D0).LT.1.D-06)GOTO 200

	
	IF(E.LT.1.D0)GOTO 100
	IF(E.GT.1.D0)GOTO 300
C
C-- CASO ELIPTICO
C
C P=PARAMETER OF CONIC SECTION
  100   A=Q/(1.D0-E)
	P=A*(1.D0-E**2)
C N = MEAN MOTION, M = MEAN ANOMALY
	N=DSQRT((GM0+alpm)/A**3)
	M=N*(T-T0)
      M=M-IDINT(M/DPI)*DPI   

	CALL SOLKEP(E,M,EX)  

C V = TRUE ANOMALY
	V=2*DATAN(DSQRT((1.D0+E)/(1.D0-E))*(DSIN(EX/2)/DCOS(EX/2)))
	R=A*(1.D0-E*DCOS(EX))
	BETA=DSQRT((GM0+alpm)/P)
	X1=R*DCOS(V)
	X2=R*DSIN(V)
	XP1=-BETA*DSIN(V)
	XP2= BETA*(E+DCOS(V))
	GOTO 1000
C
C--CASO PARABOLICO
 200    N=DSQRT( (GM0+alpm)/ (2.D0 * Q**3 ) ) 
	M=N*(T-T0)  
C	WRITE(*,*)'CASO paraBOLICO'

	XX=(3.D0 * M + DSQRT(9.D0 * M**2 + 4.D0))*0.5D0
	Z=1.D0/3.D0
	Z=XX**Z
	Z=Z-1.D0/Z
	X1= Q * (1.D0 - Z**2)
	X2= 2.D0 * Q * Z
	R=Q*(Z**2+1.D0)
	ZP=DSQRT((GM0+alpm)/(2.D0 * Q))/R
	XP1=-2.D0 * Q * Z * ZP
	XP2= 2.D0 * Q * ZP
	GOTO 1000
C
C-- CASO HIPERBOLICO
300     A=Q/(E-1.D0)
C
C	WRITE(*,*)'CASO HIPERBOLICO'

	N=DSQRT((GM0+alpm)/A**3)
	M=N*(T-T0)  

C
C F, F1: ECCENTRIC ANOMALY
	F1=0.d0
	NITER=0
320     F=F1-(M+F1-E*DSINH(F1))/(1.D0-E*DCOSH(F1))
	IF(DABS(F-F1).LT.1.D-12)GO TO 330
	F1=F
	
	NITER=NITER+1
	IF(NITER.GT.1000)THEN
C	WRITE(*,*)'PROBLEMAS DE ITERACION HIPERBOLICA'
	GOTO 330
	ENDIF

	GO TO 320
 330    CONTINUE
	AE2= A*DSQRT(E**2-1.D0)
	CHF=DCOSH(F)
	SHF=DSINH(F)
	X1= A * (E - CHF)
	X2= AE2 * SHF
	R=A*(E*CHF-1.D0)
	FP=DSQRT((GM0+alpm)/A)/R
	XP1=-A*SHF*FP
	XP2=AE2*CHF*FP
C--
C
C SINES AND COSINES OF INCLINATION I, NODE K, PERIGEE O
 1000   CK=DCOS(ANODE)
	SK=DSIN(ANODE)
	CI=DCOS(I)
	SI=DSIN(I)
	CO=DCOS(PERI)
	SO=DSIN(PERI)
C
C VECTORS P AND Q
	P1=CK*CO-SK*CI*SO
	P2=SK*CO+CK*CI*SO
	P3=SI*SO
C
	Q1=-CK*SO-SK*CI*CO
	Q2=-SK*SO+CK*CI*CO
	Q3= SI*CO
C
C COMPUTE POSITION AND VELOCITY
	X(1)=P1*X1+Q1*X2
	X(2)=P2*X1+Q2*X2
	X(3)=P3*X1+Q3*X2
C
	XP(1)=P1*XP1+Q1*XP2
	XP(2)=P2*XP1+Q2*XP2
	XP(3)=P3*XP1+Q3*XP2
	RETURN
	END

C ========================================================================      
C ========================================================================        
C DT en anios julianos
C SUBRUTINA DE ALTO USO
	SUBROUTINE MOVREL(GM0,DT,X,V,R,E,Q,IFAIL,ALPM)
	IMPLICIT REAL*8(A-H,O-Z)
	DIMENSION X(3),V(3)
	IFAIL=0   
C CONSIDERO MASA AL CALCULAR SEMIEJE Y N 
	GGM=GM0+ALPM
	VQ=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)
	A=1.D0/((2.D0/R)-VQ/GGM)
	IF(A.LT.0.D0)THEN
	  IFAIL=1
	RETURN
	ENDIF
	EC=1.D0-R/A
	U=X(1)*V(1)+X(2)*V(2)+X(3)*V(3)
	EN=SQRT(GGM/(A*A*A))  
	ES=U/(EN*A*A)
	E=DSQRT(EC*EC+ES*ES)
	Q=A*(1.D0-E)
C------------------------------------------------
C CASO CUASI PARABOLICO:

c atenti: permito valores de exc mas cerca de 1
c	IF(DABS(E-1.D0).LT.1.D-02)THEN      


	IF(DABS(E-1.D0).LT.1.D-06)THEN
	T=0.D0
	CALL ELEM(GM0,T,X,V,Q,E,FI,FN,PER,T0,ALPM)
	T=T+DT
      CALL EPHEM1(GM0,Q,E,FI,FN,PER,T0,T,X,V,ALPM)
C	R=DSQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
C ESTAMOS ELIMINANDO LA PARTICULA
        IFAIL=1
	RETURN
	ENDIF

C-----------------------------------------------             
	CALL KEPREL(ES,EC,EN,DT,DX,C,S,FP) 

	F=(A/R)*(C-1.D0)+1.D0
	G=DT+(S-DX)/EN
	FDOT=-A*EN*S/(R*FP)
	GDOT=1.D0+(C-1.D0)/FP
	X1=X(1)*F+V(1)*G
	X2=X(2)*F+V(2)*G
	X3=X(3)*F+V(3)*G
	V1=X(1)*FDOT+V(1)*GDOT
	V2=X(2)*FDOT+V(2)*GDOT
	V3=X(3)*FDOT+V(3)*GDOT
	X(1)=X1        
	X(2)=X2
	X(3)=X3
	V(1)=V1
	V(2)=V2
	V(3)=V3
	R=DSQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))   
	RETURN
	END
	
C ========================================================================      
C ========================================================================        
C RESUELVE LA ECUACION DE KEPLER ELIPTICA PARA UN ARCO DT:        
	SUBROUTINE KEPREL(ES,EC,EN,DT,X,C,S,FP)
	IMPLICIT REAL*8(A-H,O-Z)
	ICOD=0
	N=0
C	TOL=1.D-14
	TOL=1.D-13
	ENDT=EN*DT
	X=ENDT
100     S=DSIN(X)
	C=DCOS(X)
	F=X-EC*S+ES*(1.D0-C)-ENDT
	FP=1.D0-EC*C+ES*S
	FPP=EC*S+ES*C
	FPPP=EC*C-ES*S
	DX=-F/FP
	DX=-F/(FP+DX*FPP/2.D0)
	DX=-F/(FP+DX*FPP/2.D0+DX*DX*FPPP/6.D0)
	X=X+DX
	N=N+1
	IF(N.GT.10)THEN
C	 WRITE(*,*)'VOY AL DICOTOMICO'
	  GOTO 200
	ENDIF
	IF(DABS(DX).GT.TOL)GOTO 100
	RETURN
C METODO DICOTOMICO PARA ENCONTRAR UN MEJOR PUNTO INICIAL:
200      CONTINUE 
	PI=4.D0*DATAN(1.D0)
	 DPI=2.D0*PI       
	 NDIC=0
	 X0=-DPI
	 DX0=DPI/10.D0
	 S=DSIN(X0)
	 C=DCOS(X0)
	 F0=X0-EC*S+ES*(1.D0-C)-ENDT  

400      DX=DX0/(10.D0**NDIC)           
	 NITER=0
	
300      X=X0+DX
	 NITER=NITER+1
	 
	 IF(NITER.GT.100)THEN
C	 WRITE(*,*)'ERROR EN LA SOLUCION DE LA ECUACION DE KEPLER'
	 RETURN
	 ENDIF

	 S=DSIN(X)
	 C=DCOS(X)
	 F=X-EC*S+ES*(1.D0-C)-ENDT    
	 IF(F*F0.GT.0.D0)THEN
	 X0=X
	 F0=F
	 GOTO 300
	 ELSE
	 NDIC=NDIC+1
	 IF(NDIC.LT.3)GOTO 400
	 IF(ICOD.EQ.1)RETURN
	N=0      
101     S=DSIN(X)
	C=DCOS(X)
	F=X-EC*S+ES*(1.D0-C)-ENDT
	FP=1.D0-EC*C+ES*S
	FPP=EC*S+ES*C
	FPPP=EC*C-ES*S
	DX=-F/FP
	DX=-F/(FP+DX*FPP/2.D0)
	DX=-F/(FP+DX*FPP/2.D0+DX*DX*FPPP/6.D0)
	X=X+DX
	N=N+1
	IF(N.GT.200)THEN
C	 WRITE(*,*)'ERROR EN LA SOLUCION DE LA ECUACION DE KEPLER    2'
	ICOD=1
	GOTO 200
	ENDIF
	IF(DABS(DX).GT.TOL)GOTO 101
	RETURN
	 
	 ENDIF
	 END

C ========================================================================      
C ========================================================================
	SUBROUTINE LLEVA(XPLA,VPLA,RPLA,TAST0,T0PLA)
	IMPLICIT REAL*8(A-H,O-Z)
	REAL*8 MPLA
      PARAMETER (NTT=500)
C      include 'numcuerp.inc'

	DIMENSION XPLA(NTT,3),VPLA(NTT,3)
C	DIMENSION APLA(NTT),EPLA(NTT),TAST0(NTT)
	DIMENSION TAST0(NTT)
	DIMENSION MPLA(NTT),RPLA(NTT)
	DIMENSION FINX(100),FINY(100),FINZ(100)


	COMMON/CTE/GM,GM0,UAPYR
	COMMON/PMASS/MPLA
	COMMON/UNIR/NPLA,NBOD
	
c numero de pasos tal que cada paso sea aprox 1 dias
      deltajd=dabs((TAST0(NPLA+1)-T0PLA)/1.D0)   
c le agrego 10 pasos para evitar que pueda ser cero
	n1d=int(deltajd)+10
	dn1d=DBLE(n1d)
	DT=(TAST0(NPLA+1)-T0PLA)/365.25D0
	HPLA=DT/dn1d
	
	DO 4200 NSTEP=1,n1d
	
	DO 4100 I=1,NPLA
	R0=RPLA(I)
	V02=VPLA(I,1)*VPLA(I,1)+VPLA(I,2)*VPLA(I,2)+VPLA(I,3)*VPLA(I,3)
c	A0=GM*R0/(2.D0*GM-R0*V02)                   ! calcula el semieje      !
c mejor lo calculo considerando la masa
      ua0=2.d0/r0-v02/(gm0+mpla(i))
      a0=1.d0/ua0



	CALL FYG(XPLA,VPLA,HPLA,R0,A0,E0,I,ICOD,mpla)
	RPLA(I)=R0
C	APLA(I)=A0
C	EPLA(I)=E0
4100    CONTINUE
	
	 DO 5000 I=1,NPLA                       ! terminos indirectos      !
	 RP3=RPLA(I)*RPLA(I)*RPLA(I)            !                          !
	 FINX(I)=-XPLA(I,1)/RP3                 !                          !
	 FINY(I)=-XPLA(I,2)/RP3                 !                          !
	 FINZ(I)=-XPLA(I,3)/RP3                 !                          !
5000     CONTINUE
	
C CALCULO DE LOS TERMINOS DIRECTOS DE LA PERTURBACION
	 DO 2200 I=1,NPLA
	 PPLAX=0.D0
	 PPLAY=0.D0
	 PPLAZ=0.D0
c-------------------------------------------------------------------------
	 DO 2400 K=1,NPLA                    !                             !
	 IF(I.EQ.K)GOTO 2400                 !      no se autoperturba     !
c-------------------------------------------------------------------------
	 DX1=XPLA(I,1)-XPLA(K,1)
	 DX2=XPLA(I,2)-XPLA(K,2)
	 DX3=XPLA(I,3)-XPLA(K,3)
	 RRP2=DX1*DX1+DX2*DX2+DX3*DX3
	 RRP3=RRP2**1.5D0
	 FXP = -DX1/RRP3
	 FYP = -DX2/RRP3
	 FZP = -DX3/RRP3
	 PPLAX=PPLAX+MPLA(K)*(FINX(K)+FXP)
	 PPLAY=PPLAY+MPLA(K)*(FINY(K)+FYP)
	 PPLAZ=PPLAZ+MPLA(K)*(FINZ(K)+FZP)
 2400    CONTINUE

C CALCULO DE LAS PERTURBACIONES A LA VELOCIDAD:
	 VPLA(I,1)=VPLA(I,1)+HPLA*PPLAX
	 VPLA(I,2)=VPLA(I,2)+HPLA*PPLAY
	 VPLA(I,3)=VPLA(I,3)+HPLA*PPLAZ
2200     CONTINUE
4200     CONTINUE
	
	 RETURN
	 END


C ========================================================================      
C ========================================================================        
       SUBROUTINE SEXTRA (M,HA,TOL,TA,YA,D12)
CC    PROPOSITO : REALIZAR UN PASO DE INTEGRACION DEL SISTEMA DE
CC                ECUACIONES.
CC    IN /
CC       HA  :      PASO DE INTEGRACION
CC       TOL :      PRECISION CON QUE SE QUIEREN LOS RESULTADOS
CC       T0  :      INSTANTE EN QUE SE DAN LAS CONDICIONES INICIALES
CC       YA  :      COND INICIALES DEL COMETA
CC   OUT /
CC    EN YA VUELVE A SALIR LA SOLUCION DEL SISTEMA
CC       TA  :      INSTANTE EN QUE SE DA LA SOLUCION (TA=T0+HA)
CC       YS  :      POS Y VEL. DE LA ESTRELLA
CC
CC   SBR :SMIDPT
CC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION YA(12),Y(M),V(20,M)
      COMMON / PARAM / NP,NS,MAXIT,MAXEXT
      IFAIL=0
      NSTEP=2**NS
      NP2=2**NP
      H=HA
      N=1
      DO 100 I0=0,MAXIT
      I=I0+1
      ERR=0.D0
      CALL SMIDPT (M,N,H,TA,YA,Y,D12)
      IF (I.GT.1) GOTO 600
      DO 700 K=1,M
700   V(1,K)=Y(K)
      GOTO 300
600   DO 200 K=1,M
      YY=Y(K)
      V(I,K)=YY
      IP2=NP2
      JF=I0
      IF (JF.GT.MAXEXT) JF=MAXEXT
      DO 400 JB=1,JF
      J=I-JB
      DIV=IP2-1
      DY=(YY-V(J,K))/DIV
      YY=YY+DY
      V(J,K)=YY
400   IP2=NSTEP*IP2
      DERR=DABS(DY)
      IF (DERR.GT.ERR) ERR=DERR
200   CONTINUE
      IF (ERR.LE.TOL) GOTO 500
300   H=H/2.D0
100   N=2*N
      IFAIL=1
500   TA=TA+HA
      DO 800 K=1,M
800   YA(K)=V(J,K)
      RETURN
      END

C ========================================================================      
C ========================================================================
      SUBROUTINE SMIDPT (M,N,H,TA,YA,Y,D12)
C
C  PROPOSITO : INTEGRACION DE UN SISTEMA DE ECUACIONES
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C      DIMENSION YA(12),Y(12),Y0(1100),DY(1100)
      DIMENSION YA(12),Y(12),Y0(M),DY(M)
      H2=H/2.D0
      DO 1 K=1,M
1     Y0(K)=YA(K)
      T=TA+H2
      CALL DERIV (T,Y0,DY,D12)
      DO 2 K=1,M
2     Y(K)=Y0(K)+H2*DY(K)
      IF (N.EQ.1) GOTO 20
      DO 10 J=2,N
      T=T+H2
      CALL DERIV (T,Y,DY,D12)
      DO 3 K=1,M
3     Y0(K)=Y0(K)+H*DY(K)
      T=T+H2
      CALL DERIV (T,Y0,DY,D12)
      DO 4 K=1,M
4     Y(K)=Y(K)+H*DY(K)
10    CONTINUE
20    T=T+H2
      CALL DERIV (T,Y,DY,D12)
      DO 5 K=1,M
5     Y0(K)=Y0(K)+H*DY(K)
      CALL DERIV (T,Y0,DY,D12)
      DO 6 K=1,M
6     Y(K)=(Y0(K)+Y(K)+H2*DY(K))/2.D0
      RETURN
      END

C ========================================================================      
C ========================================================================
C
C        DEFINICION DE LOS PARAMETROS DE LA EXTRAPOLACION.
C
      BLOCK DATA
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON / PARAM / NP,NS,MAXIT,MAXEXT
      DATA NP,NS,MAXIT,MAXEXT /2,2,10,6/
      END

      
C ========================================================================      
C ========================================================================      
      SUBROUTINE DERIV (T,Y,DY,D12)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION Y(12),DY(12)
	COMMON/CTE/GM,GM0,UAPYR
	COMMON/A/EME1
	
c
C ***           Calcula distancias
c
      D1=DSQRT(Y(1)**2+Y(2)**2+Y(3)**2)**3
      D2=DSQRT(Y(7)**2+Y(8)**2+Y(9)**2)**3   
      Y17=Y(1)-Y(7)
      Y28=Y(2)-Y(8)
      Y39=Y(3)-Y(9)
      D12=DSQRT(Y17**2+Y28**2+Y39**2)
      D123=D12**3
      GM1=GM0/D1
      GM2=GM0/D2
      EME1D123=EME1/D123
      EME1D1=EME1/D1
      
c
c ***           Evalua las fuerzas
c
      DY(1)=Y(4)
      DY(2)=Y(5)
      DY(3)=Y(6)
      DY(4)=-GM1*Y(1)
      DY(5)=-GM1*Y(2)
      DY(6)=-GM1*Y(3)
      DY(7)=Y(10)
      DY(8)=Y(11)
      DY(9)=Y(12)
      DY(10)=-GM2*Y(7)+EME1D123*Y17-EME1D1*Y(1)
      DY(11)=-GM2*Y(8)+EME1D123*Y28-EME1D1*Y(2)
      DY(12)=-GM2*Y(9)+EME1D123*Y39-EME1D1*Y(3)

      RETURN
      END

C ========================================================================      
C ========================================================================
c calculo de constante de jacobi
	SUBROUTINE JACOBO(X11,X12,X13,X21,X22,X23,
     *V11,V12,V13,V21,V22,V23,C,uj,aj,CMASS)
	IMPLICIT REAL*8(A-H,O-Z)
	PI=4.D0*DATAN(1.D0)
      gau=0.01720209895d0
C MASA SOL  Y MU
      SO=CMASS     
      UM=UJ/(SO+UJ)
C PERIODO ORBITAL PLANETA EN DIAS
      P=2.D0*PI*AJ**1.5D0/GAU/DSQRT(SO+UJ)   
C UNIDAD DE TIEMPO ANTERIOR
      ANIO=365.25D0              
C COORD HELIOC DE PLANETA        
      XJ=X11/AJ
      YJ=X12/AJ
      ZJ=X13/AJ  
C COORD HELIOC DE PARTICULA        
      XP=X21/AJ
      YP=X22/AJ
      ZP=X23/AJ  
C BARICENTRO
      XB=XJ*UJ/(UJ+SO)
      YB=YJ*UJ/(UJ+SO)
      ZB=ZJ*UJ/(UJ+SO)
      
      XX=(XP-XB)**2
      YY=(YP-YB)**2
C DISTANCIAS AL SOL Y A PLANETA      
      R1=DSQRT(XP**2+YP**2+ZP**2)
      X2=XP-XJ
      Y2=YP-YJ
      Z2=ZP-ZJ
      R2=DSQRT(X2**2+Y2**2+Z2**2)
C VELOCIDAD BARICENTRO
      XPB=V11*UJ/(UJ+SO)
      YPB=V12*UJ/(UJ+SO)
      ZPB=V13*UJ/(UJ+SO)
C VELOCIDAD BARIC EN NUEVA UNIDAD DE TIEMPO Y DISTANCIA
      XPB=XPB/AJ*P/ANIO/2.d0/pi
      YPB=YPB/AJ*P/ANIO/2.d0/pi
      ZPB=ZPB/AJ*P/ANIO/2.d0/pi
C VELOCIDAD PARTICULA EN NUEVA UNIDAD DE TIEMPO Y DISTANCIA
      XPP=V21/AJ*P/ANIO/2.d0/pi
      YPP=V22/AJ*P/ANIO/2.d0/pi
      ZPP=V23/AJ*P/ANIO/2.d0/pi
C VELOCIDAD EN EL SISTEMA ROTANTE
      VX=XPP-XPB+(YP-YB)
      VY=YPP-YPB-(XP-XB)
      VZ=ZPP-ZPB
      VV=VX**2+VY**2+VZ**2
C JACOBI      
      C=XX+YY+2.D0*(1.D0-UM)/R1+2.D0*UM/R2-VV
      RETURN
      END
C ========================================================================
C ========================================================================
C CALCULO EFECTOS RELATIVISTAS DEBIDOS AL CUERPO CENTRAL
C version tabare 17 agosto 2008
	SUBROUTINE EFREL(XPLA,VPLA,RPLA,APLA,EPLA,NBOD,IEFREL,FREL)
	IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NTT=500)
C      include 'numcuerp.inc'
	DIMENSION FREL(NTT,3)
	DIMENSION XPLA(NTT,3),VPLA(NTT,3),RPLA(NTT)
	DIMENSION APLA(NTT),EPLA(NTT)
	COMMON/CTE/GM,GM0,UAPYR
	COMMON/RELATIV/CLUZ2,H,CMASS

 	DO 220 J=1,NBOD

      FREL(J,1)=0.D0
      FREL(J,2)=0.D0
      FREL(J,3)=0.D0

C CORRECCIONES PROMEDIADAS  a la Gallardo   apsidal
      IF(IEFREL.EQ.5) THEN
C SEMIEJE
	A=APLA(J)
C MOMENTO ANGULAR
	        HX=XPLA(J,2)*VPLA(J,3)-XPLA(J,3)*VPLA(J,2)
	        HY=XPLA(J,3)*VPLA(J,1)-XPLA(J,1)*VPLA(J,3)
	        HZ=XPLA(J,1)*VPLA(J,2)-XPLA(J,2)*VPLA(J,1)
C	        H3=(HX*HX+HY*HY+HZ*HZ)**1.5D0
C VECTOR EXC
      EX=(VPLA(J,2)*HZ-VPLA(J,3)*HY)/GM0 - XPLA(J,1)/RPLA(J)
      EY=(VPLA(J,3)*HX-VPLA(J,1)*HZ)/GM0 - XPLA(J,2)/RPLA(J)
      EZ=(VPLA(J,1)*HY-VPLA(J,2)*HX)/GM0 - XPLA(J,3)/RPLA(J)
      E2=EX*EX+EY*EY+EZ*EZ
C AHORA DECIDIMOS SI VALE LA PENA CORREGIR, VEMOS SI EL EFECTO EN ARGUMENTO
C DEL PERIHELIO ES MAYOR que 0.01 arcsec/yr  (DECI=0.259)
C      DECI=DSQRT(CMASS**3/A**5)/(1.D0-E2)
C      IF(DECI.GT.0.25) THEN
 	  FAC1=-2.d0*GM0**2/CLUZ2/A**3/(1.D0-E2)**1.5D0
          FREL(J,1)=FAC1*EX
          FREL(J,2)=FAC1*EY
          FREL(J,3)=FAC1*EZ
C      ENDIF
      ENDIF


C CORRECCION ORIGINAL DE QUINN
C DA ERRORES CUANDO LA DIST PERIHELICA ES MUY CHICA, NO RECOMENDABLE
      IF(IEFREL.EQ.2) THEN
         RV=XPLA(J,1)*VPLA(J,1)+XPLA(J,2)*VPLA(J,2)+XPLA(J,3)*VPLA(J,3)
 	  R1=RPLA(J)
 	  R3=RPLA(J)**3
 	  V2=VPLA(J,1)**2+VPLA(J,2)**2+VPLA(J,3)**2
 	  FAC1=4.D0*GM0/CLUZ2/R1-V2/CLUZ2
      FREL(J,1)=GM0/R3*(FAC1*XPLA(J,1)+4.D0/CLUZ2*VPLA(J,1)*RV)
      FREL(J,2)=GM0/R3*(FAC1*XPLA(J,2)+4.D0/CLUZ2*VPLA(J,2)*RV)
      FREL(J,3)=GM0/R3*(FAC1*XPLA(J,3)+4.D0/CLUZ2*VPLA(J,3)*RV)
      ENDIF

C CORRECCION saha tremaine
      IF(IEFREL.EQ.3) THEN
C SEMIEJE
	A=APLA(J)
      E2=EPLA(J)*EPLA(J)
 	  R1=RPLA(J)
 	  R3=RPLA(J)**3
      FAC1=3.d0*GM0**2/CLUZ2/R3*((4.d0/dsqrt(1.d0-e2)-1.d0)/A-2.D0/R1)
      FREL(J,1)=FAC1*XPLA(J,1)
      FREL(J,2)=FAC1*XPLA(J,2)
      FREL(J,3)=FAC1*XPLA(J,3)
      ENDIF

C CORRECCION nobili roxbourg
      IF(IEFREL.EQ.4) THEN
 	  R1=RPLA(J)
 	  R3=RPLA(J)**3
      FAC1=3.d0*GM0**2/CLUZ2/R3*(-2.D0/R1)
      FREL(J,1)=FAC1*XPLA(J,1)
      FREL(J,2)=FAC1*XPLA(J,2)
      FREL(J,3)=FAC1*XPLA(J,3)
      ENDIF

C CORRECCIONES PROMEDIADAS  a la Gallardo RADIAL
      IF(IEFREL.EQ.1) THEN
      E2=EPLA(J)*EPLA(J)
C AHORA DECIDIMOS SI VALE LA PENA CORREGIR, VEMOS SI EL EFECTO EN ARGUMENTO
C DEL PERIHELIO ES MAYOR que 0.01 arcsec/yr  (DECI=0.259)
C      DECI=DSQRT(CMASS**3/Apla(j)**5)/(1.D0-E2)
C      IF(DECI.GT.0.25) THEN
 	  FAC1=3.d0*GM0**2/CLUZ2/APLA(J)**3/(1.D0-E2)**1.5D0/RPLA(J)
          FREL(J,1)=FAC1*XPLA(J,1)
          FREL(J,2)=FAC1*XPLA(J,2)
          FREL(J,3)=FAC1*XPLA(J,3)
C      ENDIF
      ENDIF




 220  continue
      END

