c program atlaspr  - gallardo, beauge, giuppone 2020
c gallardo@fisica.edu.uy
c given a planetary system this program calculates an atlas of resonances
c with a test (moving) planet
c jacobi reference system
c for very weak resonances the output could be erroneous and the equilibrium
c points wrong
      implicit real*8 (a-h,j-z)
      integer kp,kt,kmax,ntot,j,ncont
      integer npla(100000),nkp(100000),nkt(100000)
      integer ncpla(10000),nckp(10000),nckt(10000)
      dimension ares(100000)
      dimension se(10),sm(10),ex(10),yn(10),ln(10),lp(10),rh(10)
      dimension sig(10),tli(10)
      dimension rp(400),deriv1(400),deriv2(400),tlib(400),md(400),
     *irh(400),sma(400)
c this 50000 dimension is for storing data particle, can be modified
c max number of points in the integral   = 1000*max(k1,k2)
      dimension vx2(50000),vy2(50000),vz2(50000),vr2(50000),vla2(50000)
      character*9 aeq(400),ace(400)
      character*40 zcomment
      character*57 zcabe1
      character*84 zcabe
      twopi = 8.0d0*datan(1.0d0)
      cero  = 0.0d0
      uno   = 1.0d0
      pi=twopi/2.d0
      g2r=pi/180.d0
      kgaus=0.01720209895d0
      kg2=kgaus**2


      zcabe1=' k2:k1   strength   e_test  a1         width      a2     '
      zcabe=zcabe1//'   width     sigma  period '


      write(*,*)'------------------------------------------------------'
      write(*,*)'                    atlaspr      '
      write(*,*)'          atlas of planetary resonances'
      write(*,*)'              version jan 9, 2021'
      write(*,*)'         Gallardo, Beauge, Giuppone 2020'
      write(*,*)'------------------------------------------------------'
      write(*,*)'input data in atlaspr.inp'
c input
      open(1,file='atlaspr.inp',status='old')
      read (1,*)  zcomment
c mass central star
      read(1,*) mst
      read (1,*)  zcomment
      read(1,*) ntot
      read (1,*)  zcomment
c read planets a e i lonod, loper, mass
      do  i=1,ntot
        read (1,*) se(i),ex(i),yn(i),ln(i),lp(i), sm(i)
      enddo
      read (1,*)  zcomment
c test planet  with index=10
        read (1,*) yn(10),ln(10),lp(10), sm(10)
      read (1,*)  zcomment
        read (1,*) emin,emax,deltae
      close (1)


      write(*,*)'max n for resonances 1:n or n:1 ?'
      read(*,*)kmax
      write(*,*)'a_min, a_max (au) ?'
      read(*,*)amin,amax

c number of rhill admitted to calculate disturbing function
c this admit a fine tunning according to the orbital inclination
c this is for the calculation of STABLE libration region. If you want
c to locate the sepatrices then punt rhtol=0.do and good luck
      rhtol=3.d0



c calculation of the valid resonances
c counter i
      i=0
c for each planet.....
      do 200 ipla=1,ntot
c kp*n_pla = kt*n_test
c kp mutiplies lambda_planet
      do kp=1,kmax
        vkp=dfloat(kp)
c kt mutiplies lambda_ptest
        do kt=1,kmax
          vkt=dfloat(kt)
c calculate semimajor axis approximately just to know if is between limits
      sem=(vkt/vkp)**2*(mst+sm(10))/(mst+sm(ipla))
      sem=sem**(1.d0/3.d0)*se(ipla)
            if (sem.le.amax.and.sem.ge.amin)  then
              i=i+1
              npla(i)=ipla
              nkp(i)=kp
              nkt(i)=kt
              ares(i)=sem
            endif
        enddo
      enddo
 200  continue
c number of possible (redundant) resonances, must be less than 100000
      idat=i
      if(idat.ge.100000) then
        write(*,*)'error: too much resonances'
        stop
      end if
c eliminate repeated, example: 4:2 is the same as 2:1
      do i=1,idat
        do j=i+1,idat
          if(npla(j).eq.npla(i)) then
            if (ares(j).eq.ares(i).and.nkp(j).ge.nkp(i)) then
              ares(j)=0.d0
            endif
          endif
        enddo
      enddo
      ncont=0
      do i=1,idat
        if (ares(i).ne.0.d0) then
            ncont=ncont+1
            nckp(ncont)=nkp(i)
            nckt(ncont)=nkt(i)
            ncpla(ncont)=npla(i)
        endif
      enddo
c resonances repeated were eliminated
c there are a total of ncont resonances  between  amax and amin


c file output
      open(2,file="resonances.dat",status="unknown",access="append")
      write(2,*) zcabe


      write(*,*)'calculating total of ',ncont,' resonances'

      do 777 ires=1,ncont
      write(*,*)'resonance ',ires




c index of the planet ncpla(ires)
        in=ncpla(ires)
c kt for test planet, kp for fixed planet
        ikt=nckt(ires)
        ikp=nckp(ires)
c define exterior planet  (i2) and interior (i1)
        if(ikt.gt.ikp) then
c test planet (10) is exterior
          k2=dfloat(ikt)
          k1=dfloat(ikp)
          i2=10
          i1=in
          else
c test planet (10) is coorbital or interior
          k2=dfloat(ikp)
          k1=dfloat(ikt)
          i1=10
          i2=in
        endif


c maximum for calculation of the integral with enough precision
      maxfac=k2


c factos to be used later
      fa01=sm(i1)/(mst+sm(i1))

c if test planet is exterior
      if(i1.ne.10) then
c semimajor axis interior fixed planet
      a1=se(i1)
c mean motion fixed planet  ignoring other planets than the 2 resonant planets
      n1=dsqrt(kg2*(mst+sm(i1))/a1**3)
c mean motion esterior test planet
      n2=k1*n1/k2
c semimajor axis a2 exterior test planet
      a2=(kg2*(mst+sm(i1)+sm(i2))/n2**2)**(1.d0/3.d0)
      endif
c if test planet in interior
      if(i1.eq.10) then
c exterior fixed planet:
      a2=se(i2)
c mean motion exterior fixed planet ignoring other planets than the 2 resonant planets
      n2=dsqrt(kg2*(mst+sm(i1)+sm(i2))/a2**3)
c mean motion interior test planet
      n1=k2*n2/k1
c semimajor axis interior test planet
      a1=(kg2*(mst+sm(i1))/n1**2)**(1.d0/3.d0)
      endif



c  loop over eccentricity of test planet
      iemax=dint((emax-emin)/deltae)+2
      do 80 ie = 1,iemax
       ex(10) = emin + deltae*dfloat(ie-1)

c exterior planet is 2
c hill
      rh(i2)=a2*(sm(i2)/(mst+sm(i2))/3.d0)**(1.d0/3.d0)

      e2=ex(i2)
      b2=a2*dsqrt(uno-e2*e2)
      j2g=yn(i2)
      j2=j2g*g2r
      l2g=ln(i2)
      p2g=lp(i2)
      l2=l2g*g2r
      p2=p2g*g2r
      ar2=p2-l2
c constants l1,m1,n1 for planet 2 roy page  93
      l12=dcos(l2)*dcos(ar2)-dsin(l2)*dsin(ar2)*dcos(j2)
      m12=dsin(l2)*dcos(ar2)+dcos(l2)*dsin(ar2)*dcos(j2)
      n12=dsin(ar2)*dsin(j2)

c constants l2,m2,n2 for planet 2, roy page  93
      l22=-dcos(l2)*dsin(ar2)-dsin(l2)*dcos(ar2)*dcos(j2)
      m22=-dsin(l2)*dsin(ar2)+dcos(l2)*dcos(ar2)*dcos(j2)
      n22=dcos(ar2)*dsin(j2)



c the interior planet is 1
c hill
      rh(i1)=a1*(sm(i1)/(mst+sm(i1))/3.d0)**(1.d0/3.d0)
      e1=ex(i1)
      b1=a1*dsqrt(uno-e1*e1)
      j1g=yn(i1)
      l1g=ln(i1)
      p1g=lp(i1)
      j1=j1g*g2r
      l1=l1g*g2r
      p1=p1g*g2r
      ar1=p1-l1
c constants l1,m1,n1 for planet 1 roy page  93
      l11=dcos(l1)*dcos(ar1)-dsin(l1)*dsin(ar1)*dcos(j1)
      m11=dsin(l1)*dcos(ar1)+dcos(l1)*dsin(ar1)*dcos(j1)
      n11=dsin(ar1)*dsin(j1)
c constants l2,m2,n2 for planet 1  roy page  93
      l21=-dcos(l1)*dsin(ar1)-dsin(l1)*dcos(ar1)*dcos(j1)
      m21=-dsin(l1)*dsin(ar1)+dcos(l1)*dcos(ar1)*dcos(j1)
      n21=dcos(ar1)*dsin(j1)

c second derivative hii
      beta1=sm(i1)*mst/(sm(i1)+mst)
      beta2=sm(i2)*(mst+sm(i1))/(sm(i2)+sm(i1)+mst)
      hii=3.d0*(k1/a1)**2/beta1 +  3.d0*(k2/a2)**2/beta2


c mutual rhill
c      rhmass=  sm(i2)/(mst+sm(i2)) + sm(i1)/(mst+sm(i1))
      rhmass=  (sm(i2)+ sm(i1))/mst
      rhs=(a2+a1)/2.d0*(rhmass/3.d0)**(1.d0/3.d0)

c number of evaluations of r(sigma) between o and 360 degrees
      isimax=360
c r(360)=r(0)
c steps of the numerical integration from 0 to 2pi*k2 in lambda particle
c can be as small as 100*int(maxfac), at your own risk
      ipasos=1000*int(maxfac)
c ==============================================================
c first calculate the data ext planet, they will be used isimax times
        do il=1,ipasos
c all possible configurations occur after k1 revolutions of the ext pla
c lambda pla 2
          la2=dfloat(il)/dfloat(ipasos)*twopi*dabs(k1)
c store vector data lambda pla2
          vla2(il)=la2
c se podria calcular en anom excentrica
c mean anomaly pla2
          am2=la2-p2
          am2=dmod(am2,twopi)
          if(am2.lt.cero) am2=am2+twopi
c solving kepler for exterior planet
          call solkep(e2,am2,aex2,initer)
          cose=dcos(aex2)
          sine=dsin(aex2)
c heliocentric distance
          r2=a2*(uno-e2*cose)
c formulas roy page 93
          x2=a2*l12*cose+b2*l22*sine-a2*e2*l12
          y2=a2*m12*cose+b2*m22*sine-a2*e2*m12
          z2=a2*n12*cose+b2*n22*sine-a2*e2*n12


c vector data rpla2
          vr2(il)=r2

c vector data xyz
          vx2(il)=x2
          vy2(il)=y2
          vz2(il)=z2



        enddo
c i have calculated all the positions for pla2 that i will use later
c ================================================================
c now define sigma1
      do isi=1,isimax
c acrit is sigma1 from 0 to 359
c close encounter indicator
        ace(isi)='         '
c libration period
        tlib(isi)=0.d0
        acrit=dfloat(isi-1)/dfloat(isimax)*360.d0
c theta is the combination of lambdas
        teta=acrit+(k1-k2)*p1g
        teta=dmod(teta,360.d0)
          if(teta.lt.cero) teta=teta+360.d0

c to radians
        tetar=teta*g2r
c disturbing function
        rtot=cero
c number of encounters inside the hill sphere
        irh(isi)=0
c minimum distance particle-planet
        mindis=9999999.d0
c given sigma1 calculate the integral
c calculation of the integral in lambda pla2
        do ilambda2=1,ipasos
c all possible configurations occur after k2 revolutions of the planet 1
c heliocentric distance
          r2=vr2(ilambda2)
c mean long between 0 and twopi*k2
          la2=vla2(ilambda2)
c radius vector for pla ext 2
          x2=vx2(ilambda2)
          y2=vy2(ilambda2)
          z2=vz2(ilambda2)
c velocity vector for pla ext 2
c          xp2=vxp2(ilambda2)
c          yp2=vyp2(ilambda2)
c          zp2=vzp2(ilambda2)

c given lambda pla exterior 2 calculate la1  planet  int 1
          la1=(k2*la2+tetar)/k1
c mean anomaly pla 1
c mean anomaly pla 1
          am1=la1-p1
          am1=dmod(am1,twopi)
          if(am1.lt.cero) am1=am1+twopi

c solving kepler for interior planet
          call solkep(e1,am1,aex1,initer)
          cose=dcos(aex1)
          sine=dsin(aex1)
c heliocentric distance
          r1=a1*(uno-e1*cose)
c formulas roy page 93
          x1=a1*l11*cose+b1*l21*sine-a1*e1*l11
          y1=a1*m11*cose+b1*m21*sine-a1*e1*m11
          z1=a1*n11*cose+b1*n21*sine-a1*e1*n11

c vector r02
          x02=fa01*x1+x2
          y02=fa01*y1+y2
          z02=fa01*z1+z2
c vector r12
          x12=x02-x1
          y12=y02-y1
          z12=z02-z1

c mutual distance
      delta02=dsqrt(x02**2+y02**2+z02**2)
      delta12=dsqrt(x12**2+y12**2+z12**2)


c register close encounters and number
          if(delta12.lt.rhs*rhtol) then
            ace(isi)='close enc'
            irh(isi)=irh(isi)+1
          endif
c register minimum distance planet-particle in rhill
          delrh=delta12/rhs
          if(delrh.lt.mindis)then
            mindis=delrh
          endif
c calculation disturbing function

          rper=mst/delta02 + sm(i1)/delta12 - (mst+sm(i1))/r2

c sumation for a crude integral
          rtot=rtot+rper
        enddo
          rtot=rtot*kg2*sm(i2)

c end of calculation of the integral
        rtot=rtot/dfloat(ipasos)
        rp(isi)=rtot
        sma(isi)=acrit
        md(isi)=mindis
      enddo
c end of calculation for all sigmas
c calculation of mean value <r>  for calculating strength sr
      vmedio=0.d0
      do is=1,isimax
        vmedio=vmedio+rp(is)
      enddo
      vmedio=vmedio/dfloat(isimax)
c find r_max and r_min
      vrmax=-999999.d0
      prmax=-999999.d0
      vrmin=999999.d0
      do is=1,isimax
        if(rp(is).lt.vrmin) then
          vrmin=rp(is)
        endif
c absolute rmax
        if(rp(is).gt.prmax) then
             prmax=rp(is)
        endif
c for calculating rmax we discard close encounters to less than rhtol*rhill
        if(rp(is).gt.vrmax.and.md(is).gt.rhtol) then
          vrmax=rp(is)
        endif
      enddo
c ready, rmax and rmin calculated discarding close encounters
c calculation of delta r
         delr=vrmax-vrmin
c if delta R is very small the resonance is unrealistic
c we compare delta R with the central term which is of the order of kg2*mst*(m1+m2)/a1
         cterm=kg2*mst*(sm(i1)+sm(i2))/a1
         if(delr/cterm.lt.1.d-18)   then
           delr=0.d0
         endif
c if delta R is very small we cannot detect any resonance
c 1.d-11 is a very reasonable limit for relative variation
c it is really a small variation
c      detec=delr/vmedio
      detec2=dabs((vmedio-vrmin)/vmedio)
c      detec=delr/vmedio
      if(delr.eq.0.d0.or.detec2.lt.1.d-11) then
         delr=0.d0
c        write(*,*)'        '
        write(*,*)'*** undetectable resonance ***'
c        write(*,*)'        '
        goto 111
      endif
c first derivatives in radians for each sigma
      do i=1,isimax
        isima1=i+1
        isime1=i-1
        if(isima1.gt.360) isima1=1
        if(isime1.lt.1) isime1=360
        deriv1(i)=(rp(isima1)-rp(isime1))/2.d0/g2r
      enddo
c second derivatives
      do i=1,isimax
        isima1=i+1
        isime1=i-1
        if(isima1.gt.360) isima1=1
        if(isime1.lt.1) isime1=360
        deriv2(i)=(deriv1(isima1)-deriv1(isime1))/2.d0/g2r
      enddo
c defining equilibrium points, no close encounters allowed
      do i=1,isimax
        aeq(i)='         '
        isime1=i-1
        if(isime1.lt.1) isime1=360
c condition: no close encounters at less than 0.5 rhill
c i am interested in detecting the location of the stable eq points
c even if they are inside 3*rhill, but not inside 0.5rhill
c        if(md(i).gt.0.5d0) then
c if the first derivative changes its sign
          valor=deriv1(i)
          if(deriv1(isime1)*valor.lt.0.d0)then
c decide if i or i-1 is the point
            iextremo=i
            if(deriv2(i).gt.0.)  then
              if(rp(isime1).lt.rp(i)) then
                iextremo=isime1
              endif
c it is a stable equilibrium point and libration period in years
              aeq(iextremo)='e. stable'

c libration period
      tlib(iextremo)=twopi/dsqrt(hii*dabs(deriv2(iextremo
     *)))/365.25d0
              else
c unstable
              if(rp(isime1).gt.rp(i))then
                iextremo=isime1
              endif
            aeq(iextremo)=' unstable'
          endif
c        endif
      endif
      enddo

c arrange stable points and periods
      ili=0
      do i=1,isimax
        if(tlib(i).ne.0.d0) then
        ili=ili+1
          tli(ili)=tlib(i)
          sig(ili)=sma(i)
        endif
      enddo
c ili=total number of stable eq points


c ancha1 is resonance full stable width  planet interior 1
 111  ff= 2.d0*dsqrt(2.d0*delr/hii/kg2)
      ancha1=k1*dsqrt(a1*(mst+sm(i1)))/mst/sm(i1)*ff
      ancha1=2.d0*ancha1
c width planet exterior 2
      ancha2=k2*dsqrt(a2*(mst+sm(i1)+sm(i2)))/(mst+sm(i1))/sm(i2)*ff
      ancha2=2.d0*ancha2

c file output

      write(2,19)int(k2),int(k1),(vmedio-vrmin),ex(10),a1,ancha1,a2,
     *ancha2,(sig(i),tli(i),i=1,ili)

c erase values tlib and sigma
            do i=1,10
          tli(i)=0.d0
          sig(i)=0.d0
      enddo




   80 continue
      write(2,*)'          '
  777 continue

      close(2)


   19 format(2i3,1p e14.6,0p f6.3,1p 4e11.4,0p f6.1,1p e11.4, 0p
     *f6.1,1p e10.3,0p f6.1,1p e10.3)
      end

c ========================================================================
c solving kepler
      subroutine solkep(ex,m,e,niter)
      implicit real*8 (a-h,o-z)
      real*8 m,mk
      tole=1.d-12
      e=m
      niter=0
 100  e0=e
      se=dsin(e0)
      ce=dcos(e0)
      es=ex*se
      ec=1.d0-ex*ce
      mk=e0-es
      u=(mk-m)/ec
      xpri=e0-u
      xseg=e0-u/(1.d0-u*es)
      e=(xpri+xseg)/2.d0
      dex=dabs(e-e0)
      niter=niter+1
      if(dex.gt.tole)goto 100
      return
      end

