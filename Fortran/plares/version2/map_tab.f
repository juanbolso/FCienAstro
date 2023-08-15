*     mapa_tab.f  - Gallardo, Beauge & Giuppone (2020) ; version 26/10/2020
*
*     Semi-analytical calculation of 2-planet mean-motion resonances between
*     bodies identified as i1 and i2 belonging to an (N=npl)-planet system.
*
*     User must specify initial outfile file name to be used for outputs,
*     typically the name of the system. Character variable is called 'sistema'.
*
*     Resonances will be calculated within a box of initial conditions in the
*     semimajor axis and eccentricity plane given with [amin,amax]x[emin,emax].
*
*     User is required to specify:
*                  - npl = total number of planets
*                  - values of i1,i2
*                  - index (iaemod) of planet whose orbit will be varied
*                  - values of amin,amax,emin,emax
*
*     Only MMRs of order < iqmax and order+degree < iqpmax will be considered.
*     Even so, only resonances with width > (amax-amin)/500 at e = (emin+max)/2
*     will be studied and results included in the output.
*
*     User must specify values of:         
*                  - iqmax
*                  - iqpmax
*
*     One output file is created for each resonance, giving relevant dynamical
*     features as function of the eccentricity between (emin,emax).
*     Columns of each output file are:
*     
*     ik2,ik1,emod,delw,a1in,a1out,a2in,a2out,neqs,sig_stab(i),tlib(i),Torbmax
*
*     where:
*
*     ik2,ik1 = resonance index (i.e. ik2*n2 - ik1*n1 ~ 0)
*     emod = eccentricity being varied
*     delw = fixed value of difference in longitudes of pericenter
*     a1in,a1out = inner and outer value of a1 of resonance separatrix
*     a2in,a2out = same, for a2 (outr resonant body)
*     neqs = number of stable fixed points
*     sig_stab(i) = value of main resonant angle of i-th stable fixed point
*     tlib(i) = libration period of fixed point
*     Torbmax = orbital period of outer resonant body
*
*     In the case of multiple stable fixed points, one will be printed per line
*      
      implicit real*8 (a-h,j-z)
      parameter(imax=50000)
      integer npl,j1,j2,neqs,nres,iq(200),ip(200)
      real*8 m(0:3),summ(3),mu(3),beta(3),a(3),e(3),inc(3),om(3),w(3)
      real*8 rh(10),sigma(imax),Rtot(imax),sig_stab(imax),tlib(imax)
      real*8 x2(imax),y2(imax),z2(imax),r2(imax),lambda2(imax)
      character(len=100) sistema,arch,anum,aden
      common /gen/ twopi,cero,uno,uno2,uno3,pi,g2r,kgaus,kg2
      common /orb/ m,mu,beta,a,e,inc,om,w,a1,a2,e1,e2
      common /ibd/ i1,i2,iaemod,ipasos,isimax
      common /res/ k1,k2
      common /rot/ l11,m11,n11,l21,m21,n21,l12,m12,n12,l22,m22,n22
      common /par/ distmin
c      
      twopi = 8.0d0*datan(1.0d0)
      cero  = 0.0d0
      uno   = 1.0d0
      uno2  = 1.0d0/2.0d0
      uno3  = 1.0d0/3.0d0
      pi    = twopi/2.d0
      g2r   = pi/180.d0
      kgaus = 0.01720209895d0
      kg2   = kgaus**2
      mjup  = 9.54d-4
      mter  = 3.04e-6

c     default integration parameters. 
      isimax = 360             ! number of evaluations of R(sigma) in [0,360]
      ipasos = 1000            ! number of steps of integration in lambda

c     general name of output file (e.g. name of the planetary system).
      sistema = 'mmr_d_cd'

c     mass of central star.
      m(0) = 0.96

c     number of secondary bodies.
      npl = 3
      
c     mass of secondary bodies [Msol].
      m(1) = 10.827893*mter
      m(2) = 15.039861*mter
      m(3) = 13.504369*mter
      
c     default initial conditions for m(1).
      a(1)   = 0.1254           ! semimajor axis [au]
      e(1)   = 0.137            ! eccentricity
      inc(1) = cero             ! inclination [rad]
      w(1)   = 47.0*g2r         ! longitude pericenter [rad]
      om(1)  = cero             ! longitude ascending node [rad]

c     default initial conditions for m(2).
      a(2)   = 0.2664           ! semimajor axis [au]
      e(2)   = 0.030            ! eccentricity
      inc(2) = cero             ! inclination [rad]
      w(2)   = 277.0*g2r        ! longitude pericenter [rad]
      om(2)  = cero             ! longitude ascending node [rad]
      
c     default initial conditions for m(3).
      a(3)   = 0.8121           ! semimajor axis [au]
      e(3)   = 0.596            ! eccentricity
      inc(3) = cero             ! inclination [rad]
      w(3)   = 183.3*g2r        ! longitude pericenter [rad]
      om(3)  = cero             ! longitude ascending node [rad]
       
c     mass factors for Jacobi reference system.
      summ(1) = m(0) + m(1)
      beta(1) = m(1)*m(0)/(m(0)+m(1))
      do i = 2,npl
       summ(i) = summ(i-1) + m(i)
       beta(i) = m(i)*summ(i-1)/summ(i)
      end do
      mu = kg2*summ

c     orbital period of outermost planet (in years).
      Torbmax = twopi/dsqrt(mu(npl)/a(npl)**3)/365.2563

c     define box limits in the semimajor axis domain (amin,amax).
      amin = 0.76
      amax = 0.85
      
c     define box limits in the eccentricity axes (emin,emax).
      emin = 0.45
      emax = 0.75
      
c     specify bodies involved in the MMR to be modeled (i1 < i2).
      i1 = 2
      i2 = 3
       
c     rotation matrix for m1 & m2 (Roy page 93).
      call rot_matrix (i1,l11,m11,n11,l21,m21,n21)
      call rot_matrix (i2,l12,m12,n12,l22,m22,n22)

c     specify whose semimajor axis and eccentricity will be varied (i1 or i2).
      iaemod = 3

c     find most relevant MMR with inner planet.      
      iqmax  = 7                ! maximum value allowed for order or degree
      iqpmax = 39               ! maximum value allowed for order+degree
      call resonances (iqmax,iqpmax,amin,amax,emin,emax,nres,ip,iq)

cccc  loop over relevant resonances.
      do 80 ires = 1,nres
      
c     specify MMR integers: k2*n(i2) - k1*n(i1) = 0.
       k2 = ip(ires) + iq(ires)
       k1 = ip(ires)
       
c     define name and open output file.
       write (anum,*) int(k2)
       write (aden,*) int(k1)
       anum = adjustl(anum)
       aden = adjustl(aden)
       arch = trim(sistema)//'_'//trim(anum)//'_'//trim(aden)//'.dat'
       arch = adjustl(arch)
       open (3,file=arch,status='replace')
       
c     depending on which body remains fixed, set values of (a1,a2) and (n1,n2).
       if (iaemod.lt.max(i1,i2)) then
        a2 = a(i2)
        n2 = dsqrt(mu(i2)/a(i2)**3) ! mean motion of m2
        n1 = k2*n2/k1               ! resonant n1
        a1 = (mu(i1)/n1**2)**uno3   ! resonant semimajor axis
       else
        a1 = a(i1)
        n1 = dsqrt(mu(i1)/a(i1)**3) ! mean motion of m1
        n2 = k1*n1/k2               ! resonant n2
        a2 = (mu(i2)/n2**2)**uno3   ! resonant semimajor axis
       end if
       
c     Hill radii.
       rh(i1) = a1*(m(i1)/summ(i1)/3.d0)**uno3
       rh(i2) = a2*(m(i2)/summ(i2)/3.d0)**uno3 
      
c     mininum allowed distance between bodies [au].
       rms = uno2*(a1+a2)*((m(i1)+m(i2))/3.0/m(0))**uno3
       distmin = 2.0*sqrt(3.0)*rms

cccc  loop over eccentricity of m(iaemod).
       iemax = 100
       do 50 ie = 1,iemax
        emod = emin + (emax-emin)*dfloat(ie-1)/dfloat(iemax-1)
        
c     set eccentricities.
        e1 = e(i1)
        e2 = e(i2)
        if (iaemod.eq.i1) e1 = emod
        if (iaemod.eq.i2) e2 = emod
        
c     calculate Rtot (mean disturbing function for sigma \in [0,twopi]).
        call rper (sigma,Rtot,Rmed,delR,Rstr)
        
c     number (neqs) of stable points, position (sig_stab) & libration periods.
        call analysis (Rtot,sigma,neqs,sig_stab,tlib)
        
c     libration widths in semimajor axis and mean-motion ratio.
        call widths (delR,a1in,a1out,a2in,a2out,n1n2in,n1n2out)
        
c     output of results.
        ik2 = k2
        ik1 = k1
        do i = 1,neqs
         write (*,100) ik2,ik1,emod,delw/g2r,a1in,a1out,a2in,a2out,
     *        neqs,sig_stab(i)/g2r,tlib(i),Torbmax
         write (3,100) ik2,ik1,emod,delw/g2r,a1in,a1out,a2in,a2out,
     *        neqs,sig_stab(i)/g2r,tlib(i),Torbmax
        end do
        flush (3)
c     
 50    continue
 51    continue
c
       close (3)
 80   continue
c
 100  format (2i5,1p6e15.5,i5,1p8e15.5)
c     
      end

      
      subroutine resonances(iqmax,iqpmax,amin,amax,emin,emax,nres,ip,iq)
      implicit real*8 (a-h,j-z)
      parameter(imax=50000)
      integer npl,j1,j2,neqs,jp,jq,iq(200),ip(200),nres
      real*8 m(0:3),summ(3),mu(3),beta(3),a(3),e(3),inc(3),om(3),w(3)
      real*8 sigma(imax),Rtot(imax)
      common /gen/ twopi,cero,uno,uno2,uno3,pi,g2r,kgaus,kg2
      common /orb/ m,mu,beta,a,e,inc,om,w,a1,a2,e1,e2
      common /ibd/ i1,i2,iaemod,ipasos,isimax
      common /res/ k1,k2
c
c     initialization.
      iq   = 0
      ip   = 0
      nres = 0                  ! counts number of MMR to be mapped

c     save initial semiajor axes and eccentricities.
      a1 = a(i1)
      a2 = a(i2)
      e1 = e(i1)
      e2 = e(i2)

c     calculate strength of MMR for e = (emin+emax)/2.
      if (iaemod.eq.i1) e1 = uno2*(emin+emax)
      if (iaemod.eq.i2) e2 = uno2*(emin+emax)

c     minimum and maximum value of n1/n2 over (amin,amax) interval.
      if (iaemod.eq.i1) then
       n1n2min = sqrt((mu(i1)/mu(i2))*(a(i2)/amax)**3)
       n1n2max = sqrt((mu(i1)/mu(i2))*(a(i2)/amin)**3)
      else
       n1n2min = sqrt((mu(i1)/mu(i2))*(amin/a(i1))**3)
       n1n2max = sqrt((mu(i1)/mu(i2))*(amax/a(i1))**3)
      end if

c     for n1/n2 = (p+q)/p within interval, loop over possible values of p & q.
      ipmax = int(n1n2max)
      do 20 jp = 1,ipmax
       do 10 jq = jp*(int(n1n2min)-1),jp*int(n1n2max)
        p = dfloat(jp)
        q = dfloat(jq)
        rat = (p+q)/p
        a2 = a1*(rat**(2.0/3.0))
        k2 = jp+jq
        k1 = jp
        if (rat.lt.n1n2min.or.rat.gt.n1n2max) goto 10
        if (min(jp,jq).gt.iqmax) goto 10
        if (jp+jq.gt.iqpmax) goto 10
c     calculate resonance width in pixels (assuming resolution=400).
        call rper (sigma,Rtot,Rmed,delR,Rstr)
        call widths (delR,a1in,a1out,a2in,a2out,n1n2in,n1n2out)
        widmax = 400.0*max((a1out-a1in)/a1,(a2out-a2in)/a2)/(amax-amin)
        if (nres.eq.0.and.widmax.gt.uno) then ! first time
         nres = 1
         ip(nres) = jp
         iq(nres) = jq
         write (*,*) ' MMR identified: (p+q),p,width =',jp+jq,jp,widmax
        else
         iflag = 0
         do i = 1,nres
          rati = abs(dfloat(ip(i)+iq(i))/dfloat(ip(i)))
          dif  = abs(rati - rat)
          if (dif.lt.1.0e-2) iflag = 1
         end do
         if (iflag.eq.0.and.widmax.gt.0.9) then  ! new relevant MMR identified
          nres = nres + 1
          ip(nres) = jp
          iq(nres) = jq
          write (*,*) ' MMR identified: (p+q),p,width =',jp+jq,jp,widmax
         end if
        end if
 10    continue
 20   continue
c
      return
      end

      
      subroutine widths (delr,a1in,a1out,a2in,a2out,n1n2in,n1n2out)
      implicit real*8 (a-h,j-z)
      real*8 m(0:3),mu(3),beta(3),a(3),e(3),inc(3),om(3),w(3)
      common /gen/ twopi,cero,uno,uno2,uno3,pi,g2r,kgaus,kg2
      common /orb/ m,mu,beta,a,e,inc,om,w,a1,a2,e1,e2
      common /ibd/ i1,i2,iaemod,ipasos,isimax
      common /res/ k1,k2
c
c     second derivative of unperturbed Hamiltonian (Hii).      
      Hii = 3.d0*(k1/a1)**2/beta(i1) + 3.d0*(k2/a2)**2/beta(i2)
      
c     libration width around inner planet (wid1).
      wid1 = 2.0*k1/m(0)/m(i1)*dsqrt(2.d0*delR*a1*(m(i1)+m(0))/Hii/kg2)
      
c     libration width around outer planet (wid2).
      wid2 = m(0)/m(i2)*dsqrt((m(i2)+m(i1)+m(0))*a2/(m(i1)+m(0))/a1)
      wid2 = (m(i1)/(m(0)+m(i1)))*(k2/k1)*wid2*wid1
      
c     semimajor axis of inner planet corresponding to each banch of separatrix.
      a1in  = a1 - wid1
      a1out = a1 + wid1
      
c     semimajor axis of outer planet corresponding to each banch of separatrix.
      a2in  = a2 - wid2
      a2out = a2 + wid2
      
c     mean-motion ratio corresponding to each banch of separatrix.
      facn = (m(0)+m(i1))/(m(0)+m(i1)+m(i2))
      n1n2in  = sqrt(facn*((a2-wid2)/a1)**3)
      n1n2out = sqrt(facn*((a2+wid2)/a1)**3)
c     
      return
      end
      
      
      subroutine analysis (Rtot,sigma,neqs,sig_stab,tlib)
      implicit real*8 (a-h,j-z)
      parameter(imax=50000)
      integer npl,j1,j2,neqs
      real*8 m(0:3),mu(3),beta(3),a(3),e(3),inc(3),om(3),w(3),rh(10)
      real*8 sig_stab(imax),sigma(imax),Rtot(imax),deriv1(imax)
      real*8 deriv2(imax),tlib(imax),x2(imax),y2(imax),z2(imax),r2(imax)
      real*8 lambda2(imax)
      common /gen/ twopi,cero,uno,uno2,uno3,pi,g2r,kgaus,kg2
      common /orb/ m,mu,beta,a,e,inc,om,w,a1,a2,e1,e2
      common /ibd/ i1,i2,iaemod,ipasos,isimax
      common /res/ k1,k2
      common /rot/ l11,m11,n11,l21,m21,n21,l12,m12,n12,l22,m22,n22
c      
      neqs = 0                  ! number of stable equilibrium solutions
      sig_stab = cero           ! value of sigma of each solution
      tlib = cero               ! corresponding zero-amplitude libration period

c     second derivative of unperturbed Hamiltonian (Hii).      
      Hii = 3.d0*(k1/a1)**2/beta(i1) + 3.d0*(k2/a2)**2/beta(i2)

c     first and second derivatives of Rtot(sigma) for each sigma [in radians].
      dels = twopi/dfloat(isimax)
      do i = 1,isimax
       ip1 = i + 1
       im1 = i - 1
       if (ip1.gt.isimax) ip1 = 1
       if (im1.lt.1)      im1 = isimax
       deriv1(i) = uno2*(Rtot(ip1)-Rtot(im1))/dels
       deriv2(i) = (Rtot(ip1)-2.0*Rtot(i)+Rtot(im1))/dels/dels
      end do 
      
c     search for equilibrium points as local maximum/minimum of Rtot(sigma).
      do 60 i = 1,isimax
       im1 = i - 1
       if (im1.lt.1) im1 = isimax
       if (deriv1(im1)*deriv1(i).lt.cero.or.deriv1(i).eq.cero) then
c     check stability (stable if d2R/dsigma2 > 0).
        if (deriv2(i).gt.cero) then ! stable
         istab = i              ! decide if i or i-1 is the point
         if (Rtot(im1).lt.Rtot(i)) istab = im1
         neqs = neqs + 1
         sig_stab(neqs) = sigma(istab)
         tlib(neqs) = twopi/dsqrt(Hii*dabs(deriv2(istab)))/365.25d0
        end if
       end if
 60   continue
c     
      return
      end

      
      subroutine rper (sigma,Rtot,Rmed,delR,Rstr)
      implicit real*8 (a-h,j-z)
      parameter(imax=50000)
      integer npl,j1,j2
      real*8 m(0:3),mu(3),beta(3),a(3),e(3),inc(3),om(3),w(3)
      real*8 sigma(imax),Rtot(imax),delta_min(imax)
      real*8 lambda2(imax),x2(imax),y2(imax),z2(imax),r2(imax)
      common /gen/ twopi,cero,uno,uno2,uno3,pi,g2r,kgaus,kg2
      common /orb/ m,mu,beta,a,e,inc,om,w,a1,a2,e1,e2
      common /ibd/ i1,i2,iaemod,ipasos,isimax
      common /res/ k1,k2
      common /rot/ l11,m11,n11,l21,m21,n21,l12,m12,n12,l22,m22,n22
      common /par/ distmin
c
      b1   = a1*dsqrt(uno-e1*e1)
      Rtot = cero
      
c     calculate position of m(i2) at ipasos points of its orbit.
      call position_m2 (i2,ipasos,lambda2,r2,x2,y2,z2)

c     initialize minimum distance between bodies.
      delta_min = 1.0d11
      
c     loop over resonant angle sigma = k1*lambda1 - k2*lambda2 + (k2-k1)*varpi1.
      do 50 isi = 1,isimax
       sigma(isi) = twopi*dfloat(isi-1)/dfloat(isimax)
       
c     theta = k1*lambda1 - k2*lambda2.
       teta = sigma(isi) - (k2-k1)*w(i1)
       teta = mod(teta,twopi)
       if (teta.lt.cero) teta = teta + twopi
       
c     loop over mean longitud of outer body.
       do 40 ilambda2 = 1,ipasos*int(k2)
        
c     given mean longitude of m2, calculate that of m1 from value of theta.
        lambda1 = (k2*lambda2(ilambda2)+teta)/k1
        
c     mean anomaly of m1.
        anom1 = lambda1 - w(i1)
        anom1 = mod(anom1,twopi)
        if (anom1.lt.cero) anom1 = anom1 + twopi
        
c     transform to eccentric anomaly.
        call solkep (e1,anom1,aex1,initer)
        cose = dcos(aex1)
        sine = dsin(aex1)
        
c     modulus of position vector.
        r1 = a1*(uno-e1*cose)
        
c     catersian coordinates.
        x1 = a1*l11*cose + b1*l21*sine - a1*e1*l11
        y1 = a1*m11*cose + b1*m21*sine - a1*e1*m11
        z1 = a1*n11*cose + b1*n21*sine - a1*e1*n11
        
c     relative position between m0 and m2.
        x02 = x1*m(i1)/(m(0)+m(i1)) + x2(ilambda2)
        y02 = y1*m(i1)/(m(0)+m(i1)) + y2(ilambda2)
        z02 = z1*m(i1)/(m(0)+m(i1)) + z2(ilambda2)
        
c     relative position between m1 and m2.
        x12 = x02 - x1
        y12 = y02 - y1
        z12 = z02 - z1
        
c     mutual distance between bodies.
        delta02 = dsqrt(x02**2+y02**2+z02**2)
        delta12 = dsqrt(x12**2+y12**2+z12**2)

c     update minimum distance for each value of resonant angle.
        delta_min(isi) = min(delta_min(isi),delta12)
        
c     calculate disturbing function.
        r_per = m(0)/delta02 + m(i1)/delta12 - (m(0)+m(i1))/r2(ilambda2)
        
c     sumation for a crude integral.
        Rtot(isi) = Rtot(isi) + r_per
 40    continue
       Rtot(isi) = Rtot(isi)*kg2*m(i2)/dfloat(ipasos)/k2
       
 50   continue
       
c     minimum value of Rtot over sigma: Rmin = min(Rtot(sigma)).
      Rmin = 1.0d30
      do isi = 1,isimax
       if (Rtot(isi).lt.Rmin) then
        Rmin = Rtot(isi)
        isi_Rmin = isi
       end if
      end do
       
c     maximum value of Rtot(sigma), excluding region beyond collision curve.
      Rmax = Rmin
      isi = isi_Rmin + 1
      do while (isi.lt.isimax.and.delta_min(isi).gt.distmin)
       Rmax = max(Rmax,Rtot(isi))
       isi = isi + 1
      end do
      isi = isi_Rmin - 1
      do while (isi.gt.0.and.delta_min(isi).gt.distmin)
       Rmax = max(Rmax,Rtot(isi))
       isi = isi - 1
      end do
      
c     average value (not associated to close encounter) of Rtot(sigma).
      Rmed = 0.0d0
      do isi = 1,isimax
       isi0 = isi - 1
       isi1 = isi + 1
       if (isi0.lt.0) isi0 = isimax
       if (isi1.gt.isimax) isi1 = 1
       dmin_min = min(delta_min(isi),delta_min(isi0),delta_min(isi1))
       if (dmin_min.gt.distmin) Rmed = Rmed + Rtot(isi)/dfloat(isimax)
      end do
       
c     amplitude of <Rtot> (i.e. delR) and resonance strength (Rstr).
       delR = Rmax - Rmin
       delR = max(delR,cero)    ! just in case....
       Rstr = Rmed - Rmin
c      
      return
      end
      
      
      subroutine position_m2 (i2,ipasos,lambda2,r2,x2,y2,z2)
      implicit real*8 (a-h,j-z)
      parameter(imax=50000)
      integer j2
      real*8 m(0:3),mu(3),beta(3),a(3),e(3),inc(3),om(3),w(3)
      real*8 lambda2(imax),r2(imax),x2(imax),y2(imax),z2(imax)
      common /gen/ twopi,cero,uno,uno2,uno3,pi,g2r,kgaus,kg2
      common /orb/ m,mu,beta,a,e,inc,om,w,a1,a2,e1,e2
      common /res/ k1,k2
      common /rot/ l11,m11,n11,l21,m21,n21,l12,m12,n12,l22,m22,n22
c     
      b2 = a2*dsqrt(uno-e2*e2)
c
      do 20 j2 = 1,ipasos*int(k2) ! evaluates over k1 periods of m2
       lambda2(j2) = k1*twopi*dfloat(j2)/dfloat(ipasos)/k2
       
c     mean and eccentric anomalies of m2.
       anom2 = mod(lambda2(j2)-w(i2),twopi)
       if (anom2.lt.cero) anom2 = anom2 + twopi
       call solkep (e2,anom2,aex2,initer)
       cose = dcos(aex2)
       sine = dsin(aex2)
       
c     modulus of position verctor.
       r2(j2) = a2*(uno-e2*cose)
       
c     cartesian coordinates (Roy page 93).
       x2(j2) = a2*l12*cose + b2*l22*sine - a2*e2*l12
       y2(j2) = a2*m12*cose + b2*m22*sine - a2*e2*m12
       z2(j2) = a2*n12*cose + b2*n22*sine - a2*e2*n12
       
 20   continue
c     
      return
      end

      
      subroutine rot_matrix (i,l1i,m1i,n1i,l2i,m2i,n2i)
      implicit real*8 (a-h,j-z)
      real*8 m(0:3),mu(3),beta(3),a(3),e(3),inc(3),om(3),w(3)
      common /orb/ m,mu,beta,a,e,inc,om,w,a1,a2,e1,e2
c
      arg = w(i)-om(i)
      l1i = dcos(om(i))*dcos(arg) - dsin(om(i))*dsin(arg)*dcos(inc(i))
      m1i = dsin(om(i))*dcos(arg) + dcos(om(i))*dsin(arg)*dcos(inc(i))
      n1i = dsin(arg)*dsin(inc(i))
      l2i =-dcos(om(i))*dsin(arg) - dsin(om(i))*dcos(arg)*dcos(inc(i))
      m2i =-dsin(om(i))*dsin(arg) + dcos(om(i))*dcos(arg)*dcos(inc(i))
      n2i = dcos(arg)*dsin(inc(i))
c
      return
      end
      
      
      subroutine solkep (ex,m,e,niter) ! Kepler's equation
      implicit real*8 (a-h,o-z)
      real*8 m,mk
c      
      tole = 1.d-12
      twopi = 8.d0*datan(1.d0)
      m = dmod(m,twopi)
      e = m
      niter = 0
 100  e0 = e
      se = dsin(e0)
      ce = dcos(e0)
      es = ex*se
      ec = 1.d0 - ex*ce
      mk = e0 - es
      u = (mk-m)/ec
      xpri = e0 - u
      xseg = e0 - u/(1.d0-u*es)
      e = (xpri+xseg)/2.d0
      dex = dabs(e-e0)
      niter = niter+1
      if (dex.gt.tole) goto 100
c      
      return
      end

