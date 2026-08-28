c
c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
c
 
       subroutine strengthforce(grpmj,rhoij,
     $     sxxi,syyi,szzi,sxyi,sxzi,syzi,sxxj,syyj,szzj,
     $     sxyj,sxzj,syzj,dmi,dmj,dx,dy,dz,dfxi,dfyi,dfzi)
c****************************************************************
c This subroutine calculates the force from the strength module *
c Only call if there is strength between particles.             *
c****************************************************************
c
      implicit none
c--input:
c     rhoi: density
c     sxxi,syyi,sxyi,sxzi,syzi: stress tensor components for particle i
c     sxxj,syyj,sxyj,sxzj,syzj: stress tensor components for particle j
c     dx=xi-xj,dy=yi-yj,dz=zi-zj (distance of particles)
c     grpmj kernel gradient
c     redi,redj: reduction from damage
c     rhoij:  average rho
      real*4 rhoij,grpmj
      real*4 sxxi,syyi,szzi,sxyi,sxzi,syzi
      real*4 sxxj,syyj,szzj,sxyj,sxzj,syzj
      real*4 dmi,dmj
      real*4 dx,dy,dz
      real*8 inf
c
c--output:
c     fxi,fyi,fzi:  acceleration terms
c 
      real*4 dfxi,dfyi,dfzi
c
c--variables used
c    sigyyij,sigzzij,sigxyij,sigxzij,sigyzij
c    tx,ty,tz
c    redi,redj - reduction factor
c
      real*4 sigxxij,sigyyij,sigzzij,sigxyij,sigxzij,sigyzij
      real*4 tx,ty,tz
      real*4 redi,redj
c
c
c  material strength if pairs belong to same object
c
      redi=1.-dmi**3
      redj=1.-dmj**3
      sigxxij=(redi*sxxi + redj*sxxj)/rhoij
      sigyyij=(redi*syyi + redj*syyj)/rhoij
c      sigzzij=(redi*sxxi - redi*syyi + redj*sxxj - redj*syyj)/rhoij
      sigzzij=(redi*szzi + redj*szzj)/rhoij
      sigxyij=(sxyi + redj*sxyj)/rhoij
      sigxzij=(sxzi + redj*sxzj)/rhoij
      sigyzij=(syzi + redj*syzj)/rhoij
      tx=sigxxij*dx + sigxyij*dy + sigxzij*dz
      ty=sigxyij*dx + sigyyij*dy + sigyzij*dz
      tz=sigxzij*dx + sigyzij*dy + sigzzij*dz

      dfxi=grpmj*tx
      dfyi=grpmj*ty
      dfzi=grpmj*tz

      inf = HUGE (0.0)

      if (dfxi .gt. inf .or. dfxi .ne. dfxi) then
          write (*,*) "dfxi is ", dfxi
          write (*,*) "dx, dy, dz, redi, redj, sigxxij, sigxyij,",
     *     " sigxzij, tx"
          write (*,*) dx, dy, dz, redi, redj, sigxxij, sigxyij, sigxzij,
     *     tx
          write (*,*) "sxxi, sxxj, sxyi, sxyj, sxzi, sxzj"
          write (*,*) sxxi, sxxj, sxyi, sxyj, sxzi, sxzj
          write (*,*) "rhoij, dmi, dmj"
          write (*,*) rhoij, dmi, dmj
      endif

      if (dfyi .gt. inf .or. dfyi .ne. dfyi) then
          write (*,*) "dfyi is ", dfyi
          write (*,*) "dx, dy, dz, redi, redj, sigxyij, sigyyij,",
     *     " sigyzij, ty"
          write (*,*) dx, dy, dz, redi, redj, sigxyij, sigyyij, sigyzij,
     *     ty
          write (*,*) "sxyi, sxyj, syyi, syyj, syzi, syzj"
          write (*,*) sxyi, sxyj, syyi, syyj, syzi, syzj
          write (*,*) "rhoij, dmi, dmj"
          write (*,*) rhoij, dmi, dmj
      endif

      if (dfzi .gt. inf .or. dfzi .ne. dfzi) then
          write (*,*) "dfzi is ", dfzi
          write (*,*) "dx, dy, dz, redi, redj, sigxzij, sigyzij,",
     *     " sigzzij, tz"
          write (*,*) dx, dy, dz, redi, redj, sigxzij, sigyzij, sigzzij,
     *     tz
          write (*,*) "sxzi, sxzj, syzi, szyj, szzi, szzj"
          write (*,*) sxzi, sxzj, syzi, syzj, szzi, szzj
          write (*,*) "rhoij, dmi, dmj"
          write (*,*) rhoij, dmi, dmj
      endif

      return
      end

      subroutine deviator(xmui,sxxi,syyi,szzi,sxyi,sxzi,
     $     syzi,epsxxi,epsyyi,epszzi,epsxyi,epsxzi,epsyzi,
     $     rxyi,rxzi,ryzi,dsxxi,dsyyi,dszzi,dsxyi,dsxzi,dsyzi)
c***************************************************************
c                                                              *
c  subroutine to compute the derivative of the deviatoric      *
c  stress tensor.                                              *
c                                                              *
c***************************************************************
c
      implicit none
c--input:
c       sxxi,syyi,sxyi,sxzi,syzi: stress tensor components
c       epsxxi,epsyyi,epsxyi,epsxzi,epsyzi: strain rate tensor components
c       rxyi,ryzi,rxzi:  rotation tensor components
c       xmui:  shear modulus
      real*4 sxxi,syyi,szzi,sxyi,sxzi,syzi
      real*4 epsxxi,epsyyi,epszzi,epsxyi,epsxzi,epsyzi
      real*4 rxyi,ryzi,rxzi
      real*4 xmui
c
c--output
c      dsxxi,dsyyi,dsxyi,dsxzi,dsyzi: rate of change of the stress tensor
      real*4 dsxxi,dsyyi,dszzi,dsxyi,dsxzi,dsyzi
c
c--variables
c      div3:  mean of the strain rate tensor
c      xmu2: 2 times shear modulus 
      real*4 div3,xmu2
c
c--compute derivative of deviatoric stress tensor
c
       div3=(epsxxi+epsyyi+epszzi)/3.0
      xmu2=xmui*2.0
c
      dsxxi=xmu2*(epsxxi-div3)+2.*sxyi*rxyi+
     $     2.*sxzi*rxzi
      dsyyi=xmu2*(epsyyi-div3)-2.*sxyi*rxyi+
     $     2.*syzi*ryzi
      dszzi=xmu2*(epszzi-div3)-2.*sxzi*rxzi+
     $     2.*syzi*ryzi
      dsxyi=xmu2*epsxyi+(syyi-sxxi)*rxyi+
     $     sxzi*ryzi+syzi*rxzi
      dsxzi=xmu2*epsxzi+(szzi-sxxi)*rxzi-
     $     sxyi*ryzi+syzi*rxyi
      dsyzi=xmu2*epsyzi+(szzi-syyi)*ryzi-
     $     sxyi*rxzi-sxzi*rxyi
c
      return
      end

      subroutine straintensor(grpmrj,dvx,dvy,dvz,dx,dy,dz,
     $     depsxxi,depsyyi,depszzi,depsxyi,depsxzi,depsyzi,
     $     drxyi,drxzi,dryzi)
c
      implicit none
c
c--input:
c grpmrj: kernel gradient/rho
c velocity change (e.g. vxi-vxj):  dvx,dvy,dvz
c separation: dx,dy,dz
      double precision grpmrj
      double precision dvx,dvy,dvz
      double precision dx,dy,dz
c
c--output
c change in strain rate tensor: depsxxi,...
c change in rotation rate: drxyi,...
c
      double precision depsxxi,depsyyi,depszzi,
     $     depsxyi,depsxzi,depsyzi
      double precision drxyi,drxzi,dryzi
c
      depsxxi=grpmrj*dvx*dx
      depsyyi=grpmrj*dvy*dy
      depszzi=grpmrj*dvz*dz
      depsxyi=-grpmrj*0.5*(dvx*dy+dvy*dx)
      depsxzi=-grpmrj*0.5*(dvx*dz+dvz*dx)
      depsyzi=-grpmrj*0.5*(dvy*dz+dvz*dy)
c
c--rotation
c
      drxyi=-grpmrj*0.5*(dvx*dy-dvy*dx)
      drxzi=-grpmrj*0.5*(dvx*dz-dvz*dx)
      dryzi=-grpmrj*0.5*(dvy*dz-dvz*dy)
c
      return
      end

      subroutine strengthdu(rhoi,sxxi,syyi,szzi,sxyi,sxzi,syzi,
     $     epsxxi,epsyyi,epsxyi,epsxzi,epsyzi,epszzi,dmi,dudt)
c****************************************************************
c This subroutine calculates the change in energy (du/dt) from  *
c the strength module.                                          *
c****************************************************************
c
      implicit none
c--input:
c     rhoi: density
c     sxxi,syyi,sxyi,sxzi,syzi: stress tensor components
c     epsxxi,epsyyi,epsxyi,epsxzi,epsyzi: strain rate tensor components
c     redi: reduction from damage
      real*4 rhoi,dmi
      real*4 sxxi,syyi,szzi,sxyi,sxzi,syzi
      real*4 epsxxi,epsyyi,epszzi,epsxyi,epsxzi,epsyzi
c
      real*4 redi
c
c--output:
c     dudt: change in energy (du/dt)
c 
      real*4 dudt
c
      redi=1.0-dmi**3
c      dudt=sxxi*epsxxi+syyi*epsyyi+(-sxxi-syyi)*epszzi+
c     $     2.0*sxyi*epsxyi+2.0*sxzi*epszzi+2.0*syzi*epszzi
c     extend to full 3D, also, is there a typo on the 2nd line?
      dudt=sxxi*epsxxi+syyi*epsyyi+szzi*epszzi+
     $     2.0*sxyi*epsxyi+2.0*sxzi*epsxzi+2.0*syzi*epsyzi
      dudt=dudt/rhoi*redi

      return
      end
      subroutine fracture(sxxi,syyi,szzi,sxyi,sxzi,syzi,pri,dmi,nflawi,
     $     ifrac,youngi,epsmini,xmi,acoefi,ddmi,stmax)
c************************************************************************
c This code determines the initiation of fraction and evolves           *
c the damage.                                                           *
c  citations here                                                       *
c************************************************************************
c--input:
c       sxxi,syyi,sxyi,sxzi,syzi: stress tensor components
c       dmi:  damage
c       pri:  pressure
c       ifrac: fracture model...  only 1 viable at the moment
c
      real*4 sxxi,syyi,szzi,sxyi,sxzi,syzi
      real*4 dmi,pri
      integer ifrac
c
c-- material data
c     youngi: Young's modulus
c     epsmini:  weakest flaw activation stress
c     xmi: I think this is the 1/m in the activation thresholds
c     nflawi: number of flaws in particle
c     acoefi:  growth=cg/(2.*h(i))
      real*4 youngi,epsmini,xmi,acoefi
      integer nflawi
c
c variables used
c     sxxip,syyip,szzip,sxyip,sxzip,syzi: modified values
c     redi: reduction from damage
c     sig1,sig2,sig3:   stress eigenvalues
c     stmax,stmin:  max,min stress
c     pop:  damage factor
c     ddmi: change in damage
      real*4 redi
      real*4 sxxip,syyip,szzip,sxyip,sxzip,syzip
      real*4 stmax,stmin,sig1,sig2,sig3
      real*4 angle
c
c-variable set (change in dmi)
c
      real*4 ddmi
c
c--constants
      real*4 tiny,epsbig,pi2
      tiny=1.0d-20
      epsbig=1.0d10
      pi2=3.14159265359/2.0
c
cxxxx need grav(i) which is some pressure term
c
      ddmi=0.d0
      if (dmi.lt.1.0) then
         redi=1.0-dmi**3
         sxxip=redi*sxxi-pri
         syyip=redi*syyi-pri
         szzip=redi*szzi-pri
         sxyip=redi*sxyi
         syzip=redi*syzi
         sxzip=redi*sxzi

         call paxis(sxxip,syyip,szzip,sxyip,sxzip,syzip,
     $        sig1,sig2,sig3)
      
         stmax=max(sig1,sig2,sig3)
         stmin=min(sig1,sig2,sig3)

c         print *, sxxi,syyi,sxyi,syzi,sxzi,pri,redi
c         print *, sxxip,syyip,szzip,sxyip,sxzip,syzip
c         print *, "stress axis: ", sig1,sig2,sig3
c         print *, ifrac,youngi,stmax,epsmini,xmi
c
c--tensile failure (Weibull)
c  -------------------------
         if (ifrac.eq.1.or.ifrac.eq.3) then
            youngi=max(youngi*(1.0-dmi**3),tiny)
            if (stmax/youngi.gt.epsmini.and.epsmini.lt.epsbig) then
               pop=min((stmax/youngi/epsmini)**xmi,real(nflawi))
               ddmi=pop**(1.d0/3.d0)*acoefi
            end if
         end if
c         print *, stmax/youngi/epsmini,xmi,
c     $        (stmax/youngi/epsmini)**xmi,
c     $        stmax/youngi,epsmini,real(nflawi)
c         print *, 'ddmi set',ddmi,pop,acoefi
c         stop
c
c--shear failure (Mohr Coulomb)
c
         if (ifrac.eq.2.or.ifrac.eq.3) then
            if (stmax/youngi.gt.epsmini) then
c
c--  I think there is a typo in the Benz/Asphaug code in this and 
c--  I have guessed a solution, but I need to read more.
c
               angle=pi2 - atan(1.0/frictioni)
               ratiofsi=(stmax-stmin)*sin(angle)+
     $              frictioni*(0.5*(stmax+stmin)+
     $              (stmax-stmin)*cos(angle))
               ratiofsi=ratiofsi/(sigma0i*redi+tiny)
               if (ratiosfi.ge.1.0) then
                  ddmi=acoefi
               end if
            end if
         end if
      end if
      return
      end
      subroutine paxis(sxxi,syyi,szzi,sxyi,sxzi,syzi,sig1,sig2,
     1                 sig3)
c****************************************************************
c                                                               *
c  subroutine performing the principal axis transformation of a *
c  symmetric tensor. The 3 eigenvalues are returned.            *
c                                                               *
c****************************************************************
c
c-- input:
c       sxxi,syyi,szzi,sxyi,sxzi,syzi: stress tensor components
      real*4 sxxi,syyi,szzi,sxyi,sxzi,syzi
c
c-- output:
c     sig1,sig2,sig3:   stress eigenvalues
c       sxxt,syyt,sxyt,sxzt,syzt: normalized stress tensor components
      real*4 sig1,sig2,sig3
      real*4 sxxt,syyt,sxyt,sxzt,syzt
c
c-- variables used
c     smax(max of tensor component magnitudes)
      real*4 smax
      real*4 xi1,xi2,xi3
      real*4 a,b,a1,t1,phi,phi1,phi2,phi3
c
c-- constants
      real*4 pi2,pi4
      pi2=6.283185307
      pi4=12.56637061
c
      sig1=0.
      sig2=0.
      sig3=0.
c
c--find temporary norm to avoid overflow
c
      smax=max(abs(sxxi),abs(syyi),abs(szzi),abs(sxyi),abs(sxzi),
     1         abs(syzi))
      if(smax.ne.0.)then
         sxxt=sxxi/smax
         syyt=syyi/smax
         szzt=szzi/smax
         sxyt=sxyi/smax
         sxzt=sxzi/smax
         syzt=syzi/smax
c
c--compute 3 invariants
c
         xi1=sxxt+syyt+szzt
         xi2=-(syyt*szzt+szzt*sxxt+sxxt*syyt) +
     1         syzt*syzt + sxzt*sxzt + sxyt*sxyt
         xi3=sxxt*syyt*szzt+2.*syzt*sxzt*sxyt -
     1       sxxt*syzt*syzt - syyt*sxzt*sxzt -
     2       szzt*sxyt*sxyt
c
c--find eigenvalue by soving cubic equation
c  y**3 + -xi1*y**2 + -xi2*y + -xi3 =0
c  (the stress tensor being symmetric, roots are real)
c
c
         a=(-3.*xi2-xi1*xi1)/3.0
         b=(-2.*xi1*xi1*xi1-9.*xi1*xi2-27.*xi3)/27.
c
c--make sure there are 3 real roots (round offs)
c
         a1=a*a*a/27.
         rroots=0.25*b*b + a1
         if(rroots.lt.0.) then
            t1=2.*sqrt(-a/3.d0)
            phi=cos(-0.5*b/sqrt(-a1))
            phi1=phi/3.0
            phi2=(phi+pi2)/3.0
            phi3=(phi+pi4)/3.0
            sig1=(t1*cos(phi1)+xi1/3.0)*smax
            sig2=(t1*cos(phi2)+xi1/3.0)*smax
            sig3=(t1*cos(phi3)+xi1/3.0)*smax
         end if
      end if
c
      return
      end

      subroutine plastic(sxxi,syyi,szzi,sxyi,sxzi,syzi,ui,dmi,umelti,
     $     yiei,vonmisesi)
c************************************************************************
c  This subroutine applies the von Mises criterion to limit the         *
c deviatoric stress (von Mises, R., 1913,  Mechanik der festen Korperim *
c plastisch deformablen Zustand. Gottin. Nachr. Math. Phys., vol. 1,    *
c pp. 582-592., se also Von Mises yield criterion in Wikipedia          *
c                                                                       *
c************************************************************************
c--input:
c       sxxi,syyi,sxyi,sxzi,syzi: stress tensor components (also output)
c       dmi:  damage
c       ui:  internal energy
c
c
      real*4 sxxi,syyi,szzi,sxyi,sxzi,syzi
      real*4 ui,dmi
c-- material data
c       umelti:  melt energy for material
c       yiei:  yielding
      real*4 umelti,yiei
c
c--output:
c     sxxi,syyi,sxyi,sxzi,syzi: stress tensor components (also input)
c     vonmisesi - factor decreasing the stress due to yield
c     pri: pressure
      real*4 vonmisesi
c
c
c     variables used
c
c     unorm= ratio energy to melt energy (ui/umelti)
c     roundoff, tiny: lower limits
c     redi: reduce factor
c     
      real*4 unorm
      real*4 roundouff,tiny
      real*4 redi
      real*4 txxi,tyyi,tzzi,txyi,txzi,tyzi
c
      tiny=1.d-15
      roundoff=1.d-5
      
      unorm=ui/umelti
      
      if (unorm.lt.roundoff) then
         yst=yiei
      else
         yst=max(yiei*(1.d0-unorm),0.)
         if (yst.le.tiny) then
            sxxi=0.d0
            syyi=0.d0
            sxyi=0.d0
            sxzi=0.d0
            syzi=0.d0
            redi=1.0-dmi**3
            if (pri.lt.0) pri=redi*pri
            return
         end if
      end if
c
      redi=1.0-dmi**3
      txxi=(redi*sxxi+tiny)/yst + tiny
      tyyi=(redi*syyi+tiny)/yst + tiny
      tzzi=(redi*szzi+tiny)/yst + tiny
      txyi=(redi*sxyi+tiny)/yst + tiny
      txzi=(redi*sxzi+tiny)/yst + tiny
      tyzi=(redi*syzi+tiny)/yst + tiny
c
      xj2=0.5*(txxi*txxi+tyyi*tyyi+tzzi*tzzi) +
     $     txyi*txyi+txzi*txzi+tyzi*tyzi+tiny
      vonmisesi=min(dsqrt(1.d0/(3.d0*xj2)),1.d0)
      if (vonmisesi.lt.0.9) then
         print *, xj2,vonmisesi
      end if
c      sxxi=vonmisesi*sxxi
c      syyi=vonmisesi*syyi
c      sxyi=vonmisesi*sxyi
c      sxzi=vonmisesi*sxzi
c      syzi=vonmisesi*syzi
c
      return
      end












