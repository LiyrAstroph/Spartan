      subroutine energycal(rdp,sigp,wip,wep,ird) 
      implicit none
      include 'const.f' 
      include 'diskvar.f' 
     
      real*8 rdp,sigp,wip,wep
      integer ird
      
      real*8 radheat,radcool,adiabatic_gam
      real*8 a13,c1,a21,a22,a23,c2,a31,a32,a33,c3
      real*8 Omega,Omegak1,Omegak2,lrdp,macc,csp,Lamb,Qrad
     &,      aspin,ametric,omegas,delta,Btemp,V,gamp,gamp_old,gamr,eng
     &,      mu,ai,ae,gi,ge,Gbi,Gbe,thetae,thetai,dai,dae
      real*8 tep,tip,nep,nip,heigp,rhop
      real*8 dsig,dwi,dwe,bb1,bb2,aa,tempa,tempc
     &,      x1,x2,x3,y1,y2,num,den,aa1,aa2,S
      real*8 brem,syn,comp,tot         
      macc=mdot
      aspin=astar/2.0d0
      ametric=1.0d0 + aspin**2.0d0 / rdp**2.0d0
     &              + aspin**2.0d0 / rdp**3.0d0
      delta = 1.0d0 - 1.0d0 /rdp + aspin**2.0d0 /rdp**2.0d0 
      
      omegas = aspin / rdp**3.0d0 / ametric
      csp=dsqrt((wip+wep)/sigp/beta)      
      Omegak1= 1.0d0/(dsqrt(2.0d0 * rdp**3.0d0) + aspin)
      Omegak2=-1.0d0/(dsqrt(2.0d0 * rdp**3.0d0) - aspin)
      
      Btemp= macc**2.0d0 / (4.0d0*pi*pi) 
     &                   / (delta * rdp**2.0d0 * sigp**2.0d0)
      if(Btemp.lt.1.0d-6)then
      V=dsqrt(Btemp*(1.0d0-Btemp))
      else
      V = dsqrt(Btemp / ( 1.0d0+Btemp ) )
      end if
      if(V.lt.1.0d-3)then
      gamr=dsqrt(1.0d0+V*V)
      else
      gamr = 1.0d0 / dsqrt(1.0d0 - V*V)
      end if
      
      Tep=mue*mp*wep/k/sigp*c**2.0d0
      tip=mui*mp*wip/k/sigp*c**2.0d0
      thetae=k*tep/mec2
      thetai=k*tip/mpc2
      
      gi=adiabatic_gam(thetai)
      ge=adiabatic_gam(thetae)
      ai=1.0d0/(gi-1.0d0) + 2.0d0*(1.0d0-beta)/beta
      ae=1.0d0/(ge-1.0d0) + 2.0d0*(1.0d0-beta)/beta
      call daie(tip,tep,dai,dae)
      Gbi= 1.0d0 + 1.0d0 / (ai*(1.0d0 + dai))
      Gbe= 1.0d0 + 1.0d0 / (ae*(1.0d0 + dae))
      
      mu=1.0d0 + ( (ai+1.0d0/beta)*wip + (ae+1.0d0/beta)*wep )/sigp
c-----------------------------------------------------------------------
      gamp=1.0d0
70    gamp_old=gamp 
      lrdp=lin + 2.0d0*pi*rdp**2.0d0 / macc * alpha * ametric**1.5d0
     &         * delta**0.5d0 * (wip+wep)/beta * gamp**3.0d0
      gamp=dsqrt(1.0d0 + lrdp**2.0d0 / (mu*gamr*rdp)**2.0d0 /ametric)
      if(lrdp.gt.1.0d2*angum(ird))then
      lrdp=angum(ird)
      gamp=dsqrt(1.0d0 + lrdp**2.0d0 / (mu*gamr*rdp)**2.0d0 /ametric)
      goto 100
      endif      
      if(dabs(gamp-gamp_old).gt.1.0d-6)goto 70
c------------------------------------------------------------------------      
100   Omega=omegas + delta**0.5d0*lrdp 
     &              /(mu*ametric**1.5d0*rdp**2.0d0*gamr*gamp)
     
      eng= mu*gamr*gamp*dsqrt(delta/ametric) + omegas*lrdp
      aa1=lrdp**2.0d0
      aa2=aspin**2.0d0*(eng**2.0d0-mu*mu)
      if(aa1.le.aa2)then
      heigp= csp *rdp**2.0 * dsqrt(mu)
     &      /dsqrt(dabs(aa1-aa2))
      else
      heigp= csp *rdp**2.0 * dsqrt(mu)
     &      /dsqrt(aa1- aa2)
      end if
      if(heigp/rdp.gt.2.0d0)heigp=2.0d0*rdp 
      heigp=heigp*rs
      rhop=sigp/dsqrt(2.0*pi)/heigp*medd/rs/c
      nep=rhop/mue/mp
      nip=rhop/mui/mp
      
      call flux_eng(tep,tip,rhop,heigp,csp*c,brem,syn,comp,tot)
      write(nfenergy,*)log10(rd(ird)),log10(2.0*pi*rd(ird)*rd(ird)*tot)
     &                  +2.0*log10(rs)
     &,          log10(2.0*pi*rd(ird)*rd(ird)*brem)+2.0*log10(rs)
     &,          log10(2.0*pi*rd(ird)*rd(ird)*syn)+2.0*log10(rs)
     &,          log10(2.0*pi*rd(ird)*rd(ird)*comp)+2.0*log10(rs)
      end subroutine     
      
c===================================
c subroutine flux_nu
c===================================     
      subroutine flux_eng(tep,tip,rhop,heigp,csp,brem,syn,comp,tot)
      implicit none
      include 'const.f'
      include 'diskvar.f'
      real*8 tep,tip,rhop,heigp,csp,brem,syn,comp,tot,totemp
      real*8 bessk,gammp
      real*8 thetae,ftheta,qbr,B,ne,rho
      
      real*8 temp_bnu,temp_nu,temp_gaun,temp_syn,temp_sync
     &,      bnu,gaunt,xm,fxm,chibr,chisyn,chinu,kappa,taunu
c
c for comptonization
c

      real*8 A,s,yitamax,taues,taueff,yita,jm
      
      integer inu
            
      thetae=k*Tep/mec2

      B=dsqrt(8.0d0*pi*(1.0d0-beta)*rhop)*Csp
      ne=rhop/mue/mp
c      write(*,*)"B,thetae,heig,rho,ne "
c      write(*,*)B,thetae,heig,rho,ne      
c-------------------------------------------
c cal fnu
c      
c bremsstrahlung 
      
      if(thetae.lt.1.0d0)then
      ftheta=4.0*dsqrt(2.0*thetae/pi**3)*(1.0+1.781*thetae**(1.34d0))
     &        +1.73*thetae**(1.5d0)*(1.0+1.1*thetae+thetae**2
     &                              -1.25*thetae**(2.5d0))
      else
      ftheta=4.5d0*thetae/pi*(dlog(0.48+1.123*thetae)+1.5)
     &       +2.30*thetae*(dlog(1.123*thetae)+1.28)
      end if
      qbr=1.48d-22*ne**2*ftheta

      temp_bnu=2.0*h/c**2      
      temp_gaun=h/k/tep

      temp_syn=4.0*pi*me*c/3.0/e/B/thetae**2
      
      if(bessk(2,1.0/thetae).lt.1.0d-150)then
      
      write(*,*)"bessk(2,1.0/thetae).lt.1.0d-150"
      temp_sync=1.0d-150
c      stop
      end if
      temp_sync=4.43d-30*4.0*pi*ne/bessk(2,1.0/thetae)      

      tot=0.0
      brem=0.0
      syn=0.0
      comp=0.0
      
      do 10 inu=2,nnu
      
c black body radiation
      temp_nu=h*nu(inu)/k/tep

      if(temp_nu.lt.100.0d0)then
      bnu=temp_bnu*nu(inu)**3/(dexp(temp_nu)-1.0d0)
      else
      bnu=temp_bnu*nu(inu)**3/(dexp(temp_nu/3.0)**3.0d0-1.0d0)
      end if
      
c cal the gaunt factor
      if(temp_nu.lt.1.0d0)then
      gaunt=temp_gaun*dsqrt(3.0d0)/pi*dlog(4.0d0/1.781d0/temp_nu)
      else
      gaunt=temp_gaun*dsqrt(3.0/pi/temp_nu)
      endif
      chibr=qbr*gaunt*dexp(-temp_nu)
c cal the synchrotron radiation
      
      xm=temp_syn*nu(inu)
      fxm=4.0505/xm**(1.0/6.0d0)
     &*(1.0+0.4/xm**(0.25d0)+0.5316/xm**(0.5d0))
     &*dexp(-1.8899d0*xm**(1.0d0/3.0))
      chisyn=temp_sync*nu(inu)*fxm
      
c cal the flux.      
      chinu=chibr
      if(bnu.gt.eps)then
      kappa=chinu/4.0/pi/bnu
      taunu=dsqrt(pi)/2.0*kappa*heigp
      if(taunu.gt.1.0d-6)then
      totemp=2.0*pi/dsqrt(3.0d0)*bnu*(1.0
     &             -dexp(-2.0*dsqrt(3.0d0)*taunu))
      else
      totemp=2.0*pi/dsqrt(3.0d0)*bnu*(2.0*dsqrt(3.0d0)*taunu
     &            -6.0*taunu**2+4.0*dsqrt(3.0d0)*taunu**3)
      end if 
      
      else

      totemp=dsqrt(pi)*heigp*chinu
      end if
      brem=brem+totemp*(nu(inu)-nu(inu-1))
c*************** 
      chinu=chisyn     
      if(bnu.gt.eps)then
      kappa=chinu/4.0/pi/bnu
      taunu=dsqrt(pi)/2.0*kappa*heigp
      if(taunu.gt.1.0d-6)then
      totemp=2.0*pi/dsqrt(3.0d0)*bnu*(1.0
     &             -dexp(-2.0*dsqrt(3.0d0)*taunu))
      else
      totemp=2.0*pi/dsqrt(3.0d0)*bnu*(2.0*dsqrt(3.0d0)*taunu
     &            -6.0*taunu**2+4.0*dsqrt(3.0d0)*taunu**3)
      end if 
      
      else
      totemp=dsqrt(pi)*heigp*chinu
      end if
      totemp=max(totemp,1.0d-150)
      syn=syn+totemp*(nu(inu)-nu(inu-1)) 
      
c******************      
      chinu=chibr+chisyn

c      write(80,*)kappa,bnu,chinu
      
      if(bnu.gt.eps)then
      kappa=chinu/4.0/pi/bnu
      taunu=dsqrt(pi)/2.0*kappa*heigp
      if(taunu.gt.1.0d-6)then
      totemp=2.0*pi/dsqrt(3.0d0)*bnu*(1.0
     &             -dexp(-2.0*dsqrt(3.0d0)*taunu))
      else
      totemp=2.0*pi/dsqrt(3.0d0)*bnu*(2.0*dsqrt(3.0d0)*taunu
     &            -6.0*taunu**2+4.0*dsqrt(3.0d0)*taunu**3)
      end if 
      
      else
      totemp=dsqrt(pi)*heigp*chinu
      end if

c##########
c cal energy enhancement factor due to comptonization
c

      
      A=1.0+4.0d0*thetae+16.0d0*thetae**2
      taues=dsqrt(2.0*PI)*ne*sigmat*heigp
      s=taues+taues**2.0
      yitamax=3.0d0/temp_nu
      if(yitamax.ge.1.0d0)then
      jm=dlog(yitamax)/dlog(A)
      yita=dexp(s*(A-1.0d0))*(1.0d0-gammp(jm+1,A*s))
     &     +yitamax*gammp(jm+1,s)
      yita=max(1.0d0,yita)
      else
      jm=0.0d0
      yita=1.0d0      
      end if  
      totemp=totemp*yita
      totemp=max(totemp,0.0d0)
      tot=tot+totemp*(nu(inu)-nu(inu-1))
      comp=comp+totemp*(yita-1.0) 
              
10    continue     

      return               
      end subroutine flux_eng
