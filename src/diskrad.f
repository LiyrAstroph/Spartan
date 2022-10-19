c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                 supported by Prof. Wang Jian-Min
c                          and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-01-
c======================================================================= 

c To calculate the radiaiton processes for disk structure.

c===================================
c  function radcool                
c  cal the cooling flux            
c===================================
      function radcool(tep,tip,rhop,heigp,csp)
      implicit none
      include 'diskvar.f'       
      real*8 radcool,tep,tip,rhop,heigp,csp
      real*8 temp
      integer inu
      call flux_nu(tep,tip,rhop,heigp,csp)
      temp=0.0d0
      do inu=1,nnu
      temp=temp+wleg(inu)*flux(inu)*2.0*nu(inu)*dlog(10.0d0)
c      write(99,*)"temp=",temp
      end do
      radcool=temp      
      end function radcool
c===================================
c subroutine flux_nu
c===================================     
      subroutine flux_nu(tep,tip,rhop,heigp,csp)
      implicit none
      include 'const.f'
      include 'diskvar.f'
      real*8 tep,tip,rhop,heigp,csp
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

      do 10 inu=1,nnu
      
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
      chinu=chibr+chisyn

c      write(80,*)kappa,bnu,chinu
      
      if(bnu.gt.eps)then
      kappa=chinu/4.0/pi/bnu
      taunu=dsqrt(pi)/2.0*kappa*heigp
      if(taunu.gt.1.0d-6)then
      flux(inu)=2.0*pi/dsqrt(3.0d0)*bnu*(1.0
     &             -dexp(-2.0*dsqrt(3.0d0)*taunu))
      else
      flux(inu)=2.0*pi/dsqrt(3.0d0)*bnu*(2.0*dsqrt(3.0d0)*taunu
     &            -6.0*taunu**2+4.0*dsqrt(3.0d0)*taunu**3)
      end if 
      
      else
c      kappa=0.0d0
c      taunu=dsqrt(pi)/2.0*kappa*heig
      flux(inu)=dsqrt(pi)*heigp*chinu
      end if

c##########
c cal energy enhancement factor due to comptonization
c

      
      A=1.0+4.0d0*thetae+16.0d0*thetae**2

c      taueff=dsqrt(taunu*(taunu+ne*sigmat*heig))
c      taues=2.0*ne*sigmat*heigp
       taues=dsqrt(2.0*PI)*ne*sigmat*heigp
c      *max(1.0d0,taueff)
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
c      write(99,*)"jm=",jm," s=",s," taueff=",taueff
c     &," yita=",yita
c      if(ird.eq.50)then
c      write(99,*)gammp(jm+1,A*s),gammp(jm+1,s)
c      write(99, *)yita,"   ",jm 
c      end if    
      flux(inu)=flux(inu)*yita
      flux(inu)=max(flux(inu),0.0d0)          
10    continue     

      return               
      end subroutine flux_nu
c===================================
c function radheat
c===================================      
      function radheat(te,ti,ne,ni)
      implicit none
      include 'const.f'
      real*8 ti,te,ne,ni,radheat
      real*8 bessk,bessk0,bessk1
      real*8 thetai,thetae,tempi,tempe,tempie,temp
      
      thetai=k*ti/mpc2
      thetae=k*te/mec2
      if((thetae.ge.1.0d-2).and.(thetai.ge.1.0d-2))then
      radheat=5.61d-32*ne*ni*(ti-te)
     &        /(bessk(2,1.0/thetae)*bessk(2,1.0/thetai))
     &        *((2.0*(thetae+thetai)**2+1.0)/(thetae+thetai)
     &        *bessk1((thetae+thetai)/thetae/thetai)
     &        +2.0*bessk0((thetae+thetai)/thetae/thetai))     
      return
      end if
      
      temp=thetae*thetai/(thetae+thetai) 
           
      if((thetae.ge.1.0d-2).and.(thetai.lt.1.0d-2))then
      radheat=5.61d-32*ne*ni*(ti-te)
     &/bessk(2,1.0/thetae)*dexp(-1.0/thetae)
     &*dsqrt(thetae/(thetae+thetai))
     &*((2.0*(thetae+thetai)**2+1)/(thetae+thetai)*
     &(1.0+3.0/8.0*temp-15.0/128.0*temp**2+15.0*21.0/6.0/8.0**3*temp**3) 
     &+2.0*(1.0-1.0/8.0*temp+9.0/128.0*temp**2
     &      -9.0*25.0/6.0/8.0**3*temp**3))
     & /(1.0+15.0/8.0*thetai+15.0*7.0/128.0*thetai**2
     &   -15.0*7.0*9.0/6.0/8.0**3*thetai**3)          
      return
      end if 
      
      if((thetae.lt.1.0d-2).and.(thetai.ge.1.0d-2))then
      radheat=5.61d-32*ne*ni*(ti-te)
     &/bessk(2,1.0/thetai)*dexp(-1.0/thetai)
     &*dsqrt(thetai/(thetae+thetai))
     &*((2.0*(thetae+thetai)**2+1)/(thetae+thetai)*
     &(1.0+3.0/8.0*temp-15.0/128.0*temp**2+15.0*21.0/6.0/8.0**3*temp**3) 
     &+2.0*(1.0-1.0/8.0*temp+9.0/128.0*temp**2
     &      -9.0*25.0/6.0/8.0**3*temp**3)) 
     & /(1.0+15.0/8.0*thetae+15.0*7.0/128.0*thetae**2
     &   -15.0*7.0*9.0/6.0/8.0**3*thetae**3)           
      return
      end if          
      
      if((thetae.lt.1.0d-2).and.(thetai.lt.1.0d-2))then
      radheat=5.61d-32*ne*ni*(ti-te)
     &*dsqrt(2.0/pi)/dsqrt(thetae+thetai)
     &*((2.0*(thetae+thetai)**2+1)/(thetae+thetai)*
     &(1.0+3.0/8.0*temp-15.0/128.0*temp**2+15.0*21.0/6.0/8.0**3*temp**3) 
     &+2.0*(1.0-1.0/8.0*temp+9.0/128.0*temp**2
     &      -9.0*25.0/6.0/8.0**3*temp**3))
     & /(1.0+15.0/8.0*thetai+15.0*7.0/128.0*thetai**2
     &   -15.0*7.0*9.0/6.0/8.0**3*thetai**3)  
     & /(1.0+15.0/8.0*thetae+15.0*7.0/128.0*thetae**2
     &   -15.0*7.0*9.0/6.0/8.0**3*thetae**3)             
      return
      end if    
      
      end function radheat
