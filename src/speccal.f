c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                 supported by Prof. Wang Jian-Min
c                          and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-01-
c======================================================================= 

c To calculate the spectrum.f

c***********************************
c subroutine spectrum
c***********************************     
      subroutine spectrum()
      implicit none
      include 'const.f'      
      include 'specvar.f'
      real*8 fluxr(nnu)
      integer ird,inu
      
      
      open(unit=30,file="../data/spectrum.dat",form="formatted")
      write(*,*)nu_min,nu_max

      omig_in_min=dlog10(h*nu_min/mec2)
      omig_in_max=dlog10(h*nu_max/mec2)
      call gauleg(omig_in_min,omig_in_max
     &,           xleg_omig,wleg_omig,nleg_omig)

      do inu=1,nnu
      nu_unred(inu)=nu(inu)
      end do 
c cal the comptonization
      call rate_numb()
      
      do 20 ird=nr,1,-1
      write(*,*)ird,rd(ird)
    
      call flux_nu(ird)  
      
      do inu=1,nnu
      fluxr(inu)=fnu_unred(inu)
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*fnu_unred(inu)*darea(ird)
     &                     *2.0
      end do
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c--------------------------------
c first comptonization
c--------------------------------
c source photons distribution: fnu_unred
c output photon distribution:flux_comp
      write(*,*)"first comptonization"
      call compn(ird)
      do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)
     &                     *2.0
      end do
c--------------------------------
c second comptonization
c--------------------------------
c source photons distribution: fnu_unred
c output photon distribution:flux_comp
      write(*,*)"second comptonization"
      do inu=1,nnu
      fnu_unred(inu)=flux_comp(inu)
      enddo
      call compn(ird)
      do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)      
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)
     &                     *2.0
      end do
c--------------------------------
c third comptonization
c--------------------------------
c source photons distribution: fnu_unred
c output photon distribution:flux_comp
      write(*,*)"third comptonization"
      do inu=1,nnu
      fnu_unred(inu)=flux_comp(inu)
      enddo
      call compn(ird)
      do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)      
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)
     &                     *2.0
      end do
c--------------------------------
c fourth comptonization
c--------------------------------
c source photons distribution: fnu_unred
c output photon distribution:flux_comp
      write(*,*)"fourth comptonization"
      do inu=1,nnu
      fnu_unred(inu)=flux_comp(inu)
      enddo
      call compn(ird)
      do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)      
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)
     &                     *2.0
      end do
c--------------------------------
c fifth comptonization
c--------------------------------
c source photons distribution: fnu_unred
c output photon distribution:flux_comp
      write(*,*)"fifth comptonization"
      do inu=1,nnu
      fnu_unred(inu)=flux_comp(inu)
      enddo
      call compn(ird)
      do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)      
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)
     &                     *2.0
      end do

c--------------------------------
c sixth comptonization
c--------------------------------
c source photons distribution: fnu_unred
c output photon distribution:flux_comp
      write(*,*)"sixth comptonization"
      do inu=1,nnu
      fnu_unred(inu)=flux_comp(inu)
      enddo
      call compn(ird)
      do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)      
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)
     &                     *2.0
      end do       
c--------------------------------
c seventh comptonization
c--------------------------------
c source photons distribution: fnu_unred
c output photon distribution:flux_comp
      write(*,*)"seventh comptonization"
      do inu=1,nnu
      fnu_unred(inu)=flux_comp(inu)
      enddo
      call compn(ird)
      do inu=1,nnu
      fluxr(inu)=fluxr(inu)+flux_comp(inu)      
      nulnu(inu)=nulnu(inu)+nu_unred(inu)*flux_comp(inu)*darea(ird)
     &                     *2.0
      end do                             

      open(unit=85,file="../data/spec/spec"//char(ird/100+48)//
     &char(mod(ird,100)/10+48)//char(mod(ird,10)+48)//".txt",
     &form="formatted")
      do inu=1,nnu
      write(85,"(E15.7,a,E15.7)")dlog10(nu(inu)),"	"
     &,                          dlog10(fluxr(inu)/2.0/pi+1.0d-150)

c note the factor 2.0*PI we divided.
      enddo
      close(85)    
           
20    continue

      do inu=1,nnu
      write(30,'(E15.7,a,E15.7)')dlog10(nu(inu)),"     "
     &                           ,dlog10(nulnu(inu)+1.0d-150)

      end do
      close(30)
      
      end subroutine spectrum
      
c***********************************
c subroutine flux_nu
c***********************************    
      subroutine flux_nu(nrd)
      implicit none
      include 'const.f'      
      include 'specvar.f'
      integer nrd

      real*8 bessk

      real*8 thetae,ftheta,ne
      real*8 qbr,bnu,gaunt,chibr,chisyn,chinu,g1,g2
      
      real*8 velo,B,xm,temp_nu,temp_bnu,temp_gaun,temp_syn
     &,      temp_sync,fxm,kappa,taunu
      integer i
      
      do i=1,nnu
      fnu_unred(i)=1.0d-150
      enddo
c cal bremsstrahlung
      thetae=k*te(nrd)/mec2
      ne=rho(nrd)/mp/mue
      
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
      temp_gaun=h/k/te(nrd)
      B=dsqrt(8.0*pi*(1.0-beta)*rho(nrd))*cs(nrd)
c      B=dsqrt(8.0*pi*(1.0-beta)/beta*rho(nrd))*cs(nrd)
      temp_syn=4.0*pi*me*c/3.0/e/B/thetae**2
      temp_sync=4.43d-30*4.0*pi*ne/bessk(2,1.0/thetae)
      
      do 10 i=1,nnu
      
c black body radiation
      temp_nu=h*nu_unred(i)/k/te(nrd)

      if(temp_nu.lt.100.0d0)then
      bnu=temp_bnu*nu_unred(i)**3.0/(dexp(temp_nu)-1.0d0)
      else
      bnu=temp_bnu*nu_unred(i)**3.0/(dexp(temp_nu/3.0)**3.0d0-1.0d0)
      end if
      
c cal the gaunt factor
      g1=temp_gaun*dsqrt(3.0d0)/pi*dlog(4.0d0/1.781d0/temp_nu)
      g2=temp_gaun*dsqrt(3.0/pi/temp_nu)

      if(temp_nu.lt.0.359d0)then
      gaunt=temp_gaun*dsqrt(3.0d0)/pi*dlog(4.0d0/1.8182d0/temp_nu)
      elseif(temp_nu.lt.0.955d0)then
      gaunt=temp_gaun*1.0d0
      else
      gaunt=temp_gaun*dsqrt(3.0/pi/temp_nu)
      endif

c      if(temp_nu.lt.1.0d0)then
c      gaunt=temp_gaun*dsqrt(3.0d0)/pi*dlog(4.0d0/0.67966d0/temp_nu)
c      else
c      gaunt=temp_gaun*dsqrt(3.0/pi/temp_nu)
c      endif
      
      chibr=qbr*gaunt*dexp(-temp_nu)
      
c cal the synchrotron radiation
      xm=temp_syn*nu_unred(i)
      fxm=4.0505/xm**(1.0/6.0d0)
     &*(1.0+0.4/xm**(0.25d0)+0.5316/xm**(0.5d0))
     &*dexp(-1.8899d0*xm**(1.0d0/3.0))
      chisyn=temp_sync*nu_unred(i)*fxm
      
c cal the flux.      
      chinu=chibr+chisyn
      
      if(bnu.gt.1.0d-150)then
      kappa=chinu/4.0/pi/bnu
      taunu=dsqrt(pi)/2.0*kappa*heig(nrd)
      if(taunu.gt.1.0d-6)then
      fnu_unred(i)=2.0*pi/dsqrt(3.0d0)*bnu*(1.0
     &             -dexp(-2.0*dsqrt(3.0d0)*taunu))
      else
      fnu_unred(i)=2.0*pi/dsqrt(3.0d0)*bnu*(2.0*dsqrt(3.0d0)*taunu
     &            -6.0*taunu**2+4.0*dsqrt(3.0d0)*taunu**3)
      end if 
      else
      fnu_unred(i)=dsqrt(pi)*heig(nrd)*chinu      
      end if 
10    continue
      return               
      end subroutine flux_nu
      
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c modification:
c Fri, Aug 29,2008
c > ne=rho(nrd)/mp/mue
c < ne=rho(nrd)/mp     
