c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                  supported by Prof. Wang Jian-Min
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-01-
c======================================================================= 

c To calculate the comptonization.
c
c all the formula can be found in paper: 
c     Coppi & Blandford 1990 MNRAS 245,453 (CB)
c

c\\\\\\\\\\\\\\\\\\\\\\\\\\
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\           
      subroutine rate_numb()
      implicit none
      include 'const.f'      
      include 'specvar.f'
            
      real*8 ratecal
      real*8 omig_mean,omig2_mean
      real*8 omi,omi2
      
      real*8 omig_in,gam,temp
      integer iomig,igam 
      
      write(*,*)"call rate_numb(nrd)"     
      do igam=1,nleg_gamm
      write(*,*)"igam=",igam
      gam=10.0d0**(xleg_gamm(igam))
          
      do iomig=1,nleg_omig
      omig_in=10.0**(xleg_omig(iomig))
      rate(iomig,igam)=ratecal(omig_in,gam)
      omi=omig_mean(omig_in,gam)
      omi2=omig2_mean(omig_in,gam)
      temp=dsqrt(3.0*(omi2-omi**2)) 
      omigm(iomig,igam)=omi
      domi(iomig,igam)=min(temp,omi)    
      end do
      end do
                 
      end subroutine rate_numb
c\\\\\\\\\\\\\\\\\\\\\\\\\\
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\      
      subroutine compn(nrd)
      implicit none
      include 'const.f'      
      include 'specvar.f'
      integer nrd
      real*8 distribution,elect_dist
      real*8 omig_in,omig_out,gam,theta,temp_gamm,temp_omig
     &,      taoes,numb_dens_omig_in(nleg_omig)
     
      integer inu,iomig,igam
      
      theta=k*te(nrd)/mec2
      taoes=heig(nrd)*rho(nrd)/mue*sigmat/mp
      
      do iomig=1,nleg_omig
      numb_dens_omig_in(iomig)=fnu_unred(iomig)/(h*nu_unred(iomig))
      end do
           
      do 50 inu=1,nnu
      omig_out=nu_unred(inu)*h/mec2
      temp_gamm=0.0d0
      
      do 40 igam=nleg_gamm,1,-1
      gam=10.0d0**(xleg_gamm(igam))
      numb_elec(igam)=elect_dist(gam,theta) 
c      if((gam-1.0).lt.omig_out)goto 40  
          
      temp_omig=0.0d0
      do 30 iomig=1,nleg_omig
      
      omig_in=10.0**(xleg_omig(iomig))
      
      if(omig_in.gt.omig_out)goto 30
      
      temp_omig=temp_omig+wleg_omig(iomig)
     &          *distribution(omig_out,iomig,igam)
     &          *rate(iomig,igam)**2
     &          *omig_in*numb_dens_omig_in(iomig)
30    continue
      temp_gamm=temp_gamm+wleg_gamm(igam)*numb_elec(igam)
     &                    *gam*taoes*temp_omig
40    continue

      flux_comp(inu)=temp_gamm*dlog(10.d0)**2*h*nu_unred(inu)
c     &                        *(1.0+taoes)
50    continue
      
      return
      end subroutine compn
c-------------------------------

c above is to integrate the kinetic equation
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c To cal the electron population
      function elect_dist(gam,theta)
      implicit none
      include 'const.f'
      real*8 elect_dist,gam,theta
      real*8 beta_gam,bessk
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      if(theta.lt.1.0d-2)then
      elect_dist=gam**2*beta_gam*dexp(-(gam-1.0d0)/theta) 
     &     /dsqrt(PI/2.0d0)/theta**1.5d0    
      else
      elect_dist=gam**2*beta_gam*dexp(-gam/theta)
     &          /theta/bessk(2,1.0/theta)
      end if
      return
      end function elect_dist
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c To cal the scatted photon distribution
c see equation 2.7 in CB
      function distribution(omig_out,iomig,igam)
      implicit none    
      include 'specvar.f'
      real*8 distribution,omig_out
      integer iomig,igam
      real*8 delta_omig,omi,disp
      
      omi=omigm(iomig,igam)
      disp=domi(iomig,igam)
c      write(*,*)disp,dabs(omig_out-omi)      
      if(disp.gt.dabs(omig_out-omi))then
      distribution=0.5d0/disp
      else
      distribution=0.0d0
      end if     
      end function distribution
c above is to cal the scattederd photon
c distribution
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
c################################################
c To cal the second momentum of the mean scattered
c photon energy.
c///////////////////////
c function omig2_mean 
c see equation 2.16 in CB 
c//////////////////////    
      function omig2_mean(omiga,gam)
      implicit none    
      include 'specvar.f'
      real*8 omig2_mean,omiga,gam

      real*8 beta_gam,omig2_alpha,ratecal
      integer i
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      omig2_mean=0.0d0
      do i=1,nleg/2
      omig2_mean=omig2_mean
     &+wleg(i)*omig2_alpha(omiga,gam,xleg(i))*(1.0-beta_gam*xleg(i))/2.0
     &+wleg(i)*omig2_alpha(omiga,gam,xleg(nleg-i+1))
     & *(1.0-beta_gam*xleg(nleg+1-i))/2.0
      end do
c      omig2_mean=omig2_mean/ratecal(omiga,gam)
      return
      end function omig2_mean
c///////////////////////
c function omig2_alpha 
c see equation 2.16 in CB 
c//////////////////////   
      function omig2_alpha(omiga,gam,mu)
      implicit none    
      include 'specvar.f'
      real*8 omiga,gam,omig2_alpha,mu
      real*8 beta_gam
      real*8 fun_omig2
      
      integer i

      beta_gam=dsqrt(1.0-1.0/gam/gam)
      omig2_alpha=0.0d0
      do i=1,nleg/2
      omig2_alpha=omig2_alpha
     & +wleg(i)*fun_omig2(xleg(i),omiga,gam,mu)
     & +wleg(i)*fun_omig2(xleg(nleg-i+1),omiga,gam,mu)
      end do
      
      omig2_alpha=omig2_alpha*3.0d0/8.0
      return      
      end function omig2_alpha
c///////////////////////
c function fun_omig2
c see equation 2.16 in CB 
c//////////////////////                         
      function fun_omig2(t,omiga,gam,mu)
      implicit none
      real*8 fun_omig2,t,omiga,gam,mu
      real*8 xx,xr,beta_gam,omig2averg
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      xx=gam*omiga*(1.0-beta_gam*mu)
      xr=1.0+(1.0-t)*xx
      
      omig2averg=gam**2*omiga**2*
     & (gam**2*(1.0-beta_gam*mu+beta_gam*t*(mu-beta_gam))**2
     & +0.5*beta_gam**2*(1-t**2)*(1-mu**2))/xr**2
      fun_omig2=omig2averg/xr**2*(xr+1.0/xr-1.0+t**2)
      return      
      end function fun_omig2
c above is to cal the second moment of the mean 
c scattered
c################################################

c**********************************************
c To cal the mean scattered photon energy.
c///////////////////////
c function omig_mean 
c see equation 2.9 in CB 
c//////////////////////    
      function omig_mean(omiga,gam)
      implicit none     
      include 'specvar.f'
      real*8 omig_mean,omiga,gam

      real*8 beta_gam,omig_alpha,ratecal
      integer i
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      omig_mean=0.0d0
      do i=1,nleg/2
      omig_mean=omig_mean
     & +wleg(i)*omig_alpha(omiga,gam,xleg(i))*(1.0-beta_gam*xleg(i))/2.0
     & +wleg(i)*omig_alpha(omiga,gam,xleg(nleg-i+1))
     & *(1.0-beta_gam*xleg(nleg+1-i))/2.0
      end do
c      omig_mean=omig_mean/ratecal(omiga,gam)
      return
      end function omig_mean
c///////////////////////
c function omig_alpha 
c see equation 2.9 in CB 
c//////////////////////   
      function omig_alpha(omiga,gam,mu)
      implicit none     
      include 'specvar.f'
      real*8 omiga,gam,omig_alpha,mu
      real*8 beta_gam
      real*8 fun_omig
      
      integer i

      beta_gam=dsqrt(1.0-1.0/gam/gam)
      omig_alpha=0.0d0
      do i=1,nleg/2
      omig_alpha=omig_alpha
     & +wleg(i)*fun_omig(xleg(i),omiga,gam,mu)
     & +wleg(i)*fun_omig(xleg(nleg-i+1),omiga,gam,mu)
      end do
      
      omig_alpha=omig_alpha*3.0d0/8.0
      return      
      end function omig_alpha
c///////////////////////
c function fun_omig 
c see equation 2.10 in CB
c//////////////////////                         
      function fun_omig(t,omiga,gam,mu)
      implicit none
      real*8 fun_omig,t,omiga,gam,mu
      real*8 xx,xr,beta_gam,omigaverg
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      xx=gam*omiga*(1.0-beta_gam*mu)
      xr=1.0+(1.0-t)*xx
      omigaverg=gam**2*omiga*(1.0-beta_gam*mu+beta_gam*t
     &          *(mu-beta_gam))/xr
      fun_omig=omigaverg/xr**2*(xr+1.0/xr-1.0+t**2)
      return      
      end function fun_omig
c above is to cal the mean scattered photon energy.
c*********************************************************
      
c==========================================
c To cal the angle-averaged scattering rate 
c///////////////////////
c function ratecal 
c see equation (2.3) in CB
c the unit is c*sigmat
c//////////////////////       
      function ratecal(omiga,gam)
      implicit none   
      include 'specvar.f'
      real*8 ratecal,omiga,gam
c      integer nleg
c      parameter(nleg=400)
c      real*8 xleg(nleg),wleg(nleg)
c      common/intleg/xleg,wleg,nleg
      real*8 beta_gam,section
      integer i
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      ratecal=0.0d0
      do i=1,nleg/2
      ratecal=ratecal
     & +wleg(i)*section(omiga,gam,xleg(i))*(1.0-beta_gam*xleg(i))/2.0
     & +wleg(i)*section(omiga,gam,xleg(nleg-i+1))
     & *(1.0-beta_gam*xleg(nleg+1-i))/2.0
      end do
      return
      end function ratecal
c///////////////////////
c function section  
c see equation (2.9) in CB
c//////////////////////   
      function section(omiga,gam,mu)
      implicit none      
      include 'specvar.f'
      real*8 omiga,gam,section,mu
      real*8 beta_gam
      real*8 x,funsec
c      integer nleg
c      parameter(nleg=400)
c      real*8 xleg(nleg),wleg(nleg)
c      common/intleg/xleg,wleg,nleg
      
      integer i,n
      beta_gam=dsqrt(1.0-1.0/gam/gam)
      x=omiga*gam*(1.0-beta_gam*mu)

      section=0.0d0
      do i=1,nleg/2
      section=section
     & +wleg(i)*funsec(xleg(i),x)
     & +wleg(i)*funsec(xleg(nleg-i+1),x)
      end do
      section=section*3.0d0/8.0
      return      
      end function section
c///////////////////////
c function funsec 
c see equation (2.11) in CB 
c//////////////////////                         
      function funsec(t,xx)
      implicit none
      real*8 funsec,t,xx
      real*8 xr
      xr=1.0+(1.0-t)*xx
      funsec=1.0/xr**2*(xr+1.0/xr-1.0+t**2)
      return      
      end function funsec
c above is to cal the angle-averaged scattering rate
c====================================================
