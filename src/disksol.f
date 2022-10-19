c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                  supported by Prof. Wang Jian-Min
c                           and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-03-
c=======================================================================

c To solve the disk equation using RK method.
c

c**********************************************
c subroutine autocal()
c**********************************************
      include 'diskflux.f'

      subroutine autocal()
      implicit none
      include 'diskvar.f'
      
      
      real*8 linold,Numold,Denold,rsi,numsi
      integer i,iter
      logical flag
c      flag=.true.
c      lin=1.825
c      call solve(flag)
c      return 
           
      flag=.false.
      iter=0
      
      if((dabs(linmax-linmin).lt.1.0d-4).or.(linmax.lt.linmin))then
      write(*,*)"Improper linmax and linmin!"
      stop
      end if
      
      lin=linmax
      iter=iter+1      
      write(*,'(i2,f8.4,f8.4,f8.4)')iter,linmin,lin,linmax
      call solve(flag)
      write(*,*)"---------------------------------------"
           
      if(.not.is_sonic)then
      write(*,*)"Too small lin_max"    
      endif
      
      linold=lin
      rsi=rson(2)-(0.0d0-dens(2))/(dens(1)-dens(2))*(rson(2)-rson(1))
      Numold=nums(2)+(rson(2)-rsi)/(rson(2)-rson(1))*(nums(1)-nums(2))    
      
      lin=linmin
      iter=iter+1
      write(*,'(i2,f8.4,f8.4,f8.4)')iter,linmin,lin,linmax
      call solve(flag)
      write(*,*)"--------------------------------------"
      
      if(is_sonic)then
      write(*,*)'Too large lin_min!'
      end if
      
120   lin=0.5d0*(linmax+linmin)
      iter=iter+1
      write(*,'(i2,f8.4,f8.4,f8.4)')iter,linmin,lin,linmax
      call solve(flag)

            
      if(is_sonic)then
      linmax=lin
      linold=lin
      rsi=rson(2)-(0.0d0-dens(2))/(dens(1)-dens(2))*(rson(2)-rson(1))
      Numold= nums(2)+(rson(2)-rsi)/(rson(2)-rson(1))*(nums(1)-nums(2)) 
      write(*,'(f8.2,f8.2,f8.2)')rsi,numold              
      else
      linmin=lin
      endif
      write(*,*)"--------------------------------------"      
      if(dabs(linmax-linmin).gt.0.1)goto 120
       
130   lin=0.5d0*(linmax+linmin)
      iter=iter+1
      write(*,'(i2,f8.4,f8.4,f8.4)')iter,linmin,lin,linmax
      call solve(flag) 
      rsi=rson(2)-(0.0d0-dens(2))/(dens(1)-dens(2))*(rson(2)-rson(1))
      Numsi=nums(2)+(rson(2)-rsi)/(rson(2)-rson(1))*(nums(1)-nums(2))                      
      if(.not.is_sonic)then
      linmin=lin
      else
       write(*,'(f8.2,f8.2,f8.2)')rsi,numsi       
       if((dabs(numsi).lt.dabs(Numold)).and.((numsi.gt.0.0d0) 
     &  .and.(nums(2)*dens(2).lt.0.0)))then  !  
       linmax=lin
       linold=lin
       Numold=numsi    
       else
       linmin=lin
       endif  
      endif
      write(*,*)"--------------------------------------"       
      if(dabs(linmax-linmin).gt.1.0d-4)goto  130
      
      flag=.true.
c      lin=lin+0.001d0
      call solve(flag)
           
      write(*,'(a,f10.6)')"Preferred Lin=",lin
      end subroutine autocal

c*****************************
c subroutine solve()
c*****************************
      subroutine solve(flag)
      implicit none
      include 'const.f'
      include 'diskvar.f'
            
      logical flag
      real*8 df(3),df0(3)
      real*8 macc,rdp,wip,wep,sigp
c for rk       
      real*8 k1(3),k2(3),k3(3),k4(3),k5(3),k6(3),x
     &,      y1,y2,y3,x0,y10,y20,y30,yerr1,yerr2,yerr3
     &,      errmax,htemp
     &,      a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43
     &,      b51,b52,b53,b54,b61,b62,b63,b64,b65
     &,      c1,c2,c3,c4,c5,c6
     &,      dc1,dc2,dc3,dc4,dc5,dc6
      parameter(a2=0.2d0,a3=0.3d0,a4=0.6d0,a5=1.0d0
     &,         a6=0.875d0
     &,         b21=0.2d0,b31=0.075d0,b32=0.225d0
     &,         b41=0.3d0,b42=-0.9d0, b43=1.2d0
     &,         b51=-11.0d0/54.0d0,b52=2.5d0,b53=-70.0d0/27.0d0
     &,         b54=35.0d0/27.0d0,b61=1631.0d0/55296.0d0
     &,         b62=175.0d0/512.0d0,b63=575.0d0/13824.0d0
     &,         b64=44275.0d0/110592.0d0,b65=253.0d0/4096.0d0
     &,         c1=37.0d0/378.0d0,c2=0.0d0,c3=250.0d0/621.0d0
     &,         c4=125.0d0/594.0d0,c5=0.0d0,c6=512.0d0/1771.0d0
     &,         dc1=c1-2825.0d0/27648.0d0,dc2=0.0d0
     &,         dc3=c3-18575.0d0/48384.0d0,dc4=c4-13525.0d0/55296.0d0
     &,         dc5=-277.0d0/14336.0d0,dc6=c6-0.25d0)     
      real*8 safety,pgrow,pshrnk,errcon,eacc
      parameter(safety=0.9d0,pgrow=-0.20d0,pshrnk=-0.25d0,errcon=1.89d-4
     &,         eacc=1.0d-4)
      integer ird,irdi,irdo,ird0,i,nl
      integer nout,nerg,nftot
      real*8 mach,nep,nip,Lamb,Qrad,heigp,tip,tep,rhop,hstep
     &,      factor,v2cs,num,den
      real*8 csp,radcool
      open(unit=60,file='../data/soltot.dat',status='unknown')                  
      irdi=1
      irdo=ndim
      
      is_sonic=.false.   
      macc=mdot
       
      hstep=-drstep
      
      rdp=rdout
      call output(irdi,flag)  
      irdi=irdi+1                
      do 30 ird=irdi,irdo,1

      if(rdp.le.rdin)goto 60
c      if(mod(ird,100).eq.0)then
c      write(*,'(i8,a,f8.3,a,f10.6)')ird,', rdp=',rdp,', hstep=',hstep
c      endif
      ird0=ird-1
c rk step 1
40    sigp=Sig(ird0)
      wip=Wi(ird0)
      wep=We(ird0)
      
      x0=dlog(rdp)      
      y10=dlog(sigp)
      y20=dlog(wip)
      y30=dlog(wep)      
       
      call deriva(df,rdp,sigp,wip,wep,ird0) 
c      write(*,*)"Here!"      
      do i=1,3
      k1(i)=hstep*df(i)
      end do  
c rk step 2        
      x=x0+a2*hstep
      y1=y10+b21*k1(1)
      y2=y20+b21*k1(2)
      y3=y30+b21*k1(3)
      
      rdp=dexp(x)
      sigp=dexp(y1)
      wip=dexp(y2)
      wep=dexp(y3)
      call deriva(df,rdp,sigp,wip,wep,ird0) 
      do i=1,3
      k2(i)=hstep*df(i)
      end do      
c rk step 3
      x=x0+a3*hstep
      y1=y10+b31*k1(1)+b32*k2(1)
      y2=y20+b31*k1(2)+b32*k2(2)
      y3=y30+b31*k1(3)+b32*k2(3)
      
      rdp=dexp(x)
      sigp=dexp(y1)
      wip=dexp(y2)
      wep=dexp(y3)
      call deriva(df,rdp,sigp,wip,wep,ird0) 
      do i=1,3
      k3(i)=hstep*df(i)
      end do 
                   
c rk step 4
      x=x0+a4*hstep
      y1=y10+b41*k1(1)+b42*k2(1)+b43*k3(1)
      y2=y20+b41*k1(2)+b42*k2(2)+b43*k3(2)
      y3=y30+b41*k1(3)+b42*k2(3)+b43*k3(3)
      
      rdp=dexp(x)
      sigp=dexp(y1)
      wip=dexp(y2)
      wep=dexp(y3)
      call deriva(df,rdp,sigp,wip,wep,ird0) 
      do i=1,3
      k4(i)=hstep*df(i)
      end do 
c rk step 5

      x=x0+a5*hstep
      y1=y10+b51*k1(1)+b52*k2(1)+b53*k3(1)+b54*k4(1)
      y2=y20+b51*k1(2)+b52*k2(2)+b53*k3(2)+b54*k4(2)
      y3=y30+b51*k1(3)+b52*k2(3)+b53*k3(3)+b54*k4(3)
      
      rdp=dexp(x)
      sigp=dexp(y1)
      wip=dexp(y2)
      wep=dexp(y3)
      call deriva(df,rdp,sigp,wip,wep,ird0) 
      do i=1,3
      k5(i)=hstep*df(i)
      end do 
c rk step 6

      x=x0+a6*hstep
      y1=y10+b61*k1(1)+b62*k2(1)+b63*k3(1)+b64*k4(1)+b65*k5(1)
      y2=y20+b61*k1(2)+b62*k2(2)+b63*k3(2)+b64*k4(2)+b65*k5(2)
      y3=y30+b61*k1(3)+b62*k2(3)+b63*k3(3)+b64*k4(3)+b65*k5(3)
      
      rdp=dexp(x)
      sigp=dexp(y1)
      wip=dexp(y2)
      wep=dexp(y3)
      call deriva(df,rdp,sigp,wip,wep,ird0) 
      do i=1,3
      k6(i)=hstep*df(i)
      end do 

      y1=c1*k1(1)+c3*k3(1)+c4*k4(1)+c6*k6(1)
      y2=c1*k1(2)+c3*k3(2)+c4*k4(2)+c6*k6(2)
      y3=c1*k1(3)+c3*k3(3)+c4*k4(3)+c6*k6(3)
      
      df0(1)=y1/hstep
      df0(2)=y2/hstep
      df0(3)=y3/hstep 
           
c      write(90,*)y1,y2,y3      
      yerr1=dc1*k1(1)+dc3*k3(1)+dc4*k4(1)+dc5*k5(1)+dc6*k6(1)                      
      yerr2=dc1*k1(2)+dc3*k3(2)+dc4*k4(2)+dc5*k5(2)+dc6*k6(2)
      yerr3=dc1*k1(3)+dc3*k3(3)+dc4*k4(3)+dc5*k5(3)+dc6*k6(3)
      
      errmax=0.0d0
      errmax=max(dabs(yerr1),dabs(yerr2))
      errmax=max(errmax,dabs(yerr3))

      errmax=errmax/eacc
      if(errmax.gt.1.0d0)then
      htemp=safety*hstep*(errmax**pshrnk)/10.0
      htemp=sign(max(dabs(htemp),0.1*dabs(hstep)),hstep)
         if(dabs(htemp).lt.drstep_min)then
         write(*,*)'   hstep underflow'
         htemp=-drstep_min
         nstep=ird0
         goto 60
         end if
      hstep=htemp
      goto 40 
      else
         if(errmax.gt.errcon)then
         htemp=safety*hstep*(errmax**pgrow)
         else
         htemp=5.0*hstep
         end if
         if(dabs(htemp).gt.drstep_max)then
         htemp=-drstep_max 
         end if

      rdp=rd(ird0)*dexp(hstep)
      rd(ird)=rdp     
      Sig(ird)=Sig(ird0)*dexp(y1)      
      Wi(ird)=Wi(ird0)*dexp(y2)
      We(ird)=We(ird0)*dexp(y3)
      if(flag)then
      call energycal(rd(ird),sig(ird),wi(ird),we(ird),ird) 
      endif
      call vcal(ird,v2cs,factor)
      v2cs=v2cs/factor
      
      if((v2cs.gt.0.95d0).and.(v2cs.lt.1.05d0))then
      call ndcal(rd(ird0),sig(ird0),wi(ird0),we(ird0),ird0) 
      rson(1)=rd(ird0)
      nums(1)=sonic_num
      dens(1)=sonic_den 
      write(*,'(f8.2,f8.2,f8.2)')rd(ird0),sonic_num,sonic_den     
      call crossonic(ird0,df0,-1.0d-1)
      
      call ndcal(rd(ird),sig(ird),wi(ird),we(ird),ird) 
      rson(2)=rd(ird)
      nums(2)=sonic_num
      dens(2)=sonic_den  
               
      write(*,'(f8.2,f8.2,f8.2)')rd(ird),sonic_num,sonic_den      
      
      rdp=rd(ird0)*dexp(-1.0d-1)
      write(*,*)'   Transonic!'
      is_sonic=.true.
      if(.not.flag)goto 60
      end if
      call output(ird,flag)
      if((dabs(htemp).gt.1.0d-3).and.((v2cs.gt.0.90)
     &.and.(v2cs.lt.1.10)))then
      htemp=-1.0d-3
      endif
      hstep=htemp
      if(rdp*dexp(hstep).lt.rdin)then
      hstep=dlog(rdin/rdp)
      end if
                                                    
      end if 
                
30    continue  
60    close(60)
      end subroutine solve
c*******************
c subroutine deriva
c*******************      
      subroutine deriva(df,rdp,sigp,wip,wep,ird0)
      implicit none
      include 'const.f' 
      include 'diskvar.f' 
     
      real*8 df(3),rdp,sigp,wip,wep
      integer ird0
      
      real*8 radheat,radcool,adiabatic_gam
      real*8 a13,c1,a21,a22,a23,c2,a31,a32,a33,c3
      real*8 Omega,Omegak1,Omegak2,lrdp,macc,csp,Lamb,Qrad
     &,      aspin,ametric,omegas,delta,Btemp,V,gamp,gamp_old,gamr,eng
     &,      mu,ai,ae,gi,ge,Gbi,Gbe,thetae,thetai,dai,dae
      real*8 tep,tip,nep,nip,heigp,rhop
      real*8 dsig,dwi,dwe,bb1,bb2,aa,tempa,tempc
     &,      x1,x2,x3,y1,y2,num,den,aa1,aa2,S
                 
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
      if(lrdp.gt.1.0d2*angum(ird0))then
      lrdp=angum(ird0)
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
      if(heigp/rdp.gt.2.0d0)heigp=1.0d0*rdp 
      heigp=heigp*rs
      rhop=sigp/dsqrt(2.0*pi)/heigp*medd/rs/c
      nep=rhop/mue/mp
      nip=rhop/mui/mp
      
      Lamb=dsqrt(pi)*heigp*radheat(tep,tip,nep,nip)
      Lamb=Lamb*rs**2/medd/c**2
      
      Qrad=radcool(tep,tip,rhop,heigp,csp*c)
      Qrad=Qrad*rs**2/medd/c**2
       
c------------------------------------------------------------------------
      a13= -beta * mu*gamr**2.0d0*Btemp/(1.0d0+Btemp)**2.0d0
      c1=   beta * mu*gamr**2.0d0*Btemp/(1.0d0+Btemp)**2.0d0
     &          * sigp*(rdp-0.5d0)/(rdp**2.0d0 * delta)
     &    -0.5d0*beta * mu*gamp**2.0d0*sigp*ametric/rdp**2.0/delta
     &          * (Omega-Omegak1)*(Omega-Omegak2)/Omegak1/Omegak2 
     
      tempa= - 4.0d0*(pi*alpha*gamp**3.0d0)**2.0d0/(mu*gamr*beta)
     &     *rdp**2.0d0*ametric**2.0d0*delta*(wip+wep)/beta/macc
      
      a22= (1.0d0-visf)*tempa
      aa1= -1.0d0/(Gbi-1.0d0)*macc/sigp
      a21= -aa1 + a22
      a23= - Gbi/(Gbi-1.0d0)*macc*wip/sigp**2.0d0  
c=======================================================================      
      x1= 1.0d0/rdp**2.0d0 - 2.0d0 * aspin**2.0d0 / rdp**3.0d0
      x2= delta**0.5d0 / rdp**2.0d0 / ametric**1.5d0
     & *(  x1/2.0d0/delta - 2.0d0/rdp 
     &   + 1.5d0*(2.0d0 + 3.0d0/rdp) *aspin**2.0d0/rdp**3.0d0/ametric )
      x3= aspin/rdp**4.0d0 /ametric
     & * (-3.0d0 + aspin**2.0d0/ametric/rdp**2.0d0*(2.0d0+3.0d0/rdp))
      
      if(rdp.eq.rd(ird0))then
      y1=0.0d0
      y2=0.0d0
      else
      y1=(gamp**2.0d0/mu/gamr-gamdp(ird0)**2.0d0/mud(ird0)/gamdr(ird0))
     &   /(rdp-rd(ird0))
      y2= ( 1.0d0 / (mu*gamr*gamp)
     &     -1.0d0 / (mud(ird0)*gamdr(ird0)*gamdp(ird0))
     &    ) / (rdp-rd(ird0))
      endif

      tempc=2.0d0*pi*alpha * gamp**4.0d0 * ametric**2.0d0 * rdp**2.0d0
     &     *( wip + wep ) / beta
     &     *(  2.0d0*pi*alpha*(wip+wep)/beta /macc
     &          *(gamp**2.0d0 /mu/gamr*x1 +delta*y1 )
     &       + lin*(1.0d0/mu/gamr/gamp * x2 
     &            + delta**0.5d0/rdp**2.0d0/ametric**1.5d0*y2)
     &       + x3 )
     
c=======================================================================
      c2=(1.0d0-visf)*tempc - macc*wip/sigp/rdp + 2.0d0*pi*rdp*Lamb
      
      a31=visf*tempa
      aa2= 1.0d0/(Gbe-1.0d0)*macc/sigp
      a32=  aa2 + a31
      a33= -Gbe/(Gbe-1.0d0)*macc*wep/sigp**2.0d0
      
      c3= visf*tempc - macc*wep/sigp/rdp +2.0d0*pi*rdp*(Qrad-Lamb)
c-----------------------------------------------------------------------
      num=(c3-a31*c1)*aa1 - (c2-a21*c1)*aa2
      den=(a33-a31*a13)*aa1 - (a23-a21*a13)*aa2
      S=num/den

      dsig=S
      dwi=(a22*c1-c2)/aa1 - (a22*a13-a23)/aa1 * S
      dwe=(c2-a21*c1)/aa1 - (a23-a21*a13)/aa1 * S
     
      dsig=dsig*rdp/sigp    
      dwi=dwi*rdp/Wip
      dwe=dwe*rdp/Wep
      
      df(1)=dsig
      df(2)=dwi
      df(3)=dwe
      end subroutine deriva
      
c***************************
c subroutine crossonic()
c***************************

      subroutine crossonic(ird0,df,hstep)
      implicit none
      include 'const.f'
      include 'diskvar.f'
      integer ird0
      real*8 df(3),hstep
      integer j,ird
      real*8 y1,y2,y3,rdp
      
      y1=df(1)*hstep
      y2=df(2)*hstep
      y3=df(3)*hstep
      
      ird=ird0+1
      
      rdp=rd(ird0)*dexp(hstep)
      rd(ird)=rdp     
      Sig(ird)=Sig(ird0)*dexp(y1)      
      Wi(ird)=Wi(ird0)*dexp(y2)
      We(ird)=We(ird0)*dexp(y3)
      
      end subroutine crossonic

c****************************
c subroutine daie(ti,te)
c****************************
      subroutine daie(tip,tep,dai,dae)
      implicit none
      include 'const.f'
      include 'diskvar.f'
      real*8 tip,tep,dai,dae
      
      real*8 adiabatic_gam
      real*8 thetai,thetae,dti,del
      real*8 gi,ge,ai,ae,aid,aed
      
      thetai=k*tip/mpc2
      thetae=k*tep/mec2
      del=1.0d-3
      
      gi=adiabatic_gam(thetai)
      ge=adiabatic_gam(thetae)
     
      ai=1.0d0/(gi-1.0d0) + 2.0d0*(1.0d0-beta)/beta
      ae=1.0d0/(ge-1.0d0) + 2.0d0*(1.0d0-beta)/beta
      
      thetai= thetai * (1.0d0+del)
      thetae= thetae * (1.0d0+del)
      gi=adiabatic_gam(thetai)
      ge=adiabatic_gam(thetae)     

      aid=1.0d0/(gi-1.0d0) + 2.0d0*(1.0d0-beta)/beta
      aed=1.0d0/(ge-1.0d0) + 2.0d0*(1.0d0-beta)/beta
            
      dai=dlog(aid/ai)/dlog(1.0d0+del)
      dae=dlog(aed/ae)/dlog(1.0d0+del)
      
      end subroutine daie

c******************************
c subroutine adiabatic_gam(theta)
c******************************
      function adiabatic_gam(theta)
      implicit none
      real*8 theta,adiabatic_gam
      real*8 bessk,bessk1
      if(theta.lt.2.0d-3)then
      adiabatic_gam = 1.0d0 + (8.0d0 - 9.95d0*theta) / 12.0d0
      return
      end if

      adiabatic_gam=1.0d0
     &  +theta / ( ( 3.0d0*bessk(3,1.0d0/theta) + bessk1(1.0d0/theta) )
     &                /4.0d0/bessk(2,1.0d0/theta) -1.0d0 )
      return      
      end function adiabatic_gam

c******************************
c subroutine output()
c******************************
      subroutine output(ird,flag)
      implicit none
      include 'const.f'
      include 'diskvar.f'
      
      integer ird
      logical flag
      
      real*8 adiabatic_gam
      real*8 Delta,Aspin,Ametric,Btemp,V,gamr,gamp,gamp_old,heigp
     &,      thetai,thetae,ai,ae,gi,ge,mu,lrd,rdp,macc,omegas,omega
     &,      eng,csp,lkep,aa1,aa2
      integer nftot,nout
      nftot=60
      
      macc=mdot
      rdp=rd(ird)
     
      Aspin=astar/2.0d0 
      Delta=1.0d0 - 1.0d0/rdp + Aspin**2.0d0 / rdp**2.0d0
      Ametric=1.0d0 + Aspin**2.0d0/rdp**2.0d0 
     &              + Aspin**2.0d0/rdp**3.0d0
      omegas = aspin / rdp**3.0d0 / ametric 
            
      Ti(ird)=mui*mp*Wi(ird)/k/Sig(ird)*c**2.0d0
      Te(ird)=mue*mp*We(ird)/k/Sig(ird)*c**2.0d0
c      write(*,*)"Ti=",log10(Ti(ird)),",Te=",log10(Te(ird))
            
      Csp=dsqrt((Wi(ird)+We(ird))/Sig(ird)/beta)      

      Btemp=macc**2.0d0/(4.0d0*pi*pi)/Delta/rdp**2.0d0/Sig(ird)**2.0d0      
      V=dsqrt(Btemp/(1.0+Btemp))
      
c      write(*,*)V,Btemp

      gamr=1.0d0/sqrt(1.0d0-V*V)
      thetae=k*te(ird)/mec2
      thetai=k*ti(ird)/mpc2
      gi=adiabatic_gam(thetai)
      ge=adiabatic_gam(thetae)
      ai=1.0d0/(gi-1.0d0)+2.0d0*(1.0d0-beta)/beta
      ae=1.0d0/(ge-1.0d0)+2.0d0*(1.0d0-beta)/beta
      mu=1.0d0+((ai+1.0d0/beta)*Wi(ird)+(ae+1.0d0/beta)*We(ird))
     &          /Sig(ird) 
C------------------------------------------------------------------------
C cal angular momentum L and gamma_phi by interation
c-----------------------------------------------------------------------
      gamp=1.0d0
80    gamp_old=gamp 
      lrd=lin + 2.0d0*pi*rdp**2.0d0 / macc * alpha * ametric**1.5d0
     &         * delta**0.5d0 * (wi(ird)+we(ird))/beta * gamp**3.0d0
      gamp=dsqrt(1.0d0 + lrd**2.0d0 / (mu*gamr*rdp)**2.0d0 /ametric)
      if(lrd.gt.1.0d2*angum(ird-1).and.(ird.gt.1))then
      lrd=angum(ird-1)
      gamp=dsqrt(1.0d0 + lrd**2.0d0 / (mu*gamr*rdp)**2.0d0 /ametric)
      goto 110
      endif      
      if(dabs(gamp-gamp_old).gt.1.0d-6)goto 80
c------------------------------------------------------------------------ 
110   angum(ird)=lrd
      Omega=omegas + delta**0.5d0*lrd 
     &              /(mu*ametric**1.5d0*rdp**2.0d0*gamr*gamp)

      eng= mu*gamr*gamp*dsqrt(delta/ametric) + omegas*lrd
      aa1=lrd**2.0d0
      aa2=aspin**2.0d0*(eng**2.0d0-mu*mu)
      if(aa1.le.aa2)then
      heigp= csp *rdp**2.0 * dsqrt(mu)
     &      /dsqrt(dabs(aa1-aa2))
      else
      heigp= csp *rdp**2.0 * dsqrt(mu)
     &      /dsqrt(aa1 - aa2)
      end if
      if(heigp/rdp.gt.2.0d0)heigp=1.0d0*rdp       
  
C=======================================================================
      heig(ird)=heigp
      mud(ird)=mu
      gamdr(ird)=gamr
      gamdp(ird)=gamp
      Omeg(ird)=Omega
      vel(ird)=V      
C=======================================================================
      lkep=dsqrt(0.5d0)*(rdp**2.0d0-dsqrt(2.0d0*rdp)*aspin+aspin**2.0d0)
     &/rdp**0.75d0/dsqrt(rdp**1.5d0-1.5d0*rdp**0.5d0+dsqrt(2.0d0)*aspin)

      nstep=ird
C output to file.
      if(flag)then       
      write(nftot,'(11(e15.7,a))')dlog10(rdp),"      ",dlog10(Ti(ird))
     &,"    ",dlog10(Te(ird)),"      ",Sig(ird)*medd/rs/c
     &,"    ",V,"      ",Omega,"   ",mu,"      "
     &,heigp/rdp,"      ",Csp," ",lrd,"   ",lkep
c      write(*,*)rdp**1.5d0-1.5d0*rdp**0.5d0+dsqrt(2.0d0)*aspin
      end if
      end subroutine output
       
c********************************
c     subroutine vcal()
c********************************      
      subroutine vcal(ird,v2cs,factor)
      implicit none
      include 'const.f'
      include 'diskvar.f'
      real*8 v2cs,factor
      integer ird
      
      real*8 adiabatic_gam
      real*8 Delta,Aspin,Ametric,Btemp,V,gamr,gamp,gamp_old,tip,tep
     &,      thetai,thetae,ai,ae,gi,ge,mu,lrd,rdp,macc,Gbi,Gbe
     &,      dai,dae,csp
      real*8 b1,b2,b3
     
      macc=mdot
      rdp=rd(ird)
     
      Aspin=astar/2.0d0 
      Delta=1.0d0 - 1.0d0/rdp + Aspin**2.0d0 / rdp**2.0d0
      Ametric=1.0d0 + Aspin**2.0d0/rdp**2.0d0 
     &              + Aspin**2.0d0/rdp**3.0d0
      
      csp=dsqrt((Wi(ird)+We(ird))/beta/Sig(ird))      
      Tip=mui*mp*Wi(ird)/k/Sig(ird)*c**2.0d0
      Tep=mue*mp*We(ird)/k/Sig(ird)*c**2.0d0     

      Btemp=macc**2.0d0/(4.0d0*pi*pi)/Delta/rdp**2.0d0/Sig(ird)**2.0d0      
      V=dsqrt(Btemp/(1.0+Btemp))

      gamr=1.0d0/sqrt(1.0d0-V*V)
      thetae=k*tep/mec2
      thetai=k*tip/mpc2
      gi=adiabatic_gam(thetai)
      ge=adiabatic_gam(thetae)
      ai=1.0d0/(gi-1.0d0)+2.0d0*(1.0d0-beta)/beta
      ae=1.0d0/(ge-1.0d0)+2.0d0*(1.0d0-beta)/beta
      call daie(tip,tep,dai,dae)
      Gbi= 1.0d0 + 1.0d0 / (ai*(1.0d0 + dai))
      Gbe= 1.0d0 + 1.0d0 / (ae*(1.0d0 + dae))      
      mu=1.0d0+((ai+1.0d0/beta)*Wi(ird)+(ae+1.0d0/beta)*We(ird))
     &          /Sig(ird) 
C------------------------------------------------------------------------
C cal angular momentum L and gamma_phi by interation
c-----------------------------------------------------------------------
      gamp=1.0d0
90    gamp_old=gamp 
      lrd=lin + 2.0d0*pi*rdp**2.0d0 / macc * alpha * ametric**1.5d0
     &         * delta**0.5d0 * (wi(ird)+we(ird))/beta * gamp**3.0d0
      gamp=dsqrt(1.0d0 + lrd**2.0d0 / (mu*gamr*rdp)**2.0d0 /ametric)
      
      if(dabs(gamp-gamp_old).gt.1.0d-6)goto 90
c------------------------------------------------------------------------
      b1= beta*mu*gamr**4.0d0/(1.0d0+Btemp)**2.0d0
      
      b2= beta*(Gbi*Wi(ird)+Gbe*We(ird))/(Wi(ird)+We(ird)) 
     
      b3= (visf*(Gbe-1.0d0)+(1.0d0-visf)*(Gbi-1.0d0)) 
     &   *alpha**2.0d0 *gamr * gamp**6.0d0 * Ametric**2.0d0
     &   /(1.0d0+Btemp)**2.0d0 
     
      v2cs=V/csp
      factor=dsqrt((b2+b3)/b1)
      end subroutine vcal  
      
c*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
c*-                                                                   -*
c*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
      subroutine ndcal(rdp,sigp,wip,wep,ird0)
      implicit none
      include 'const.f' 
      include 'diskvar.f' 
     
      real*8 df(3),rdp,sigp,wip,wep
      integer ird0
      
      real*8 radheat,radcool,adiabatic_gam
      real*8 a13,c1,a21,a22,a23,c2,a31,a32,a33,c3
      real*8 Omega,Omegak1,Omegak2,lrdp,macc,csp,Lamb,Qrad
     &,      aspin,ametric,omegas,delta,Btemp,V,gamp,gamp_old,gamr,eng
     &,      mu,ai,ae,gi,ge,Gbi,Gbe,thetae,thetai,dai,dae
      real*8 tep,tip,nep,nip,heigp,rhop
      real*8 dsig,dwi,dwe,bb1,bb2,aa,tempa,tempc
     &,      x1,x2,x3,y1,y2,num,den,aa1,aa2,S
                 
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
      if(lrdp.gt.1.0d2*angum(ird0))then
      lrdp=angum(ird0)
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
      if(heigp/rdp.gt.2.0d0)heigp=1.0d0*rdp 
      heigp=heigp*rs
      rhop=sigp/dsqrt(2.0*pi)/heigp*medd/rs/c
      nep=rhop/mue/mp
      nip=rhop/mui/mp
      
      Lamb=dsqrt(pi)*heigp*radheat(tep,tip,nep,nip)
      Lamb=Lamb*rs**2/medd/c**2
      
      Qrad=radcool(tep,tip,rhop,heigp,csp*c)
      Qrad=Qrad*rs**2/medd/c**2
       
c------------------------------------------------------------------------
      a13= -beta * mu*gamr**2.0d0*Btemp/(1.0d0+Btemp)**2.0d0
      c1=   beta * mu*gamr**2.0d0*Btemp/(1.0d0+Btemp)**2.0d0
     &          * sigp*(rdp-0.5d0)/(rdp**2.0d0 * delta)
     &    -0.5d0*beta * mu*gamp**2.0d0*sigp*ametric/rdp**2.0/delta
     &          * (Omega-Omegak1)*(Omega-Omegak2)/Omegak1/Omegak2 
     
      tempa= - 4.0d0*(pi*alpha*gamp**3.0d0)**2.0d0/(mu*gamr*beta)
     &     *rdp**2.0d0*ametric**2.0d0*delta*(wip+wep)/beta/macc
      
      a22= (1.0d0-visf)*tempa
      aa1= -1.0d0/(Gbi-1.0d0)*macc/sigp
      a21= -aa1 + a22
      a23= - Gbi/(Gbi-1.0d0)*macc*wip/sigp**2.0d0  
c=======================================================================      
      x1= 1.0d0/rdp**2.0d0 - 2.0d0 * aspin**2.0d0 / rdp**3.0d0
      x2= delta**0.5d0 / rdp**2.0d0 / ametric**1.5d0
     & *(  x1/2.0d0/delta - 2.0d0/rdp 
     &   + 1.5d0*(2.0d0 + 3.0d0/rdp) *aspin**2.0d0/rdp**3.0d0/ametric )
      x3= aspin/rdp**4.0d0 /ametric
     & * (-3.0d0 + aspin**2.0d0/ametric/rdp**2.0d0*(2.0d0+3.0d0/rdp))
      
      if(rdp.eq.rd(ird0))then
      y1=0.0d0
      y2=0.0d0
      else
      y1=(gamp**2.0d0/mu/gamr-gamdp(ird0)**2.0d0/mud(ird0)/gamdr(ird0))
     &   /(rdp-rd(ird0))
      y2= ( 1.0d0 / (mu*gamr*gamp)
     &     -1.0d0 / (mud(ird0)*gamdr(ird0)*gamdp(ird0))
     &    ) / (rdp-rd(ird0))
      endif

      tempc=2.0d0*pi*alpha * gamp**4.0d0 * ametric**2.0d0 * rdp**2.0d0
     &     *( wip + wep ) / beta
     &     *(  2.0d0*pi*alpha*(wip+wep)/beta /macc
     &          *(gamp**2.0d0 /mu/gamr*x1 +delta*y1 )
     &       + lin*(1.0d0/mu/gamr/gamp * x2 
     &            + delta**0.5d0/rdp**2.0d0/ametric**1.5d0*y2)
     &       + x3 )
     
c=======================================================================
      c2=(1.0d0-visf)*tempc - macc*wip/sigp/rdp + 2.0d0*pi*rdp*Lamb
      
      a31=visf*tempa
      aa2= 1.0d0/(Gbe-1.0d0)*macc/sigp
      a32=  aa2 + a31
      a33= -Gbe/(Gbe-1.0d0)*macc*wep/sigp**2.0d0
      
      c3= visf*tempc - macc*wep/sigp/rdp +2.0d0*pi*rdp*(Qrad-Lamb)
c-----------------------------------------------------------------------
      num=(c3-a31*c1)*aa1 - (c2-a21*c1)*aa2
      den=(a33-a31*a13)*aa1 - (a23-a21*a13)*aa2
     
      sonic_num=num*rdp/sigp    
      sonic_den=den
      
      end subroutine ndcal
