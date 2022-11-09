c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                  supported by Prof. Wang Jian-Min
c                           and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-03-
c=======================================================================

c To set the initial conditions and the boundary conditions.


c***********************************
c subroutine initial()
c***********************************      
      subroutine initial()
      implicit none
      include 'const.f'
      include 'diskvar.f'
      integer nfin,i
      real*8 nulog_min,nulog_max,rdlog_min,rdlog_max
     &,      drdlog
      
      nfin=10
      open(unit=nfin,file='./data/datain.txt',status='unknown')
      read(nfin,adaf)
      read(nfin,odeint)
      read(nfin,auto)
      close(nfin)
      rdin=(1.0d0+dsqrt(1.0-astar*astar))/2.0d0+0.01d0
c      write(*,*)rdin,astar
c      pause
c  note that consider the raditive efficiency for medd.       
      medd=1.39d18*m
      rs=2.95d5*m
      
      nulog_min=dlog10(numin)
      nulog_max=dlog10(numax)
      
      call gauleg(nulog_min,nulog_max,xleg,wleg,nnu)
      
      do i=1,nnu
      nu(i)=10.0d0**(xleg(i))
      end do            
      end subroutine initial 
      
c***********************************
c subroutine boundary()
c
c***********************************     
      subroutine boundary()
      implicit none
      include 'const.f'
      include 'diskvar.f'
      real*8 Tvir,Denst,lout,omegas,aspin,ametric,delta,temp1,temp2
      real*8 omegak,csp
      integer irdout
 
      irdout=1
           
      rd(irdout)=rdout
      
      mud(irdout)=1.0d0
      gamdr(irdout)=1.0d0
      gamdp(irdout)=1.0d0
           
      aspin=astar/2.0d0
      ametric=1.0d0 + aspin**2.0d0/rdout**2.0d0 
     &              + aspin**2.0d0/rdout**3.0d0
      delta=  1.0d0 - 1.0d0/rdout + aspin**2.0d0/rdout**2.0d0
      
c angular velocity
      omegas= aspin /rdout**3.0d0/ametric
      Omegak=dsqrt(0.5d0/rdout)/rdout    
      Omeg(irdout)=0.8*Omegak
      lout= rdout**2.0d0 * ametric**1.5d0 * (Omeg(irdout) - omegas)
     &                                    / delta**0.5d0
      angum(irdout)=lout
c      write(*,*)log10(omeg(irdout)),log10(lout)

c temperature      
      Tvir=(gammaa-1.0d0)*mp/k/rdout*c**2/2.0d0
      Te(irdout)=0.1*Tvir
      Ti(irdout)=0.1*Tvir
c      write(*,*)"T=",dlog10(0.1*Tvir)


c surface density
      temp1= mdot * (lout-lin) / (2.0*pi) *beta
     &      /(rdout**2.0d0 * alpha * ametric**1.5d0 *delta**0.5d0 ) 
      temp2=k/mp*(Ti(irdout)/mui+Te(irdout)/mue) 
      Denst=temp1/temp2*c**2
c      write(*,*)"Sigma=",Denst*medd/rs/c
      Sig(irdout)=Denst

c vertically integrated pressure
      Wi(irdout)=k*Ti(irdout)/mui/mp*Denst/c**2
      We(irdout)=k*Te(irdout)/mue/mp*Denst/c**2 
c      write(*,*)"Wi=",Wi(irdout),",  We=",We(irdout)

c sound speed, consider the magnetic pressure.
      Csp=dsqrt((Wi(irdout)+We(irdout))/Denst/beta)  
          
c      heig(irdout)=csp/omegak
      heig(irdout)=csp/Omeg(irdout)
      end subroutine boundary           
