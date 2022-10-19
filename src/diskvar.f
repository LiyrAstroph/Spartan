c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                  supported by Prof. Wang Jian-Min
c                           and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-03-
c=======================================================================

c To define the variables used in the code for calculate the global
c structure of ADAF.

      real*8 m,mdot,alpha,beta,gammaa,astar,visf,lin
      common/adaf/m,mdot,alpha,beta,gammaa,astar,visf,lin
      namelist/adaf/m,mdot,alpha,beta,gammaa,astar,visf
      
      real*8 medd,rs
      common/units/medd,rs
      namelist/units/medd,rs
      
      integer ndim
      parameter(ndim=5000)
      integer nnu,nrd
      real*8 xleg(ndim),wleg(ndim),rdin,rdout,numin,numax,drstep
     &,      drstep_min,drstep_max
      common/odeint/nnu,nrd,xleg,wleg,rdin,rdout,numin,numax,drstep
     &,      drstep_min,drstep_max      
      namelist/odeint/nnu,nrd,rdout,numin,numax,drstep
     &,             drstep_min,drstep_max
      
      real*8 flux(ndim),nu(ndim)
      common/flux/flux,nu
      
      real*8 Wi(ndim),We(ndim),Sig(ndim),Te(ndim),Ti(ndim)
     &,      Omeg(ndim),rd(ndim),gamdr(ndim),gamdp(ndim),mud(ndim)
     &,      heig(ndim),angum(ndim),vel(ndim)
      common/equa/Wi,We,Sig,Te,Ti,Omeg,rd,gamdr,gamdp,mud,heig
     &,           angum,vel
     
      logical is_sonic
      integer nstep
      real*8 nums(2),dens(2),rson(2),sonic_num,sonic_den,linmax,linmin
      common/auto/nums,dens,rson,sonic_num,sonic_den,linmax,linmin,nstep
     &,            is_sonic
      namelist/auto/linmax,linmin
      
       
