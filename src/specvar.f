c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                 supported by Prof. Wang Jian-Min
c                          and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-01-
c======================================================================= 

c To define the variables for calculation of spectrum.

      integer nnu,nnr
      parameter(nnu=200,nnr=200) 
      real*8 medd,rs
      parameter(medd=1.39d18,rs=2.95d5)
      
      real*8 rdin,rdout
      common/rbound/rdin,rdout
      
      integer nleg_gamm
      parameter(nleg_gamm=200)
      real*8 xleg_gamm(nleg_gamm),wleg_gamm(nleg_gamm)
      common/intleg_gamm/xleg_gamm,wleg_gamm   
      
      integer nleg_omig
      parameter(nleg_omig=nnu)
      real*8 xleg_omig(nleg_omig),wleg_omig(nleg_omig)
      common/intleg_omig/xleg_omig,wleg_omig 
      
      integer nleg
      parameter(nleg=100)
      real*8 xleg(nleg),wleg(nleg)
      common/intleg/xleg,wleg
      
      real*8 nu(nnu),nulnu(nnu),nu_unred(nnu),fnu_unred(nnu)
     &,      flux_comp(nnu),darea(nnr)
      common/lumn/nu,nulnu,nu_unred,fnu_unred,flux_comp,darea
    
      integer nr
      real*8 rd(nnr),mach(nnr),te(nnr),ti(nnr),tao(nnr)
     &, heig(nnr),cs(nnr),rho(nnr)
      common/sol/rd,mach,te,ti,tao,heig,cs,rho,nr
      
      real*8 rate(nleg_omig,nleg_gamm),numb_elec(nleg_gamm)
     &,domi(nleg_omig,nleg_gamm),omigm(nleg_omig,nleg_gamm)
      common/rate/rate,numb_elec,domi,omigm
      
      real*8 m,mdot,alpha,beta,gammaa,astar,visf
      common/adaf/m,mdot,alpha,beta,gammaa,astar,visf
      namelist/adaf/m,mdot,alpha,beta,gammaa,astar,visf
      
      real*8 omig_in_min,omig_in_max,nu_min,nu_max
      real*8 gam_min,gam_max
      common/omig/omig_in_min,omig_in_max,nu_min,nu_max
      common/gam/gam_min,gam_max            
