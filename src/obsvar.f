c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                 supported by Prof. Wang Jian-Min
c                          and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-01-
c======================================================================= 
c variables definition

      integer nobs,nrd
      parameter(nobs=1000)
      real*8 rd(nobs),gamr(nobs),gamp(nobs),omega(nobs),v(nobs)
      real*8 log_rd(nobs),log_gamr(nobs),log_gamp(nobs)
     &,      log_v(nobs)
      common/sol/rd,gamr,gamp,omega,v,nrd
      common/logsol/log_rd,log_gamr,log_gamp,log_v
      

      integer nnu
      parameter(nnu=200)
      real*8 nu(nnu),flux(nobs,nnu),log_nu(nnu),log_flux(nobs,nnu)
      common/flux/nu,flux,log_nu,log_flux
            
      real*8 astar
      common/spec/astar
      
      
