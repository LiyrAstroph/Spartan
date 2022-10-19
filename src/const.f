c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                  supported by Prof. Wang Jian-Min
c                           and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-03-
c=======================================================================

c To define constants.

      real*8 h,k,c,mp,me,G,msun,sigmat,pi,e,kappaes
      parameter(h=6.6262d-27,k=1.3807d-16,c=2.9979d10,mp=1.6726d-24
     &,         G=6.672d-8,msun=1.989d33,sigmat=6.652d-25
     &,         me=9.1095d-28,pi=3.14159265d0,e=4.8032d-10
     &,         kappaes=0.34d0)
      real*8 mpc2,mec2
      parameter(mpc2=mp*c**2,mec2=me*c**2)
    
      real*8 mui,mue
      parameter(mui=1.23d0,mue=1.14d0)
      
      real*8 eps
      parameter(eps=1.0d-100) 

c------------------------------------------------------------------------
c 1. check up on Tuesday March 11 2008, OK
c
      
