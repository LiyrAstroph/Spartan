c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                  supported by Prof. Wang Jian-Min
c                           and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-03-
c=======================================================================

      subroutine datagen
      implicit none     
      include 'const.f'
      include 'diskvar.f'

      real*8 csp,amach,rmin,rmax
      integer ni,i,nend
      namelist/nrows/nend
                  
      open(unit=12,file='./data/adaf.dat',status='unknown')
      open(unit=14,file='./data/nrows.txt',status='unknown')
      open(unit=18,file='./data/sol_for_spec.dat',status='unknown')
      open(unit=20,file='./data/rdisk.dat',status='unknown')

      rmax=1.0d3
      rmin=(1.0d0+dsqrt(1.0d0-astar*astar))/2.0d0+0.20d0
      
30    format(E15.7,a,E15.7,a,E15.7,a,E15.7,a,E15.7,a,E15.7)
31    format(E15.7,a,E15.7,a,E15.7,a,E15.7,a,E15.7)      
      write(*,*)"nstep=",nstep
      
      nend=0
      ni=nstep/200
      
      do i=1,nstep,ni
      
      if((rd(i).ge.rmin).and.(rd(i).le.rmax))then
      
      csp=dsqrt((Wi(i)+We(i))/Sig(i)/beta) 
      amach=vel(i)/csp
      
      write(12,30)dlog10(rd(i)),"     ",amach,"      ",dlog10(te(i))
     &,"    ",dlog10(ti(i)),"   ",dlog10(sig(i)*medd/rs/c)
     &,"  ",dlog10(heig(i))
     
      write(18,30)rd(i),"	",gamdr(i),"	"
     &,gamdp(i),"	",omeg(i),"	",vel(i)
     
      write(20,31)rd(i)*2.0d0,"   ",gamdr(i),"	",gamdp(i),"	"
     &,0.5*omeg(i),"	",vel(i)
     
      nend=nend+1
      end if      
      enddo  

c output the last step 
c      if((i-ni).lt.nstep)then
c      i=nstep
      
c      csp=dsqrt((Wi(i)+We(i))/Sig(i)/beta) 
c      amach=vel(i)/csp
      
c      write(12,30)dlog10(rd(i)),"     ",amach,"      ",dlog10(te(i))
c     &,"    ",dlog10(ti(i)),"   ",dlog10(sig(i)*medd/rs/c)
c     &,"  ",dlog10(heig(i))
     
c      write(18,30)rd(i),"	",gamdr(i),"	"
c     &,gamdp(i),"	",omeg(i),"	",vel(i)
     
c      write(20,31)rd(i)*2.0d0,"   ",gamdr(i),"	",gamdp(i),"	"
c     &,0.5*omeg(i),"	",vel(i)
     
c      nend=nend+1
c      endif            
      
      write(14,nrows)  
      write(*,*)"nend=",nend  
      close(12)
      close(14)
      close(18)
      close(20)           
      end subroutine datagen
