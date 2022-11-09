c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                 supported by Prof. Wang Jian-Min
c                          and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-01-
c======================================================================= 

c To calculate the spectrum.

c-----------------------------------------------------------------------
c References: 
c 1. Coppi & Blandford, 1990, MNRAS, 245,453-469
c 2. Narayan & Yi, 1995, ApJ, 452,710-735
c 3. Ozel, Psaltis & Narayan, 2000, ApJ, 541, 234-249
c 4. Manmoto, Mineshige & Kusunose, 1997, ApJ, 489, 791-803
c 5. Press, et al, Numerical Recipes in Fortran, Cambrige Univ. Press, 1992
c-----------------------------------------------------------------------

      include 'specfun.f'
      include 'speccomp.f'
      include 'speccal.f'

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      program main
      
      implicit none
      include 'specvar.f'
      
      call system('rm ./data/spec/spec*')
                  
      call initial()
      call spectrum()
      
      write(*,adaf)
         
      end program main
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
c***********************************
c subroutine initial
c***********************************   
      subroutine initial()
      implicit none
      include 'const.f'
      include 'specvar.f'
     
      integer i,j,nend
      
      namelist/nrows/nend 
      
      real*8 temp_cs,velo,zmax,area
      
      nu_min=1.0d9
      nu_max=1.0d22 
      
      omig_in_min=dlog10(h*nu_min/mec2)
      omig_in_max=dlog10(h*nu_max/mec2)
      
      gam_min=dlog10(1.001d0)
      gam_max=dlog10(50.0d0)
      
      call gauleg(omig_in_min,omig_in_max
     &,           xleg_omig,wleg_omig,nleg_omig)
     
      call gauleg(gam_min,gam_max
     &,           xleg_gamm,wleg_gamm,nleg_gamm)
     
      call gauleg(-1.0d0,1.0d0,xleg,wleg,nleg)

      open(unit=50,file='./data/datain.txt',status='unknown')
      read(50,adaf)
      close(50) 
      write(*,adaf)
      open(unit=50,file='./data/nrows.txt',status='unknown')
      read(50,nrows) 
      close(50)
      write(*,nrows)
      
      nr=nend         
      rdin=(1.0d0+dsqrt(1.0d0-astar*astar))/2.0d0+0.2d0
c      rdin=0.5d3
      rdout=1.0d3
      do i=1,nleg_omig
      nu(i)=10.0**(xleg_omig(i))*mec2/h
      nulnu(i)=0.0d0
      fnu_unred(i)=0.0d0
      flux_comp(i)=0.0d0
      end do
   
      open(unit=10,file='./data/adaf.dat',form='formatted'
     &,status='unknown')
      
      j=0
      do i=1,nend
      j=j+1
      read(10,*)rd(j),mach(j),te(j),ti(j),tao(j),heig(j)
      rd(j)=10.0d0**(rd(j))
 
      if(rd(j).lt.rdin)then
      nr=j-1
      goto 11
      endif
      
      if(rd(j).gt.rdout)then
      j=j-1
      nr=j
      endif
           
      te(j)=10.0d0**(te(j))
      ti(j)=10.0d0**(ti(j))
      tao(j)=10.0d0**(tao(j))
      heig(j)=10.0d0**(heig(j))*rs*m
      nr=j
      end do
      
11    close(10)     

      temp_cs=k/mp/beta
      do i=1,nr
      cs(i)=dsqrt(temp_cs*(te(i)/mue+ti(i)/mui))
      rho(i)=tao(i)/heig(i)/dsqrt(2.0*pi)
      end do
      
      darea(1)=pi*((rd(1)+rdout)**2/4.0-(rd(1)+rd(2))**2/4.0)
     &           *m**2*rs**2      
      do i=2,nr-1
      darea(i)=pi*((rd(i-1)+rd(i))**2/4.0-(rd(i)+rd(i+1))**2/4.0)
     &           *m**2*rs**2
      end do
      darea(nr)=pi*((rd(nr)+rd(nr-1))**2/4.0-(rdin+rd(nr))**2/4.0)
     &           *m**2*rs**2 
     
c cal redshift at radius rd
c      zmax=1.0d10
c      do i=1,nr
c      velo=mach(i)*cs(i)
c      if(velo.gt.c)then
c      redshift(i)=1.0d-3
c      write(*,*)"vel.gt.c!"
c      else
c      redshift(i)=dsqrt(1.0-rdin/rd(i))*(1.0-velo**2/c**2)
c      end if
c      zmax=min(zmax,redshift(i))
c      end do
c      
c      do i=1,nr
c      redshift(i)=1.0d0
c      enddo  
c      write(*,*)zmax  
      end subroutine initial      

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c modification:
c Fri, Aug 29,2008
c > rho(i)=tao(i)/heig(i)/dsqrt(2.0*pi)
c < rho(i)=tao(i)/heig(i)/2.0
