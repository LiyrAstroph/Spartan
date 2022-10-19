c=======================================================================
c************************     Spartan      *****************************
c                        
c                         By Li Yan-Rong
c                  supported by Prof. Wang Jian-Min
c                           and Prof. Yuan Ye-Fei
c                    Email: liyanrong@ihep.ac.cn
c                        2008-01-16--2008-03-
c=======================================================================
c
c The main file for disk structure
c
c To compile this code, following files are needed:
c 
c diskset.f   ####   set up the code, boundary and initial conditions
c disksol.f   ####   solve the equations set
c diskrad.f   ####   calculate the radiaiton term
c specfun.f   ####   special function used in the code
c datagen.f   ####   data generation for spectrum and obs
c diskvar.f   ####   variables definition
c constat.f   ####   some constants
c


      include 'diskset.f'
      include 'disksol.f'
      include 'diskrad.f'
      include 'specfun.f'
      include 'datagen.f'
      
      program main
      implicit none

      write(*,*)"======================================================"
      write(*,*)"==================     Spartan      =================="
      write(*,*)"======================================================"           
      call initial()
      call boundary()
      call autocal()
      call datagen()
      
      end program main

