************************************************************************
                            Spartan
************************************************************************

**Sp**\ ectrum of **a**\ dvection-dominated acc\ **r**\ e\ **t**\ ion flow in optic\ **a**\ lly thi\ **n**\ case.

This name is to memorialize the film

  300 Spartans

in which 300 Spartans profoundly impressed upon me with their 
bravery and persistence.

This code aims to calculate the spectrum of optically thin ADAF by 
solving the disk equation numerically. It adopts the notable geometrically
thin assumption so as to apply the integration over the height, i.e. to
obtain the height-integrated equation sets.


How to use this code
====================
under the command shell, type

.. code:: bash
   
   FC=gfortran make           # compile the code 
   ./disk                     # solve disk structure
   ./spec                     # cal intrinsic spectrum
   ./obs                      # cal observed spectrum

Note that your system must have installed Fortran 77/90 compiler. Change 
"**gfortran**" to the corresponding compiler in your system. 

The code will read the input option in **data/datain.txt** and calculate the structure and
spectrum of ADAF.


Input Option 
=============
The input option sees in **data/datain.txt**, which looks like::

  &ADAF 
  M = 5.2d7, 
  MDOT =7.0d-3,
  ALPHA =  0.1d0, 
  BETA =  0.5d0, 
  GAMMAA =  1.6666666666667d0,
  ASTAR = 0.0d0,
  VISF = 1.0d-3/
  
 &ODEINT
  NNU =  400,
  NRD  =  9001,
  RDOUT  =  1.0d4,
  NUMIN  = 1.0d9,
  NUMAX  = 1.0d22,
  DRSTEP  = 1.0d-3,
  DRSTEP_MIN=1.0d-4,
  DRSTEP_MAX=1.0d-2/
  
 &AUTO
  LINMAX = 2.5d0,
  LINMIN = 1.0d0/
 &OBS
  THETA_OBS = 30.0d0/


Here, the namelist **ADAF** includes ADAF parameters::

  M: black hole mass, in a unit of solar mass  
  MDOT: dimensionless accretion rate (accretion rate in a unit of Medd=1.39e18*M g/cm**2)
  ALPHA: viscosity parameter 
  BETA: the ratio of gas pressure to total pressure 
  GAMMAA: the adiabatic index
  ASTAR: the black hole spin parameter (-1~1)
  VISF: the fraction of viscous dissipation energy deposited into electrons.

The namelist **ODEINT** includes settings for solving ODE equations of ADAF::

  NNU: number of fequencies grid
  NRD: number of radial grid 
  RDOUT: outer radius in unit of Rs=2Rg
  NUMIN: minimum frequency
  NUMAX: maximum frequency 
  DRSTEP: initial radial step size 
  DRSTEP_MIN: minimum radial step size 
  DRSTEP_MAX: maximum radial step size 

The namelist **AUTO** includes the minimum and maximum specific angular momentum accreted into the black hole .

The namelist **OBS** includes the inclination angle (in a unit of degree) of the ADAF to the observer.


Outputs 
=======
The output data sees in data/, including

.. code:: bash 

    spectrum.dat               # intrinsic spectrum
    specobs.dat                # observed spectrum
    soltot.dat                 # disk solution
  
    sol_for_spec.dat           # disk solution used for cal spectrum
                               # with coarser radius grid
  
    spec/specxxx.txt           # spectrum at each radius
                               # see radius at sol_for_spec.dat

Plotting
========
See the Juypter notebook **plot_spartan.ipynb** in the folder **data/** for how to use and visualize the outputs.

Reference
=========
If you use this code, please cite our paper 
`Li et al. 2009, ApJ, 699, 513 <https://ui.adsabs.harvard.edu/abs/2009ApJ...699..513L/abstract>`_.

log
========

* **Wed, Mar 12, 2008**

  tag: auto_me  

  Try to make the code operate automatically, i.e. self-consistently determine
  the proper the lin once the boundary are given.

* **Thu, Mar 13, 2008**

  tag: autos  auto + sonic

  use transonic point to automate the calculation. I find that for larger lin, 
  the slope of surface density change from positive to negative, while for
  smaller lin, it is just reverse. Also, the correct solution must be transonic. 

* **Tue, May 6, 2008**

  tag: theta_obs set theta_obs as an input parameter in datain.txt.

* **Wed, Aug 27, 2008**

  A serious error foud. The factor 0.5 was missed in normalization (By Prof.Yuan).

* **Wed, Sep 3rd, 2008**

  1. Rewrite the spectrum code to improve the precision and computational speed.
  2. Rewrite the obs code using Prof. Yuan's subroutines.
