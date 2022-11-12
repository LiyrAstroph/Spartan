************************************************************************
                            Spartan
************************************************************************

**Sp**\ ectrum of **a**\ dvection-dominated acc\ **r**\ e\ **t**\ ion flow in optic\ **a**\ lly thi\ **n**\ case.

This name is to memorialize the film

  300 Spartans

in which 300 Spartans profoundly impressed upon me with their 
bravery and persistence.

This code aims to calculate the spectrum of optically thin ADAF by 
solving the disk equation numerically. It adopt the notable geometrically
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

The input option sees in data/datain.txt.

The output data sees in data/, including

.. code:: bash 

    spectrum.dat               # intrinsic spectrum
    specobs.dat                # observed spectrum
    soltot.dat                 # disk solution
  
    sol_for_spec.dat           # disk solution used for cal spectrum
                               # with coarser radius grid
  
    spec/specxxx.txt           # spectrum at each radius
                               # see radius at sol_for_spec.dat

Reference
=========
If you use this code, please cite our paper 
`Li et al. 2008, ApJ, 699, 513 <https://ui.adsabs.harvard.edu/abs/2009ApJ...699..513L/abstract>`_

log
========

Wed, Mar 12, 2008
tag: auto_me

Try to make the code operate automatically, i.e. self-consistently determine
the proper the lin once the boundary are given.

Thu, Mar 13, 2008
tag: autos  auto + sonic

use transonic point to automate the calculation. I find that for larger lin, 
the slope of surface density change from positive to negative, while for
smaller lin, it is just reverse. Also, the correct solution must be transonic. 

Tue, May 6, 2008
tag: theta_obs set theta_obs as an input parameter in datain.txt.

Wed, Aug 27, 2008
A serious error foud. The factor 0.5 was missed in normalization (By Prof.Yuan).

Wed, Sep 3rd, 2008
1. Rewrite the spectrum code to improve the precision and computational speed.
2. Rewrite the obs code using Prof. Yuan's subroutines.
