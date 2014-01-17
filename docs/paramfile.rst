.. _parfiles:

Working with Parameter Files
================================


Parameter files can be read into ``Paramfile`` objects, manipulated, and written back out. Many parameters are accessible as object attributes::

   >>> par = mcfost.Paramfile('examples/ref2.19.para')
   >>> par.version
   2.17
   >>> par.grid_n_rad
   100
   >>> par.RT_imin
   45
 
Attribute names generally correspond to variable names as stated in the MCFOST Manual. 
The attributes provide read/write access to all the settings in that parameter file.::

    >>> print "Distance:", par.distance, 'pc'
    >>> print "Wavelength range:", par.lambda_min, " to", par.lambda_max, "microns"
    >>> print "Changing distance to 200 pc"
    >>> par.distance = 200
    >>> print "Distance:", par.distance, 'pc'

    Distance: 140.0 pc
    Wavelength range: 0.1  to 3000.0 microns
    Changing distance to 200 pc
    Distance: 200 pc

     
Derived Quantities
----------------------

There are also some additional readable quantities which are computed from the parameters, not directly stated in the parameter file::

   >>> print "Wavelengths", par.wavelengths
   >>> print "Inclinations in RT mode: ", par.inclinations
   Wavelengths [  1.10859065e-01   1.36242824e-01   1.67438785e-01   2.05777785e-01
       2.52895390e-01   3.10801666e-01   3.81966929e-01   4.69427133e-01
       ...
       2.20195083e+03   2.70613864e+03]
   Inclinations in RT mode:  [ 45.]

Right now these are calculated only when the parameter file is created and don't automatically update if you change any parameters.
Likewise if you set their values directly nothing is updated for consistency and they're not used when writing files, so just treat these as read only.


Repeatable parameters for zones and dust properties
------------------------------------------------------

Parameters that can be repeated multiple times are stored as lists of dicts. Here's a model with two density zones::

   >>> print len(par.density_zones)
   2
   >>> print par.density_zones[0]['dust_mass']
   0.001

Since dust properties are different for every zone, and each zone can have multiple species of dust,
they're stored within the density zone dicts as another list of dicts. This can lead to some pretty complicated indexing::
 
   >>> print par.density_zones[0]['dust'][0]['aexp']   
   3.5  # aexp for the first dust species in the first zone
   >>> print par.density_zones[2]['dust'][1]['aexp']   
   2.0  # aexp for the second dust species in the third zone. 


As a convenience for the common case of there being only one of these, some shortcuts are defined:

 * **paramfile.density**: The first density_zone in a given model, equivalent to ``paramfile.density_zones[0]``
 * **paramfile.dust**: The first dust grain properties in a given model, equivalent to ``paramfile.density_zones[0]['dust'][0]``
 * **paramfile.star**: The first star in a given model, equivalent to ``paramfile.stars[0]``


Parameter File Input & Output
--------------------------------

The parameter file reading code attempts to handle several recent versions of the parameter file format, 
but this is not fully tested. 

.. warning::
        Versions of MCFOST parameter files prior to 2.15 are unlikely to work.     

For output, the code should always output the most recent version it's been updated to handle. Currently this is 2.19. ::

   >>> paramfile.writeto('desired_output_filename.par')  # the output file will be updated to 2.19
   >>>                                                   # regardless of the original version read in


Creating Grids of Parameter Files
--------------------------------------

The ``grid_generator`` function provides a flexible interface for creating a regular grid of parameter files. 
First, define a dictionary that lists the parameters to vary and the desired values. Then pass that and a
reference parameter file to the grid_generator function::

   >>> pars_to_vary = {'dust_amax': [100, 500., 1000], 'dust_mass': [1e-4, 3e-5, 1e-5], 'm_star': [2.0, 3.0]}
   >>> mcfost.grid_generator('ref2.19.para', pars_to_vary, filename_prefix='gridtest', start_counter=6)
   
The parameter files will be written to disk using filenames starting with the specified ``filename_prefix`` followed by an
integer counter starting at ``start_counter``. 


Paramfile class detailed documentation
-----------------------------------------

.. autoclass:: mcfost.Paramfile
   :members:

Function detailed documentation
-------------------------------------------

.. autofunction:: mcfost.grid_generator
