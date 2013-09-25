
Overview
===========



Working with Parameter Files
-------------------------------


Parameter files can be read into objects, manipulated, and written back out. Many parameters are accessible as object attributes::

   >>> par = mcfostpy.Paramfile('/path/to/some_parameters.para')
   >>> par.version
   2.17
   >>> par.grid_nrad
   100
   >>> 

The resulting object's attributes provide read/write access to all the settings in that parameter file.::

    >>> print "Distance:", par.distance, 'pc'
    >>> print "Wavelength range:", par.wavelengths_min, " to", par.wavelengths_max, "microns"
    >>> print "Changing distance to 200 pc"
    >>> par.distance = 200
    >>> print "Distance:", par.distance, 'pc'

    Distance: 140.0 pc
    Wavelength range: 0.1  to 3000.0 microns
    Changing distance to 200 pc
    Distance: 200 pc

     
This also includes some quantities which are computed from the parameters, not directly stated in the parameter file::

   >>> print "Wavelengths", par.wavelengths
   >>> print "Inclinations in RT mode: ", par.inclinations
   Wavelengths [  1.10859065e-01   1.36242824e-01   1.67438785e-01   2.05777785e-01
       2.52895390e-01   3.10801666e-01   3.81966929e-01   4.69427133e-01
       ...
       2.20195083e+03   2.70613864e+03]
   Inclinations in RT mode:  [ 45.]


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


**Input & Output**

The parameter file reading code attempts to handle several recent versions of the parameter file format, 
but this is not fully tested. 

.. warning::
    Currently not all parameters are implemented to be read in - this needs more work.

For output, the code should always output the most recent version it's been updated to handle. Currently this is 2.17. ::

   >>> paramfile.writeto('desired_output_filename.par')



.. warning::
    I just re-did most of the output formatting to use Python's new-style string formatting instead of
    old C-style % formatting. Needs more work still. In particular Boolean output as T/F instead of True/False needs to be 
    fixed. 

    2.17 code not fully tested. 


Working with Model Results
----------------------------





Detailed Function Reference
-----------------------------
.. autoclass:: mcfostpy.Paramfile
    :members:



