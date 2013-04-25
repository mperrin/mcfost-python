
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



Detailed Function Reference
-----------------------------
.. autoclass:: mcfostpy.Paramfile
    :members:



