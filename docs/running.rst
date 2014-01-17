.. _running:

Running MCFOST
======================


Once you have created one or more parameter files, you will want to actually run MCFOST on them. This can be
commanded from within Python::

  >>> mcfost.run_one_file('myparfile.par', wavelengths=[0.6])
  INFO:mcfost:Running MCFOST for: ./myparfile.par
  INFO:mcfost:Relocating myparfile.par to a subdirectory of the same name.
  INFO:mcfost:Computing SED for myparfile/myparfile.par
  INFO:mcfost:Computing image at 0.6 microns for myparfile/myparfile.par
  INFO:mcfost:Calculation complete.

The above line will run MCFOST to compute the SED and an 0.6 micron image for the disk parameters specified in that file. By default it will first
create a subdirectory named ``myparfile``, and then move ``myparfile.par`` into that subdirectory, before running MCFOST itself to create 
the ``data_*`` output subdirectories. Thus you'll end up with a subdirectory called ``myparfile`` containing all the computed results for that 
parameter file. If you don't want this relocation to occur then just set ``move_to_subdir=False`` in the call to ``run_one_file``.

The MCFOST process's standard output will be saved in the file ``output.log`` in that directory. 

Note that there are individual functions ``run_sed`` and ``run_image`` that compute the SED and image respectively.


Running many files
--------------------

To repeat that process for many files in a directory, use ``run_all_files``::

  >>> mcfost.run_all_files('/some/path/to/files')
  INFO:mcfost:Found 18 files to process: example1.par, example2.par, [...] example18.par
  INFO:mcfost:Running MCFOST for: ./example1.par
  [and so on]



Function detailed documentation
-------------------------------------------

.. autofunction:: mcfost.run_one_file

.. autofunction:: mcfost.run_sed

.. autofunction:: mcfost.run_image

.. autofunction:: mcfost.run_all_files
