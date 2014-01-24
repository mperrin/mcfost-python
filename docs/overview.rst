
Overview
===========


Parameter files can be read in to ``Paramfile`` objects, manipulated, and written back out.  See :ref:`parfiles` for details.

.. code-block:: python

   >>> par = mcfost.Paramfile('ref2.19.para')
   >>> par.distance = 200
   >>> par.writeto('example.par')

MCFOST calculations may then be executed for a given parameter file, or for many parameter files at once (see :ref:`running`):

.. code-block:: python

   >>> par = mcfost.run_one_file('example.par')
   INFO:mcfost:Running MCFOST for: ./example.par
   [...]
   INFO:mcfost:Calculation complete.

 
After which the results may be examined and plotted (see :ref:`modelresults`):

.. code-block:: python

  >>> res = mcfost.ModelResults('example')
  >>> res.describe()
  Model results in examples for ref2.19.para
    Model has 15 inclinations from 60.0 to 90.0
    SED computed from 0.1 - 3500.0 microns using 30 wavelengths
    Images computed for 1 wavelengths: [ 0.8] micron 
  >>> res.sed.plot(title='Example SED plot', nlabels=4)

.. plot::
  :height: 300

  import mcfost
  res = mcfost.ModelResults('example')
  res.sed.plot(title='Example SED plot', nlabels=4)


You presumably also have observed data you wish to compare to, and
those can be loaded and displayed as well. 
See :ref:`observations` for details, including the required file formats.
The interfaces for working with 
model results and observations are similar:

.. code-block:: python

  >>> obs = mcfost.observations('example_obs')
  >>> obs.describe()
  >>> obs.sed.plot()


There are included functions for calculating :math:`\chi^2` values for model fitting, currently supporting SEDs only. Image fitting is a work in progress. 
The long term goal is to provide an interface for MCMC model fitting via `emcee <http://dan.iel.fm/emcee/current>`_.  
 




