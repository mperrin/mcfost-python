
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


To compare observed SED's and images with the MCFOST model outputs, 
there are included functions for calculating :math:`\chi^2` values.
See :ref:`chisqr` for details. Once you're read in the observations 
and model data into the structures defined above, compute :math:`\chi^2` 
values:

.. code-block:: python

  >>> imagechi = image_chisqr(model,obs,wavelength=0.8)
  >>> sedchi = sed_chisqr(model, obs)

The larger goal is to provide an interface for MCMC model fitting 
via `emcee <http://dan.iel.fm/emcee/current>`_.  
An example of this is provided in mcmccall.py. In addition to the MCMC 
framework described in detail at `emcee <http://dan.iel.fm/emcee/current>`_, 
this includes three functions to calculate the likelihood and prior 
distributions for each combination of model parameters. 

``lnprior`` defines the prior distribution used for each parameter. This is 
left up to the user. In the example provided, we use constrain the acceptable 
parameter value ranges using step functions.

``lnprobab`` calculates the likelihood of the model given the data using the 
model uncertainties and the calculated :math:`\chi^2` values. It is left up 
to the user how to combine likelihoods for images at different wavelengths, 
SEDs, and spectra. In this example, we have left the weighting between the 
image and SED as a free parameter in the MCMC fit. This function obtains the 
uncertainties and :math:`\chi^2` values from the third and final function. 

``mcmcwrapper`` is responsible for returning the uncertainties and 
:math:`\chi^2` values given a set of model parameter values to be used for 
the likelihood calculation. The required inputs to this function are a 
parameter file, the parameter values for a step in the MCMC chain, and a 
directory. It then creates the model image and SED, reads in the observables, 
and returns the uncertainties :math:`\chi^2` values. The ``mcmcwrapper`` 
function employs all of the tools for interacting with the observations, model 
data, parameter files, and MCFOST. 
