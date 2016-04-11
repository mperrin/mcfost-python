
Chi Squared Fitting
-----------------------

We've included functions for comparing the observations and models for both images and SEDs via :math:`\chi^2` values. To use these functions, the user must first read the observations and MCFOST model into the correct data structure as defined in :ref:`observations` and :ref:`modelresults`. The :math:`\chi^2` values can then be computed as described here:

.. code-block:: python

  >>> from mcfost.chisqr import image_chisqr
  >>> from mcfost.chisqr import sed_chisqr
  >>> imagechi = image_chisqr(model,obs,wavelength=0.8)
  >>> sedchi = sed_chisqr(model, obs)

The image :math:`\chi^2` function provides several options for the user. The model image can be convolved with an instrumental PSF. The model and observed images can be normalized to the 'total' or max of the observed data. Image registration can be done using either 'sub_pixel' or 'integer_pixel' accuracy.  

The SED :math:`\chi^2` calculator allows the user to vary the distance to the central star and the extinction properties. This was adapted from the fit_dist_extinct3.pro MCRE file designed to allow distance and extinction to vary when computing the SED chsqr. For a given inclination, this function returns the best fit distance, extinction, Rv and chisqrd value.

Note: MCMC works best with continuous parameter distributions. While MCFOST is capable or producing a range of inclination values for every model SED and image file created, we recommend generating only one inclination per MCMC walker step. 


Function Detailed documentation
--------------------------------------------

.. autofunction:: mcfost.chisqr.sed_chisqr

.. autofunction:: mcfost.chisqr.image_chisqr
