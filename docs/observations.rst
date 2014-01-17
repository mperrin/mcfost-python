.. _observations:

Working with Observations
============================

Examples and more docs to be added. 

Required File Formats and Organization
----------------------------------------

These are by design precisely compatible with the organization scheme adopted by ``MCRE``, the IDL Monte Carlo Results Explorer package by Perrin et al. 

The observations summary index
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All observations should be stored in one subdirectory. That directory must contain a file named ``summary.txt`` providing an index of the available files. This might look like::

        hktau_sed.txt           SED
        example1.fits           image   1.1
        example2.fits           image   10
        psf_image1.fits         psf     10
        tinytim_psf_F110W.fits  psf     1.1

In other words, it lists file names, what those files contain (image/psf/SED/visibility) and (for all types other than SED and SPECTRUM) the wavelength. There can be images and PSFs at multiple wavelengths, but can be only one SED file.   This file should be created by hand when you assemble the SED and images. 

Allowed input file types are as follows:
 
 * **SED**:                    SED as a text file (see below)
 * **SPECTRUM**:               Spectrum as text file. Same format as SED
 * **IMAGE**:          Image as a 2D FITS file. 
 * **IMAGE_UNCERT**:   Uncertainty image, 2D FITS file, must match image in size
 * **MASK**:                   Good/bad pixel mask; 2D FITS file, must match image in size

The following options from MCRE are not yet supported in the ``mcfost`` python package.
 * **POLIMAGE**:               Polarization percentage images, 2D FITS file, must match image in size
 * **PSF**:                    PSF. Must match image in pixel scale but can be smaller in size
 * **V2DATA**:         OIR V^2 interferometry data
 * **MMUVDATA**:       Millimeter V^2 interferometry data


The SED file
^^^^^^^^^^^^^^^^^


The name of this file can be anything, as set in the ``summary.txt`` file. Fluxes and uncertainties must be given in Janskys. The file format should be as follows::
 
   wavelength     flux in Jy     uncertainty in Jy     comment. 

Lines starting with '#' are interpreted as comments. This file is created manually when you assemble the SED from data or the literature. For example::

    # WAVELEN FLUX   FLUXERR         COMMENTS/SOURCE
    0.550000 0.0082 0.002           ACS_F555W     
    0.795000 0.0124 0.005           ACS_F814W    
    1.23500 0.019   0.005           2MASS_J     
    1.66200 0.029   0.005           2MASS_H    
    2.15900 0.050   0.015           2MASS_Ks  
    3.550   0.188   0.010           IRAC_1   
    4.439   0.182   0.010           IRAC_2
    100     0.00    5.2             IRAS_100
    <etc.>


To specify a point with an upper limit only, set the Flux column to zero, and the Uncertainy to the 1 sigma upper limit value, as in the last line of the above example. 

Images
^^^^^^^^

These should be FITS files, of the exact same file dimensions in X and Y as the MCFOST outputs. Image uncertainty files give the 1-sigma uncertainty for the corresponding image. 
Mask files are binary masks indicating which pixels should be included when calculating :math:`\chi^2` statistics. 


Observation classes detailed descriptions
---------------------------------------------

.. autoclass:: mcfost.Observations
   :members:

.. autoclass:: mcfost.models.ObservedSED
   :members:
   :undoc-members:

.. autoclass:: mcfost.models.ObservedImage
   :members:

