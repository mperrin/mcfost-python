import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy.interpolate

from . import models

import logging
_log = logging.getLogger('mcfost')

def sed_chisqr(modelresults, observations, dof=1,
        write=True,
    plot=True, save=True,
    vary_distance=False, distance_range=None,
    vary_AV=False, AV_range=[0,10],
    **kwargs):

    """ Compute chi^2 for a given model 

    Not written yet - this is just a placeholder

    Parameters
    -----------
    dof : int
        Number of degrees of freedom to use in computing the
        reduced chi^2. Will be incremented appropriately if 
        vary_distance or vary_AV are set.
    vary_distance : bool
        Allow the distance to the target to vary
    distance_range : 2 element iterable
        Min and max allowable distance, in pc
    vary_AV : bool
        Allow optical extinction (A_V) to vary
    AV_range : 2 element iterable
        Min and max allowable A_V
    RV : float
        Reddening law color parameter for extinction
    write : bool
        Write output to disk as chi^2 results FITS file in the model directory?
        Default is True
    save : bool
        Save results to disk as FITS bintable? Default is True.
    plot : bool
        Display plot of chi^2 fit?

    """

#    if not isinstance(modelresults, models.ModelResults):
#        print type(modelresults)
#        raise ValueError("First argument to sed_chisqr must be a ModelResults object")
#
#    if not isinstance(observations, models.Observations):
#        raise ValueError("Second argument to sed_chisqr must be an Observations object")
    my_dof = dof
    if vary_distance:
        raise NotImplementedError("Varying distance not yet implemented")
        my_dof += 1
    if vary_AV:
        _log.info("Computing chi^2 while allowing A_V to vary between {0} and {1} with R_V={2}".format(AV_range[0], AV_range[1], RV) )
        my_dof += 1


    # observed wavelengths and fluxes 
    obs_wavelengths = observations.sed.wavelength
    obs_nufnu = observations.sed.nu_fnu
    obs_nufnu_uncert = observations.sed.nu_fnu_uncert
    # model wavelengths and fluxes. Fluxes will be a 2D array: [ninclinations, nwavelengths]
    mod_wavelengths = modelresults.sed.wavelength
    mod_nufnu= modelresults.sed.nu_fnu
    mod_inclinations = modelresults.parameters.inclinations

    if plot:
        observations.sed.plot()

    ninc = mod_inclinations.size
    chi2s = np.zeros(ninc)
    avs = np.zeros(ninc)
    rvs = np.zeros(ninc)
    distances = np.zeros(ninc) + modelresults.parameters.distance
    # iterate over inclinations, computing estimated observations and fitting
    # allow reddening to vary
    for i in range(ninc):

        interpolator = scipy.interpolate.interp1d(mod_wavelengths, mod_nufnu[i], kind='linear', copy=False)
        est_mod_nufnu = interpolator(obs_wavelengths)

        if plot:
            modelresults.sed.plot(inclination=mod_inclinations[i], overplot=True)
            ax = plt.gca()
            color = ax.lines[-1].get_color()
            ax.plot(obs_wavelengths, est_mod_nufnu, color=color, marker='s', linestyle='none')

        if vary_AV:
            raise NotImplementedError('not yet')
        else:
            chi2 = ((obs_nufnu.value - est_mod_nufnu)**2 / obs_nufnu_uncert.value**2).sum()

        _log.info( "inclination {0} : {1:4.1f} deg has chi2 = {2:5g}".format(i, mod_inclinations[i], chi2))
        #print "obs = ", obs_nufnu.value
        #print "mod = ", est_mod_nufnu


        chi2s[i] = chi2

    if save:
        import astropy.table
        import os
        tabledict = {'INCLINATION': mod_inclinations, 'SED_CHISQ': chi2s, 'AV':avs, 'RV':rvs, 'DISTANCE':distances}
        meta = {'SOURCE': 'Python mcfost.chisqr.sed_chisqr()', 
                'VARY_AV': vary_AV,
                'VARYDIST': vary_distance,
                'AVRANGE': str(AV_range),
                'DISTRNGE': str(distance_range) }
        tab = astropy.table.Table(tabledict, meta=meta)

        tab.write(os.path.join(modelresults.directory, 'observables', 'py_sed_chisq.fits'))

    return chi2s


def image_chisqr(modelresults, observations, wavelength=None, write=True,
        normalization='total', registration='integer_pixel_shift'):
    """
    Not written yet - this is just a placeholder

    Parameters
    ------------
    wavelength : float
        Wavelength of the image to compute the chi squared of. 
    write : bool
        If set, write output to a file, in addition to displaying
        on screen. 

    """
    raise NotImplementedError("Not implemented yet")




# Required output format for Christophe's GA tools: 
#    
#
#

