from . import models

#placeholder for eventually having chi sqr fitting functions here

def sed_chisqr(modelresults, observations, write=True,
    clobber=False, plot=True,
    vary_distance=False, distance_range=None,
    log_chisqr = False):

    """ Compute chi^2 for a given model 

    Not written yet - this is just a placeholder

    Parameters
    -----------
    log_chisqr : bool
        Use log of Chi^2 instead?
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
    clobber : bool
        Overwrite existing results? Default is False
    plot : bool
        Display plot of chi^2 fit?

    """
    raise NotImplementedError("Not implemented yet")

    if not isinstance(modelresults, models.ModelResults):
        raise ValueError("First argument to sed_chisqr must be a ModelResults object")

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

