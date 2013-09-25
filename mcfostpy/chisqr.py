
#placeholder for eventually having chi sqr fitting functions here

def sed_chisqr(modelresults, observations, write=True,
    clobber=False, plot=True,
    vary_distance=False, distance_range=None,
    log_chisqr = False):

    """ Compute chi^2 for a given model 

    **Not written yet - this is just a placeholder **

    Parameters
    -----------
    write : bool
        Write output to disk as chi^2 results FITS file in the model directory?
        Default is True
    clobber : bool
        Overwrite existing results? Default is False
    plot : bool
        Display plot of chi^2 fit? 
    log_chisqr : bool
        Use log of 
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
    """
    raise NotImplementedError("Not implemented yet")


def image_chisqr(modelresults, observations, wavelength=None, write=True):
    """
    **Not written yet - this is just a placeholder **
    """
    raise NotImplementedError("Not implemented yet")



