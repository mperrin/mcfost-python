import astropy

# Functions for actually running MCFOST models under Python control


def run_all(filename_or_Paramfile, wavelengths=[]):
    run_sed(filename_or_Paramfile)

    for wl in wavelengths:
        run_image(filename_or_Paramfile, wl)

def run_sed(filename_or_Paramfile,raytrace=True,  *args ):
    """ Run a MCFOST calculation of the SED"""
    _log.info("Computing SED for {0}".format(filename_or_Paramfile))

    if isinstance(filename_or_Paramfile, str):
        filename = filename_or_Paramfile
    else: # assume it's an instance of a Paramfile structure
        filename = filename_or_Paramfile.filename

    directory = os.path.dirname(os.path.abspath(filename))




def run_image(filename_or_Paramfile, wavelength, raytrace=True, *args):
    """ Run a MCFOST calculation of the image """
    _log.info("Computing image at {1} microns for {0}".format(filename_or_Paramfile, wavelength))
    pass



#=====================================================================
# Functions for MCMC interface


def mcmc_evaluate_params(modelparams, base_paramfile, paramnames, observations):
    """ Given a set of parameters, set up and evaluate the model parameter file,
    and return the log of the probability relative to the observations 


    Parameters
    ------------
    modelparams : tuple of floats
        Free parameters to evaluate for this particular instance
    base_paramfile : mcfost.Paramfile instance
        The base parameter file to use for all other parameters
    paramnames : tuple of strings 


    """


def grid_generator(base_paramfile, paramsdict, filename_prefix=None, start_counter=1):
    """Generate a grid of parameter files by varying specified MCFOST parameters

    Parameters
    -------------
    base_paramfile : filename
        MCFOST parameter file to use as the basis for this grid.
    paramsdict : dictionary of lists
        Dictionary specifying the parameters to be varied. See notes below
        for example. 
    filename_prefix : string, optional
        Filename prefix to use for the output parameter files. If not set, the
        name of the base_paramfile is used. 
    start_counter : integer, optional
        Starting value for the integer counter that distinguishes each output filename.
        Default value is 1. 

    Returns
    ----------
    filename_list : list
        List of output filenames

    Notes
    --------

    The parameters to vary are specified using a dict, the keys of which are parameter names
    and the values of which give the set of values for each parameters.  For instance:

        pars = {'dust_amax': [100, 500., 1000], 'dust_mass': [1e-4, 3e-5, 1e-5]}

    would create a grid of 9 models iterating over those values for total disk dust mass
    and maximum grain size. 

    For the list of allowable parameter name strings, any attribute names of the 
    Paramfile object are allowed, plus the shortcuts defined in the Paramfile.set_parameter()
    function.
    """
    import itertools

    # Read param file from disk into an object, if necessary
    if isinstance(base_paramfile, basestring):
        from .paramfiles import Paramfile
        base_paramfile = Paramfile(base_paramfile)
    if filename_prefix is None: 
        filename_prefix = os.path.splitext(os.path.basename(base_paramfile.filename))[0] 

    counter = start_counter

    keys = paramsdict.keys()

    # Of course the Python library has exactly the function we want here!
    params_generator = itertools.product(*paramsdict.values())

    output_fns = []

    for params in params_generator:
        for parname, parvalue in zip(keys, params):
            base_paramfile.set_parameter(parname, parvalue)

        output_fn = "{0}_{1:05d}.par".format(filename_prefix, counter)
        base_paramfile.writeto(output_fn)
        outout_fns.append(output_fn)
        counter+=1

    return output_fns


