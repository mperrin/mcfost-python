import os
import subprocess
import logging
_log = logging.getLogger('mcfost')


import astropy

# Functions for actually running MCFOST models under Python control


def run_all_files(directory, loop_forever=False, **kwargs):
    """ Run all parameter files in a directory.

    Parameters
    -----------
    loop_forever : bool
        If set, once all files are complete, wait indefinitely for
        more parameter files to be written and then process those too.

    Notes
    -------
    See run_one_file for other allowable parameters.
    """

    import glob

    keepgoing = True
    while keepgoing:
        parfiles = glob.glob(os.path.join(directory, "*.par"))
        _log.info("Found {0} files to process: {1}".format(len(parfiles), ", ".join(parfiles)))
        for filename in parfiles:
            run_one_file(filename, **kwargs)
        keepgoing = loop_forever


def run_one_file(filename, wavelengths=[], move_to_subdir=True, delete_previous=True):
    """ Run a given parameter file.

    Parameters
    --------------
    wavelengths : iterable of floats
        Wavelengths to compute images for. Leave empty to just run the SED.
    move_to_subdir : bool
        Should we create a subdirectory for this parameter file, and move
        the parameter file into that subdirectory, before running it?
        This is useful for handling the output of grid_generator()
    delete_previous : bool
        Should we delete any previous calculation results if they already
        exist? Default is true.
    """
    if not isinstance(filename, str):
        raise TypeError("First argument to run_all must be a filename.")
    if not os.path.exists(filename):
        raise IOError("Nonexistent file: "+filename)

    _log.info("Running MCFOST for: "+filename)

    if move_to_subdir:
        _log.info("Relocating {0} to a subdirectory of the same name.".format(os.path.basename(filename)))
        modelname = os.path.splitext(os.path.basename(filename))[0]
        modeldir = os.path.splitext(os.path.dirname(filename))[0]
        #os.chdir(modeldir)
        #_log.info( "model directory"+modeldir)
        os.mkdir(modelname)
        subprocess.call('chmod -R g+w '+modelname,shell=True)
        subprocess.call('mv '+filename+' '+modelname,shell=True)

        filename = os.path.join(modelname, os.path.basename(filename))


    run_sed(filename, delete_previous=delete_previous)

    for wl in wavelengths:
        run_image(filename, wl, delete_previous=delete_previous)

    _log.info("Calculation complete.")


def run_sed(filename,raytrace=True, delete_previous=True, *args ):
    """ Run a MCFOST calculation of the SED

    Parameters
    ------------
    raytrace : bool
        Should raytrace mode be used? Default is true.
    """

    if not isinstance(filename, str):
        raise TypeError("First argument to run_sed must be a filename.")
    if not os.path.exists(filename):
        raise IOError("Nonexistent file: "+filename)

    _log.info("Computing SED for {0}".format(filename))

    directory = os.path.dirname(os.path.abspath(filename))

    # Possibly clean up previous outputs
    # MCFOST itself will do this once only and after that
    # will fail, so it's useful to do this cleaning here.
    outputpath = os.path.join(directory, 'data_th')
    if delete_previous and os.path.isdir(outputpath):
        _log.debug("Removing previous calculation results")
        import shutil
        shutil.rmtree(outputpath, ignore_errors=True)

    cmdstr = 'mcfost '+os.path.basename(filename)
    if raytrace: cmdstr += " -rt"

    subprocess.call("echo  '>> "+ cmdstr+"' >> output.log",shell=True, cwd=directory)
    result = subprocess.call(cmdstr+' >> output.log',shell=True, cwd=directory)
    subprocess.call('chmod -R g+w *',shell=True, cwd=directory)
    _log.info("SED results written to {0}".format(outputpath))
    _log.debug("Result: {0}".format(result))


def run_image(filename, wavelength, raytrace=True,  delete_previous=True, *args):
    """ Run a MCFOST calculation of the image

    Parameters
    ------------
    wavelength : float
        Wavelength in microns for the image
    raytrace : bool
        Should raytrace mode be used? Default is true.

    """
    if not isinstance(filename, str):
        raise TypeError("First argument to run_image must be a filename.")
    if not os.path.exists(filename):
        raise IOError("Nonexistent file: "+filename)

    _log.info("Computing image at {1} microns for {0}".format(filename, wavelength))

    directory = os.path.dirname(os.path.abspath(filename))

    # Possibly clean up previous outputs
    # MCFOST itself will do this once only and after that
    # will fail, so it's useful to do this cleaning here.
    outputpath = os.path.join(directory, 'data_'+str(wavelength))
    if delete_previous and os.path.isdir(outputpath):
        _log.debug("Removing previous calculation results")
        import shutil
        shutil.rmtree(outputpath, ignore_errors=True)



    cmdstr = 'mcfost '+os.path.basename(filename)+' -img '+str(wavelength)
    if raytrace: cmdstr += " -rt"

    subprocess.call("echo  '>> "+ cmdstr+"' >> output.log",shell=True, cwd=directory)
    subprocess.call(cmdstr +' >> output.log',shell=True, cwd=directory)
    subprocess.call('chmod -R g+w *',shell=True, cwd=directory)

def run_dust_prop(filename,  delete_previous=True,  sed=False, *args):
    """ Run a MCFOST calculation of the dust properties

    Parameters
    ------------
    filename : string
        Filename
    delete_previous : bool
    sed : bool
        if True, it computes the SED and the dust properties (option +dust_prop),
        else it computes only the dust properties (option -dust_prop)

    """
    if not isinstance(filename, str):
        raise TypeError("First argument to run_dust_prop must be a filename.")
    if not os.path.exists(filename):
        raise IOError("Nonexistent file: "+filename)

    _log.info("Computing dust properties for {0}".format(filename))

    directory = os.path.dirname(os.path.abspath(filename))

    # Possibly clean up previous outputs
    # MCFOST itself will do this once only and after that
    # will fail, so it's useful to do this cleaning here.
    outputpath = os.path.join(directory, 'data_dust')
    if delete_previous and os.path.isdir(outputpath):
        _log.debug("Removing previous calculation results")
        import shutil
        shutil.rmtree(outputpath, ignore_errors=True)
    if sed:
        cmdstr = 'mcfost '+os.path.basename(filename)+' +dust_prop'
    else:
        cmdstr = 'mcfost '+os.path.basename(filename)+' -dust_prop'
    subprocess.call("echo  '>> "+ cmdstr+"' >> output.log",shell=True, cwd=directory)
    subprocess.call(cmdstr +' >> output.log',shell=True, cwd=directory)
    subprocess.call('chmod -R g+w *',shell=True, cwd=directory)


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
    if isinstance(base_paramfile, str):
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
