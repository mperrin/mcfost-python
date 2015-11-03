import numpy as np
import matplotlib
import emcee


def mcmccall(params,paramrange, directory, paramfile):
    """
    PARAMETERS
    ----------
    params - numpy array
         Parameters to be used in this call
         to emcee.

         For the purposes of ESO Halpha 569,
         these are inclination, scale height,
         dust mass, amax, beta, alpha, and
         rstar.

    paramrange - numpy array
         Array defining the parameter ranges
         to be tested in the emcee call.

    directory - string
         Directory to output all of the data for the
         MCMC calculations.

    paramfile - string
         Parameter file to begin with.

    USAGE
    -----
    Run emcee for a given parameter file and parameter
    ranges.

    """




#Below is from emcee website example
    ndim, nwalkers = 7, 100
    ivar = 1. / np.random.rand(ndim)
    p0 = [np.random.rand(ndim) for i in range(nwalkers)]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ivar])
    sampler.run_mcmc(p0, 1000)





def lnprob(params, directory, paramfile):

    """
    PARAMETERS
    ----------
    Params - numpy array
         Parameters to be used in this call
         to emcee.

         For the purposes of ESO Halpha 569,
         these are inclination, scale height,
         dust mass, amax, beta, alpha, and
         rstar.

    directory - string
         Directory to output all of the data for the
         MCMC calculations.

    paramfile - string
         Parameter file to begin with.


    USAGE
    -----
    Computes and returns the log of the likelihood
    distribution for a given model.

    """

    imageuncert, imagechi, seduncert, sedchi = mcmcwrapper(params, directory,paramfile)


    lnpimage = -0.5*np.log(2*np.pi)*imageuncert.size-0.5*imagechi-np.sum(-np.log(imageuncert))

    lnpsed = - 0.5*np.log(2*np.pi)*seduncert.size-0.5*sedchi-np.sum(-np.log(seduncert))

    return lnpimage + lnpsed




def mcmcwrapper(params, directory, paramfile):

    """
    PARAMETERS
    ----------
    Params - numpy array
         Parameters to be used in this call
         to emcee.

         For the purposes of ESO Halpha 569,
         these are inclination, scale height,
         dust mass, amax, beta, alpha, and
         rstar.

    directory - string
         Directory to output all of the data for the
         MCMC calculations.

    paramfile - string
         Parameter file to begin with.


    USAGE
    -----
    Takes a parameter file, the variable parameters
    and a directory. Creates the model image and
    SED, computes the chisqr, reads in the
    observables, and returns the uncertainties and
    Chi^2 values. This is called by the function
    that calculates the likelihood function.

    """

    return imageuncert, imagechi, seduncert, sedchi
