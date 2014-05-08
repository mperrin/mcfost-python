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


def lnprob_wrapper_esoha569(params, observations, )





def lnprob(params, directory, paramfile, weights):

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
    
    weights - numpy array
         Two element array giving the weights to 
         apply to the sed [0] and image [1] chi2.


    USAGE
    -----
    Computes and returns the log of the likelihood 
    distribution for a given model.

    """

    imageuncert, imagechi, seduncert, sedchi = mcmcwrapper(params, directory,paramfile)


    lnpimage = -0.5*np.log(2*np.pi)*imageuncert.size-0.5*imagechi-np.sum(-np.log(imageuncert))

    lnpsed = - 0.5*np.log(2*np.pi)*seduncert.size-0.5*sedchi-np.sum(-np.log(seduncert))

    return weights[0]*lnpimage + weights[1]*lnpsed




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

    """
    Step 1: Get parameters for this step for emcee.
    - Right a standalone function that will take in the list of parameters to vary and write the parameter files. This is unique to each object.
    Step 2: Write new parameter files with these values.
    Step 3: Run mcfost and store images and SEDs for this file. (Switch for running seds or images separately. )
    Step 4: Calculate the chi2 values for the given model and SED.
    Step 5: Pass the values to the log probability function. 
    """


    return, imageuncert, imagechi, seduncert, sedchi
