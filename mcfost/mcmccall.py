
#import numpy as np
#import matplotlib
#import emcee


#def mcmccall(params,paramrange, directory, paramfile):
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

"""
#Below is from emcee website example
    ndim, nwalkers = 7, 100
    ivar = 1. / np.random.rand(ndim)
    p0 = [np.random.rand(ndim) for i in range(nwalkers)]

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[ivar])
    sampler.run_mcmc(p0, 1000)
"""

import numpy as np
import matplotlib.pyplot as plt
#import emcee
import os
import subprocess 
#import models
import acor
#import utils
from mcfost.chisqr import image_chisqr
from mcfost.chisqr import sed_chisqr
from mcfost.paramfiles import Paramfile
from mcfost.models import ModelResults
from mcfost.models import Observations
#from mcfost import 

from emcee import PTSampler

import logging
_log = logging.getLogger('mcfost')


# Define the log likelihood 
def lnprobab(theta):

    """
    PARAMETERS
    ----------

    Params - numpy array
         Parameters to be used in this call
         to emcee.

         For the purposes of the  ESO Halpha 569,
         example, these are:
         theta[0] = inclination
         theta[1] = scale_height
         theta[2] = disk_mass
         theta[3] = max_grain_size
         theta[4] = alpha
         theta[5] = beta
         theta[6] = weights 

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

    imageuncert, imagechi, seduncert, sedchi = mcmcwrapper(theta)


    lnpimage = -0.5*np.log(2*np.pi)*imageuncert.size-0.5*imagechi-np.sum(-np.log(imageuncert))

    lnpsed = -0.5*np.log(2*np.pi)*seduncert.size-0.5*sedchi-np.sum(-np.log(seduncert))

    return theta[6]*lnpimage + (1.0-theta[6])*lnpsed

# Define the log Priors 
def lnprior(theta):
    
    inc    = theta[0]
    ho     = theta[1]
    mass   = theta[2]
    amax   = theta[3]
    alpha  = theta[4]
    beta   = theta[5]
    weight = theta[6]
    
    # include priors here
    if (inc < 65.0 or inc > 90.0):
        return -np.inf    
    
    if (ho < 5.0 or ho > 25.0):
        return -np.inf    

    if (np.log10(mass) < -5.0 or np.log10(mass) > -3.0):
        return -np.inf

    if (np.log10(amax) < 2.0 or np.log10(amax) > 4.0):
        return -np.inf

    if (alpha < -2.0 or alpha > 0.0):
        return -np.inf
    
    if (beta < 1.0 or beta > 1.5):
        return -np.inf

    if (weight < 0.3 or weight > 0.7):
        return -np.inf

    # otherwise ...
    return 0.0


def mcmcwrapper(theta):

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

    Step 1: Get parameters for this step for emcee.
    - Right a standalone function that will take in the list of parameters to vary and write the parameter files. This is unique to each object.
    Step 2: Write new parameter files with these values.
    Step 3: Run mcfost and store images and SEDs for this file. (Switch for running seds or images separately. )
    Step 4: Calculate the chi2 values for the given model and SED.
    Step 5: Pass the values to the log probability function. 
    """
    # STEP 1: This is passed via theta 
    # STEP 2:
    olddir=os.path.abspath(os.curdir)
    maindir = '/Users/swolff/Disks/mcfost-python/esoha569/'
    #par = Paramfile(maindir+'esoha569.para')
    par = Paramfile('/Users/swolff/Disks/mcfost-python/esoha569/esoha569.para')
    par.RT_imax = theta[0]
    par.RT_imin = theta[0]
    par.RT_n_incl = 1
    par.set_parameter('scale_height',theta[1])
    par.set_parameter('dust_mass',theta[2])
    par.set_parameter('dust_amax',theta[3])
    par.set_parameter('flaring',theta[4])
    par.set_parameter('surface_density',theta[5])
    # Do I need to write these, or just do this in memory?
    # write the parameter file in the default directory
    fnstring = "{0:0.4g}".format(theta[0])+'_'+"{0:0.4g}".format(theta[1])+'_'+"{0:0.4g}".format(theta[2])+'_'+"{0:0.4g}".format(theta[3])+'_'+"{0:0.4g}".format(theta[4])+'_'+"{0:0.4g}".format(theta[5])+'_'+"{0:0.4g}".format(theta[6])
    par.writeto(fnstring+'.para')
    modeldir = olddir+'/'+fnstring
    os.mkdir(modeldir)
    subprocess.call('chmod -R g+w '+modeldir,shell=True)
    subprocess.call('mv '+fnstring+'.para '+modeldir,shell=True)
    os.chdir(modeldir)

    #STEP 3:
    # run mcfost in the given directory
    subprocess.call('mcfost '+fnstring+'.para -rt',shell=True)
    subprocess.call('mcfost '+fnstring+'.para -img 0.8 -rt',shell=True)

    #STEP 4: 
    model = ModelResults('/Users/swolff/Disks/mcfost-python/'+fnstring)
    obs = Observations('/Users/swolff/Dropbox (GPI)/MCFOST_Testing/data')
    imagechi = image_chisqr(model,obs,wavelength=0.8)
    sedchi = sed_chisqr(model, obs, dof=1, write=True, save=False, vary_distance=False,vary_AV=True, AV_range=[0.0])#AV_range=np.arange(0.0,10.25,0.25))
    sedchi=sedchi[0]
    imagechi=imagechi[0]
    sedstring = 'SED {0}'.format(sedchi)+'/n'
    imagestring = 'Image {0}'.format(imagechi)
    
    f = open('chisqrd.txt','w')
    f.write(sedstring)
    f.write(imagestring)
    f.close()
    
    seduncertainty = obs.sed.nu_fnu_uncert
    imageuncertainty = obs.images[0.8].uncertainty
    #seduncert = np.sum(seduncertainty.value)/len(seduncertainty)
    
    #imageuncert = np.sum(imageuncertainty)/len(imageuncertainty)
    # remove model image and sed
    #subprocess.call('pwd',shell=True)
    #subprocess.call('rm -r data_th',shell=True)
    #subprocess.call('rm -r data_0.8',shell=True)
    os.chdir(olddir)

    #STEP 5:
    return imageuncertainty, imagechi, seduncertainty.value, sedchi




######################################################## 

ntemps = 20
nwalkers = 4096
ndim = 7

sampler=PTSampler(ntemps, nwalkers, ndim, lnprobab, lnprior, threads=4)


#############################################################
# Initialize the walkers. The best technique seems to be
# to start in a small ball around the a priori preferred position.
# Dont worry, the walkers quickly branch out and explore the
# rest of the space.

# inclination: w0
# scale height: w1
# dust mass: w2
# max grain size: w3
# alpha: w4
# beta: w5
# weights: w6

w0 = np.random.uniform(73.0,80.0,size=(ntemps,nwalkers))
w1 = np.random.uniform(10,20,size=(ntemps,nwalkers))
w2 = np.random.uniform(0.0001,0.001,size=(ntemps,nwalkers))
w3 = np.random.uniform(100.0,3000.0,size=(ntemps,nwalkers))
w4 = np.random.uniform(-2.0,0.0,size=(ntemps,nwalkers))
w5 = np.random.uniform(1.2,1.4,size=(ntemps,nwalkers))
w6 = np.random.uniform(0.4,0.6,size=(ntemps,nwalkers))

p0 = np.dstack((w0,w1,w2,w3,w4,w5,w6))

niter = 400
nburn = np.int(0.2*niter)


for p, lnprob, lnlike in sampler.sample(p0, iterations=nburn):
    pass

sampler.reset()

print 'Burn in complete'


"""
#Parallelize after this step. Look into whether it's better to 
#vary number of walkers or number of steps across different machines.
"""


for p, lnprob, lnlike in sampler.sample(p, lnprob0=lnprob,
                                        lnlike0=lnlike,
                                        iterations=niter, thin=10):
    pass

"""
diskname='esoha569'
inclres = np.ndarray.flatten(sampler.flatchain[0,:,:,0])
hores = np.ndarray.flatten(sampler.flatchain[0,:,:,1])
massres = np.ndarray.flatten(sampler.flatchain[0,:,:,2])
amaxres = np.ndarray.flatten(sampler.flatchain[0,:,:,3])
alphares = np.ndarray.flatten(sampler.flatchain[0,:,:,4])
betares = np.ndarray.flatten(sampler.flatchain[0,:,:,5])
weightres = np.ndarray.flatten(sampler.flatchain[0,:,:,5])
filename = diskname + '_inclinations.txt'
np.savetxt(filename, inclres)
filename = diskname + '_scale_heights.txt'
np.savetxt(filename, hores)
filename = diskname + '_dust_mass.txt'
np.savetxt(filename, massres)
filename = diskname + '_amax.txt'
np.savetxt(filename, amaxres)
filename = diskname + '_alpha.txt'
np.savetxt(filename, alphares)
filename = diskname + '_beta.txt'
np.savetxt(filename, betares)
filename = diskname + '_weights.txt'
np.savetxt(filename, weightres)
ChainStats = np.zeros((7,3))
ha = np.ndarray.flatten(sampler.flatchain[0,:,:,0])
ChainStats[0] = acor.acor(ha)
ha = np.ndarray.flatten(sampler.flatchain[0,:,:,1])
ChainStats[1] = acor.acor(ha)
ha = np.ndarray.flatten(sampler.flatchain[0,:,:,2])
ChainStats[2] = acor.acor(ha)
ha = np.ndarray.flatten(sampler.flatchain[0,:,:,3])
ChainStats[3] = acor.acor(ha)
ha = np.ndarray.flatten(sampler.flatchain[0,:,:,4])
ChainStats[4] = acor.acor(ha)
ha = np.ndarray.flatten(sampler.flatchain[0,:,:,5])
ChainStats[5] = acor.acor(ha)
ha = np.ndarray.flatten(sampler.flatchain[0,:,:,6])
ChainStats[6] = acor.acor(ha)
filename = diskname + '_ChainStats.txt'
np.savetxt(filename, ChainStats)


"""

