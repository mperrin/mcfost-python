#!/usr/bin/env python

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
import triangle
# Necessary for mpipool:
# import emcee
#from emcee.utils import MPIPool
#import sys

#from mcfost import 

from emcee import PTSampler

import logging
_log = logging.getLogger('mcfost')



# Set a few Parameters
model_code = 'grater'
belts = 1
maindirectory = '/astro/4/mperrin/eods/Grater_EMCEE/'
parameters = ['g','radius','inn_slope','out_slope','inclination','xoffset','yoffset','pa','flux']
ndim = 9
imwavelength = 0.002

# Define the log likelihood 
def lnprobab(theta):



    """
    PARAMETERS
    ----------
    Theta:
    Parameters to vary. 
         theta[0] = g; henyey-grenstein parameter
         theta[1] = radius of the belt (AU)
         theta[2] = inner slope (power law exp)
         theta[3] = outer slope (power law exp)
         theta[4] = Inclination (degrees, 90 = edge-on)
         theta[5] = x offset of center (N/S)
         theta[6] = y offset of center (E/W)
         theta[7] = Position Angle (Degrees CCW from N)
         theta[8] = Flux scaling

         ** If there are multiple belts, parameters 0-8 repeat. **
         

    USAGE
    -----
    Computes and returns the log of the likelihood 
    distribution for a given model.

    """


    imageuncert, imagechi = mcmcwrapper(theta)

    lnpimage = -0.5*np.log(2*np.pi)*imageuncert.size-0.5*imagechi-np.sum(-np.log(imageuncert))

    return lnpimage





# Define the log Priors 
def lnprior(theta):
    
    g =           theta[0]
    radius =      theta[1]
    inn_slope =   theta[2]
    out_slope =   theta[3]
    inclination = theta[4]
    xoffset =     theta[5]
    yoffset =     theta[6]
    pa =          theta[7]
    flux =        theta[8]
    
    # include priors here
    if (g < -1.0 or g > 1.0):
        return -np.inf    
    
    if (radius < 150.0 or radius > 300.0):
        return -np.inf    

    if (inn_slope < 5.0 or inn_slope > 50.0):
        return -np.inf

    if (out_slope < -30.0 or out_slope > -2.0):
        return -np.inf

    if (inclination < 40 or inclination > 70):
        return -np.inf
    
    if (xoffset < -10.0 or xoffset > 10.0):
        return -np.inf

    if (yoffset < -10.0 or yoffset > 10.0):
        return -np.inf

    if (pa < 350.0 or pa > 360.0):
        return -np.inf

    if (flux < 1.5e-7 or flux > 1.5e-5):
        return -np.inf

######################################################
###################### Outer Ring ####################
######################################################

"""
    gb =           theta[9]
    radiusb =      theta[10]
    inn_slopeb =   theta[11]
    out_slopeb =   theta[12]
    inclinationb = theta[13]
    xoffsetb =     theta[14]
    yoffsetb =     theta[15]
    pab =          theta[16]
    fluxb =        theta[17]
    
    # include priors here
    if (gb < -1.0 or gb > 1.0):
        return -np.inf    
    
    if (radiusb < 300.0 or radiusb > 500.0):
        return -np.inf    

    if (inn_slopeb < 2.0 or inn_slopeb > 50.0):
        return -np.inf

    if (out_slopeb < -30.0 or out_slopeb > -2.0):
        return -np.inf

    if (inclinationb < 40 or inclinationb > 70):
        return -np.inf
    
    if (xoffsetb < -50.0 or xoffsetb > 50.0):
        return -np.inf

    if (yoffsetb < -50.0 or yoffsetb > 50.0):
        return -np.inf

    if (pab < 350.0 or pab > 360.0):
        return -np.inf

    if (fluxb < 1.0e-8 or fluxb > 1.0e-6):
        return -np.inf
 

"""


    # otherwise ...
    return 0.0




def mcmcwrapper(theta):

    """
    PARAMETERS
    ----------
    Params - numpy array
         Parameters to be used in this call
         to emcee.

    USAGE
    -----
    Takes a parameter file, the variable parameters
    and a directory. Creates the model image and
    computes the chisqr, reads in the 
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
    maindir = maindirectory
 
    # Do I need to write these, or just do this in memory?
    # write the parameter file in the default directory
    fnstring = "{0:0.4g}".format(theta[0])
    for i in np.arange(ndim-1):
        fnstring += '_'+"{0:0.4g}".format(theta[i+1])
    modeldir = olddir+'/'+fnstring
    

    try:
        os.mkdir(modeldir)
    except:
        subprocess.call('rm -r '+modeldir,shell=True)
        print 'removed model directory: ', modeldir 
        os.mkdir(modeldir)

    subprocess.call('chmod -R g+w '+modeldir,shell=True)
    subprocess.call('mv '+fnstring+'.para '+modeldir,shell=True)
    os.chdir(modeldir)

    # Now write the parameter file from the theta values.
    f = open(fnstring+'.para','w')
    f.write(str(1)+"\n")
    for i in np.arange(ndim):
       f.write(str(theta[i])+"\n")
    f.close() 

    os.mkdir('data_0.002')

    #modeldir 3:
    # run grater in the given directory
    script_file_name = 'grater_mcmc, '+olddir+', '+modeldir
    subp = subprocess.Popen("idl -quiet -e" + script_file_name + "'",
            stderr = subprocess.PIPE, stdout = subprocess.PIPE, shell = True)

    (idl_stdout, _) = subp.communicate()

    #STEP 4: 
    try:
        model = ModelResults(maindir+fnstring)
    except IOError:
        #raise IOError('Model Failed')
        os.chdir(olddir)
        return np.asarray([1000000.0,1000000.0]), 1000000.0
        
    try:
        obs
    except NameError:
        #print "well, it WASN'T defined after all!"
        obs = Observations(maindir+'data')

    imagechi = image_chisqr(model,obs,wavelength=imwavelength, inclinationflag=False, convolvepsf=False)
    imagechi = imagechi[0]
    imagestring = 'Image {0}'.format(imagechi)
    
    f = open('chisqrd.txt','w')
    f.write(imagestring)
    f.close()
    
    imageuncertainty = obs.images[imwavelength].uncertainty
    

    subprocess.call('pwd',shell=True)
    os.chdir(olddir)


    # force close all model images:
    model.images.closeimage()  

    #STEP 5:
    return imageuncertainty, imagechi



######################################################## 

ntemps = 1
nwalkers = 10
ndim = 9

sampler=PTSampler(ntemps, nwalkers, ndim, lnprobab, lnprior)
# Use the MPIPool instead
"""
pool = MPIPool()
if not pool.is_master():
    pool.wait()
    sys.exit(0)
"""
#sampler = PTSampler(ntemps, nwalkers, ndim, lnprobab, lnprior, threads=15)

#############################################################
# Initialize the walkers. The best technique seems to be
# to start in a small ball around the a priori preferred position.
# Dont worry, the walkers quickly branch out and explore the
# rest of the space.


w0 = np.random.uniform(-1.0,1.0,size=(ntemps,nwalkers))
w1 = np.random.uniform(150.0,300.0,size=(ntemps,nwalkers))
w2 = np.random.uniform(5.0,50.0,size=(ntemps,nwalkers))
w3 = np.random.uniform(-30.0,-2.0,size=(ntemps,nwalkers))
w4 = np.random.uniform(40.0,70.0,size=(ntemps,nwalkers))
w5 = np.random.uniform(-10.0,10.0,size=(ntemps,nwalkers))
w6 = np.random.uniform(-10.0,10.0,size=(ntemps,nwalkers))
w7 = np.random.uniform(350.0,360.0,size=(ntemps,nwalkers))
w8 = np.random.uniform(1.5e-7,1.5e-5,size=(ntemps,nwalkers))

p0 = np.dstack((w0,w1,w2,w3,w4,w5,w6,w7,w8))
#print p0.shape
niter = 2000
nburn = np.int(0.01*niter)

diskname='HD141569A'
#f = open("burnchain.dat", "w")
#f.close()

######################################################
##################  BURN IN STAGE ####################
######################################################



for p, lnprob, lnlike in sampler.sample(p0, iterations=nburn):
    #samples = sampler.chain[:,:,:].reshape((-1,ndim))
    #fig = triangle.corner(samples, labels=["$i$","$Ho$","$M$","$amax$","$alpha$","$beta$","$w$"])#,truths=[i_true,Ho_true,M_true,amax_true,alpha_true,beta_true,w_true])
    #fig.savefig("burntriangle.png")
    pass
 #   for k in range(nwalkers):
 #       for j in range(ntemps):
 #           f.write('{0} \n'.format(sampler.acceptance_fraction[j,k]))
    

#f.close()

    
print 'Burn in complete'

#pool.close()

samples = sampler.chain[:,:,:].reshape((-1,ndim))
fig = triangle.corner(samples, labels=parameters)
#,truths=[i_true,Ho_true,M_true,amax_true,alpha_true,beta_true,w_true])
fig.savefig("burntriangle.png")

sampler.reset()


#########################################################
################## RESTART FROM FILE ####################
#########################################################
"""
olddiskname = '/Users/swolff/EMCEE/emceetests/141202chaind/esoha569'
filename = olddiskname + '_inclinations.txt'
inclres = np.loadtxt(filename)
filename = olddiskname + '_scale_heights.txt'
hores = np.loadtxt(filename)
filename = olddiskname + '_dust_mass.txt'
massres = np.loadtxt(filename)
filename = olddiskname + '_amax.txt'
amaxres = np.loadtxt(filename)
filename = olddiskname + '_alpha.txt'
alphares = np.loadtxt(filename)
filename = olddiskname + '_beta.txt'
betares = np.loadtxt(filename)
filename = olddiskname + '_weights.txt'
weightres = np.loadtxt(filename)

nchain = 50
inclresd = inclres.reshape((2,50,nchain),order='C')
horesd = hores.reshape((2,50,nchain),order='C')
massresd = massres.reshape((2,50,nchain),order='C')
amaxresd = amaxres.reshape((2,50,nchain),order='C')
alpharesd = alphares.reshape((2,50,nchain),order='C')
betaresd = betares.reshape((2,50,nchain),order='C')
weightresd = weightres.reshape((2,50,nchain),order='C')


fullchain = np.zeros((2,50,nchain,7))
fullchain[:,:,:,0] = inclresd
fullchain[:,:,:,1] = horesd
fullchain[:,:,:,2] = massresd
fullchain[:,:,:,3] = amaxresd
fullchain[:,:,:,4] = alpharesd
fullchain[:,:,:,5] = betaresd
fullchain[:,:,:,6] = weightresd
p = fullchain[:,:,49,:]

#Parallelize after this step. Look into whether it's better to 
#vary number of walkers or number of steps across different machines.

#pos = open("chainpositions.dat", "w")
#probs = open("chainlnprobs.dat", "w")
#likes = open("chainlikes.dat", "w")

for p, lnprob, lnlike in sampler.sample(p, iterations=50,storechain=True):
#        for k in range(nwalkers):
#            for j in range(ntemps):
#                f.write('{0} \n'.format(sampler.acceptance_fraction[j,k]))

 #   pass
#f.close()

"""

gres = np.ndarray.flatten(sampler.chain[:,:,:,0])
radiusres = np.ndarray.flatten(sampler.chain[:,:,:,1])
innsloperes = np.ndarray.flatten(sampler.chain[:,:,:,2])
outsloperes = np.ndarray.flatten(sampler.chain[:,:,:,3])
inclinationres = np.ndarray.flatten(sampler.chain[:,:,:,4])
xoffsetres = np.ndarray.flatten(sampler.chain[:,:,:,5])
yoffsetres = np.ndarray.flatten(sampler.chain[:,:,:,6])
pares = np.ndarray.flatten(sampler.chain[:,:,:,7])
fluxres = np.ndarray.flatten(sampler.chain[:,:,:,8])

arrays = [gres, radiusres, innsloperes, outsloperes, inclinationres, xoffsetres, yoffsetres, pares, fluxres]

for i in arange(ndim):
    filename = diskname + '_'+parameters[i]+'.txt'
    np.savetxt(filename, arrays[i])
   
"""
    ChainStats = np.zeros((7,3))
    ha = np.ndarray.flatten(sampler.chain[:,:,:,0])
    ChainStats[0] = acor.acor(ha)
    ha = np.ndarray.flatten(sampler.chain[:,:,:,1])
    ChainStats[1] = acor.acor(ha)
    ha = np.ndarray.flatten(sampler.chain[:,:,:,2])
    ChainStats[2] = acor.acor(ha)
    ha = np.ndarray.flatten(sampler.chain[:,:,:,3])
    ChainStats[3] = acor.acor(ha)
    ha = np.ndarray.flatten(sampler.chain[:,:,:,4])
    ChainStats[4] = acor.acor(ha)
    ha = np.ndarray.flatten(sampler.chain[:,:,:,5])
    ChainStats[5] = acor.acor(ha)
    ha = np.ndarray.flatten(sampler.chain[:,:,:,6])
    ChainStats[6] = acor.acor(ha)
    filename = diskname + '_ChainStats.txt'
    np.savetxt(filename, ChainStats)

samples = sampler.chain[0,:,:,:].reshape((-1,ndim))
fig = triangle.corner(samples, labels=["$i$","$Ho$","$M$","$amax$","$alpha$","$beta$","$w$"])#,truths=[i_true,Ho_true,M_true,amax_true,alpha_true,beta_true,w_true])
fig.savefig("triangle.png")

"""

import pdb #@@@
pdb.set_trace() #@@@
print 'stop here' #@@@

