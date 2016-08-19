import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import astropy.io.fits as fits
import astropy.io.ascii as ioascii
import os
import logging
_log = logging.getLogger('mcfost')

from . import paramfiles
from . import utils
from . import models


def plot_seds(dir="./",
        observed_sed_filename=None, **kwargs):
    """
        Plot the SEDs for the current directory.

    This is just a convenience wrapper for the plotting code in ModelResults.

    Parameters
    -----------
    parfilename, dir: strings
        the usual
    overplot : bool
        Should this overplot or make a new plot?
    observed_sed_filename : string or None
        Filename for observed SED to overplot. If None, will automatically check
        several plausible filenames + paths.
    inclination : string 'all' or float
        'all' to plot all inclinations, else just the
        desired inclination (specifically the SED whose
        inclination is closest to the specified floating point value)

    """

    mr = models.ModelResults(dir)
    mr.sed.plot( **kwargs)



def plot_lir_lstar(parfilename=None,dir="./", inclination=0):
    """ Estimate L_IR/L_star for a model.

    The MCFOST model must have been computed using separate components = T in the parameter
    file for this to work.
    """
    par = paramfiles.Paramfile(parfilename,dir)

    # Try RT SED if it exists. If not, fall back to SED2 file instead
    if os.path.exists( os.path.join(dir,'data_th/sed_rt.fits.gz')):
        RayTraceModeSED = True
        sed = fits.getdata( os.path.join(dir, 'data_th/sed_rt.fits.gz'))
    elif os.path.exists( os.path.join(dir,'data_th/sed2.fits.gz')):
        RayTraceModeSED = False
        sed = fits.getdata(os.path.join(dir,'data_th/sed2.fits.gz'))
    elif ~ os.path.isdir(os.path.join(dir,'data_th/')):
        raise IOError('Cannot find a data_th subdirectory inside '+dir)
    else:
        raise IOError('Cannot find either a sed_rt or sed2 file in data_th inside '+dir)

    lambd = par.wavelengths
    phi=1

    try:
        flux_nufnu = sed[6,phi-1,inclination,:]   #  see if the file has components separated
    except:
        raise IOError("The supplied MCFOST SED model was not computed with separate components. Please recompute to enable distinguishing between star and disk flux.")


    plt.clf()
    for i in range(9):
        ls = "--" if i > 4 else "-"
        flux_nufnu = sed[i,phi-1,inclination,:]   # this is in nu F_nu units.
        plt.loglog(lambd, flux_nufnu, label="channel %d" % i, linestyle=ls)
    plt.legend()

    # channels appear to be:
    # 0 total
    # 1 ??
    # 2 ??
    # 3 ??
    # 4 star
    # 5,6 disk

    plt.draw()



    def integrate_channel(i):
        flux_nufnu = sed[i,phi-1,inclination,:]   # this is in nu F_nu units.
        flux_flambda = flux_nufnu/lambd
        integrated_flux = np.trapz(flux_flambda, lambd)
        return integrated_flux

    ltot = integrate_channel(0)
    lstar = integrate_channel(4)
    ldisk = integrate_channel(6) +  integrate_channel(5)


    print("Integrated disk + star flux is %.2e W m^-2" % ltot)
    print("Integrated star flux is %.2e W m^-2" % lstar)
    print("Integrated disk flux is %.2e W m^-2" % ldisk)

    print("")
    #print "L_IR / L_star  =  %.2e "   % (ldisk/lstar)
    print("L_IR / L_star  =  %.2e "   % ((ltot-lstar)/lstar))



def plot_images(parfilename=None,dir="./", overplot=False, psf_fwhm=None, **kwargs):
    """
        Plot all available images for the current model

    keywords:


    """
    import glob
    par = paramfiles.Paramfile(parfilename,dir)

    if not overplot: plt.cla()

    ims = glob.glob(os.path.join(dir,"data_*/RT.fits.gz"))

    wavelens =  [i[i.find('_')+1:i.find("/RT.fits.gz")] for i in ims]
    _log.debug("     Wavelengths found: "+str( wavelens))


    for w,ct in zip(wavelens, np.arange(len(wavelens))):
        plt.subplot(len(wavelens), len(wavelens), ct+1)
        plot_image(w, par=par, dir=dir, **kwargs)


def plot_image(wavelength, parfilename=None,par=None, dir="./", overplot=False, inclination=80, cmap=None, ax=None,
        polarization=False, polfrac=False,
        psf_fwhm=None,
        vmin=None, vmax=None, dynamic_range=1e6):
    """ Show one image from an MCFOST model

    Parameters
    ----------
    wavelength : string
        note this is a string!! in microns. TBD make it take strings or floats
    inclination : float
        Which inclination to plot, in the case where there are multiple images? For RT computations
        with just one, then this parameter is essentially ignored
    overplot : bool
        Should we overplot on current axes (True) or display a new figure? Default is False
    cmap : matplotlib Colormap instance
        Color map to use for display.
    ax : matplotlib Axes instance
        Axis to display into.
    vmin, vmax :    scalars
        Min and max values for image display. Always shows in log stretch
    dynamic_range : float
        default vmin is vmax/dynamic range. Default dynamic range is 1e6.
    psf_fwhm : float or None
        convolve with PSF of this FWHM? Default is None for no convolution
    parfilename : string, optional
    par : dict, optional
    dir :  string

    """
    if not overplot:
        if ax is None: plt.clf()
        else: plt.cla()


    if ax is None: ax = plt.gca()
    if par is None:
        par = paramfiles.Paramfile(parfilename,dir)
    if cmap is None:
        cmap = plt.cm.gist_heat
        cmap = plt.cm.gist_gray
        cmap.set_over('white')
        cmap.set_under('black')
        cmap.set_bad('black')

    rt_im = fits.getdata(os.path.join(dir,"data_"+wavelength,"RT.fits.gz"))
    inclin_index = utils.find_closest(par.inclinations, inclination)
    #print "using image %s, inc=%f" % (  str(inclin_index), par['im:inclinations'][inclin_index]  )

    #ax = plt.subplot(151)
    image = rt_im[0,0,inclin_index,:,:]
    if vmax is None:
    #image.shape = image.shape[2:] # drop leading zero-length dimensions...
        vmax = image.max()
    if vmin is None:
        vmin=vmax/dynamic_range
    norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax, clip=True)

    if psf_fwhm is not None:
        raise NotImplementedError("This code not yet implemented")
        import optics
        #psf = optics.airy_2d(aperture=10.0, wavelength=1.6e-6, pixelscale=0.020, shape=(99,99))
        psf = pyfits.getdata('/Users/mperrin/Dropbox/AAS2013/models_pds144n/manual/keck_psf_tmp4.fits')
        from astropy.nddata import convolve
        convolved = convolve(image, psf, normalize_kernel=True)

        image = convolved

    if polfrac:
        ax = plt.subplot(121)

    ax = utils.imshow_with_mouseover(ax, image,  norm=norm, cmap=cmap)
    ax.set_xlabel("Offset [pix]")
    ax.set_ylabel("Offset [pix]")
    ax.set_title(wavelength+" $\mu$m")


    if polarization==True:

        # dont' show every pixel's vectors:
        showevery = 5

        imQ = rt_im[1,0,inclin_index,:,:] * -1  #MCFOST sign convention requires flip for +Q = up
        imU = rt_im[2,0,inclin_index,:,:] * -1  #MCFOST sign convention

        imQn = imQ/image
        imUn = imU/image
        polfrac_ar = np.sqrt(imQn**2 + imUn**2)
        theta = 0.5* np.arctan2(imUn, imQn) + np.pi/2
        vecX = polfrac_ar * np.cos(theta) *-1
        vecY = polfrac_ar * np.sin(theta)

        Y, X = np.indices(image.shape)
        Q = plt.quiver(X[::showevery, ::showevery], Y[::showevery, ::showevery], vecX[::showevery, ::showevery], vecY[::showevery, ::showevery],
                headwidth=0, headlength=0, headaxislength=0.0, pivot='middle', color='white')
        plt.quiverkey(Q, 0.85, 0.95, 0.5, "50% pol.", coordinates='axes', labelcolor='white', labelpos='E', fontproperties={'size':'small'}) #, width=0.003)
#        ax = plt.subplot(152)
#        ax.imshow(imQn, vmin=-1, vmax=1, cmap=matplotlib.cm.RdBu)
#        ax.set_title('Q')
#        ax = plt.subplot(153)
#        ax.imshow(imUn, vmin=-1, vmax=1, cmap=matplotlib.cm.RdBu)
#        ax.set_title('U')
#        ax = plt.subplot(154)
#        ax.imshow(vecX, vmin=-1, vmax=1, cmap=matplotlib.cm.RdBu)
#        ax.set_title('vecX')
#        ax = plt.subplot(155)
#        ax.imshow(vecY, vmin=-1, vmax=1, cmap=matplotlib.cm.RdBu)
#        ax.set_title('vecY')
#        ax.colorbar()

    #plt.colorbar()

        if polfrac:
            ax2 = plt.subplot(122)

            cmap2 = plt.cm.gist_heat
            cmap2 = plt.cm.jet
            cmap2.set_over('red')
            cmap2.set_under('black')
            cmap2.set_bad('black')

            ax2 = imshow_with_mouseover(ax2, polfrac_ar, vmin=0, vmax=1,  cmap=cmap2)
            ax2.set_title("Polarization fraction")

            cax = plt.axes([0.92, 0.25, 0.02, 0.5])
            plt.colorbar(ax2.images[0], cax=cax, orientation='vertical')




def plot_dust(directory='./', noerase=False, parameters=None):
    """ plot_dust

    Plot the dust scattering properties as calculated by MCFOST
    for some given collection of grains.

    Parameters
    -----------
    directory : string
        a MCFOST results directory
    noerase : bool
        don't erase current plots in active window.

    History:
        Python port ~2011 of IDL plot_dust.py circa 2009 by Marshall Perrin
    """


    #if not os.path.isdir(dir):
       #print "ERROR - invalid directory"


    if not os.path.isfile(os.path.join(directory,'kappa.fits.gz')):
        # if no dust files in that directory, check also a data_dust subdirectory
        directory = os.path.join(directory, 'data_dust')
        # and then try again:
        if not os.path.isfile(os.path.join(directory,'kappa.fits.gz')):
            # and if we still can't find any then give up
            raise IOError("No dust properties files exist in that directory! You need to compute some using MCFOST...")

    if parameters is None:
        parameters = paramfiles.find_paramfile(directory=directory)

    lambd = fits.getdata(os.path.join(directory, "lambda.fits.gz"))
    g = fits.getdata(os.path.join(directory, "g.fits.gz"))
    albedo = fits.getdata(os.path.join(directory, "albedo.fits.gz"))
    kappa = fits.getdata(os.path.join(directory, "kappa.fits.gz"))
    if os.path.exists(os.path.join(directory, "polar.fits.gz")):
        has_polar= True
        polar = fits.getdata(os.path.join(directory, "polar.fits.gz"))
    elif os.path.exists(os.path.join(directory, "polarizability.fits.gz")):
        has_polar= True
        polar = fits.getdata(os.path.join(directory, "polarizability.fits.gz"))
    else:
        has_polar = False

    if os.path.exists(os.path.join(directory, "phase_function.fits.gz")):
        has_phase= True
        phase = fits.getdata(os.path.join(directory, "phase_function.fits.gz"))
    else:
        has_phase = False



    if noerase is False:  plt.clf()
    plt.subplots_adjust(top=0.98, hspace=0.3, wspace=0.4)
    ax = plt.subplot(321)
    plt.semilogx(lambd, kappa)
    plt.ylabel("Opacity $\kappa$ [cm$^2$/g]")
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(5) )

    #----
    ax = plt.subplot(323)
    plt.semilogx(lambd, albedo)
    plt.ylabel("Albedo")
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(5) )

    #----
    ax = plt.subplot(325)
    plt.semilogx(lambd, g)
    plt.ylabel("g= <cos $\\theta$>")
    plt.xlabel("Wavelength [$\mu$m]")
    ax.set_ybound([0.0, 1.0])
    ax.yaxis.set_major_locator( matplotlib.ticker.MaxNLocator(5) )
    #ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])

    #----- now switch to show the phase function versus angle ----

    # these variables are used for both phase function and polarization
    pwaves = np.asarray([1.0, 2.5, 10.0, 100.0])
    colors = ['blue', 'green', 'orange', 'red']
    theta = np.r_[0:180]


    ax = plt.subplot(222)
    plt.xlabel("Scattering Angle [degrees]")
    plt.ylabel("Phase Function")
    plt.xticks(np.arange(7)*30)

    if has_phase:
        for i in np.arange(pwaves.size):
           c = utils.find_closest(lambd, pwaves[i])
           plt.semilogy(phase[:,c].ravel(), label=str(pwaves[i])+" $\mu$m", color=colors[i] )

        prop = matplotlib.font_manager.FontProperties(size=10)
        plt.legend(prop=prop)
    else:
        plt.text(90,0.5, 'Phase function information\n not present', horizontalalignment='center', verticalalignment='center')
    #----- now switch to show the polarization versus angle ----

    ax = plt.subplot(224)
    # move this 4th plot down slightly for visual offset
    #pos = ax.get_position()
    #pa = pos.get_points()
    #pa[0,1] *= 0.5
    #pos.set_points(pa)
    #ax.set_position(pos)

    theta = np.r_[0:180]
    plt.xlabel("Scattering Angle [degrees]")
    plt.ylabel("Polarization $[-S_{12}/S_{11}]$")
    ymin = -1 if polar.min() < 0 else 1

    plt.axis(xmin=0,xmax=180, ymin=ymin,ymax=1)
    plt.xticks(np.arange(7)*30)
    if has_polar:
        for i in np.arange(pwaves.size):
           c = utils.find_closest(lambd, pwaves[i])
           plt.plot(polar[:,c].ravel(), label=str(pwaves[i])+" $\mu$m" , color=colors[i])

        prop = matplotlib.font_manager.FontProperties(size=10)
        plt.legend(prop=prop)
    else:
        plt.text(90,0.5, 'Polarization information\n not present', horizontalalignment='center', verticalalignment='center')
