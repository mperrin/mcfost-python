import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits


def plot_seds(parfilename=None,dir="./", overplot=False, nlabels=None, alpha=0.75,
        inclination='all'):
    """
        Plot the SEDs for the current directory.

    This assumes regular seds from MC, rather than the RT calculated seds with the -rt1 option.

    Parameters
    -----------
    parfilename, dir: string
        the usual
    overplot : bool

    nlabels : int or None
        limit the number of inclination labels in the legend? set to None to show all (default)
    alpha : float
        Transparency for plots

    inclination : string 'all' or float
        'all' to plot all inclinations, else just the
        desired inclination.


    """
    parfilename = get_current_parfile(parfile=parfilename,dir=dir)
    par = readparfile(parfilename)

    if dir == './': title = parfilename
    else: title = dir

    if not overplot:
        plt.cla()


    # Try RT SED if it exists. If not, fall back to SED2 file instead
    if os.path.exists(dir+os.sep+'data_th/sed_rt.fits.gz'):
        RayTraceModeSED = True
        sed = pyfits.getdata(dir+os.sep+'data_th/sed_rt.fits.gz')
    else:
        RayTraceModeSED = False
        sed = pyfits.getdata(dir+os.sep+'data_th/sed2.fits.gz')

    phi=1

    lambd = par['lambda']


    # plot model SED for all inclinations
    ninc = sed.shape[2]
    labelstep = 1 if nlabels is None else ninc/nlabels

    if str(inclination)== 'all':
        for inc in range(ninc):
            if RayTraceModeSED:
                label = "%4.1f$^o$" % (par['im:inclinations'][inc])
            else:
                incmin=nplt.arccos(1.-((inc)*1.0)/ninc)*180./3.1415926
                incmax=nplt.arccos(1.-((inc+1  )*1.0)/ninc)*180./3.1415926
                label = "%4.1f - %4.1f$^o$" % (incmin,incmax)


            flux = sed[0,phi-1,inc,:]
            if nplt.mod(inc,labelstep) !=0: # allow skipping some labels if many are present
                label=None
            plt.loglog(lambd, flux, color=((ninc-inc)*1.0/ninc, 0, 0), label=label, alpha=alpha)
    else:
        wmin = nplt.argmin( abs(par['im:inclinations'] - inclination))
        print "Closest inclination found to %f is %f. " % (inclination, par['im:inclinations'][wmin])
        label = "%4.1f$^o$" % (par['im:inclinations'][wmin])
        flux = sed[0,phi-1,wmin,:]
        plt.loglog(lambd, flux, color=((ninc-wmin)*1.0/ninc, 0, 0), label=label, alpha=alpha)
           

    plt.xlabel("$\lambda$ ($\mu$m)")
    plt.ylabel("$\\nu F_\\nu$ (W m$^{-2}$)")
    plt.title("SED for "+title)
    plt.gca().xaxis.set_major_formatter(NicerLogFormatter())


    # observed SED. Try looking in a couple different places for this file.
    possible_sed_files = [ os.path.join(dir, "observed_sed.txt"), 
            os.path.join(dir, "data"+os.sep+"observed_sed.txt"), 
            os.path.join("data", "observed_sed.txt"),
            os.path.join("..", "data"+os.sep+"observed_sed.txt"),
            os.path.join("..", "..", "data"+os.sep+"observed_sed.txt")]
    sedfile=None
    for possible_name in possible_sed_files:
        if _VERBOSE: _log.info("Checking for observed SED at "+possible_name)
        if os.path.exists(possible_name):
            sedfile=possible_name
            if _VERBOSE: _log.info("Found observed SED at "+sedfile)
            break

    # if an observed SED is found, plot it.
    if sedfile is not None and os.path.exists(sedfile):
        #raise NotImplemented('Rewrite this with asciitable or atpy?')
        #seddata = nplt.asarray(asciidata.open(sedfile)[0:3])
        #obswavelen = seddata[0]
        #flux = seddata[1]*(1.e-26)*(3.e14/obswavelen) # convert from Jy to nu Fnu in W m ^-2
        #errflux=seddata[2]*(1.e-26)*(3.e14/obswavelen) #  convert from Jy to nu Fnu in W m ^-2

        observed = asciitable.read(sedfile, Reader=asciitable.Tab)
        flux = observed.Flux*(1.e-26)*(3.e14/observed.Wavelen) # convert from Jy to nu Fnu in W m ^-2
        errflux=observed.Uncert*(1.e-26)*(3.e14/observed.Wavelen) #  convert from Jy to nu Fnu in W m ^-2

        w_meas = nplt.where(flux > 0)
        w_upper = nplt.where(  (np.isnan(flux) | (flux == 0)) & np.isfinite(errflux) )

        plt.errorbar(observed.Wavelen[w_meas], flux[w_meas], yerr=errflux[w_meas], label="Observed", color='blue', fmt='o')
        plt.plot(observed.Wavelen[w_upper], errflux[w_upper], 'rv')
        #print(flux)

    plt.legend(prop = {'size':10})
    plt.draw()

def plot_lir_lstar(parfilename=None,dir="./", inclination=0):
    """ Estimate L_IR/L_star for a model.

    The MCFOST model must have been computed using separate components = T in the parameter
    file for this to work.
    """
    parfilename = get_current_parfile(parfile=parfilename,dir=dir)
    par = readparfile(parfilename)

    sed = pyfits.getdata('data_th/sed2.fits.gz')

    phi=1

    lambd = par['lambda']

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
        integrated_flux = nplt.trapz(flux_flambda, lambd)
        return integrated_flux

    ltot = integrate_channel(0)
    lstar = integrate_channel(4)
    ldisk = integrate_channel(6) +  integrate_channel(5)


    print "Integrated disk + star flux is %.2e W m^-2" % ltot
    print "Integrated star flux is %.2e W m^-2" % lstar
    print "Integrated disk flux is %.2e W m^-2" % ldisk

    print ""
    #print "L_IR / L_star  =  %.2e "   % (ldisk/lstar)
    print "L_IR / L_star  =  %.2e "   % ((ltot-lstar)/lstar)



def show_images(parfilename=None,dir="./", overplot=False, verbose=True, psf_fwhm=None, **kwargs):
    """
        Plot all available images for the current model

    keywords:


    """
    par = readparfile(get_current_parfile(parfile=parfilename,dir=dir))

    if not overplot: plt.cla()

    ims = glob.glob(dir+os.sep+"data_*/RT.fits.gz")

    wavelens =  [i[i.find('_')+1:i.find("/RT.fits.gz")] for i in ims]
    if verbose or _VERBOSE: _log.info("     Wavelengths found: "+str( wavelens))


    for w,ct in zip(wavelens, nplt.arange(len(wavelens))):
        plt.subplot(len(wavelens), len(wavelens), ct+1)
        image(w, par=par, dir=dir, **kwargs)


def show_image(wavelength, parfilename=None,par=None, dir="./", overplot=False, inclination=80, cmap=None, ax=None, 
        polarization=False, polfrac=False,
        psf_fwhm=None, 
        vmin=None, vmax=None):
    """ Show one image from an MCFOST model

    Parameters
    ----------
    wavelength : string
        note this is a string!! in microns. TBD make it take strings or floats
    inclination : float
        Which inclination to plot, in the case where there are multiple images? For RT computations
        with just one, then this parameter is essentially ignored
    overplot : bool
    cmap : matplotlib Colormap instance
    ax : matplotlib Axes instance
    vmin, vmax :    scalars
        for image display. Always shows in log stretch
    psf_fwhm : bool 
        convolve with PSF of this FWHM
    parfilename : string, optional
    par : dict, optional
    dir :  string


    """
    if not overplot:
        if ax is None: plt.clf()
        else: plt.cla()


    if ax is None: ax = plt.gca()
    if par is None:
        par = readparfile(get_current_parfile(parfile=parfilename,dir=dir))
    if cmap is None:
        cmap = plt.cm.gist_heat
        cmap = plt.cm.gist_gray
        cmaplt.set_over('white')
        cmaplt.set_under('black')
        cmaplt.set_bad('black')

    rt_im = pyfits.getdata(dir+os.sep+"data_"+wavelength+os.sep+"RT.fits.gz")
    inclin_index = _find_closest(par['im:inclinations'], inclination)
    #print "using image %s, inc=%f" % (  str(inclin_index), par['im:inclinations'][inclin_index]  )

    #ax = plt.subplot(151)
    image = rt_im[0,0,inclin_index,:,:]
    if vmax is None:
    #image.shape = image.shape[2:] # drop leading zero-length dimensions...
        vmax = image.max()
    if vmin is None:
        vmin=vmax/1e8
    norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax, clip=True)

    if psf_fwhm is not None:
        _log.info("Convolving with hard-coded Keck 1.6 um PSF - FIXME!!!!")
        import optics
        #psf = optics.airy_2d(aperture=10.0, wavelength=1.6e-6, pixelscale=0.020, shape=(99,99))
        psf = pyfits.getdata('/Users/mperrin/Dropbox/AAS2013/models_pds144n/manual/keck_psf_tmp4.fits')
        from astropy.nddata import convolve
        convolved = convolve(image, psf, normalize_kernel=True)

        image = convolved

    #ax.imshow(image,  norm=norm, cmap=cmap)
    if polfrac:
        ax = plt.subplot(121)

    ax = imshow_with_mouseover(ax, image,  norm=norm, cmap=cmap)
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
        polfrac_ar = nplt.sqrt(imQn**2 + imUn**2)
        theta = 0.5* nplt.arctan2(imUn, imQn) + np.pi/2
        vecX = polfrac_ar * nplt.cos(theta) *-1
        vecY = polfrac_ar * nplt.sin(theta)

        Y, X = nplt.indices(image.shape)
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



    plt.draw()

def plot_dust(dir='./', noerase=False):
    """ plot_dust

    Plot the dust scattering properties as calculated by MCFOST
    for some given collection of grains.

    Parameters
    -----------
    dir : string
        a MCFOST results directory
    noerase : bool
        don't erase current plots in active window.

    History:
        Python port ~2011 of IDL plot_dust.py circa 2009 by Marshall Perrin
    """


    if not os.path.isdir(dir):
       print "ERROR - invalid directory"

    if not os.path.isfile(dir + 'data_dust/kappa.fits.gz'):
       print("No dust properties files exist in that directory! need to compute some using MCFOST...")
       return False

    if noerase is False:  plt.clf()

    lambd = fits.getdata(dir + "data_dust/lambda.fits.gz")
    g = fits.getdata(dir + "data_dust/g.fits.gz")
    albedo = fits.getdata(dir + "data_dust/albedo.fits.gz")
    kappa = fits.getdata(dir + "data_dust/kappa.fits.gz")
    if os.path.exists(dir + "data_dust/polar.fits.gz"):
        has_polar= True
        polar = fits.getdata(dir + "data_dust/polar.fits.gz")
    else:
        has_polar = False

    plt.subplots_adjust(top=0.98, hspace=0.3)
    plt.subplot(411)
    plt.semilogx(lambd, kappa)
    plt.ylabel("Opacity $\kappa$")

    #----
    plt.subplot(412)
    plt.semilogx(lambd, albedo)
    plt.ylabel("Albedo")

    #----
    ax = plt.subplot(413)
    plt.semilogx(lambd, g)
    plt.ylabel("g= <cos $\\theta$>")
    plt.xlabel("Wavelength")
    ax.set_ybound([0.0, 1.0])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])


    #----- now switch to show the polarization versus angle ----

    ax = plt.subplot(414)
    pos = ax.get_position()  # move it down slightly for visual offset
    pa = pos.get_points()
    pa[0,1] *= 0.5
    pos.set_points(pa)
    ax.set_position(pos)

    polwaves = nplt.asarray([1.0, 2.5, 10.0, 100.0])
    theta = nplt.r_[0:180]
    plt.xlabel("Scattering Angle [degrees]")
    plt.ylabel("Polarization $[-S_{12}/S_{11}]$")
    plt.axis(xmin=0,xmax=180, ymin=0,ymax=1)
    plt.xticks(np.arange(7)*30)
    if has_polar:
        for i in nplt.arange(polwaves.size):
           c = _find_closest(lambd, polwaves[i])
           plt.plot(polar[:,c].ravel(), label=str(polwaves[i])+" $\mu$m" )

    prop = matplotlib.font_manager.FontProperties(size=10)
    plt.legend(prop=prop)

    #plt.draw()

