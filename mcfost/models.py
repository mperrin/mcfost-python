
# Model class - main class for interacting with an entire model

import warnings
import os
import glob
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import astropy.units as units
import logging
_log = logging.getLogger('mcfostpy')

#from . import plotting
from .paramfiles import Paramfile, find_paramfile

from . import utils


class ModelImageLoader(object):
    """ Helper class to implement on-demand loading of 
    image data via a dict interface. The dict looks like it
    always contains all available images, and can be indexed via
    either floats or strings interchangably. 

    This returns a PrimaryHDU object, not a HDUlist, under the
    assumption that MCFOST output files don't ever have more than one 
    extension.

    Most users will not need to instantiate this directly. Just
    use it transparently via the modelresults objectL

    Example
    ---------
    mod = mcfostpy.ModelResults()
    im_1_micron = mod.image_data[1.0]  # looks like array access, but is actually a file read

    """
    def __init__(self, modelresults):
        self._modelresults = modelresults

    def keys(self):
        return self._modelresults.image_wavelengths

    def __getitem__(self, key):

        return fits.open( self._getpath(key))[0]

    def _getpath(self, key):
        canonical_wavelength = self._modelresults._standardize_wavelength(key)
        return os.path.join(self._modelresults.directory,"data_"+canonical_wavelength,"RT.fits.gz")

    @property   
    def filenames(self):
        return [self._getpath(k) for k in self.keys()]



class ModelResults(object):
    """ One MCFOST model result set, possibly containing multiple inclinations 
    
    This class provides an object oriented interface to MCFOST model results, including SEDs 
    and images.  It also casts results so that the astropy.units framework can be used
    
    For now this class assumes all your model results have been precomputed. 
    I.e. it will not call mcfost itself to run anything new. 
    """
    def __init__(self, directory="./"):
        """ Initialize model. 
        
        This does read in parameters, but does NOT automatically load
        SED or image information. Those will be loaded automatically into memory when needed. 
        This optimization is to aid in loading large model sets at once.

        Parameters
        -----------
        directory : str
            Directory path to the root dir of the model (i.e. the parent directoy containing 
            data_th, data_1.0 etc)
        """

        self.directory = directory
        self._paramfilename = find_paramfile(self.directory)

        if self._paramfilename is None:
            raise IOError("Could not find a parameter file in that directory")
        self.parameters = Paramfile(self._paramfilename)


        # check for existence of SEDs but don't load yet?
        if not os.path.exists(os.path.join(self.directory, 'data_th', 'sed_rt.fits.gz')):
            raise IOError("There does not appear to be an SED model output (data_th/sed_rt.fits.gz) in the model directory {0}".format(self.directory))
        self._sed_data = None # see the property sed_data for the actual loading below

        # check for available images, 
        # and keep a list of available wavelengths, both as a s
        image_waves = glob.glob(os.path.join(self.directory, 'data_[0123456789]*'))
        if image_waves is None:
            warnings.warn("No model image results were found in that directory")
        self._wavelengths_lookup = {}
        for wl in image_waves:
            #print wl
            wlstring = os.path.basename(wl).split('_')[1]
            self._wavelengths_lookup[float(wlstring)] = wlstring

        self.image_data = ModelImageLoader(self)

    @property
    def image_wavelengths(self):
        #return self._wavelengths_lookup.values()
        return np.asarray(self._wavelengths_lookup.values(), dtype=float) * units.micron



    def __repr__(self):
        return "<MCFOST ModelResults in directory '{self.directory}'>".format(self=self)
   
    @property 
    def sed_data(self):
        """ Attribute that provides access to SED results data (from sed_rt.fits.gz files)

        Implement lazy loading for SED data - no files are read until the 
        user tries to access this attribute, at which point the SED data is
        automaatically loaded """
        # if it's already loaded then just return it
        if self._sed_data is not None: 
            return self._sed_data * (units.W / units.m**2)
        else:

            if os.path.exists( os.path.join(self.directory,'data_th/sed_rt.fits.gz')):
                self._RayTraceModeSED = True
                self._sed_data = fits.getdata( os.path.join(self.directory, 'data_th/sed_rt.fits.gz'))
                _log.debug("loading SED data from RT mode SED")
                return self._sed_data * (units.W / units.m**2)
            elif os.path.exists( os.path.join(self.directory,'data_th/sed2.fits.gz')):
                self._RayTraceModeSED = False
                self._sed_data = fits.getdata(os.path.join(self.directory,'data_th/sed2.fits.gz'))
                _log.debug("loading SED data from MC mode SED") 
                return self._sed_data * (units.W / units.m**2)
            else:
                _log.error("No SED data present in "+self.directory+"!")



    @property
    def sed_wavelengths(self):
        """ Wavelengths used for SED calculation"""
        return self.parameters.wavelengths * units.micron



    def _standardize_wavelength(self, wavelength):
        """ Utility function to return a "standardized" representation
        of a wavelength, such that e.g. "2" microns, "2.0" microns, 
        "2.00e0" microns etc are all the same thing, and that thing is
        consistent with the directory names on disk produced by MCFOST"""
        wl = float(wavelength)
        try:
            return self._wavelengths_lookup[wl]
        except:
            raise ValueError("This model does not have images for wavelength="+str(wavelength)+" microns. Please check your data and try again.")

    def plot_SED(self, inclination='all', title=None, overplot=False, 
            nlabels=None, alpha=0.75, **kwargs):
        """ Plot SED(s). Set inclination=value to plot just one, or leave blank for all
    
        This will plot either an RT or MC mode SED. By default RT is tried first.

        Parameters
        -----------
        overplot : bool
            Should this overplot or make a new plot?
        nlabels : int or None
            limit the number of inclination labels in the legend? set to None to show all (default)
        alpha : float
            Transparency for plots
        inclination : string 'all' or float
            'all' to plot all inclinations, else just the
            desired inclination (specifically the SED whose 
            inclination is closest to the specified floating point value)


        """


        if title is None:
            title=self.directory

        dummy = self.sed_data # force auto loading of SED data if that hasn't already happened...

        if not overplot: plt.cla()
        # plot model SED for all inclinations
        ninc = len(self.parameters.inclinations) # or self.parameters.im_rt_ninc  should be equivalent here...
        phi = 1

        if str(inclination)== 'all':
            labelstep = 1 if nlabels is None else ninc/nlabels
            for inc in range(ninc):
                if self._RayTraceModeSED:
                    label = "%4.1f$^o$" % (self.parameters.inclinations[inc])
                else:
                    incmin=np.arccos(1.-((inc)*1.0)/ninc)*180./3.1415926
                    incmax=np.arccos(1.-((inc+1  )*1.0)/ninc)*180./3.1415926
                    label = "%4.1f - %4.1f$^o$" % (incmin,incmax)


                flux = self.sed_data[0,phi-1,inc,:]
                if np.mod(inc,labelstep) !=0: # allow skipping some labels if many are present
                    label=None
                plt.loglog(self.parameters.wavelengths, flux, color=((ninc-inc)*1.0/ninc, 0, 0), label=label, alpha=alpha)
        else:
            wmin = np.argmin( abs(self.parameters.im_inclinations) - inclination)
            print "Closest inclination found to %f is %f. " % (inclination, self.parameters.im_inclinations[wmin])
            label = "%4.1f$^o$" % (self.parameters.im_inclinations[wmin])
            flux = self.sed_data[0,phi-1,wmin,:]
            plt.loglog(self.parameters.wavelengths, flux, color=((ninc-wmin)*1.0/ninc, 0, 0), label=label, alpha=alpha)

        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("$\\nu F_\\nu$ (W m$^{-2}$)")
        plt.title("SED for "+title)
        plt.gca().xaxis.set_major_formatter(utils.NicerLogFormatter())


    def plot_image(self, wavelength0, overplot=False, inclination=None, cmap=None, ax=None, 
            axes_units='AU',
            polarization=False, polfrac=False,
            psf_fwhm=None, 
            vmin=None, vmax=None, dynamic_range=1e6):
        """ Show one image from an MCFOST model

        Parameters
        ----------
        wavelength : string
            note this is a string!! in microns. TBD make it take strings or floats
        inclination : float
            Which inclination to plot, in the case where there are multiple images?
            If not specified, defaults to the median inclination of whichever RT mode
            inclinations were computed.
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
        axes_units : str
            Units to label the axes in. May be 'AU', 'arcsec', 'deg', or 'pixels'
        psf_fwhm : float or None 
            convolve with PSF of this FWHM? Default is None for no convolution

        

        To Do:
        * Add convolution with PSFs
        * Add coronagraphic occulters

        """

        wavelength = self._standardize_wavelength(wavelength0)
        if inclination is None: 
            inclination = self.parameters.inclinations[len(self.parameters.inclinations)/2]

        if not overplot:
            if ax is None: plt.clf()
            else: plt.cla()


        if ax is None: ax = plt.gca()

        # Read in the data and select the image of 

        imHDU = self.image_data[wavelength]


        inclin_index = utils.find_closest(self.parameters.inclinations, inclination)
        image = imHDU.data[0,0,inclin_index,:,:]


        # Set up color mapping
        if cmap is None:
            from copy import deepcopy
            #cmap = plt.cm.gist_heat
            cmap = deepcopy(plt.cm.gist_gray)
            cmap.set_over('white')
            cmap.set_under('black')
            cmap.set_bad('black')
        if vmax is None: vmax = image.max()
        if vmin is None: vmin=vmax/dynamic_range
        norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax, clip=True)


        # calculate the display extent of the image

        # check if the image header has CDELT keywords, if so use those in preference to
        # the plate scale specified in the parameter file.
        # This will gracefully handle the case of the user overriding the default pixel
        # scale on the command line to MCFOST.
        cd1 = imHDU.header['CDELT1']
        cd2 = imHDU.header['CDELT2']
        pixelscale = np.asarray([cd1,cd2])  # in degrees

        if axes_units.lower() == 'deg':
            pass
        elif axes_units.lower() == 'arcsec':
            pixelscale *= 3600
        elif axes_units.lower() == 'au':
            pixelscale *= 3600 * self.parameters.distance
        elif axes_units.lower() == 'pixels' or axes_units.lower() == 'pixel':
            pixelscale = np.ones(2)
        else:
            raise ValueError("Unknown unit for axes_units: "+axes_units)

        halfsize = (np.asarray(image.shape))/2 * pixelscale
        extent = [-halfsize[0], halfsize[0], -halfsize[1], halfsize[1]]

        
        



        # Optional: convolve with PSF
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

        # Display the image!
        ax = utils.imshow_with_mouseover(ax, image,  norm=norm, cmap=cmap, extent=extent)
        ax.set_xlabel("Offset [{unit}]".format(unit=axes_units))
        ax.set_ylabel("Offset [{unit}]".format(unit=axes_units))
        ax.set_title("Image for "+wavelength+" $\mu$m")


        if polarization==True:
            # overplot polarization vectors

            # don't show every pixel's vectors:
            showevery = 5

            imQ = imHDU.data[1,0,inclin_index,:,:] * -1  #MCFOST sign convention requires flip for +Q = up
            imU = imHDU.data[2,0,inclin_index,:,:] * -1  #MCFOST sign convention, ditto

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
                # Display a 2nd panel showing the polarization fraction
                ax2 = plt.subplot(122)

                cmap2 = plt.cm.gist_heat
                cmap2 = plt.cm.jet
                cmap2.set_over('red')
                cmap2.set_under('black')
                cmap2.set_bad('black')

                ax2 = utils.imshow_with_mouseover(ax2, polfrac_ar, vmin=0, vmax=1,  cmap=cmap2)
                ax2.set_title("Polarization fraction")

                cax = plt.axes([0.92, 0.25, 0.02, 0.5])
                plt.colorbar(ax2.images[0], cax=cax, orientation='vertical')

    def describe(self):
        """ Return a descriptive brief paragraph on what results are present """
        print "Model results for {self._paramfilename}".format(self=self)
        if self.sed_data is None:
            print "    No SED data present"
        else:
            print "    SED computed from {par.wavelengths_min} - {par.wavelengths_max} microns using {par.nwavelengths} wavelengths".format(par=self.parameters)
        print "    Images computed for {0} wavelengths: {1}".format(len(self.image_wavelengths), self.image_wavelengths)



class Observations(object):
    """ Class for observational data on a target 
    
    To the extent possible, the API is consistent with that of the ModelResults class.
    That is, the names and functionality of attributes and methods should be similar
    to the extent possible given the differences in the underlying data.

    """

    sed = None
    """ Observed SED"""

    images = None
    """ Container of images?  TBD what type of container. Indexed by wavelength?  """

    image_wavelengths = []
    """ List of wavelengths for which we have images """


    def __init__(self, directory=None):
        """ Initialize based on a directory containing several results files.
        For now this follows the conventions used by the IDL results explorer code, but
        this is likely to change as the code evolves...
        """

        if directory is not None:
            if not os.path.isdir(directory): raise ValueError("Not a valid directory: "+directory)
            self.directory = directory
            summaryfile = os.path.join(self.directory, 'summary.txt')
            if not os.path.exists(summaryfile):
                raise ValueError("Cannot find index file of observed data: summary.txt")

            summary = open(summaryfile).readlines()

            #self.filenames = np.zeros(len(summary), dtype=str)
            #self.types= np.zeros(len(summary), dtype=str)
            #self.wavelengths= np.zeros(len(summary), dtype=str)
            filenames = []
            types = []
            wavelengths = []

            for i, line in enumerate(summary):
                # filename, type, wavelength (optional)
                parts = line.split()
                print i, parts
                filenames.append(parts[0])
                types.append(parts[1].lower())
                wavelengths.append(parts[2] if len(parts) >= 3 else None)
                if parts[1].lower() == 'image': self.image_wavelengths.append(parts[2])
                elif parts[1].lower() == 'sed': self.sed = ObservedSED(parts[0])
            self.filenames = np.asarray(filenames)
            self.types = np.asarray(types)
            self.wavelengths = np.asarray(wavelengths)

     
            
    
        


class ObservedSED(object):
    """ Observed SED class. 
    Reads observations from disk; returns them as as astropy Units objects 
    """

    def __init__(self, filename=None, uncertainty=None, mask=None, format='ascii.no_header'):

        # temporary hard coded default for development
        if filename is None: 
            filename = '/Users/mperrin/data/mcfost/models_esoha569/data/observed_sed.txt'


        self._sed_table = ascii.read(filename, format=format, 
                names=['wavelength','flux','uncertainty','source'])

    @property
    def flux(self):
        return self._sed_table['flux']  * (units.W / units.m**2)

    @property
    def wavelength(self):
        return self._sed_table['wavelength'] *units.micron

    @property
    def uncertainty(self):
        return self._sed_table['uncertainty'] * (units.W / units.m**2)


    def plot(self, title=None, overplot=False, 
             alpha=0.75, **kwargs):
        """ Plot observed SED
    
        Parameters
        -----------
        overplot : bool
            Should this overplot or make a new plot?
        alpha : float
            Transparency for plots

        """

        if title is None:
            title='Observed SED from '+self.directory
        label = 'Observed SED'

        dummy = self.sed_data # force auto loading of SED data if that hasn't already happened...

        if not overplot: plt.cla()
        flux = self.sed_data
        plt.loglog(self.sed_wavelengths, flux, label=label, alpha=alpha)

        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("$\\nu F_\\nu$ (W m$^{-2}$)")
        plt.title(title)
        plt.gca().xaxis.set_major_formatter(utils.NicerLogFormatter())



class ObservedImage(object):
    """ Container class for an observed image, 
    optionally along with associated uncertainty image and pixel mask.

    """
    def __init__(self, filename, uncertainty=None, mask=None):
        pass



    def plot(self, wavelength, overplot=False, cmap=None, ax=None, 
            vmin=None, vmax=None, which='image'):
        """
        Display one image

        Parameters
        ----------
        wavelength : float
           desired wavelength to plot into 
        overplot : bool:q


        cmap : matplotlib color map
            desired comor map
        ax : matplotlib axis
            axis to plot into
        vmin, vmax : floats
            min and max for image display range. 

        """
        wm = np.where( (self.types == which) and (self.wavelengths == wavelength0) )

        imagefilename = self.filenames[wm]
        print "Filename for image: "+imagefilename
        raise NotImplementedError("Not implemented yet!")


