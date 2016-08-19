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
import collections
import logging
_log = logging.getLogger('mcfost')

#from . import plotting
from .paramfiles import Paramfile, find_paramfile

from . import utils


class ModelImageCollection(collections.Mapping):

    """ 
    Helper collection class to implement on-demand loading of
    image data via a dict interface. The dict looks like it
    always contains all available images, and can be indexed via
    either floats or strings interchangably.

    This returns a PrimaryHDU object, not a HDUlist, under the
    assumption that MCFOST output files don't ever have more than one
    extension.

    Most users will not need to instantiate this directly. Just
    use it transparently via a ModelResults object

    Example
    ---------
    mod = mcfostpy.ModelResults()
    im_1_micron = mod.image_data[1.0]  # looks like array access, but is
                    # actually an automatic file read first, then access

    """
    def __init__(self, modelresults):
        self._modelresults = modelresults
        self._loaded_fits = dict()

    def _getpath(self, key):
        canonical_wavelength = self._modelresults._standardize_wavelength(key)
        return os.path.join(self._modelresults.directory,"data_"+canonical_wavelength,"RT.fits.gz")

    @property
    def wavelengths(self):

        return np.asarray(list(self._modelresults._wavelengths_lookup.values()), dtype=float) * units.micron
        #return self._modelresults.image_wavelengths

    @property
    def filenames(self):
        return [self._getpath(k) for k in self.keys()]

    # Implement magic methods to make this compliant with the collections.Mapping abstract base class
    # that way it will behave just like a python Dict
    def keys(self):
        # cast these as a list
        return self._modelresults._wavelengths_lookup.values()
        #return list(self._modelresults.image_wavelengths.value)

    def closeimage(self):
        for k in self.keys():
            canonical_wavelength = self._modelresults._standardize_wavelength(k)
            try:
                self.loaded_fits[canonical_wavelength].close()
            except:
                print( "Failed to close model image files.")


    def __len__(self):
        return len(self.keys())

    def __iter__(self):
        for i in self.keys():
            yield self[i]

    def __getitem__(self, key):
        canonical_wavelength = self._modelresults._standardize_wavelength(key)

        if canonical_wavelength not in self._loaded_fits.keys():
            self._loaded_fits[canonical_wavelength] = fits.open( self._getpath(canonical_wavelength))[0]
        
        return self._loaded_fits[canonical_wavelength]


class MCFOST_Dataset(object):
    """ Base class - can either be model results or
    observations. Implements common behavior to each.
    """
    def __init__(self, directory="./"):
        pass

    @property
    def image_wavelengths(self):
        """ Wavelengths for which images are present, in microns """
        return np.asarray(list(self._wavelengths_lookup.values()), dtype=float) * units.micron

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

#----------------------------------------------------

class ModelResults(MCFOST_Dataset):
    """ One MCFOST model result set, possibly containing multiple inclinations

    This class provides an object oriented interface to MCFOST model results, including SEDs
    and images.  It also casts results so that the astropy.units framework can be used

    For now this class assumes all your model results have been precomputed.
    I.e. it will not call mcfost itself to run anything new.

    Because of how MCFOST works, in practice this class may actually represent multiple
    inclinations all calculated for the same set of physical parameters.
    """

    parameters = None
    """ Paramfile object for this model's parameters"""

    sed = None
    """ ModelSED object containing the computed SED """

    def __init__(self, directory="./", parameters=None):
        """ Initialize model.

        This does read in parameters, but does NOT automatically load
        SED or image information. Those will be loaded automatically into memory when needed.
        This optimization is to aid in loading large model sets at once.

        Parameters
        -----------
        directory : str
            Directory path to the root dir of the model (i.e. the parent directoy containing
            data_th, data_1.0 etc)
        parameters : filename
            Filename of MCFOST parameter file. Can be omitted if there is one unambiguous
            parameter file in the specified directory.
        """

        self.directory = os.path.abspath(directory)

        if parameters is None:
            self._paramfilename = find_paramfile(self.directory)
        else:
            self._paramfilename = parameters

        if self._paramfilename is None:
            raise IOError("Could not find a parameter file in that directory")
        self.parameters = Paramfile(self._paramfilename)


        # check for available images,
        # and keep a list of available wavelengths, both as a s
        image_waves = glob.glob(os.path.join(self.directory, 'data_[0123456789]*'))
        if image_waves is None:
            warnings.warn("No model image results were found in that directory.")
        self._wavelengths_lookup = {}
        for wl in image_waves:
            wlstring = os.path.basename(wl).split('_')[1]
            self._wavelengths_lookup[float(wlstring)] = wlstring

        self.images = ModelImageCollection(self)


        # Check for SED
        try:
            self.sed = ModelSED(directory=os.path.join(self.directory, 'data_th') )
        except:
            self.sed = None
            warnings.warn('No ray-traced SED results were found in that directory.')

    def __repr__(self):
        return "<MCFOST ModelResults in directory '{self.directory}'>".format(self=self)

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

        raise DeprecationWarning("ModelResults.plot_sed is deprecated; use Modelresults.sed.plot() instead.")


        self.sed.plot(inclination=inclination, title=title, overplot=overplot,
                nlabels=nlabels, alpha=alpha, **kwargs)

        return


    def plot_image(self, wavelength0, overplot=False, inclination=None, cmap=None, ax=None,
            axes_units='AU',
            polarization=False, polfrac=False,
            psf_fwhm=None,
            vmin=None, vmax=None, dynamic_range=1e6, colorbar=False):
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
        colorbar : bool
            Also draw a colorbar


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

        imHDU = self.images[wavelength]


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
        ax = utils.imshow_with_mouseover(ax, image,  norm=norm, cmap=cmap, extent=extent, origin='lower')
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

        if colorbar==True:
            cb = plt.colorbar(mappable=ax.images[0])
            cb.set_label("Intensity [$W/m^2/pixel$]")


    def plot_dust_properties(self, *args, **kwargs):
        """ Plot dust scattering properties


        Plot the dust scattering properties as calculated by MCFOST
        for some given collection of grains.

        Note: For this to work you must have first saved dust properties
        to disk using the ``mcfost somefile.par -dust_prop``
        """

        from . import plotting
        plotting.plot_dust(directory=self.directory, parameters = self.parameters, *args, **kwargs)

    def describe(self):
        """ Return a descriptive brief paragraph on what results are present """
        print("Model results in {self.directory} for {parfn}".format(self=self, parfn = os.path.basename(self._paramfilename)))
        print("    Model has {0} inclinations from {1} to {2}".format(len(self.parameters.inclinations), min(self.parameters.inclinations), max(self.parameters.inclinations)))
        if self.sed is None:
            print("    No SED data present")
        else:
            print("    SED computed from {par.lambda_min} - {par.lambda_max} microns using {par.nwavelengths} wavelengths".format(par=self.parameters))
        print("    Images computed for {0} wavelengths: {1}".format(len(self.image_wavelengths), self.image_wavelengths))

        if os.path.exists(os.path.join(self.directory, 'data_dust/kappa.fits.gz')):
            print("    Dust scattering properties have been saved to disk for that model grain population.")

    def calc_chisqr(self, observed, weights=None, **kwargs):
        """ Calculate chi^2 statistic based on comparison with observations


        Parameters
        -----------
        observed : MCFOST.Observations instance
            Available observational data to be compared with models
        weights : list of floats
            Weights for computing the combined chi^2 merging all sub-
            chi^2 values.  Provide these in the order [SED, shortest wavelength image,
            second shortest wavelength image, ... longest wavelength image].
        save : bool
            Save chi^2 results to disk

        Notes
        -------

        Many other parameters are accepted and passed to the individual chi^2 functions.
        See mcfost.chisqr.sed_chisqr and mcfost.chisqr.image_chisqr for details.

        Note that the chi^2 values calculated will typically be an array, with one chi^2 value
        for each inclination present in the MCFOST models.

        """
        from . import chisqr

        if weights is None:
            weights = [1.0]* (self.image_wavelengths.size + 1)

        if self.sed is None:
            raise IOError("No SED data present; cannot compute chi^2.")

        _log.info("Computing SED chi^2")
        sed_chi2 = chisqr.sed_chisqr(self, observed, **kwargs)  * weights[0]

        return sed_chi2
        im_chi2 = []
        for i in range(self.image_wavelengths.size):
            _log.info("Computing image chi^2 for {0} microns".format(self.image_wavelengths[i].value) )
            im_chi2.append( chisqr.image_chisqr(self, observed, wavelength=self.image_wavelengths[i],
                **kwargs) *  weights[i+1] )


        #_log.info("Resulting chi^2s: {0} SED, {1} images".format(sed_chi2,  ",".join(im_chi2)) )


#----------------------------------------------------

class Observations(MCFOST_Dataset):
    """ Class for observational data on a target

    To the extent possible, the API is consistent with that of the ModelResults class.
    That is, the names and functionality of attributes and methods should be similar
    to the extent possible given the differences in the underlying data.

    """

    sed = None
    """ Observed SED"""

    images = None
    """ Container for model images.  """

    image_wavelengths = []
    """ List of wavelengths for which we have images """


    def __init__(self, directory='.'):
        """ Initialize based on a directory containing several results files.
        For now this follows the conventions used by the IDL results explorer code, but
        this is likely to change as the code evolves...

        Parameters
        -------------
        directory : string
            Directory path for where to find observations. Default is current working directory.


        """

        if directory is None or not os.path.isdir(directory):
            raise ValueError("Not a valid directory: "+directory)

        self.directory = directory
        summaryfile = os.path.join(self.directory, 'summary.txt')
        if not os.path.exists(summaryfile):
            raise ValueError("Cannot find index file of observed data: summary.txt")

        summaryf = open(summaryfile)
        summary = summaryf.readlines()
        #summary = open(summaryfile).readlines()

        filenames = []
        types = []
        wavelengths = []

        for i, line in enumerate(summary):
            # filename, type, wavelength (optional)
            parts = line.split()
            _log.info("Found observations: {0} is {1}".format(parts[0], " at ".join(parts[1:])) )
            thisfile = os.path.join(self.directory, parts[0])
            filenames.append(thisfile)
            types.append(parts[1].lower())
            wavelengths.append(parts[2] if len(parts) >= 3 else None)
            if parts[1].lower() == 'image': self.image_wavelengths.append(parts[2])
            elif parts[1].lower() == 'sed': self.sed = ObservedSED(thisfile)
        self.file_names = np.asarray(filenames)
        self.file_types = np.asarray(types)
        self.file_wavelengths = np.asarray(wavelengths)

        self.images = OBSImageCollection(self)
        summaryf.close()

    def __repr__(self):
        return "<MCFOST Observations in directory '{self.directory}'>".format(self=self)

    def describe(self):
        """ Return a descriptive brief paragraph on what results are present """
        print("Observations in {self.directory}".format(self=self))
        if self.sed is None:
            print("    No SED data present")
        else:
            print("    SED from {0} - {1} microns at {2} wavelengths".format(self.sed.wavelength.min(), self.sed.wavelength.max(), self.sed.wavelength.size))
        print("    Images available for {0} wavelengths: {1}".format(len(self.image_wavelengths), self.image_wavelengths))



#----------------------------------------------------

class MCFOST_SED_Base(object):
    """ Base SED class, implements common functionality

    Attributes include:
        .flux, .uncertainty: in Jy
        .nu_fnu, .nu_fnu_uncert: in W/m^2
        .wavelength:  in microns
        .frequency: in Hz
    """
    def __init__(self, ):
        pass

    @property
    def frequency(self):
        """ Frequency in Hertz, as an astropy.Quantity"""
        import astropy.constants as const
        return (const.c/self.wavelength).to(units.Hz)

    @property
    def nu_fnu(self):
        """ Spectral energy distribution in W/m^2, as an astropy.Quantity"""
        return ( self.flux * self.frequency ).to(units.W/units.m**2)

    @property
    def nu_fnu_uncert(self):
        """ Uncertainty in spectral energy distribution in W/m^2, as an astropy.Quantity"""
        return ( self.uncertainty * self.frequency ).to(units.W/units.m**2)


class ModelSED(MCFOST_SED_Base):
    """ Model SED class
    Reads model SED from disk; returns them as as astropy Units objects
    """

    directory = None
    """ Directory containing the model results."""

    filename = None
    """ Filename for the model results SED FITS file."""

    def __init__(self, directory="./",filename=None, parameters=None):
        self._sed_type = 'Model'
        self.directory = os.path.abspath(directory)
        if filename is None:
            filename = os.path.join(self.directory, 'sed_rt.fits.gz')
        if not os.path.exists(filename): raise IOError("SED files does not exist: "+filename)
        self.filename = filename

        if parameters is None:
            parameters= Paramfile(directory=self.directory)
        self.parameters=parameters

        self._sed_data = None # lazy loading

    @property
    def inclinations(self):
        """Inclinations for which the model was computed, in degrees"""
        return self.parameters.inclinations

    @property
    def wavelength(self):
        """ Wavelength in microns, as an astropy.Quantity"""
        return self.parameters.wavelengths * units.micron

    @property
    def nu_fnu(self):
        """ Spectral energy distribution in W/m^2, as an astropy.Quantity

        Note this will be a 2D array containing the SED for all inclinations that were calculated
        in that model.

        This implements lazy loading for SED data; no files are read until the
        user tries to access this attribute, at which point the SED data is
        automatically loaded."""
        # The axes of the SED data from MCFOST are [stokes, azimuth, inclination, wavelength] in Python axis order

        # if it's already loaded then just return it
        if self._sed_data is not None:
            # return the [inclination, wavelength] array for total intensity at first azimuth
            return self._sed_data[0,0] * (units.W / units.m**2)
        else:

            if os.path.exists( os.path.join(self.directory,'sed_rt.fits.gz')):
                self._RayTraceModeSED = True
                self._sed_data = fits.getdata( os.path.join(self.directory, 'sed_rt.fits.gz'))
                _log.debug("loading SED data from RT mode SED")
            elif os.path.exists( os.path.join(self.directory,'sed2.fits.gz')):
                self._RayTraceModeSED = False
                self._sed_data = fits.getdata(os.path.join(self.directory,'sed2.fits.gz'))
                _log.debug("loading SED data from MC mode SED")
            # Also look one data_th deeper for back compatibility with the case where
            # this class was created pointing at the model root directory
            elif os.path.exists( os.path.join(self.directory,'data_th/sed_rt.fits.gz')):
                self._RayTraceModeSED = True
                self._sed_data = fits.getdata( os.path.join(self.directory, 'data_th/sed_rt.fits.gz'))
                _log.debug("loading SED data from RT mode SED")
            elif os.path.exists( os.path.join(self.directory,'data_th/sed2.fits.gz')):
                self._RayTraceModeSED = False
                self._sed_data = fits.getdata(os.path.join(self.directory,'data_th/sed2.fits.gz'))
                _log.debug("loading SED data from MC mode SED")
            else:
                _log.error("No SED data present in "+self.directory+"!")


            # return the [inclination, wavelength] array for total intensity at first azimuth
            return self._sed_data[0,0] * (units.W / units.m**2)


    @property
    def uncertainty(self):
        """ Uncertainty in spectral energy distribution in W/m^2, as an astropy.Quantity"""
        return np.zeros_like(self.nu_fnu.value) * (units.Jy)

    @property
    def flux(self):
        """ Flux in Janskys, as an astropy.Quantity"""
        return (self.nu_fnu / self.frequency).to(units.Jy)


    def plot(self, title=None, inclination='all', overplot=False,
             alpha=0.75, marker='None', linestyle='-',
             nlabels=None, legend=True,
             color='red',color_imin='blue', color_imax='red',
             **kwargs):
        """ Plot observed SED

        Parameters
        -----------
        inclination : string or float
            Either the string 'all' to plot all inclinations or
            a floating point value to plot the closest available
            inclination to that value.
        nlabels : int or None
            limit the number of inclination labels in the legend? set to None to show all (default)
        legend : bool
            Draw a legend for the different lines?
        overplot : bool
            Should this overplot or make a new plot?
        alpha : float
            Matplotlib alpha channel setting for plot transparency.
        marker : str
            Matplotlib plot marker specification
        color : str or tuple
            Matplotlib color specification, for the case of plotting a single inclination
        color_imin, color_imax : str or tuple
            Matplotlib color specifications for min and max inclinations, for the case
            of plotting all inclinations in the model.
        linestyle : str
            Matplotlib line style specification

        """

        if title is None:
            title=self._sed_type+' SED from '+self.filename
        label = self._sed_type +' SED'

        if alpha is None:
            alpha = 'ramp' if str(inclination) == 'all' else 0.75

        if not overplot: plt.cla()

        if str(inclination) == 'all':
            # Plot all inclinations
            ninc = self.inclinations.size

            if nlabels is None:
                label_multiple=1
            else:
                label_multiple = np.ceil(float(ninc)/nlabels)

            c_imin = np.asarray(matplotlib.colors.colorConverter.to_rgba(color_imin))
            c_imax = np.asarray(matplotlib.colors.colorConverter.to_rgba(color_imax))
            for i in range(ninc):
                # Scale colors from color_min to color_max
                relative_pos = float(i)*(ninc+1)/ninc**2
                mycolor = tuple( c_imin*(1-relative_pos) + c_imax*relative_pos)
                label = '$i={inc:.1f}^\circ$'.format(inc=self.inclinations[i]) if np.mod(i,label_multiple) == 0 else None


                plt.loglog(self.wavelength.to(units.micron).value, self.nu_fnu[i].to(units.W/units.m**2).value,
                        label=label, linestyle=linestyle, marker=marker, color=mycolor, alpha=alpha )
        else:
            # Plot one inclination
            iclosest = np.abs(self.inclinations - float(inclination)).argmin()
            label = '$i={inc:.1f}^\circ$'.format(inc=self.inclinations[iclosest])

            plt.loglog(self.wavelength.to(units.micron).value, self.nu_fnu[iclosest].to(units.W/units.m**2).value,
                label=label, linestyle=linestyle, marker=marker, color=color, alpha=alpha )



        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("$\\nu F_\\nu$ (W m$^{-2}$)")
        plt.title(title)
        ax = plt.gca()
        ax.xaxis.set_major_formatter(utils.NicerLogFormatter())

        if legend:
            plt.legend(loc='upper right')
        return ax


class ObservedSED(MCFOST_SED_Base):
    """ Observed SED class.
    Reads observations from disk; returns them as as astropy Units objects
    """

    def __init__(self, filename=None, uncertainty=None, mask=None, format='no_header'):
        self._sed_type = 'Observed'

        # temporary hard coded default for development

        if filename is None:
            filename = './data/observed_sed.txt'
        self.filename = filename

        self._sed_table = ascii.read(filename, format=format,
                names=['wavelength','flux','uncertainty','source'])

    @property
    def flux(self):
        """ Flux in Janskys, as an astropy.Quantity"""
        return self._sed_table['flux']  * (units.Jy) # (units.W / units.m**2)

    @property
    def wavelength(self):
        """ Wavelength in microns, as an astropy.Quantity"""
        return self._sed_table['wavelength'] *units.micron

    @property
    def uncertainty(self):
        """ Uncertainty in flux in Janskys, as an astropy.Quantity"""
        return self._sed_table['uncertainty'] * (units.Jy)

    def plot(self, title=None, overplot=False,
             alpha=0.75, marker='o', color='blue', linestyle='None', **kwargs):
        """ Plot the observed SED.

        Parameters
        -----------
        overplot : bool
            Should this overplot or make a new plot?
        alpha : float
            Transparency for plots
        title : string
            Title for plot.

        Matplotlib keywords such as marker, linestyle, color, etc are also supported.

        """

        if title is None:
            title='Observed SED from '+self.filename
        label = 'Observed SED'


        if not overplot: plt.cla()

        nu_fnu = ( self.flux * self.frequency ).to(units.W/units.m**2)
        nu_fnu_uncert = ( self.uncertainty * self.frequency ).to(units.W/units.m**2)



        plt.loglog(self.wavelength.value, nu_fnu.value, label=label, linestyle=linestyle, marker=marker, color=color, alpha=alpha )
        # for some reason the errorbar function hangs if fed Quantities so extract the values first:
        plt.errorbar( self.wavelength.value, nu_fnu.value, yerr =nu_fnu_uncert.value, linestyle=linestyle, marker=marker, color=color, alpha=alpha)

        plt.xlabel("Wavelength ($\mu$m)")
        plt.ylabel("$\\nu F_\\nu$ (W m$^{-2}$)")
        plt.title(title)
        ax = plt.gca()
        ax.xaxis.set_major_formatter(utils.NicerLogFormatter())
        return ax

#----------------------------------------------------


class ObservedImage(object):
    """ Container class for an observed image,
    optionally along with associated uncertainty image and pixel mask.

    """
    def __init__(self, filename, uncertaintyfn=None, maskfn=None, wavelength=None, psffn=None):
        self.filename = filename
        self.uncertainty_filename = uncertaintyfn
        self.mask_filename = maskfn
        self.wavelength = wavelength
        self.psf_filename = psffn


#Make this more robust by assigning appropriate versions of each below when None
        self.image = fits.getdata(self.filename)
        if uncertaintyfn != None:
            self.uncertainty = fits.getdata(self.uncertainty_filename)
        self.mask = fits.getdata(self.mask_filename)
        self.psf = fits.getdata(self.psf_filename)

    

    def __repr__(self):
        return "<Observed image at {0} microns>".format(self.wavelength)



    def show(self, overplot=False, cmap=None, ax=None,
            vmin=None, vmax=None, which='image'):
        """
        Display one image

        Parameters
        ----------
        overplot : bool
            Overplot into existing axes or create new axes?
        cmap : matplotlib color map
            desired color map
        ax : matplotlib axis
            axis to plot into
        vmin, vmax : floats
            min and max for image display range.

        """
#        wm = np.where( (self.types == which) and (self.wavelengths == wavelength0) )

 #       imagefilename = self.filenames[wm]
 #       print "Filename for image: "+imagefilename
 #       raise NotImplementedError("Not implemented yet!")

class OBSImageCollection(collections.Mapping):
    """ Helper collectin class to implement on-demand loading of 
    image data via a dict interface. The dict looks like it
    always contains all available images, and can be indexed via
    either floats or strings interchangably. 

    This returns a PrimaryHDU object, not a HDUlist, under the
    assumption that MCFOST output files don't ever have more than one 
    extension.

    Most users will not need to instantiate this directly. Just
    use it transparently via a ModelResults object

    Example
    ---------
    mod = mcfostpy.ModelResults()
    im_1_micron = mod.image_data[1.0]  # looks like array access, but is 
                    # actually an automatic file read first, then access

    """
    def __init__(self, observations):
        #self._observations = observations
        #self._loaded_fits = dict()

        self.observed_images = dict()
        
        for n in observations.image_wavelengths:
            ind = np.where(observations.file_wavelengths == n)  
            file_types = observations.file_types[ind]
            filenames = observations.file_names[ind]
            for i in np.arange(len(ind[0])):
                if file_types[i].lower() == 'psf':
                    psffn = filenames[i]
                elif file_types[i].lower() == 'mask':
                    maskfn = filenames[i]
                elif file_types[i].lower() == 'image':
                    filename = filenames[i]
                elif file_types[i].upper() == 'IMAGE_UNCERT':
                    uncertfn = filenames[i]

            self.observed_images[float(n)]=ObservedImage(filename, uncertaintyfn=uncertfn, maskfn=maskfn, wavelength=n, psffn=psffn)



    @property
    def wavelengths(self):
        return np.asarray(self.observations.image_wavelengths, dtype=float) * units.micron
        #return self._modelresults.image_wavelengths

    @property   
    def filenames(self):
        return [self._getpath(k) for k in self.keys()]

    # Implement magic methods to make this compliant with the collections.Mapping abstract base class
    # that way it will behave just like a python Dict
    def keys(self):
        # cast these as a list
        return observations.image_wavelengths
        #return list(self._modelresults.image_wavelengths.value)

    def __len__(self):
        return len(self.keys())

    def __iter__(self):
        for i in self.keys():
            yield self[i]

    def __getitem__(self, key):
        canonical_wavelength = float(key)

        return self.observed_images[canonical_wavelength]


        imagefilename = self.filenames[wm]
        print("Filename for image: "+imagefilename)
        raise NotImplementedError("Not implemented yet!")

class ModelImage(object):
    """ Class for a model image, at a single wavelength
    but potentially multiple inclinations.
    """
    def __init__(self, directory, wavelength):
        self.directory = directory
        self.wavelength = wavelength

    def show(self):
        pass
