
# Model class - main class for interacting with an entire model

from . import paramfiles
from . import plotting

import warnings

class MCFOSTmodel(object):
    """ One MCFOST model result set, possibly containing multiple inclinations 
    
    This class provides an object oriented interface to MCFOST model results, including SEDs 
    and images. 
    
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

        self._directory = directory
        self._paramfilename = paramfiles.find_paramfile(self._directory)

        if self._paramfilename is None:
            raise IOError("Could not find a parameter file in that directory")
        self.parameters = paramfiles.Paramfile(self._paramfilename)


        # check for existence of SEDs but don't load yet?
        self._loaded_SEDs = False
        if not os.path.isdir(os.path.join(self._directory, 'data_th', 'sed_rt.fits.gz')):
            raise IOError("There does not appear to be an SED model output (data_th/sed_rt.fits.gz) in the model directory {0}".format(self._directory))

        # check for available images
        image_waves = glob.glob(os.path.join(self._directory, 'data_[0123456789]*'))
        if image_waves is None:
            warnings.warn("No model image results were found in that directory")
        print image_waves


   
    def _check_and_load_SEDs(self):
        """ check if SED data has been read from disk, and
        read it if necessary"""
        if self._loaded_SEDs: return
        # load SEDs here
        self._loaded_SEDs = True

    def _standardize_wavelength(self, wavelength):
        """ Utility function to return a "standardized" representation
        of a wavelength, such that e.g. "2" microns, "2.0" microns, 
        "2.00e0" microns etc are all the same thing, and that thing is
        consistent with the directory names on disk produced by MCFOST"""
        wl = float(wavelength)
        try:
            return self._wavelengths_lookup[wl]
        except:
            raise ValueError("This model does not have images for wavelength "+string(wavelength)+". Please check your data and try again.")

    def _check_and_load_images(self, wavelength):
        """ check if SED data has been read from disk, and
        read it if necessary"""
        wavelen = self._standardize_wavelength(wavelength)
        if self._loaded_images[wavelen]: return
        # load images here
        self._loaded_images[wavelen] = True

    def plot_SED(self, inclination=None, *args, **kwargs):
        """ Plot SED(s). Set inclination=value to plot just one, or leave blank for all"""
        plotting.plot_sed(parfilename=self._paramfilename, dir=self._directory)


    def plot_image(self, wavelength, inclination=None):
        pass



