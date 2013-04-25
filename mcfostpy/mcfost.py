# public code
import numpy as np
import os
import os.path
import matplotlib.pyplot as p
import pyfits
import glob
from scipy import ndimage
#import asciidata
import matplotlib.font_manager
import pywcs
import coords
import asciitable

from util import NicerLogFormatter
import logging
_log = logging.getLogger('mcfost')


_VERBOSE=True

try:
    from IPython.Debugger import Tracer; stop = Tracer()
except:
    pass

__doc__ = """
 Python routines for interacting with MCFOST output:


     - readparfile, writeparfile    functions for interacting with parameter files, represented as nested dicts
     - parfile_set_pixscale         Utility function for adjusting computation grid size to achieve a certain
                                    desired pixel scale in arcsec/pix.


    - seds                          Plot SED for all inclinations of a model
    - images                        Plot images for all available wavelengths of a model
    - image                         Plot image for one particular wavelength
    - plot_dust                     Plot dust properties (requires a 'data_dust' output dir produced by
                                    running mcfost -dust_prop in the model dir)

    - cut_psf_to_match, cut_image_to_match      Functions for cropping science images to match MCFOST models



 Marshall Perrin, mperrin@stsci.edu

"""
####################### Internal Utility routines #####################
def _find_closest(list_, item):
    """ return the index for the entry closest to a desired value """
    minelt = min(list_,key=lambda x:abs(x-item))
    return np.where(np.asarray(list_) == minelt)[0][0]

def imshow_with_mouseover(ax, *args, **kwargs):
    """ Wrapper for pyplot.imshow that sets up a custom mouseover display formatter
    so that mouse motions over the image are labeled in the status bar with
    pixel numerical value as well as X and Y coords.

    Why this behavior isn't the matplotlib default, I have no idea...
    """
    myax = ax.imshow(*args, **kwargs)
    aximage = ax.images[0].properties()['array']
    # need to account for half pixel offset of array coordinates for mouseover relative to pixel center,
    # so that the whole pixel from e.g. ( 1.5, 1.5) to (2.5, 2.5) is labeled with the coordinates of pixel (2,2)
    report_pixel = lambda x, y : "(%5.1f, %5.1f)   %g" % \
        (x,y, aximage[np.floor(y+0.5),np.floor(x+0.5)])
    ax.format_coord = report_pixel

    return ax


####################### Parameter Files #####################

class Parfile(object):
    """ Object class interface to MCFOST parameter files 

    This wraps and subsumes the previous functional interface
    in a higher level interface. Mostly still TBD...


    """
    def __init__(self, filename=None, dir="./", **kwargs):
        self._dict = readparfile(filename, dir=dir, **kwargs)
        self._filename= filename
        self._dir=dir

    def __getattr__(self, key):
        return self._dict[key]

    def write(self, outfilename):
        writeparfile(outfilename, self._dict)


def readparfile(filename=None,dir="./",silent=False,return_text=False,verbose=False):
    """
        Read in an MCFOST par file into a dictionary.

        Note that some elements of that dict are in turn smaller sub-dicts. This nested structure
        is required because some of the parameter structures can repeat multiple times in complex models.


        WARNING: Not all parameters are read in now, just the most useful ones. TBD fix later.

    ----------
    filename, dir: string
        input file and dir. Leave filename blank to search in dir.
    return_text : bool
        Return the entire parfile text as an entry in the dict
    """

    if filename is not None:
        if os.path.isdir(filename): #allow flexible input arguments...
            dir=filename
            filename=None
        dir = os.path.dirname(filename)

    if filename is None:
        if not dir.endswith(os.path.sep): dir+=os.path.sep

        fs = glob.glob(os.path.join(dir,"*.par"))
        filename=fs[0]
    text = open(filename,'r').readlines()

    def getparts(text,linenum,nparts):
        return text[linenum].split(None,nparts)[0:nparts]
    def getparts2(linenum,nparts):
        return text[linenum].split(None,nparts)[0:nparts]
    def getfloats(*args):
        return np.asarray(getparts2(*args), dtype=float)

    def get1part(linenum, partnum):
        return text[linenum].split(None,partnum+1)[partnum-1]
    def set1part(par, key,linenum, partnum, typ):
        if verbose: _log.info( str( (key,linenum, partnum, typ,)))
        try:
            par[key] = typ(text[linenum].split(None,partnum+1)[partnum-1])
            if verbose: print par[key]
        except:
            raise IOError("Could not parse line %d:'%s' item %d for '%s'" % (linenum, text[linenum][:-1], partnum, key))


    par={}
    par['filename'] = filename
    if return_text:
        par['text'] = text

    par['version'] = get1part(0,1)


    if verbose or _VERBOSE: _log.info("Parsing MCFOST parameter file version "+par['version'])

    # figure out offsets of the fixed lookup table
    if float(par['version']) >= 2.15:
        d1=-1   # removal of nbr_parallel_loop
        d2=-4   # and also 'consider disk emissions' flag
        d3=-5   # and also gas-to-dust parameter
    elif float(par['version']) >= 2.11:
        d1=0
        d2= -2
        d3=d2
    elif float(par['version']) >= 2.10:
        d1=0
        d2= -1
        d3=d2
    else:
        d1=0
        d2=0
        d3=d2
    #print "offset is ", d

    # TODO add distance and other image parameters



    # read in fixed parameters
    fixed_keys=(
    ('nbr_photons_eq_th', 4+d1, 1, float),
    ('nbr_photons_lambda',5+d1, 1, float),
    ('nbr_photons_image', 6+d1, 1, float),
    ('nlambda', 9+d1, 1, int),
    ('lambda_min', 9+d1, 2, float),
    ('lambda_max', 9+d1, 2, float),
    ('flag_ltemp', 10+d1, 1, str),
    ('flag_lsed', 10+d1, 2, str),
    ('flag_lcomplete', 10+d1, 3, str),
    ('grid:type', 18+d2, 1, int),
    ('grid:nrad', 19+d2, 1, int),
    ('grid:nz', 19+d2, 2, int),
    ('grid:n_az', 19+d2, 3, int),
    ('grid:n_rad_in', 19+d2, 4, int),
    ('grid:ntheta', 22+d2, 1, int),         # for SED
    ('grid:nphi',22+d2 , 2, int),
    ('im:nx',22+d2 , 3 if float(par['version']) < 2.15 else 1, int),
    ('im:ny',22+d2 , 4 if float(par['version']) < 2.15 else 2, int),
    ('im:inc_min', 24+d2, 1, float),
    ('im:inc_max', 24+d2, 2, float),
    ('im:ninc',24+d2 , 3, int),           # for RT
    ('im:centered',24+d2 , 4, str),
    ('distance', 25+d2, 1, float),
    ('disk_pa', 26+d2, 1, float),  # only for > v2.09
    ('settling_flag', 39+d3, 1, str),
    ('settling_exp' ,39+d3, 2,float),
    ('settling_a' ,39+d3, 3,float),
    ('sublimate_flag' ,40+d3, 1, str))


    for key, line, item, typ in fixed_keys:
        set1part(par,key, line, item, typ)

    if float(par['version']) >=2.15:
        set1part(par, 'im:size_au', 22+d2, 3, float),
        set1part(par, 'im:zoom', 22+d2, 4, float),
        i=39  # points to nzones parameter # note: advance this counter ** after** you have read something in...
    else:
        set1part(par, 'nbr_parallel_loop', 3, 1, int),
        set1part(par, 'gas_to_dust', 38+d2, 1, float),
        set1part(par, 't_start' ,42+d3, 1,float),
        set1part(par, 'sublimtemp' ,42+d3, 2,float),
        i=43  # points to nzones parameter # note: advance this counter ** after** you have read something in...

   

    # compute or look up wavelength solution
    if par['flag_lcomplete'] == 'T':
        # use log sampled wavelength range
        lambda_inc=np.exp( np.log(par['lambda_max']/par['lambda_min'])/(par['nlambda']))
        lambd = par['lambda_min'] * lambda_inc**(np.arange(par['nlambda'])+0.5)
    else:
        # load user-specified wavelength range from file
        par['lambda_file'] = get1part(11,1)
        possible_lambda_files = [ dir+os.sep+par['lambda_file'], 
             dir+os.sep+"data_th"+os.sep+par['lambda_file'], 
             os.getenv('MCFOST_UTILS')+os.sep+'Lambda'+os.sep+par['lambda_file'] ]
        lambda_file = None
        for possible_name in possible_lambda_files:
            if _VERBOSE: _log.info("Checking for wavelengths file at "+possible_name)
            if os.path.exists(possible_name):
                lambda_file=possible_name
                if _VERBOSE: _log.info("Found wavelengths file at "+lambda_file)
                break

        if not os.path.exists(lambda_file):
            raise IOError('Cannot find requested wavelength file: '+lambda_file)
        lambd= asciitable.read (lambda_file, 
                data_start=0,names=['wavelength'])['wavelength']
        
    par['lambda'] = lambd


    # compute inclinations used for SED
    inc = np.arange(par['grid:ntheta'])+1
    n_thet = par['grid:ntheta']
    incmin    = np.arccos(1.-((inc-1  )*1.0)/n_thet)*180./3.1415926
    incmax    = np.arccos(1.-((inc    )*1.0)/n_thet)*180./3.1415926
    incmiddle = np.arccos(1.-((inc-0.5)*1.0)/n_thet)*180./3.1415926
    par['sed:inclinations'] = incmiddle


    # compute inclinations used for RT
    ##    ;How to calculate the inclinations depends on whether you're using the values centered in those bins or not.
    dtor = np.pi/180
    if par['im:centered'] == "T":
        # find the average cos(i) in each bin (averaged linearly in cosine space)
        a1 = np.linspace(np.cos(par['im:inc_max']*dtor), np.cos(par['im:inc_min']*dtor), par['im:ninc']+1)
        a2 = (a1[0:-1]+a1[1:])/2
        par['im:inclinations'] = (np.arccos(a2)/dtor)[::-1]
    else:
        # just use the bottom value for each bin
        a1 = np.linspace(np.cos(par['im:inc_max']*dtor), np.cos(par['im:inc_min']*dtor), par['im:ninc'])
        par['im:inclinations'] = (np.arccos(a1)/dtor)[::-1]



    # read in the different zones
    set1part(par, 'nzones' ,i, 1,int)
    i+=3

    for idensity in range(par['nzones']):
        density={}
        density_keys =(('type', i+0, 1, int),
                ('dust_mass', i+1, 1, float),
                ('scale_height', i+2, 1, float),
                ('ref_radius', i+2, 2, float),
                ('r_in', i+3, 1, float),
                ('r_out', i+3, 2, float),
                ('size_neb', i+3, 3, float),
                ('edge', i+3, 4 if float(par['version']) <2.15 else 3, float),
                ('flaring_exp', i+4, 1, float),
                ('surfdens_exp', i+5, 1, float))
        for key, line, item, typ in density_keys:
            set1part(density,key, line, item, typ)
        i+=6
        if idensity == 0:
            par['geometry']=density
        else:
            par[('geometry',idensity)]=density

    i+=2
    cavity_keys=(('cavity_flag', i+0, 1, str),
                ('cavity_height', i+1, 1, float),
                ('cavity_ref_radius', i+1, 2, float),
                ('cavity_flaring_exp', i+2, 1, float))
    for key, line, item, typ in cavity_keys:
        set1part(par,key, line, item, typ)


    # read in the dust geometry. One set of dust props **per zone**
    i+=5
    for idensity in range(par['nzones']):
        if idensity ==0: zonekey='geometry'
        else: zonekey=('geometry',idensity)
        dust={}
        set1part(dust,'nspecies', i, 1,int)
        i+=1
        for idust in range(dust['nspecies']):

            if float(par['version']) >= 2.12:
                # version 2.12 or higher, allow for multi-component grains
                dust_keys = (('ncomponents', i+0, 1, int),
                        ('mixing_rule', i+0, 2, int),
                        ('porosity', i+0, 2, float),
                        ('mass fraction', i+0, 3, float))
                for key, line, item, typ in dust_keys:
                    set1part(dust,key, line, item, typ)
                i+=1
                if dust['ncomponents'] >1: raise NotImplementedError("Need multi-component parsing code!")
                for icomponent in range(dust['ncomponents']):
                    dust_keys = (('filename', i+1, 1, str),
                        ('vol_fraction', i+1, 2, str))
                    for key, line, item, typ in dust_keys:
                        set1part(dust,'filename', i, 1, str)

                    i+=1
                # now the heating and grain size properties
                dust_keys = ( ('heating', i+0, 1, int),
                        ('amin', i+1, 1, float),
                        ('amax', i+1, 2, float),
                        ('aexp', i+1, 3, float),
                        ('ngrains', i+1, 4, int))
                for key, line, item, typ in dust_keys:
                    set1part(dust,key, line, item, typ)
                i+=2

            else:
                # earlier versions than 2.12, so only one component allowed. 
                dust_keys = (('filename', i+0, 1, str),
                        ('porosity', i+0, 2, float),
                        ('mass fraction', i+0, 3, float),
                        ('heating', i+1, 1, int),
                        ('amin', i+2, 1, float),
                        ('amax', i+2, 2, float),
                        ('aexp', i+2, 3, float),
                        ('n grains', i+2, 4, int))
                for key, line, item, typ in dust_keys:
                    set1part(dust,key, line, item, typ)
                i+=3
            if idust == 0:
                par[zonekey]['dust']=dust
                if idensity == 0: par['dust'] = dust # special case for only one kind of dust globally
            else:
                par[zonekey][('dust',idust)]=dust


    #
    #
    #    #skip all but dust for now...
    #    # and assume only 1 dust type
    #
    #    i=62
    #    par[('dust',0)]={}
    #    t = getparts(text,i,3); i+=1
    #    par[('dust',0)]['filename']=t[0]
    #    par[('dust',0)]['porosity']=float(t[1])
    #    par[('dust',0)]['mass fraction']=float(t[2])
    #    t = getparts(text,i,1); i+=1
    #    par[('dust',0)]['heating']=float(t[0])
    #    p = getparts(text,i,4); i+=1
    #    par[('dust',0)]['amin']=float(p[0])
    #    par[('dust',0)]['amax']=float(p[1])
    #    par[('dust',0)]['aexp']=float(p[2])
    #    par[('dust',0)]['nbr_grains']=float(p[3])
    #


    i+= 12 # skip molecular RT settings

    set1part(par,'nstar', i, 1,int)
    i+=1
    for istar in range(par['nstar']):
        star_keys = [('temp', i+0, 1, float) ,
                ('radius', i+0, 2, float),
                ('mass', i+0, 3, float),
                ('x', i+0, 4, float),
                ('y', i+0, 5, float),
                ('z', i+0, 6, float),
                ('spectrum', i+1, 1, str)]
        if float(par['version']) >= 2.11:
            star_keys.append( ('fUV', i+2, 1, float) )
            star_keys.append( ('slope_fUV', i+2, 2, float) )
            nlines_star=3
        else:
            nlinse_2tar=2
        star={}
        for key, line, item, typ in star_keys:
                set1part(star,key, line, item, typ)
        i+=nlines_star

        if istar == 0:
            par['star']=star
        else:
            par[('star',istar)]=star


    # now try to read in command line options which can override
    # some of the settings.
    # set some defaults first.
    par['im:zoom'] =1.0
    par['im:raytraced'] = False
    i+=1
    if len(text) >= i:
        optionsline = text[i]
        options = optionsline.split()
        for j in np.arange(len(options)):
            if options[j] == "-img":
                par['im:wavelength'] = options[j+1]
            elif options[j] == "-resol":
                par['im:nx'] = np.int(options[j+1])
                par['im:ny'] = np.int(options[j+2])
            elif  options[j] == "-zoom":
                par['im:zoom'] = np.int(options[j+1])
            elif  options[j] == "-rt":
                par['im:raytraced'] = True
    else:
        if verbose: print "could not read in command line options from MCFOST; using default grid settings"


    # compute a few extra things
    xpix_au = (2*par['geometry']['size_neb'] / par['im:zoom']) / par['im:nx']  # au per pixel
    xpix_arcsec = xpix_au / par['distance']# /zoom_factor       # arcsec per pixel
    ypix_au = (2*par['geometry']['size_neb'] / par['im:zoom']) / par['im:ny']  # au per pixel
    ypix_arcsec = ypix_au / par['distance']# /zoom_factor       # arcsec per pixel

    par['xscale au/pix'] = xpix_au
    par['xscale arcsec/pix'] = xpix_arcsec
    par['yscale au/pix'] = ypix_au
    par['yscale arcsec/pix'] = ypix_arcsec
    par['im:xsize au'] =  (2*par['geometry']['size_neb'] / par['im:zoom'])
    par['im:ysize au'] =  (2*par['geometry']['size_neb'] / par['im:zoom'])

    par['im:xsize arcsec'] = par['im:xsize au'] / par['distance']
    par['im:ysize arcsec'] = par['im:ysize au'] / par['distance']



    if return_text: return (text,par)
    else: return par


def writeparfile(outname='sample.par', pardict=None):
    """ Write an MCFOST parameter file to disk


    WARNING: Not all parameters are written yet, just the most useful ones. TBD fix later.


    Parameters
    ----------
    outname : string
        Output file name
    pardict : dict
        MCFOST parameter information as a nested set of dictionaries.


    HISTORY
    --------
    2013-01-05 updated to version 2.15

    """
    par = pardict

    getkeys = lambda l: tuple([par[li] for li in l])

    # if we read in a <2.14 version file, we have a nbr_parallel_loop factor that
    # was multiplicatively merged in with the other nbrs for later versions, so
    # do that multiplication here before writing out a >= 2.15 parfile
    try:
        scale_nbr = par['nbr_parallel_loop']
    except:
        scale_nbr = 1

    template = """2.15                      mcfost version

#Number of photon packages
  %-10d              nbr_photons_eq_th  : T computation
  %-10d              nbr_photons_lambda : SED computation
  %-10d              nbr_photons_image : images computation
""" %  tuple( scale_nbr * np.asarray(getkeys(['nbr_photons_eq_th', 'nbr_photons_lambda', 'nbr_photons_image'])))

    template+= """
#Wavelength
  %-3d %5.1f %6d        n_lambda
  T T T                   compute temperature?, compute sed?, use default wavelength grid ?
  IMLup.lambda            wavelength file (if previous parameter is F)
  F T                     separation of different contributions, stokes parameters

#Grid geometry and size
  %-10d              1 = cylindrical, 2 = spherical
  %-3d %3d %2d %3d          n_rad, nz (or n_theta), n_az, n_rad_in
  """ % getkeys(['nlambda', 'lambda_min', 'lambda_max', 'grid:type', 'grid:nrad', 'grid:nz', 'grid:n_az', 'grid:n_rad_in'])


# FIXME - remove ntheta, nphi, add in size-AU (which is the 'size_neb') variable now)
# update MC nbin line in reader
    template+= """
#Maps
  %3d %3d %4d %4d  1.0   grid (nx,ny), size [AU], zoom factor
  10 1   1                MC : N_bin_incl, N_bin_az, # of bin where MC is converged
  %-4.1f  %4.1f  %d   %s  RT: imin, imax, n_incl, centered ?
  %-6.1f                  distance (pc)
  %-6.1f                  disk PA
   """ % getkeys(['grid:ntheta', 'grid:nphi', 'im:nx', 'im:ny', 'im:inc_min', 'im:inc_min', 'im:ninc', 'im:centered', 'distance', 'disk_pa'])

    template+="""
#Scattering method
  0                      0=auto, 1=grain prop, 2=cell prop
  1                      1=Mie, 2=hg (2 implies the loss of polarizarion) 

#Symetries
  T                      image symmetry
  T                      central symmetry
  T                      axial symmetry (important only if N_phi > 1)

#Dust global properties
  %s  %6.2f  %5.1f       dust settling, exp_strat, a_strat
  F                      sublimate dust
  F  0.0                 viscous heating, viscosity

#Number of zones : 1 zone = 1 density structure + corresponding grain properties
  %-d
  """ % getkeys(['settling_flag', 'setting_exp', 'setting_a', 'nzones'])

    if par['nzones'] > 1:
        raise NotImplemented("write par file code does not yet support multiple zones!")

    getsubkeys = lambda k, l: tuple([par[k][li] for li in l])
    template+="""
#Density structure
  %-1d                       zone type : 1 = disk, 2 = envelope
  %-10.3e %4d             disk dust mass,  gas-to-dust mass ratio
  %-5.1f  %6.1f           scale height, reference radius (AU), unused for envelope
  %-6.1f  %6.1f %6.1f   Rin, Rout (or Rc), edge falloff (AU)
  %-6.1f                  flaring exponent, unused for envelope
  %-6.1f                  surface density exponent (or -gamma for tappered-edge disk), usually < 0
    """  % getsubkeys('geometry', ['type', 'dust_mass', 'scale_height', 'ref_radius', 'r_in', 'r_out', 'size_neb', 'edge', 'flaring_exp', 'surfdens_exp'])

    template+="""
#Cavity : everything is empty above the surface
 %s                        cavity ?
 %5.1f %5.1f              height, reference radius (AU)
 %5.1f                    flaring exponent
    """ % getkeys(['cavity_flag', 'cavity_height', 'cavity_ref_radius', 'cavity_flaring_exp'])

    if par['dust']['ndust'] > 1: raise NotImplemented("No support for multi dust types yet")
    template+="""
#Grain properties
  %-2d                      Number of species
  1  2 %5.2f  1.0   N_components, mixing rule (1 = EMT or 2 = coating),  porosity, mass fraction
  %s  %4.1f   optical indices file , porosity, mass fraction
  %-2d                      heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE
  %-5.1f %6.1f  %4.1f %3d   amin, amax, aexp, nbr_grains

#Molecular RT settings
  T T T 15.              lpop, laccurate_pop, LTE, profile width
  0.2                    v_turb (delta)
  1                      nmol
  co@xpol.dat 6          molecular data filename, level_max
  1.0 50                 vmax (m.s-1), n_speed
  T 1.e-6 abundance.fits.gz   cst molecule abundance ?, abundance, abundance file
  T  3                   ray tracing ?,  number of lines in ray-tracing
  1 2 3                  transition numbers
  """ % getsubkeys('dust', ['ndust', 'filename', 'porosity','mass fraction','heating','amin','amax','aexp','n grains'])

    if par['nstar'] >1 : raise NotImplemented("Multiple stars not supported")
    starkeys = lambda l: tuple([par['star'][li] for li in l])
    template+="""
#Star properties
  1 Number of stars
  %-7.1f %5.1f %5.1f %5.1f %5.1f %5.1f F       Temp, radius (solar radius),M (solar mass),x,y,z (AU), is a blackbody?
  %s
  0.0    2.2  fUV, slope_fUV
    """  % starkeys(['temp', 'radius', 'mass', 'x', 'y', 'z', 'spectrum'])

    outfile = open(outname, 'w')
    outfile.write(template)
    outfile.close()
    print "  ==>> "+outname




def get_current_parfile(parfile=None, dir="./",verbose=False, wavelength=None):
    """
        Find a MCFOST par file in a specified directory.

        By default, look in the current (hopefully top level) directory of a model,
        but if a wavelength is specified, look inside the appropriate data_## directory

    KEYWORDS:
        wavelength=     string with wavelength name
        verbose         Boolean.
        dir=            directory to look in. Default is current directory

    """

    if parfile is not None:
        # TODO validate its existence, check if path name relative to dir?
        output = parfile
    else:
        if _VERBOSE: _log.info("Looking for par files in dir: "+dir)
        dir = os.path.expanduser(dir)
        if wavelength is not None:
            dir = dir+"/data_%s/" % wavelength
            _log.info("Since wavelength is %s, looking in %s subdir" % (wavelength,dir))
        l = glob.glob(dir+os.sep+"*.par")
        if len(l) == 1:
            output = l[0]
        elif len(l) > 1:
            _log.warning("Multiple par files found... returning first")
            output = l[0]
        else:
            _log.error("No par files found!")
            output = None

    if _VERBOSE: 
        _log.info('Par file found: '+str(output))
    return output


def parfile_set_pixscale(pardict, desired_pixelscale):
    """ Adjust sampling to make the pixel scale as desired

    desired_pixelscale : in arcsec per AU
    """
    # pixelscale = ( (2 * size_neb) / distance ) / Npix
    # therefore:
    # npix = ( (2 * size_neb) / distance ) / pixelscale
    desired_npix = ( (2 * pardict['geometry']['size_neb']) / pardict['distance'] ) / desired_pixelscale
    pardict['im:nx'] = desired_npix
    pardict['im:ny'] = desired_npix
    print "Image size set to (%d, %d) to get pixel size of %f for %f AU at %d pc" % (pardict['im:nx'], pardict['im:ny'], desired_pixelscale, 2 * pardict['geometry']['size_neb'], pardict['distance'])

####################### Utility Routines & Fit Preparation #####################

def cut_psf_to_match(filename, outname=None, cutsize=79):
    """
        make a cutout from some PSF (presumably from Tiny Tim?)
        to use with MCFOST
    """
    input = pyfits.open(filename)

    if outname is None:
        outname="mcfost_"+filename

    # find the location of the PSF peak
    mx, my = ndimage.maximum_position(input[0].data)

    cut = input[0].data[mx-cutsize/2:mx-cutsize/2+cutsize, my-cutsize/2:my-cutsize/2+cutsize]
    input[0].data=cut

    mx2, my2 = ndimage.maximum_position(cut)
    input[0].header.update('CENTER_X', mx2, comment="PSF Center")
    input[0].header.update('CENTER_Y', my2, comment="PSF Center")
    input[0].header.add_history("  Created by mcfost.cut_psf_to_match (Python) ")
    input[0].header.add_history("  Created from "+filename)

    input.writeto(outname)

    print "Cutout "+str(cutsize)+" pixels across writen to "+outname


def getpixscale(header):
    """  Get pixel scale from header

        usage: getpixscale(FITSheader)

        sort of like IDL's getrot.pro:
        Given a FITS header with astrometry, figure out what the
        pixel scale is.

        This actually uses the full PyWCS+Coords machinery, to handle the
        general case of all the differnt ways to specify this (CD vs CDELT etc).

        There ought to be a better way to do this?

        NOTE: This will **fail** for WCS headers that pywcs doesn't understand, like polarization.

    """
    wcs = pywcs.WCS(header)
    ans = wcs.wcs_pix2sky( np.array([0,1]), np.array([0,0]), 0)
    p1 =  coords.Position( wcs.wcs_pix2sky( np.array([0]),  np.array([0]), 0) )
    p2 =  coords.Position( wcs.wcs_pix2sky( np.array([0]),  np.array([1]), 0) )
    sep = p1.angsep(p2)

    return sep.arcsec()


def cut_image_to_match(imagefile,wavelength, parfile=None, dir="./",verbose=True, outputdir="./", returndata=False, pixelscale=None):
    print dir
    par = readparfile(get_current_parfile(parfile=parfile,dir=dir,wavelength=wavelength))

    # figure out how big the MCFOST images are
    imsize = np.array([ par['im:nx'], par['im:ny'] ]) # full size in pixels
    imsize_AU = par['im:xsize au']  # size in AU
    imsize_arcsec = par['im:xsize arcsec']  # size in AU

    if verbose: print "MCFOST model image is %f x % f pixels" % (imsize[0], imsize[1])
    if verbose: print "MCFOST model image is %f x % f ''" % (imsize_arcsec, imsize_arcsec)

    # grab the input FITS image, find the center, and cut out the appropriate subimage.
    # **note** this code assumes a bunch of stuff about the input files, mostly that
    # they are from my HST NICMOS IDL code.

    im = pyfits.open(imagefile)
    cenx  = im[0].header['NEWCEN_X']
    ceny  = im[0].header['NEWCEN_Y']

    if pixelscale is None:
        scale = getpixscale(im[0].header)
    else:
        scale = pixelscale

    cutsize = round(imsize_arcsec/scale)
    sx = cenx - cutsize/2
    sy = ceny - cutsize/2
    if verbose: print "given the pixel scale of %f arcsec/pixel" % scale
    if verbose: print "I need to cut out a %d x %d chunk of the image to match that." % (cutsize, cutsize)


    if im[0].data.ndim == 3:
        # polarized mode data in cubes
        newregion = im[0].data[:, sy:sy+cutsize, sx:sx+cutsize]
        newregion_rescaled = np.zeros((newregion.shape[0],imsize[0],imsize[1]))
        for i in range(newregion.shape[0]):
            #if np.any(np.isnan(newregion[i,:,:])): method='nearest'
            #else:
            method='spline'
            # TODO fix NaNs
            newregion_rescaled[i,:,:] = congrid(newregion[i,:,:], imsize, method=method, center=True)

    else:
        newregion = im[0].data[sy:sy+cutsize, sx:sx+cutsize]
        newregion_rescaled = congrid(newregion, imsize, method='spline', center=True)

        newregion_err = im[1].data[sy:sy+cutsize, sx:sx+cutsize]
        newregion_err_rescaled = congrid(newregion_err, imsize, method='spline', center=True)




    # resample that image to match the scale of the MCFOST image.
    hlist = pyfits.HDUList()
    hlist.append( pyfits.PrimaryHDU(newregion_rescaled) )
    hlist.append( pyfits.ImageHDU(newregion_err_rescaled, name="ERR") )
    #outhdulist = pyfits.HDUList(output)
    targname =  im[0].header['TARGNAME']
    outfile = "cropped_%s_wave=%s.fits" % (targname, wavelength)
    hlist.writeto(outfile, clobber=True)
    # also write uncertainty to a separate file
    outfile2 = "uncert_%s_wave=%s.fits" % (targname, wavelength)
    pyfits.PrimaryHDU(newregion_err_rescaled).writeto(outfile2)

    if verbose:
            print "  ===>>> %s  [sci]" % outfile
            print "  ===>>> %s  [err]" % outfile2
    if returndata: return (im, newregion, newregion_rescaled)
    #return par

####################### Plotting & Analyzing Model Outputs #####################

def seds(parfilename=None,dir="./", overplot=False, nlabels=None, alpha=0.75,
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
        p.cla()


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
                incmin=np.arccos(1.-((inc)*1.0)/ninc)*180./3.1415926
                incmax=np.arccos(1.-((inc+1  )*1.0)/ninc)*180./3.1415926
                label = "%4.1f - %4.1f$^o$" % (incmin,incmax)


            flux = sed[0,phi-1,inc,:]
            if np.mod(inc,labelstep) !=0: # allow skipping some labels if many are present
                label=None
            p.loglog(lambd, flux, color=((ninc-inc)*1.0/ninc, 0, 0), label=label, alpha=alpha)
    else:
        wmin = np.argmin( abs(par['im:inclinations'] - inclination))
        print "Closest inclination found to %f is %f. " % (inclination, par['im:inclinations'][wmin])
        label = "%4.1f$^o$" % (par['im:inclinations'][wmin])
        flux = sed[0,phi-1,wmin,:]
        p.loglog(lambd, flux, color=((ninc-wmin)*1.0/ninc, 0, 0), label=label, alpha=alpha)
           

    p.xlabel("$\lambda$ ($\mu$m)")
    p.ylabel("$\\nu F_\\nu$ (W m$^{-2}$)")
    p.title("SED for "+title)
    p.gca().xaxis.set_major_formatter(NicerLogFormatter())


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
        #seddata = np.asarray(asciidata.open(sedfile)[0:3])
        #obswavelen = seddata[0]
        #flux = seddata[1]*(1.e-26)*(3.e14/obswavelen) # convert from Jy to nu Fnu in W m ^-2
        #errflux=seddata[2]*(1.e-26)*(3.e14/obswavelen) #  convert from Jy to nu Fnu in W m ^-2

        observed = asciitable.read(sedfile, Reader=asciitable.Tab)
        flux = observed.Flux*(1.e-26)*(3.e14/observed.Wavelen) # convert from Jy to nu Fnu in W m ^-2
        errflux=observed.Uncert*(1.e-26)*(3.e14/observed.Wavelen) #  convert from Jy to nu Fnu in W m ^-2

        w_meas = np.where(flux > 0)
        w_upper = np.where(  (np.isnan(flux) | (flux == 0)) & np.isfinite(errflux) )

        p.errorbar(observed.Wavelen[w_meas], flux[w_meas], yerr=errflux[w_meas], label="Observed", color='blue', fmt='o')
        p.plot(observed.Wavelen[w_upper], errflux[w_upper], 'rv')
        #print(flux)

    p.legend(prop = {'size':10})
    p.draw()


def lir_lstar(parfilename=None,dir="./", inclination=0):
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


    p.clf()
    for i in range(9):
        ls = "--" if i > 4 else "-"
        flux_nufnu = sed[i,phi-1,inclination,:]   # this is in nu F_nu units.
        p.loglog(lambd, flux_nufnu, label="channel %d" % i, linestyle=ls)
    p.legend()

    # channels appear to be: 
    # 0 total
    # 1 ??
    # 2 ??
    # 3 ??
    # 4 star
    # 5,6 disk

    p.draw()



    def integrate_channel(i):
        flux_nufnu = sed[i,phi-1,inclination,:]   # this is in nu F_nu units.
        flux_flambda = flux_nufnu/lambd
        integrated_flux = np.trapz(flux_flambda, lambd)
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




    








def images(parfilename=None,dir="./", overplot=False, verbose=True, psf_fwhm=None, **kwargs):
    """
        Plot all available images for the current model

    keywords:


    """
    par = readparfile(get_current_parfile(parfile=parfilename,dir=dir))

    if not overplot: p.cla()

    ims = glob.glob(dir+os.sep+"data_*/RT.fits.gz")

    wavelens =  [i[i.find('_')+1:i.find("/RT.fits.gz")] for i in ims]
    if verbose or _VERBOSE: _log.info("     Wavelengths found: "+str( wavelens))


    for w,ct in zip(wavelens, np.arange(len(wavelens))):
        p.subplot(len(wavelens), len(wavelens), ct+1)
        image(w, par=par, dir=dir, **kwargs)


def image(wavelength, parfilename=None,par=None, dir="./", overplot=False, inclination=80, cmap=None, ax=None, 
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
        if ax is None: p.clf()
        else: p.cla()


    if ax is None: ax = p.gca()
    if par is None:
        par = readparfile(get_current_parfile(parfile=parfilename,dir=dir))
    if cmap is None:
        cmap = p.cm.gist_heat
        cmap = p.cm.gist_gray
        cmap.set_over('white')
        cmap.set_under('black')
        cmap.set_bad('black')

    rt_im = pyfits.getdata(dir+os.sep+"data_"+wavelength+os.sep+"RT.fits.gz")
    inclin_index = _find_closest(par['im:inclinations'], inclination)
    #print "using image %s, inc=%f" % (  str(inclin_index), par['im:inclinations'][inclin_index]  )

    #ax = p.subplot(151)
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
        ax = p.subplot(121)

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
        polfrac_ar = np.sqrt(imQn**2 + imUn**2)
        theta = 0.5* np.arctan2(imUn, imQn) + np.pi/2
        vecX = polfrac_ar * np.cos(theta) *-1
        vecY = polfrac_ar * np.sin(theta)

        Y, X = np.indices(image.shape)
        Q = p.quiver(X[::showevery, ::showevery], Y[::showevery, ::showevery], vecX[::showevery, ::showevery], vecY[::showevery, ::showevery],
                headwidth=0, headlength=0, headaxislength=0.0, pivot='middle', color='white')
        p.quiverkey(Q, 0.85, 0.95, 0.5, "50% pol.", coordinates='axes', labelcolor='white', labelpos='E', fontproperties={'size':'small'}) #, width=0.003)
#        ax = p.subplot(152)
#        ax.imshow(imQn, vmin=-1, vmax=1, cmap=matplotlib.cm.RdBu)
#        ax.set_title('Q')
#        ax = p.subplot(153)
#        ax.imshow(imUn, vmin=-1, vmax=1, cmap=matplotlib.cm.RdBu)
#        ax.set_title('U')
#        ax = p.subplot(154)
#        ax.imshow(vecX, vmin=-1, vmax=1, cmap=matplotlib.cm.RdBu)
#        ax.set_title('vecX')
#        ax = p.subplot(155)
#        ax.imshow(vecY, vmin=-1, vmax=1, cmap=matplotlib.cm.RdBu)
#        ax.set_title('vecY')
#        ax.colorbar()

    #p.colorbar()

        if polfrac:
            ax2 = p.subplot(122)

            cmap2 = p.cm.gist_heat
            cmap2 = p.cm.jet
            cmap2.set_over('red')
            cmap2.set_under('black')
            cmap2.set_bad('black')

            ax2 = imshow_with_mouseover(ax2, polfrac_ar, vmin=0, vmax=1,  cmap=cmap2)
            ax2.set_title("Polarization fraction")

            cax = p.axes([0.92, 0.25, 0.02, 0.5])
            p.colorbar(ax2.images[0], cax=cax, orientation='vertical')



    p.draw()

def plot_dust(dir='./', noerase=False):
    """
      NAME:  mc_dust

         Plot the dust scattering properties as calculated by MCFOST
         for some given collection of grains.

      INPUTS:
      KEYWORDS:
         dir=        a MCFOST results directory

         /noerase    don't erase
         window=        window to output to
         charsize=    text size for labels.


      OUTPUTS:

         draws four plots on the screen.

      HISTORY:
        Began 2009-05-22 17:33:13 by Marshall Perrin
    """


    if not os.path.isdir(dir):
       print "ERROR - invalid directory"

    if not os.path.isfile(dir + 'data_dust/kappa.fits.gz'):
       print("No dust properties files exist in that directory! need to compute some using MCFOST...")
       return False

    if noerase is False:  p.clf()

    lambd = pyfits.getdata(dir + "data_dust/lambda.fits.gz")
    g = pyfits.getdata(dir + "data_dust/g.fits.gz")
    albedo = pyfits.getdata(dir + "data_dust/albedo.fits.gz")
    kappa = pyfits.getdata(dir + "data_dust/kappa.fits.gz")
    polar = pyfits.getdata(dir + "data_dust/polar.fits.gz")


    #par = mc_readparfile((file_search(dir + "data_dust/*.par"))[0])

    #_sys_p.multi = concatenate([0, 1, 5])
    #----
    #multiplot() ; _sys_p.position += concatenate([0, -0.1, 0, -0.1])
    p.subplots_adjust(top=0.98, hspace=0.3)
    p.subplot(411)
    p.semilogx(lambd, kappa)
    p.ylabel("Opacity $\kappa$")


    #----
    p.subplot(412)
    p.semilogx(lambd, albedo)
    p.ylabel("Albedo")


    #----
    ax = p.subplot(413)
    p.semilogx(lambd, g)
    p.ylabel("g= <cos $\\theta$>")
    p.xlabel("Wavelength")
    ax.set_ybound([0.0, 1.0])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])




    #----- now switch to show the polarization versus angle ----

    ax = p.subplot(414)
    pos = ax.get_position()  # move it down slightly for visual offset
    pa = pos.get_points()
    pa[0,1] *= 0.5
    pos.set_points(pa)
    ax.set_position(pos)

    polwaves = np.asarray([1.0, 2.5, 10.0, 100.0])
    theta = np.r_[0:180]
    p.xlabel("Scattering Angle [degrees]")
    p.ylabel("Polarization $[-S_{12}/S_{11}]$")
    p.axis(xmin=0,xmax=180, ymin=0,ymax=1)
    p.xticks(np.arange(7)*30)
    for i in np.arange(polwaves.size):
       c = _find_closest(lambd, polwaves[i])
       p.plot(polar[:,c].ravel(), label=str(polwaves[i])+" $\mu$m" )

    prop = matplotlib.font_manager.FontProperties(size=10)
    p.legend(prop=prop)

    p.draw()

####################### Main #####################
logging.basicConfig(level=logging.INFO)

if __name__=='__main__':
    #par = readparfile('fomalhaut.para')
    #par['distance']  *= 20
    #parfile_set_pixscale(par, 0.014)
    #writeparfile('fom_distant.par', par)
    logging.basicConfig(level=logging.INFO)
    pass



