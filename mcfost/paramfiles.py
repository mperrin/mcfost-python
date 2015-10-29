import os
import numpy as np
import matplotlib
import matplotlib.pyplot as pl
import logging
import glob
import astropy, astropy.io.ascii
_log = logging.getLogger('mcfost')

# this lets you put "stop()" in your code to have a debugger breakpoint
from IPython.core.debugger import Tracer; stop = Tracer()

# some extremely simple classes to serve as structs.


class Paramfile(object):
    """ Object class interface to MCFOST parameter files 


    Example: 

    par = Parfile('somefile.par')


    """

    _minimum_version = 2.15  # minimum MCFOST version for this code to run 
    def __init__(self, filename=None, directory="./", **kwargs):

        # we jump through some hurdles here to allow the
        # user to provide either a filename OR a directory
        # as the first positional argument, and attempt to do
        # the right thing in either case. 
        if filename is not None:
            if os.path.isdir(filename): 
                dir=filename
                filename=None
            dir = os.path.dirname(filename)
        if filename is None and directory is not None:
            filename = find_paramfile(directory=directory)
            if filename is None:
                raise IOError("Could not find any MCFOST parameter file in directory={0}".format(directory))

        self.filename= filename
        self._directory=directory
        self._readparfile(**kwargs)

#    def __getattr__(self, key):
#        stop()
#        return self._dict[key]

    def __getitem__(self, key):
        # enable dict-like access as well, for convenience
        return self.__dict__[key]

    def _readparfile(self, silent=False,return_text=False,verbose=False):
        """
        Read in an MCFOST par file into a dictionary.

        Note that some elements of that dict are in turn smaller sub-dicts. This nested structure
        is required because some of the parameter structures can repeat multiple times in complex models.


        WARNING: Not all parameters are read in now, just the most useful ones. Work in progress.

        Variable names are typically consistent with the names given in the documentation or the
        comments in the reference parameter files.

        Parameters
        ----------
        return_text : bool
            Return the entire parfile text as an entry in the dict
        """

        #################
        #
        # Warning - the following is pretty ugly code in many places, due to the
        #   need to accomodate the changes of the input file format over time.
        #   Code has lots of hacks and workarounds added ad hoc as needed. 
        #   Could use additional refactoring and cleanup...
        #
        #################

   
        text = open(self.filename,'r').readlines()

        #--- First we define a variety of utility functions to be used below --
        # utility functions for pulling out specific items from a line

#        def set1part(key,linenum, partnum, typ):
#            """Parse out one item and set it as an attribute on this object. """
#            if verbose: _log.info( "Looking for {key} on line {linenum} item {partnum} in '{text}'".format(key=key,linenum=linenum, partnum=partnum, type=typ, text=text[linenum].rstrip()))
#            try:
#                if typ == bool:
#                    # Python doesn't automatically cast T and F to True and False so implement
#                    # that manually here. 
#                    val = str(text[linenum].split(None,partnum+1)[partnum-1])
#                    self.__dict__[key] = val == 'T'
#                else:
#                    self.__dict__[key] = typ(text[linenum].split(None,partnum+1)[partnum-1])
#                if verbose: _log.info("    found: {0} = {1}".format(key, self.__dict__[key]))
#            except:
#                raise IOError("Could not parse line %d:'%s' item %d for '%s'" % (linenum, text[linenum][:-1], partnum, key))
#
        def set1partOfDict(some_dict, key,linenum, partnum, typ):
            """ Parse out one item from the text file, cast its type, and
            save it into a dict"""
            if verbose: _log.info( "Looking for {key} on line {linenum} item {partnum} in '{text}'".format(key=key,linenum=linenum, partnum=partnum, type=typ, text=text[linenum].rstrip()))
            #some_dict[key] = typ(text[linenum].split(None,partnum+1)[partnum-1])
            try:
                if typ == bool:
                    # Python doesn't automatically cast T and F to True and False so implement
                    # that manually here. 
                    val = str(text[linenum].split(None,partnum+1)[partnum-1])
                    some_dict[key] = val == 'T'
                    #_log.info(" for {0} found {1} which is {2}".format(key, val, some_dict[key]))
                else:
                    some_dict[key] = typ(text[linenum].split(None,partnum+1)[partnum-1])
                if verbose: _log.info("    found: {0} = {1}".format(key, some_dict[key]))
            except:
                raise IOError("Could not parse line %d:'%s' item %d for '%s'" % (linenum, text[linenum][:-1], partnum, key))
        def set1part(key,linenum, partnum, typ):
            """ Same as above but a shorthand for using this own object's __dict__"""
            return set1partOfDict(self.__dict__, key, linenum, partnum, typ)

        def dict2recarray(somedict):
            """ Convert dict to recarray while preserving types and names """
            dt = []
            vals = []    
            for item in somedict.items():
                deftype =  type(item[1])
                # force longer string field because the default behavior chops it to 1 character?!?
                if deftype == str: deftype='S80' 
                dt.append( (item[0], deftype) )

                vals.append(item[1])
            mydtype = np.dtype(dt)
            return np.rec.array(tuple(vals),  dtype=mydtype)



        self.fulltext = text

        set1part('version',0,1,float)

        _log.debug("Parsing MCFOST parameter file version %.2f " % self.version)

        if self.version < self._minimum_version: raise Exception('Parameter file version must be at least {ver:.2f}'.format(ver=self._minimum_version))

        #-- Number of photon packages --
        lineptr = 3 
        if (float(self.version) < 2.15): 
            set1part( 'nbr_parallel_loop', lineptr, 1, int) ; lineptr+=1
        else:
            self.nbr_parallel_loop = 1
        set1part('nbr_photons_eq_th', lineptr, 1, float) ; lineptr+=1
        set1part('nbr_photons_lambda',lineptr, 1, float) ; lineptr+=1
        set1part('nbr_photons_image', lineptr, 1, float) ; lineptr+=3

        #--  Wavelength --
        set1part('nwavelengths', lineptr, 1, int)
        set1part('lambda_min', lineptr, 2, float)
        set1part('lambda_max', lineptr, 3, float) ; lineptr+=1
        set1part('l_temp', lineptr, 1, bool) # should be bool type?
        set1part('l_sed',  lineptr, 2, bool)
        set1part('l_complete', lineptr, 3, bool) ; lineptr+=1
        set1part('wavelengths_file', lineptr, 1, str) ; lineptr+=1
        set1part('l_separate', lineptr, 1, bool) 
        set1part('l_stokes', lineptr, 2, bool)  ; lineptr+=3

        # compute or look up wavelength solution
        if self.l_complete:
            # use log sampled wavelength range
            wavelengths_inc=np.exp( np.log(self.lambda_max/self.lambda_min)/(self.nwavelengths) )
            self.wavelengths = self.lambda_min * wavelengths_inc**(np.arange(self.nwavelengths)+0.5)
        else:
            # load user-specified wavelength range from file on disk
            possible_wavelengths_files = [ 
                os.path.join(self._directory, self.wavelengths_file), 
                os.path.join(self._directory,"data_th",self.wavelengths_file), 
                #os.path.join(os.getenv('MY_MCFOST_UTILS'),'Lambda',self.wavelengths_file) ,
                os.path.join(os.getenv('MCFOST_UTILS'),'Lambda',self.wavelengths_file) ]
            wavelengths_file = None
            for possible_name in possible_wavelengths_files:
                _log.debug("Checking for wavelengths file at "+possible_name)
                if os.path.exists(possible_name):
                    wavelengths_file=possible_name
                    _log.debug("Found wavelengths file at "+wavelengths_file)
                    break

            if not os.path.exists(wavelengths_file):
                raise IOError('Cannot find requested wavelength file: '+wavelengths_file)
            self.wavelengths = astropy.io.ascii.read(wavelengths_file, data_start=0,names=['wavelength'])['wavelength'] 
 
        #-- Grid geometry and size --
        set1part('grid_type', lineptr, 1, int) ; lineptr+=1
        set1part('grid_n_rad', lineptr, 1, int)
        set1part('grid_nz', lineptr, 2, int) # should be n_theta if type=spherical
        set1part('grid_n_az', lineptr, 3, int)
        set1part('grid_n_rad_in', lineptr, 4, int) ; lineptr+=3

        #--  Maps (Images) --
        set1part('im_nx',lineptr , 1 if float(self.version) >= 2.15 else 3, int),
        set1part('im_ny',lineptr , 2 if float(self.version) >= 2.15 else 4, int),
        if float(self.version) >= 2.15:
            set1part('im_map_size', lineptr , 3, float)
        lineptr+=1
        set1part('MC_n_incl',lineptr , 1, int if self.version >= 2.17 else float),
        set1part('MC_n_az',lineptr , 2, int) ; lineptr+=1
        set1part('RT_imin', lineptr, 1, float)
        set1part('RT_imax', lineptr, 2, float)
        set1part('RT_n_incl',lineptr , 3, int)           
        set1part('RT_centered',lineptr , 4, bool) ; lineptr+=1
        set1part('distance', lineptr, 1, float) ; lineptr+=1
        if float(self.version) > 2.09:
            set1part('disk_pa', lineptr, 1, float) ; lineptr+=1  # only for > v2.09
        else:
            self.disk_pa = 0
        lineptr+=2


        # compute inclinations used for SED
        # FIXME this needs updating for SED RT mode
        # Ignore the priore complication here - assume we always have RT SEDs so the inclinations match?

        #n_thet = self.grid_ntheta
        #inc = np.arange(n_thet)+1
        #incmin    = np.arccos(1.-((inc-1  )*1.0)/n_thet)*180./3.1415926
        #incmax    = np.arccos(1.-((inc    )*1.0)/n_thet)*180./3.1415926
        #incmiddle = np.arccos(1.-((inc-0.5)*1.0)/n_thet)*180./3.1415926
        #self.sed_inclinations = incmiddle


        # compute inclinations used for RT
        ##    ;How to calculate the inclinations depends on whether you're using the values centered in those bins or not.
        dtor = np.pi/180
        if self.RT_centered:
            # find the average cos(i) in each bin (averaged linearly in cosine space)
            a1 = np.linspace(np.cos(self.RT_imax*dtor), np.cos(self.RT_imin*dtor), self.RT_n_incl+1)
            a2 = (a1[0:-1]+a1[1:])/2
            self.inclinations = (np.arccos(a2)/dtor)[::-1]
        else:
            # just use the bottom value for each bin
            a1 = np.linspace(np.cos(self.RT_imax*dtor), np.cos(self.RT_imin*dtor), self.RT_n_incl)
            self.inclinations = (np.arccos(a1)/dtor)[::-1]


        #-- Scattering Method --
        set1part('scattering_method', lineptr, 1, int) ; lineptr+=1
        set1part('scattering_mie_hg', lineptr, 1, int) ; lineptr+=3
 
        #-- Symmetries --
        set1part('l_image_symmetry', lineptr, 1, bool) ; lineptr+=1
        set1part('l_central_symmetry', lineptr, 1, bool) ; lineptr+=1
        set1part('l_axial_symmetry', lineptr, 1, bool) ; lineptr+=3

        #-- Disk physics / dust global properties --
        if float(self.version) >2.15:
            set1part('dust_settling', lineptr, 1, int) 
        else:
            # convert old T/F setting to current integer setting
            set1part('dust_settling', lineptr, 1, bool)
            self.dust_settling= 1 if self.dust_settling else 2

        set1part('settling_exp_strat' ,lineptr, 2,float) 
        set1part('settling_a_strat' ,lineptr, 3,float) ; lineptr +=1
        if float(self.version) >=2.19:
            set1part('l_radial_migration', lineptr, 1, bool) ; lineptr+=1
        else:
            self.l_radial_migration = False
        set1part('l_sublimate_dust' ,lineptr, 1, bool) ; lineptr+=1
        if float(self.version) >=2.19:
            set1part('l_hydrostatic_equilibrium', lineptr, 1, bool) ; lineptr+=1
        else:
            self.l_hydrostatic_equilibrium = False
        set1part('l_viscous_heating' ,lineptr, 1, bool) 
        set1part('alpha_viscosity' ,lineptr, 2, float) ; lineptr+=1
        lineptr+=2


        #-- Number of Zones --
        # read in the different zones
        set1part( 'nzones' ,lineptr, 1,int) ; lineptr+=3


#        if float(self.version) >=2.15:
#            set1part( 'im_size_au', 22+d2, 3, float),
#            if float(self.version) < 2.19: 
#                set1part( 'im_zoom', 22+d2, 4, float)
#            else:
#                self.__dict__['im_zoom'] = 1.0 
#            i=39  # points to nzones parameter # note: advance this counter ** after** you have read something in...
#        else:
#            set1part( 'gas_to_dust', 38+d2, 1, float),
#            set1part( 't_start' ,42+d3, 1,float),
#            set1part( 'sublimtemp' ,42+d3, 2,float),
#            i=43  # points to nzones parameter # note: advance this counter ** after** you have read something in...
#


        #-- Density Structure (1 per zone) --
        self.density_zones=[]
        for idensity in range(self.nzones):
            density={}
            set1partOfDict(density, 'zone_type', lineptr, 1, int) ; lineptr+=1
            set1partOfDict(density, 'dust_mass', lineptr, 1, float)
            if self.version >= 2.15:
                set1partOfDict(density, 'gas_to_dust_ratio', lineptr, 2, float)
            lineptr+=1

            set1partOfDict(density, 'scale_height', lineptr, 1, float)
            set1partOfDict(density, 'reference_radius', lineptr, 2, float)
            if self.version >= 2.19:
                set1partOfDict(density, 'debris_disk_vertical_profile_exponent', lineptr, 3, float)
            else:
                density['debris_disk_vertical_profile_exponent'] = 2
            lineptr+=1
            set1partOfDict(density, 'r_in', lineptr, 1, float),
            set1partOfDict(density, 'edge', lineptr, 2 if self.version >=2.19 else (4 if self.version <2.15 else 3), float),
            set1partOfDict(density, 'r_out', lineptr, 3 if self.version >=2.19 else 2, float)
            if self.version < 2.15:
                set1partOfDict(density, 'im_map_size', lineptr, 3, float), # used to be size_neb here
            if self.version >= 2.19:
                set1partOfDict(density, 'r_critical', lineptr, 4, float) # only used for tapered edge
            else:
                density['r_critical'] = density['r_out']
            lineptr+=1

            set1partOfDict(density, 'flaring_exp', lineptr, 1, float) ; lineptr+=1
            set1partOfDict(density, 'surface_density_exp', lineptr, 1, float) 
            if self.version > 2.15:
                set1partOfDict(density, 'gamma_exp', lineptr, 2, float)
            else:
                density['gamma_exp'] = 0.0
            lineptr+=3

            self.density_zones.append(density)

        #-- Cavity --
        set1part('cavity_flag', lineptr, 1, str) ; lineptr+=1
        set1part('cavity_height', lineptr, 1, float)
        set1part('cavity_ref_radius', lineptr, 2, float) ; lineptr+=1
        set1part('cavity_flaring_exp', lineptr, 1, float) ; lineptr+=3


        #-- Grain properties --
        # read in the dust grain properties. One set of dust props **per zone**, each of which can contain multiple species
        # These are each stored as a list of dicts per each zone.
        for idensity in range(self.nzones):
            #self.density_zones['dust_nspecies']=[]
            set1partOfDict(self.density_zones[idensity], 'dust_nspecies', lineptr, 1,int)
            lineptr+=1
            self.density_zones[idensity]['dust']=[]
            for idust in range(self.density_zones[idensity]['dust_nspecies']):
                dust={}
                if float(self.version) >= 2.12:
                    # version 2.12 or higher, allow for multi-component grains
                    if self.version >= 2.17:
                        set1partOfDict(self.density_zones[idensity],'grain_type', lineptr, 1, str)

                    set1partOfDict(dust, 'ncomponents',  lineptr, 2 if self.version >=2.17 else 1, int),
                    set1partOfDict(dust, 'mixing_rule',  lineptr, 3 if self.version >=2.17 else 2, int),
                    set1partOfDict(dust, 'porosity',     lineptr, 4 if self.version >=2.17 else 3, float),
                    set1partOfDict(dust, 'mass_fraction',lineptr, 5 if self.version >=2.17 else 4, float)
                    lineptr+=1
                    #if dust['ncomponents'] >1: raise NotImplementedError("Need multi-component parsing code!")
                    for icomponent in range(dust['ncomponents']):
                        dust_keys = (('filename', lineptr, 1, str),
                            ('volume_fraction', lineptr, 2, float))
                        for key, line, item, typ in dust_keys:
                            set1partOfDict(dust,key, line, item, typ)

                        lineptr+=1
                    # now the heating and grain size properties
                    dust_keys = ( ('heating', lineptr+0, 1, int),
                            ('amin', lineptr+1, 1, float),
                            ('amax', lineptr+1, 2, float),
                            ('aexp', lineptr+1, 3, float),
                            ('ngrains', lineptr+1, 4, int))
                    for key, line, item, typ in dust_keys:
                        set1partOfDict(dust,key, line, item, typ)
                    lineptr+=2

                else:
                    # earlier versions than 2.12, so only one component allowed. 
                    dust_keys = [('filename', i+0, 1, str),
                            ('porosity', i+0, 2, float),
                            ('mass_fraction', i+0, 3, float),
                            ('heating', i+1, 1, int),
                            ('amin', i+2, 1, float),
                            ('amax', i+2, 2, float),
                            ('aexp', i+2, 3, float),
                            ('n_grains', i+2, 4, int)]
                    for key, line, item, typ in dust_keys:
                        set1partOfDict(dust,key, line, item, typ)
                    i+=3
                #add any missing keywords using defaults. This is to handle
                # parsing earlier versions of the parameter file that lacked some settings
                defaults = (('volume_fraction',1.0), ('grain_type', 'Mie'))
                for defkey, defval in defaults:
                    if defkey not in dust.keys() : dust[defkey] = defval

                self.density_zones[idensity]['dust'].append(dust)
        lineptr+=2

        # molecular RT settings
        set1part('l_pop', lineptr, 1, bool) 
        set1part('l_accurate_pop', lineptr, 2, bool)
        set1part('l_LTE', lineptr, 3, bool) 
        set1part('molecular_profile_width', lineptr, 4, float)  ; lineptr+=1
        set1part('molecular_v_turb', lineptr, 1, float)  ; lineptr+=1
        set1part('molecular_nmol', lineptr, 1, float)  ; lineptr+=1
        set1part('molecular_filename', lineptr, 1, str)  
        set1part('molecular_level_max', lineptr, 2, float)  ; lineptr+=1
        set1part('molecular_vmax', lineptr, 1, float)  
        set1part('molecular_n_speed', lineptr, 2, float)  ; lineptr+=1
        set1part('l_molecular_abundance', lineptr, 1, bool)  
        set1part('molecular_abundance', lineptr, 2, float)  
        set1part('molecular_abundance_file', lineptr, 3, str) ; lineptr+=1
        set1part('l_molecular_raytrace', lineptr, 1, bool)  
        set1part('molecular_raytrace_n_lines', lineptr, 2, float); lineptr+=1
        lineptr+=1 # skip transition numbers for now

        lineptr+=2

        #-- Star propertoes --
        set1part('nstar', lineptr, 1,int) ; lineptr +=1
        self.stars = []
        for istar in range(self.nstar):
            star_dict={}
            set1partOfDict(star_dict, 'temp', lineptr, 1, float) 
            set1partOfDict(star_dict, 'radius', lineptr, 2, float)
            set1partOfDict(star_dict, 'mass', lineptr, 3, float)
            set1partOfDict(star_dict, 'x', lineptr, 4, float)
            set1partOfDict(star_dict, 'y', lineptr, 5, float)
            set1partOfDict(star_dict, 'z', lineptr, 6, float) 
            set1partOfDict(star_dict, 'l_is_blackbody', lineptr, 7, bool) ; lineptr +=1
            set1partOfDict(star_dict, 'spectrum', lineptr, 1, str) ; lineptr +=1
            if float(self.version) >= 2.11:
                set1partOfDict(star_dict, 'fUV', lineptr, 1, float) 
                set1partOfDict(star_dict, 'slope_fUV', lineptr, 2, float)  ; lineptr+=1
            # convert star from dict into recarray
            self.stars.append(dict2recarray(star_dict))

        #---- Command line options ----
        # now try to read in command line options which can override
        # some of the settings.
        # set some defaults first.
        self.im_zoom =1.0
        self.im_raytraced = False
        if len(text) > lineptr:
            optionsline = text[lineptr]
            options = optionsline.split()
            for j in np.arange(len(options)):
                if options[j] == "-img":
                    self.im_wavelength = options[j+1]
                elif options[j] == "-resol":
                    self.im_nx = np.int(options[j+1])
                    self.im_ny = np.int(options[j+2])
                elif  options[j] == "-zoom":
                    self.im_zoom = np.int(options[j+1])
                elif  options[j] == "-rt":
                    self.im_raytraced = True
        else:
            if verbose: print "could not read in command line options from MCFOST; using default grid settings"

        #--- Derived quantities ---

        # compute a few extra things for convenience
        if self.im_nx == self.im_ny:
            self.im_xsize_au =  (2*self.im_map_size)
            self.im_ysize_au =  (2*self.im_map_size) 
        else: 
            # FIXME update if X and Y pixel sizes differ
            raise NotImplementedError('Need to implement size calc for unequal image axes')


        #self.im_yscale_au =  (2*self.density_zones[0]['size_neb'] / self.im_zoom) / self.im_ny  # au per pixel
        self.im_xscale_au = self.im_xsize_au / self.im_nx  # au per pixel
        self.im_yscale_au = self.im_ysize_au / self.im_ny  # au per pixel

        self.im_xsize_arcsec  = self.im_xsize_au  / self.distance
        self.im_ysize_arcsec  = self.im_ysize_au  / self.distance
        self.im_xscale_arcsec = self.im_xscale_au / self.distance
        self.im_yscale_arcsec = self.im_yscale_au / self.distance

        # done reading in parameter file


    # convenience access properties for the case of there only being one of something
    @property
    def star(self):
        return self.stars[0]

    @property
    def density(self):
        return self.density_zones[0]

    @property
    def dust(self):
        return self.density_zones[0]['dust'][0]

    def __str__(self):
        """ Return a nicely formatted text parameter file. Currently returns v2.17 format

        HISTORY
        --------
        2013-01-05 updated to version 2.15
        2013-04-24 Substantial code rewrites & cleanup. Updated to version 2.17
        2014-01-12 Updated to version 2.19

        """
        #par = self._dict

        #getkeys = lambda l: tuple([par[li] for li in l])

        # if we read in a <2.14 version file, we have a nbr_parallel_loop factor that
        # was multiplicatively merged in with the other nbrs for later versions, so
        # do that multiplication here before writing out a >= 2.15 parfile
        try:
            scale_nbr = self.nbr_parallel_loop
        except:
            scale_nbr = 1

        # Utility function to work around issue with formatting numpy recarrays
        # see http://mail.scipy.org/pipermail/numpy-discussion/2013-June/066796.html
        def recarray2dict(somerecarray):
            mydict = {}
            #print somerecarray.dtype.names
            for name, typecode in somerecarray.dtype.descr:
                cast = float if 'f' in typecode else str
                mydict[name] = cast(somerecarray[name])
            return mydict

        # string formatting helper
        def bstr(somebool):
            return "T" if somebool else "F"


        template = """2.19                      mcfost version

#Number of photon packages
  {self.nbr_photons_eq_th:<10.5g}              nbr_photons_eq_th  : T computation
  {self.nbr_photons_lambda:<10.5g}              nbr_photons_lambda : SED computation
  {self.nbr_photons_image:<10.5g}              nbr_photons_image : images computation

#Wavelength
  {self.nwavelengths:<3d} {self.lambda_min:<5.1f} {self.lambda_max:<7g}       n_lambda, lambda_min, lambda_max [microns]
  {str_l_temp:1s} {str_l_sed:1s} {str_l_complete:1s}                   compute temperature?, compute sed?, use default wavelength grid ?
  {self.wavelengths_file}           wavelength file (if previous parameter is F)
  {str_l_separate:1s} {str_l_stokes:1s}                     separation of different contributions?, stokes parameters?

#Grid geometry and size
  {self.grid_type:>1d}                       1 = cylindrical, 2 = spherical, 3 = Voronoi tesselation (this is in beta, please ask Christophe)
  {self.grid_n_rad:3d} {self.grid_nz:3d} {self.grid_n_az:2d} {self.grid_n_rad_in}           n_rad (log distribution), nz (or n_theta), n_az, n_rad_in

#Maps
  {self.im_nx:<3d} {self.im_ny:3d} {self.im_map_size:5.1f}        grid (nx,ny), size [AU]
  {self.MC_n_incl} {self.MC_n_az}               MC : N_bin_incl, N_bin_az
  {self.RT_imin:<4.1f}  {self.RT_imax:<4.1f}  {self.RT_n_incl:>2d} {str_RT_centered}    RT: imin, imax, n_incl, centered ?
  {self.distance:<6.2f}                 distance (pc)
  {self.disk_pa:<6.2f}                  disk PA
  """.format(self=self, 
          str_l_temp=bstr(self.l_temp), # is there a better way to get a Bool to print as a single letter?
          str_l_sed=bstr(self.l_sed),
          str_l_complete=bstr(self.l_complete),
          str_l_separate=bstr(self.l_separate),
          str_l_stokes=bstr(self.l_stokes),
          str_RT_centered=bstr(self.RT_centered),
          )


    # FIXME - remove ntheta, nphi, add in size-AU (which is the 'size_neb') variable now)
    # update MC nbin line in reader
        template+="""
#Scattering method
  {self.scattering_method:1d}                      0=auto, 1=grain prop, 2=cell prop
  {self.scattering_mie_hg:1d}                      1=Mie, 2=hg (2 implies the loss of polarizarion) 

#Symmetries
  {str_image_symmetry:1s}                      image symmetry
  {str_central_symmetry:1s}                      central symmetry
  {str_axial_symmetry:1s}                      axial symmetry (important only if N_phi > 1)

#Disk physics
  {self.dust_settling:1d}  {self.settling_exp_strat:<6.2f}  {self.settling_a_strat:<6.2f}      dust_settling (0=no settling, 1=parametric, 2=Dubrulle, 3=Fromang), exp_strat, a_strat (for parametric settling)
  {str_radial_migration}                       dust radial migration
  {str_sublimate_dust}                      sublimate dust
  {str_hydrostatic_equilibrium}                       hydrostatic equilibrium
  {str_viscous_heating}  {self.alpha_viscosity:4.1g}                 viscous heating, alpha_viscosity

#Number of zones : 1 zone = 1 density structure + corresponding grain properties
  {self.nzones}
      """.format(self=self,
          str_image_symmetry=bstr(self.l_image_symmetry),
          str_central_symmetry=bstr(self.l_central_symmetry),
          str_axial_symmetry=bstr(self.l_axial_symmetry),
          str_sublimate_dust=bstr(self.l_sublimate_dust),
          str_radial_migration=bstr(self.l_radial_migration),
          str_hydrostatic_equilibrium=bstr(self.l_hydrostatic_equilibrium),
          str_viscous_heating=bstr(self.l_viscous_heating),
          )

        if self.nzones > 1:
            raise NotImplemented("write par file code does not yet support multiple zones!")

        getsubkeys = lambda k, l: tuple([par[k][li] for li in l])
        template+="""
#Density structure
  {zone[zone_type]:1d}                       zone type : 1 = disk, 2 = tapered-edge disk, 3 = envelope, 4 = debris disk, 5 = wall
  {zone[dust_mass]:<10.2e} {zone[gas_to_dust_ratio]:<5.1f}        dust mass,  gas-to-dust mass ratio
  {zone[scale_height]:<5.1f}  {zone[reference_radius]:<6.1f} {zone[debris_disk_vertical_profile_exponent]:<6.1f}          scale height, reference radius (AU), unused for envelope, vertical profile exponent (only for debris disk)
  {zone[r_in]:<6.1f}  {zone[edge]:<6.1f} {zone[r_out]:<6.1f} {zone[r_critical]:<6.1f}  Rin, edge, Rout, Rc (AU) Rc is only used for tappered-edge & debris disks (Rout set to 8*Rc if Rout==0)
  {zone[flaring_exp]:<8.3f}                  flaring exponent, unused for envelope
  {zone[surface_density_exp]:<6.3f} {zone[gamma_exp]:<6.3f}                surface density exponent (or -gamma for tappered-edge disk or volume density for envelope), usually < 0, -gamma_exp (or alpha_in & alpha_out for debris disk)
        """.format(zone = self.density_zones[0])

        template+="""
#Cavity : everything is empty above the surface
  {self.cavity_flag}                        cavity ?
  {self.cavity_height:<5.1f}  {self.cavity_ref_radius:<5.1f}             height, reference radius (AU)
  {self.cavity_flaring_exp:<5.1f}                    flaring exponent
        """.format(self=self) #% getkeys(['cavity_flag', 'cavity_height', 'cavity_ref_radius', 'cavity_flaring_exp'])

        if len(self.density_zones[0]['dust']) > 1: raise NotImplemented("No support for multi dust types yet")
        template+="""
#Grain properties
  {zone[dust_nspecies]:<2d}                      Number of species""".format(zone=self.density_zones[0])

        template+="""
  Mie {dust[ncomponents]:<2d} {dust[mixing_rule]:<1d} {dust[porosity]:<5.2f} {dust[mass_fraction]:<5.2f} 0.9     Grain type (Mie or DHS), N_components, mixing rule (1 = EMT or 2 = coating),  porosity, mass fraction, Vmax (for DHS)
  {dust[filename]:s}  {dust[volume_fraction]:<4.1f}   Optical indices file, volume fraction
  {dust[heating]:<2d}                      Heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE
  {dust[amin]:<6.3f} {dust[amax]:6.1f}  {dust[aexp]:4.1f} {dust[ngrains]:3d}   amin, amax, aexp, nbr_grains
      """.format(dust = self.density_zones[0]['dust'][0])
      
        template+=""" 
#Molecular RT settings
  T T T 15.              lpop, laccurate_pop, LTE, profile width
  0.2                    v_turb (delta)
  1                      nmol
  co@xpol.dat 6          molecular data filename, level_max
  1.0 50                 vmax (km.s-1), n_speed
  T 1.e-6 abundance.fits.gz   cst molecule abundance ?, abundance, abundance file
  T  3                   ray tracing ?,  number of lines in ray-tracing
  1 2 3                  transition numbers
      """ 

        template+="""
#Star properties
  1 Number of stars"""
        for star in self.stars:
            template+="""
  {star[temp]:<7.1f} {star[radius]:5.1f} {star[mass]:5.1f} {star[x]:5.1f} {star[y]:5.1f} {star[z]:5.1f} {str_is_blackbody:1s}       Temp, radius (solar radius),M (solar mass),x,y,z (AU), is a blackbody?
  {star[spectrum]:s}
  {star[fUV]:<5.1f}  {star[slope_fUV]:<5.1f}     fUV, slope_fUV """.format(star=recarray2dict(star), str_is_blackbody = 'T' if star['l_is_blackbody'] else 'F')
        
        # % starkeys(['temp', 'radius', 'mass', 'x', 'y', 'z', 'spectrum'])
        return template


    def writeto(self, outname='sample.par'):
        """ Write an MCFOST parameter file to disk. Currently outputs v2.17 param files


        WARNING: Not all parameters are allowed to vary, just the most useful ones. 
        Lots of the basics are still hard coded. 
        TBD fix later.


        Parameters
        ----------
        outname : string
            Output file name

        """


        outfile = open(outname, 'w')
        outfile.write( str(self))
        outfile.close()
        print "  ==>> "+outname

    def set_parameter(self, paramname, value):
        """ Helper function for parameter setting. 

        For many parameters this is a trivial one-line wrapper. 
        But for some of the nested structures inside e.g. the density or
        dust properties structs, it's convenient to define some
        shortcuts for use in iterating over grids or MCMC/GA ensembles. 
        This implements those shortcuts. 

        Note: Shortcut names are consistent with those defined in the
        IDL MCRE re_translationtable.txt configuration file.

        """
        
        if paramname == 'm_star':
            self.star['mass'] = value
        elif paramname == 't_star':
            self.star['temp'] = value
        elif paramname == 'r_star':
            self.star['radius'] = value
        elif paramname == 'dust_settling':
            self['dust_settling'] = value
        # When iterating over properties of disk geometry, 
        # generally we mean to do so over the first density zone.
        elif paramname == 'dustmass' or paramname == 'dust_mass':
            self.density_zones[0]['dust_mass'] = value
        elif paramname == 'r_in':
            self.density_zones[0]['r_in'] = value
        elif paramname == 'r_out':
            self.density_zones[0]['r_out'] = value
        elif paramname == 'flaring':
            self.density_zones[0]['flaring_exp'] = value
        elif paramname == 'surface_density':
            self.density_zones[0]['surface_density_exp'] = value
        elif paramname == 'gamma_exp':
            self.density_zones[0]['gamma_exp']=value
        elif paramname == 'scaleheight' or paramname == 'scale_height':
            self.density_zones[0]['scale_height'] = value
        elif paramname == 'zone_type':
            self.density_zones[0]['zone_type'] = value
        # likewise, almost always when iterating over dust properties we want to
        # iterate over the first dust component of the first zone. 
        # If you want to do something more complicated, you'll have to
        # modify this or write your own code...
        elif paramname == 'dust_exponent':
            self.density_zones[0]['dust'][0]['aexp'] = value
        elif paramname == 'dust_amax':
            self.density_zones[0]['dust'][0]['amax'] = value
        elif paramname == 'dust_amin':
            self.density_zones[0]['dust'][0]['amin'] = value
        elif paramname == 'dust_species':
            self.density_zones[0]['dust'][0]['filename'] = value
        elif paramname == 'dust_porosity':
            self.density_zones[0]['dust'][0]['porosity'] = value
        else:
            try: 
                self[paramname] = value
            except:
                raise ValueError("Don't know how to set a parameter named '{0}'".format(paramname))

def find_paramfile(directory="./",  parfile=None, verbose=False, wavelength=None):
    """ Find a MCFOST par file in a specified directory

    By default, look in the current directory of a model,
    but if a wavelength is specified, look inside the appropriate data_## directory.

    This looks for any file named 'something.par' or 'something.para'.

    KEYWORDS:
    wavelength : string
        string with wavelength name
    verbose :Boolean.
        whether to be more verbose in output
    directory : string
        directory to look in. Default is current directory

    """

    if parfile is not None:
        # TODO validate its existence, check if path name relative to dir?
        output = parfile
    else:
        _log.debug("Looking for par files in dir: "+directory)
        directory = os.path.expanduser(directory)
        if wavelength is not None:
            directory = os.path.join(directory, "data_%s" % wavelength)
            _log.info("Since wavelength is %s, looking in %s subdir" % (wavelength,dir))
        l = glob.glob(os.path.join(directory, "*.par"))
        l+= glob.glob(os.path.join(directory, "*.para"))
        if len(l) == 1:
            output = l[0]
        elif len(l) > 1:
            _log.warning("Multiple par files found... returning first: "+l[0])
            output = l[0]
        else:
            _log.error("No par files found!")
            output = None

    _log.debug('Par file found: '+str(output))
    return output


