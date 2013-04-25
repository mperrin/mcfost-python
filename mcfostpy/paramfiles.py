import os
import numpy as np
import matplotlib
import matplotlib.pyplot as pl
import logging
_log = logging.getLogger('mcfost')
logging.basicConfig(level=logging.INFO)

# this lets you put "stop()" in your code to have a debugger breakpoint
from IPython.core.debugger import Tracer; stop = Tracer()

_VERBOSE = False


# some extremely simple classes to serve as structs.

class Star():
    temp = 4000
    radius = 1.0
    mass = 1.0
    x = 0.0
    y = 0.0
    z = 0.0
    spectrum = ''
    fUV = 0.0
    flope_fuv = 0.0


class Paramfile(object):
    """ Object class interface to MCFOST parameter files 


    Example: 

    par = Parfile('somefile.par')




    """
    def __init__(self, filename=None, dir="./", **kwargs):

        # we jump through some hurdles here to allow the
        # user to provide either a filename OR a directory
        # as the first positional argument, and attempt to do
        # the right thing in either case. 
        if filename is not None:
            if os.path.isdir(filename): 
                dir=filename
                filename=None
            dir = os.path.dirname(filename)
        if filename is None and dir is not None:
            filename = get_current_paramfile(dir=dir)

        self.filename= filename
        self._readparfile(**kwargs)
        self._dir=dir

#    def __getattr__(self, key):
#        stop()
#        return self._dict[key]

    def __getitem__(self, key):
        # enable dict-like access as well, for convenience
        return self.__dict__[key]

    def _readparfile(self, silent=False,return_text=False,verbose=True):
        """
        Read in an MCFOST par file into a dictionary.

        Note that some elements of that dict are in turn smaller sub-dicts. This nested structure
        is required because some of the parameter structures can repeat multiple times in complex models.


        WARNING: Not all parameters are read in now, just the most useful ones. TBD fix later.

        Parameters
        ----------
        filename, dir: string
            input file and dir. Leave filename blank to search in dir.
        return_text : bool
            Return the entire parfile text as an entry in the dict
        """

        #################
        #
        # Warning - the following is pretty ugly code in many places, due to the
        #   need to accomodate the changes of the input file format over time.
        #   Code has lots of hacks and workarounds added ad hoc as needed. 
        #   Should probably be rewritten entirely with a cleaner version...
        #
        #################

   
        text = open(self.filename,'r').readlines()

        # utility functions for pulling out specific items from a line
        def getparts(text,linenum,nparts):
            return text[linenum].split(None,nparts)[0:nparts]
        def getfloats(*args):
            return np.asarray(getparts(*args), dtype=float)

        def get1part(linenum, partnum):
            return text[linenum].split(None,partnum+1)[partnum-1]
        # Parse out one item and set it as an attribute on this object.
        def set1part(key,linenum, partnum, typ):
            if verbose: _log.info( str( (key,linenum, partnum, typ,)))
            try:
                self.__dict__[key] = typ(text[linenum].split(None,partnum+1)[partnum-1])
                if verbose: print self.__dict__[key]
            except:
                raise IOError("Could not parse line %d:'%s' item %d for '%s'" % (linenum, text[linenum][:-1], partnum, key))


        # Ditto, but set it on some other dict
        def set1partOfDict(some_dict, key,linenum, partnum, typ):
            if verbose: _log.info( str( (key,linenum, partnum, typ,)))
            try:
                some_dict[key] = typ(text[linenum].split(None,partnum+1)[partnum-1])
                if verbose: print some_dict[key]
            except:
                raise IOError("Could not parse line %d:'%s' item %d for '%s'" % (linenum, text[linenum][:-1], partnum, key))


        self.fulltext = text

        set1part('version',0,1,float)

        if verbose or _VERBOSE: _log.info("Parsing MCFOST parameter file version %.2f " % self.version)

        # figure out line offsets of the fixed lookup table implemented below
        if float(self.version) >= 2.15:
            d1=-1   # removal of nbr_parallel_loop
            d2=-4   # and also 'consider disk emissions' flag
            d3=-5   # and also gas-to-dust parameter
        elif float(self.version) >= 2.11:
            d1=0
            d2= -2
            d3=d2
        elif float(self.version) >= 2.10:
            d1=0
            d2= -1
            d3=d2
        else:
            d1=0
            d2=0
            d3=d2
        #print "offset is ", d

        # TODO add distance and other image parameters


        # read in fixed-in-position parameters
        #  In the following table, the columns are:
        # Parameter_name   row   item_in_row   dtype
        fixed_keys=(
        ('nbr_photons_eq_th', 4+d1, 1, float),
        ('nbr_photons_lambda',5+d1, 1, float),
        ('nbr_photons_image', 6+d1, 1, float),
        ('nwavelengths', 9+d1, 1, int),
        ('wavelengths_min', 9+d1, 2, float),
        ('wavelengths_max', 9+d1, 3, float),
        ('ltemp_flag', 10+d1, 1, str),
        ('lsed_flag', 10+d1, 2, str),
        ('lcomplete_flag', 10+d1, 3, str),
        ('wavelengths_file', 11+d1, 1, str),
        ('grid_type', 18+d2, 1, int),
        ('grid_nrad', 19+d2, 1, int),
        ('grid_nz', 19+d2, 2, int),
        ('grid_n_az', 19+d2, 3, int),
        ('grid_n_rad_in', 19+d2, 4, int),
        #('im_ntheta', 22+d2, 1, int),         # for SED
        #('im_nphi',22+d2 , 2, int),
        ('im_nx',22+d2 , 3 if float(self.version) < 2.15 else 1, int),
        ('im_ny',22+d2 , 4 if float(self.version) < 2.15 else 2, int),
        ('im_mc_n_bin_incl',23+d2 , 1, int if self.version >= 2.17 else float),
        ('im_mc_n_bin_az',23+d2 , 2, int),
        ('im_rt_inc_min', 24+d2, 1, float),
        ('im_rt_inc_max', 24+d2, 2, float),
        ('im_rt_ninc',24+d2 , 3, int),           # for RT
        ('im_rt_centered',24+d2 , 4, str),
        ('distance', 25+d2, 1, float),
        ('disk_pa', 26+d2, 1, float),  # only for > v2.09
        ('scattering_method', 29+d2, 1, int), 
        ('scattering_mie_hg', 30+d2, 1, int), 
        ('settling_flag', 39+d3, 1, str),
        ('settling_exp' ,39+d3, 2,float),
        ('settling_a' ,39+d3, 3,float),
        ('sublimate_flag' ,40+d3, 1, str))


        for key, line, item, typ in fixed_keys:
            set1part(key, line, item, typ)

        if float(self.version) >=2.15:
            set1part( 'im_size_au', 22+d2, 3, float),
            set1part( 'im_zoom', 22+d2, 4, float),
            i=39  # points to nzones parameter # note: advance this counter ** after** you have read something in...
        else:
            set1part( 'nbr_parallel_loop', 3, 1, int),
            set1part( 'gas_to_dust', 38+d2, 1, float),
            set1part( 't_start' ,42+d3, 1,float),
            set1part( 'sublimtemp' ,42+d3, 2,float),
            i=43  # points to nzones parameter # note: advance this counter ** after** you have read something in...

        # convert the boolean keys from string to bool
        bool_keys = ['ltemp_flag', 'lsed_flag', 'lcomplete_flag', 'im_rt_centered']
        for k in bool_keys:
            self.__dict__[k] = (self.__dict__[k].upper() == 'T')

        # compute or look up wavelength solution
        if self.lcomplete_flag:
            # use log sampled wavelength range
            wavelengths_inc=np.exp( np.log(self.wavelengths_max)/self.wavelengths_min)/(self.nwavelengths)
            self.wavelengths = self.wavelengths_min * wavelengths_inc**(np.arange(self.nwavelengths)+0.5)
        else:
            # load user-specified wavelength range from file
            self.wavelengths_file = get1part(11,1)
            possible_wavelengths_files = [ 
                os.path.join(dir, self.wavelengths_file), 
                os.path.join(dir,"data_th",self.wavelengths_file), 
                os.path.join(os.getenv('MCFOST_UTILS'),'Lambda',self.wavelengths_file) ]
            wavelengths_file = None
            for possible_name in possible_wavelengths_files:
                if _VERBOSE: _log.info("Checking for wavelengths file at "+possible_name)
                if os.path.exists(possible_name):
                    wavelengths_file=possible_name
                    if _VERBOSE: _log.info("Found wavelengths file at "+wavelengths_file)
                    break

            if not os.path.exists(wavelengths_file):
                raise IOError('Cannot find requested wavelength file: '+wavelengths_file)
            self.wavelengths = asciitable.read(wavelengths_file, data_start=0,names=['wavelength'])['wavelength'] 
            

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
        if self.im_rt_centered:
            # find the average cos(i) in each bin (averaged linearly in cosine space)
            a1 = np.linspace(np.cos(self.im_inc_max*dtor), np.cos(self.im_inc_min*dtor), self.im_ninc+1)
            a2 = (a1[0:-1]+a1[1:])/2
            self.inclinations = (np.arccos(a2)/dtor)[::-1]
        else:
            # just use the bottom value for each bin
            a1 = np.linspace(np.cos(self.im_rt_inc_max*dtor), np.cos(self.im_rt_inc_min*dtor), self.im_rt_ninc)
            self.inclinations = (np.arccos(a1)/dtor)[::-1]


        # read in the different zones
        set1part( 'nzones' ,i, 1,int)
        i+=3

        self.density_zones=[]
        for idensity in range(self.nzones):
            density={}
            density_keys =[('zone_type', i+0, 1, int),
                    ('dust_mass', i+1, 1, float),
                    ('scale_height', i+2, 1, float),
                    ('ref_radius', i+2, 2, float),
                    ('r_in', i+3, 1, float),
                    ('r_out', i+3, 2, float),
                    ('size_neb', i+3, 3, float),
                    ('edge', i+3, 4 if self.version <2.15 else 3, float),
                    ('flaring_exp', i+4, 1, float),
                    ('surfdens_exp', i+5, 1, float)]
            if self.version >= 2.15:
                density_keys.append( ('gas_to_dust', i+1, 2, float))
            for key, line, item, typ in density_keys:
                set1partOfDict(density,key, line, item, typ)
            i+=6
            #add any missing keywords using defaults. This is to handle
            # parsing earlier versions of the parameter file that lacked some settings
            defaults = (('gamma_exp',0.0), )
            for defkey, defval in defaults:
                if defkey not in density.keys() : density[defkey] = defval

            self.density_zones.append(density)

        i+=2
        cavity_keys=(('cavity_flag', i+0, 1, str),
                    ('cavity_height', i+1, 1, float),
                    ('cavity_ref_radius', i+1, 2, float),
                    ('cavity_flaring_exp', i+2, 1, float))
        for key, line, item, typ in cavity_keys:
            set1part(key, line, item, typ)


        # read in the dust geometry. One set of dust props **per zone**, each of which can contain multiple species
        # These are each stored as a list of dicts per each zone.
        i+=5
        for idensity in range(self.nzones):
            #self.density_zones['dust_nspecies']=[]
            set1partOfDict(self.density_zones[idensity], 'dust_nspecies', i, 1,int)
            i+=1
            self.density_zones[idensity]['dust']=[]
            for idust in range(self.density_zones[idensity]['dust_nspecies']):
                dust={}
                if float(self.version) >= 2.12:
                    # version 2.12 or higher, allow for multi-component grains
                    dust_keys = [('ncomponents',  i+0, 1 if self.version < 2.17 else 2, int),
                                 ('mixing_rule',  i+0, 2 if self.version < 2.17 else 3, int),
                                 ('porosity',     i+0, 3 if self.version < 2.17 else 4, float),
                                 ('mass_fraction',i+0, 4 if self.version < 2.17 else 5, float)]
                    if self.version >= 2.17: dust_keys.append(('grain_type',i+0, 1, str))
                    for key, line, item, typ in dust_keys:
                        set1partOfDict(dust,key, line, item, typ)
                    i+=1
                    if dust['ncomponents'] >1: raise NotImplementedError("Need multi-component parsing code!")
                    for icomponent in range(dust['ncomponents']):
                        dust_keys = (('filename', i, 1, str),
                            ('volume_fraction', i, 2, float))
                        for key, line, item, typ in dust_keys:
                            set1partOfDict(dust,key, line, item, typ)

                        i+=1
                    # now the heating and grain size properties
                    dust_keys = ( ('heating', i+0, 1, int),
                            ('amin', i+1, 1, float),
                            ('amax', i+1, 2, float),
                            ('aexp', i+1, 3, float),
                            ('ngrains', i+1, 4, int))
                    for key, line, item, typ in dust_keys:
                        set1partOfDict(dust,key, line, item, typ)
                    i+=2

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


        i+= 12 # skip molecular RT settings

        set1part('nstar', i, 1,int)
        i+=1
        self.star = []
        for istar in range(self.nstar):
            nlines_star=3
            star_keys = [('temp', i+0, 1, float) ,
                    ('radius', i+0, 2, float),
                    ('mass', i+0, 3, float),
                    ('x', i+0, 4, float),
                    ('y', i+0, 5, float),
                    ('z', i+0, 6, float),
                    ('spectrum', i+1, 1, str)]
            if float(self.version) >= 2.11:
                star_keys.append( ('fUV', i+2, 1, float) )
                star_keys.append( ('slope_fUV', i+2, 2, float) )
                nlines_star+=1
            star={}
            for key, line, item, typ in star_keys:
                    set1partOfDict(star,key, line, item, typ)
            i+=nlines_star
            self.star.append(star)


        # now try to read in command line options which can override
        # some of the settings.
        # set some defaults first.
        self.im_zoom =1.0
        self.im_raytraced = False
        i+=1
        if len(text) >= i:
            optionsline = text[i]
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


        # compute a few extra things for convenience
        self.im_xsize_au =  (2*self.density_zones[0]['size_neb'] / self.im_zoom)
        self.im_ysize_au =  (2*self.density_zones[0]['size_neb'] / self.im_zoom)
        self.im_xsize_arcsec = self.im_xsize_au / self.distance
        self.im_ysize_arcsec = self.im_ysize_au / self.distance

        self.im_xscale_au =  (2*self.density_zones[0]['size_neb'] / self.im_zoom) / self.im_nx  # au per pixel
        self.im_yscale_au =  (2*self.density_zones[0]['size_neb'] / self.im_zoom) / self.im_ny  # au per pixel
        self.im_xscale_arcsec = self.im_xscale_au / self.distance
        self.im_yscale_arcsec = self.im_yscale_au / self.distance

        # done reading in parameter file



    def tostring(self):
        """ Return a nicely formatted text parameter file """
        #par = self._dict

        #getkeys = lambda l: tuple([par[li] for li in l])

        # if we read in a <2.14 version file, we have a nbr_parallel_loop factor that
        # was multiplicatively merged in with the other nbrs for later versions, so
        # do that multiplication here before writing out a >= 2.15 parfile
        try:
            scale_nbr = self.nbr_parallel_loop
        except:
            scale_nbr = 1

        # FIXME rewrite the following using new style string formatting!


        template = """2.17                      mcfost version

#Number of photon packages
  {self.nbr_photons_eq_th:<10.5g}              nbr_photons_eq_th  : T computation
  {self.nbr_photons_lambda:<10.5g}              nbr_photons_lambda : SED computation
  {self.nbr_photons_image:<10.5g}              nbr_photons_image : images computation

#Wavelength
  {self.nwavelengths:<3d} {self.wavelengths_min:<5.1f} {self.wavelengths_max:<7g}       n_lambda lambda_min lambda_max
  T T T                   compute temperature?, compute sed?, use default wavelength grid ?
  {self.wavelengths_file}           wavelength file (if previous parameter is F)
  F T                     separation of different contributions, stokes parameters

#Grid geometry and size
  {self.grid_type:>1d}                       1 = cylindrical, 2 = spherical
  {self.grid_nrad:3d} {self.grid_nz:3d} {self.grid_n_az:2d} {self.grid_n_rad_in}           n_rad, nz (or n_theta), n_az, n_rad_in

#Maps
  {self.im_nx:<3d} {self.im_ny:3d} {self.im_ny:<3d} {self.im_zoom:<4.1f}   grid (nx,ny), size [AU], zoom factor
  {self.im_mc_n_bin_incl} {self.im_mc_n_bin_az}               MC : N_bin_incl, N_bin_az, # of bin where MC is converged
  {self.im_rt_inc_min:<4.1f}  {self.im_rt_inc_max:<4.1f}  {self.im_rt_ninc:>2d} {self.im_rt_centered}  RT: imin, imax, n_incl, centered ?
  {self.distance:<6.2f}                 distance (pc)
  {self.disk_pa:<6.2f}                  disk PA

  """.format(self=self)


    # FIXME - remove ntheta, nphi, add in size-AU (which is the 'size_neb') variable now)
    # update MC nbin line in reader
        template+="""
#Scattering method
  {self.scattering_method:1d}                      0=auto, 1=grain prop, 2=cell prop
  {self.scattering_mie_hg:1d}                      1=Mie, 2=hg (2 implies the loss of polarizarion) 

#Symetries
  T                      image symmetry
  T                      central symmetry
  T                      axial symmetry (important only if N_phi > 1)

#Dust global properties
  {self.settling_flag}  {self.settling_exp:<6.2f}  {self.settling_a:<6.2f}      dust settling, exp_strat, a_strat
  {self.sublimate_flag}                      sublimate dust
  F  0.0                 viscous heating, alpha_viscosity

#Number of zones : 1 zone = 1 density structure + corresponding grain properties
  {self.nzones}
      """.format(self=self)

        if self.nzones > 1:
            raise NotImplemented("write par file code does not yet support multiple zones!")

        getsubkeys = lambda k, l: tuple([par[k][li] for li in l])
        template+="""
#Density structure
  {zone[zone_type]:1d}                       zone type : 1 = disk, 2 = tapered-edge disk, 3 = envelope
  {zone[dust_mass]:<10.2e} {zone[gas_to_dust]:<5.1f}        disk dust mass,  gas-to-dust mass ratio
  {zone[scale_height]:<5.1f}  {zone[ref_radius]:<6.1f}           scale height, reference radius (AU), unused for envelope
  {zone[r_in]:<6.1f}  {zone[r_out]:<6.1f} {zone[edge]:<6.1f}   Rin, Rout (or Rc), edge falloff (AU)
  {zone[flaring_exp]:<8.3f}                  flaring exponent, unused for envelope
  {zone[surfdens_exp]:<6.3f} {zone[gamma_exp]:<6.3f}                surface density exponent (or -gamma for tappered-edge disk), usually < 0; -gamma_exp
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
      1.0 50                 vmax (m.s-1), n_speed
      T 1.e-6 abundance.fits.gz   cst molecule abundance ?, abundance, abundance file
      T  3                   ray tracing ?,  number of lines in ray-tracing
      1 2 3                  transition numbers
      """ 

        #starkeys = lambda l: tuple([self.star][li] for li in l])
        template+="""
#Star properties
  1 Number of stars"""
        for star in self.star:
            template+="""
  {star[temp]:<7.1f} {star[radius]:5.1f} {star[mass]:5.1f} {star[x]:5.1f} {star[y]:5.1f} {star[z]:5.1f} F       Temp, radius (solar radius),M (solar mass),x,y,z (AU), is a blackbody?
  {star[spectrum]:s}
  {star[fUV]:<5.1f}  {star[slope_fUV]:<5.1f}  fUV, slope_fUV
        """.format(star=star)
        
        # % starkeys(['temp', 'radius', 'mass', 'x', 'y', 'z', 'spectrum'])
        return template


    def writeto(self, outname='sample.par'):
        """ Write an MCFOST parameter file to disk. Currently outputs v2.15 param files


        WARNING: Not all parameters are allowed to vary, just the most useful ones. 
        Lots of the basics are still hard coded. 
        TBD fix later.


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


        outfile = open(outname, 'w')
        outfile.write(self.tostring())
        outfile.close()
        print "  ==>> "+outname




def get_current_paramfile(parfile=None, dir="./",verbose=False, wavelength=None):
    """ Find a MCFOST par file in a specified directory

    By default, look in the current directory of a model,
    but if a wavelength is specified, look inside the appropriate data_## directory

    KEYWORDS:
    wavelength : string
        string with wavelength name
    verbose :Boolean.
        whether to be more verbose in output
    dir : string
        directory to look in. Default is current directory

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


