# Various utility functions 
import matplotlib
import numpy as np
import logging

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


# purely cosmetic: print e.g. '100' instead of '10^2' for axis labels with small exponents
class NicerLogFormatter(matplotlib.ticker.LogFormatter):
    """ A LogFormatter subclass to print nicer axis labels
        e.g. '100' instead of '10^2'

        Parameters                                                                                                   ----------
        threshhold : int
            Absolute value of log base 10 above which values are printed in exponential notation.
            By default this is 3, which means you'll get labels like                                                     10^-4, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10^4 ...
                                                                                                                     usage:
          ax = gca()
          ax.set_yscale("log")
          ax.set_xscale("log")
          ax.xaxis.set_major_formatter(NiceLogFormatter())
    """
    def __init__(self, threshhold=3):
        self.threshhold = threshhold
    def __call__(self,val,pos=None):
        if abs(np.log10(val)) > self.threshhold:
            return "$10^{%d}$" % np.log10(val)
        elif val >= 1:
            return "%d"%val
        else:
            return "%.1f"%val



def setup_logging(level='INFO',  filename=None, verbose=False):
    """ Simple wrapper function to set up convenient log defaults, for
    users not familiar with Python's logging system.

    """
    import logging
    _log = logging.getLogger('mcfost')

    lognames=['mcfost']

    if level.upper() =='NONE':
        # disable logging
        lev = logging.CRITICAL  # we don't generate any CRITICAL flagged log items, so
                                # setting the level to this is effectively the same as ignoring
                                # all log events. FIXME there's likely a cleaner way to do this.
        if verbose: print "No log messages will be shown."
    else:
        lev = logging.__dict__[level.upper()] # obtain one of the DEBUG, INFO, WARN, or ERROR constants
        if verbose: print "Log messages of level {0} and above will be shown.".format(level)

    for name in lognames:
        logging.getLogger(name).setLevel(lev)
        _log.debug("Set log level of {name} to {lev}".format(name=name, lev=level.upper()))

    # set up screen logging
    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')
    if verbose: print("Setup_logging is adjusting Python logging settings.")


    if str(filename).strip().lower() != 'none':
        hdlr = logging.FileHandler(filename)

        formatter = logging.Formatter('%(asctime)s %(name)-10s: %(levelname)-8s %(message)s')
        hdlr.setFormatter(formatter)

        for name in lognames:
            logging.getLogger(name).addHandler(hdlr)

        if verbose: print("Log outputs will also be saved to file "+filename)
        _log.debug("Log outputs will also be saved to file "+filename)

def find_closest(list_, item):
    """ return the index for the entry closest to a desired value """
    minelt = min(list_,key=lambda x:abs(x-item))
    return np.where(np.asarray(list_) == minelt)[0][0]

def ccm_extinction(Rv, lambda_ang):

    """ 
    Python implementation of the idl_lib extinction correction function
    to be called by the SED chisqrd fitting method in the python version of
    MCRE. Accepts Rv, the reddening index at V (Default = 3.1) and the 
    wavelength in Angstroms. Extinction curve A(lambda)/A(V) is returned.
    """
    lambda_ang = np.asarray(lambda_ang)
    inv_lam = 1.0/lambda_ang
    #print 'inv_lam',inv_lam
    s = len(lambda_ang)
    a = np.zeros((s))
    b = np.zeros((s)) # confirm proper syntax

    # Range that CCM restrict it to.
    ir = inv_lam<=1.1
    #print 'ir',ir
    c_ir = len(ir)
    #flags = choose(greater(inv_lam,1.1),(-1,1))
    
    a[ir] = 0.574*inv_lam[ir]**1.61
    b[ir] = -0.527*inv_lam[ir]**1.61

    #opt = where inv_lam > 1.1 and inv_lam <= 3.3 then c_opt
    opt = ((inv_lam > 1.1) & (inv_lam <= 3.3))
    c_opt = len(opt)
    y = np.asarray(inv_lam[opt] - 1.82)
    a[opt] = 1+ 0.17699*y-0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5-0.77530*y**6 + 0.32999*y**7
    b[opt] = 1.41338*y + 2.28306*y**2 +1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7

    #uv = where inv_lam > 3.3 and inv_lam <= 8.0 then c_uv
    uv = ((inv_lam > 3.3) & (inv_lam <= 8.0))
    c_uv = len(uv)
    x = inv_lam[uv]
    xm = x - 5.9
    fa = -0.04473*(xm)**2 - 0.009779*(xm)**3
    fb = 0.2130*(xm)**2 + 0.1207*(xm)**3
    nulls = xm <= 0
    fa[nulls] = 0.0
    fb[nulls] = 0.0
            
    a[uv] = 1.752 - 0.316*x - 0.104/( (x-4.67)**2 +0.341) + fa
    b[uv] = -3.090 + 1.825*x + 1.206/( (x-4.62)**2 + 0.263) + fb

    # Compute the extintion at each wavelength and return
    A_lambda = np.asarray((a+b)/Rv)
    
    return A_lambda



