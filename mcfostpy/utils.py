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
    # set up screen logging
    logging.basicConfig(level=logging.INFO,format='%(name)-10s: %(levelname)-8s %(message)s')
    if verbose: print("Setup_logging is adjusting Python logging settings:")

    name='mcfostpy'

    if level.upper() =='NONE':
        # disable logging
        lev = logging.CRITICAL  # we don't generate any CRITICAL flagged log items, so
                                # setting the level to this is effectively the same as ignoring
                                # all log events. FIXME there's likely a cleaner way to do this.
        if verbose: print "No log messages will be shown."
    else:
        lev = logging.__dict__[level.upper()] # obtain one of the DEBUG, INFO, WARN, or ERROR constants
        if verbose: print "Log messages of level {0} and above will be shown.".format(level)

    logging.getLogger(name).setLevel(lev)


    if str(filename).strip().lower() != 'none':
        hdlr = logging.FileHandler(filename)

        formatter = logging.Formatter('%(asctime)s %(name)-10s: %(levelname)-8s %(message)s')
        hdlr.setFormatter(formatter)

        for name in lognames:
            logging.getLogger(name).addHandler(hdlr)

        if verbose: print("Log outputs will also be saved to file "+filename)

def find_closest(list_, item):
    """ return the index for the entry closest to a desired value """
    minelt = min(list_,key=lambda x:abs(x-item))
    return np.where(np.asarray(list_) == minelt)[0][0]


