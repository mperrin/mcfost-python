# Various utility functions 


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


