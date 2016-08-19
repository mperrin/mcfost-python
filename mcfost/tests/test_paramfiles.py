import mcfost

import os
import glob
import numpy as np

TEST_DATA_PATH = '/Users/mperrin/Dropbox (Personal)/Documents/software/mcfost'

LATEST_VER = '3.0'

def do_test_read_one_par_file(filename, verbose=True):
    """
    Test the reading of parameter files.

    This relies on a copy of the 'refX.XX.para' reference example parameter
    files being present in this directory. This is not automatically updated,
    so right now somebody needs to manually update this whenever the MCFOST
    version number changes.
    """


    #check that we can read the file in
    pf = mcfost.Paramfile(filename, verbose=verbose)

    #check at least one attribute, from near the end of the file. 
    # if we've read this in correctly then it's at least plausible the whole
    # file is read OK. 
    assert isinstance(pf.version,float)

    assert pf.star.temp == 4000, "Did not get the expected stellar temperature from param file"

    # try a couple other parameters just for kicks. 
    assert pf.lambda_max == 3000, "Did not get the expected lambda_max from param file"
    assert pf.distance == 140, "Did not get the expected distance from param file"
    assert pf.l_pop == True, "Did not get the expected l_pop from param file"


def test_read_all_par_files():
    filenames = glob.glob(os.path.join(os.path.dirname(__file__), 'ref*.para'))
    for fn in filenames:
        do_test_read_one_par_file(fn)


def test_write_par_file(delete=False):
    import tempfile
    tmp = tempfile.NamedTemporaryFile(delete=delete)
    outname = tmp.name

    filename = os.path.join(os.path.dirname(__file__), 'ref{}.para'.format(LATEST_VER))
    pf = mcfost.Paramfile(filename)

    #write it out
    pf.writeto(outname)
    # and see if we can read it back in
    pf2 = mcfost.Paramfile(outname)

    # and see if the values in the output file match those in the original file.
    for key in pf2.__dict__.keys():
        # some particular keys we don't expect to round trip exactly identically:
        if key in ['filename', 'fulltext']:
            continue
        # but all the others should:
        assert pf[key] == pf2[key], "Did not round-trip the value for {}: {} vs {}".format(key, pf[key], pf2[key])


def test_properties(filename=None):
    """ Test the parameter file virtual attributes which are
    implemented via function calls via @properties
    """
    if filename is None:
        filename = os.path.join(os.path.dirname(__file__), 'ref{}.para'.format(LATEST_VER))
    pf = mcfost.Paramfile(filename)


    wavelens = pf.wavelengths
    assert isinstance(wavelens, np.ndarray)
    if pf.l_sed_complete:
        # we should only do this test if we're using default wavelengths
        # and not reading the wvelengths from an input FITS file.
        assert wavelens.size == pf.nwavelengths


    # test the abbreviations:
    assert pf.star is pf.stars[0]
    assert pf.density is pf.density_zones[0]
    assert pf.dust is pf.density_zones[0]['dust'][0]
