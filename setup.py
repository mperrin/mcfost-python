### See http://packages.python.org/distribute/setuptools.html

# The following trick copied from astropy setup
import sys
import imp
try:
    # This incantation forces distribute to be used (over setuptools) if it is
    # available on the path; otherwise distribute will be downloaded.
    import pkg_resources
    distribute = pkg_resources.get_distribution('distribute')
    if pkg_resources.get_distribution('setuptools') != distribute:
        sys.path.insert(1, distribute.location)
        distribute.activate()
        imp.reload(pkg_resources)
except:  # There are several types of exceptions that can occur here
    from distribute_setup import use_setuptools
    use_setuptools()

from setuptools import setup, find_packages


# set up package metadata for eventual use in PyPI
setupargs = {
    'name'          :       'mcfostpy',
    'version'       :      	"0.01", 
    'description'   :       'Python utilities for MCFOST radiative transfer simulations',
    'fullname'      :       'MCFOST Python Tools',
    'author'        :     	"Christophe Pinte & Marshall Perrin",
    'author_email'  :      	"",
    'url'           :  		"http://",
    'download_url'           :  		"http://",  # will be replaced below
    'requires'      :       ['astropy','numpy', 'matplotlib'],
    'classifiers'   :   [
        "Programming Language :: Python",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Development Status :: 4 - Beta"
        ],
    'long_description': """
Python tools for working with MCFOST radiative transfer models.
"""
    }


# Now call the standard setup function from distutils
setup(**setupargs)



