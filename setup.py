from setuptools import setup, find_packages

setup(
    name = 'mcfostpy',
    fullname = 'MCFOST Python Tools',
    version = "0.01",
    description = 'Python utilities for MCFOST radiative transfer simulations',
    long_description = """
    Python tools for working with MCFOST radiative transfer models.
    """,

    author = "Christophe Pinte & Marshall Perrin",
    author_email = "",
    url = "https://github.com/cpinte/mcfost-python",

    requires = ['astropy','numpy', 'matplotlib'],
    packages = ['mcfost'],
    classifiers = [
        "Programming Language :: Python",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Development Status :: 4 - Beta"
        ]
)
