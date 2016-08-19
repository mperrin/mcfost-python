
Documentation for MCFOST Python tools
=================================================


.. _intro:


MCFOST is a Monte Carlo and ray tracing radiative transfer model code for simulations of
circumstellar disks and other astrophysical phenomena, by `C. Pinte et al. <http://adsabs.harvard.edu/abs/2006A%26A...459..797P>`_.

The ``mcfost`` Python package is not that code. Rather it is a suite of associated utilities for
running simulations using MCFOST, displaying outputs, and performing model fits to observations. 
While the names are the same, context should generally make it clear which is which!


This documentation assumes readers are already familiar with MCFOST itself and have read the MCFOST Manual. 
Some references and comparisons are made to the MCFOST Results Explorer tool (in IDL), but familiarity with that isn't required.

.. warning::
        For the time being, this documentation is still very incomplete. As is the software itself!


The ``mcfost`` Python package is open source, developed on github in the repository ``mcfost-python``. Contributions and pull requests welcomed.

Contents
--------

.. toctree::
   :maxdepth: 1

   overview.rst
   paramfile.rst
   running.rst
   modelresults.rst
   observations.rst
   plots.rst
   chisqr.rst
   testing.rst

Indices and tables
------------------

* :ref:`genindex`
* :ref:`search`

Documentation last updated on |today|

