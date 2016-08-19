.. _testing:

Automated Testing
================================

This package contains a few unit tests which can be used to verify basic functionality and correct operation. 
The test suite is still in early development and is very incomplete for now. 


These tests are designed to be run using the py.test framework. See http://docs.pytest.org

To run, from the source code directory simply execute py.test::

    unix> py.test
    ======================== test session starts ==========================
    platform darwin -- Python 2.7.11 -- py-1.4.20 -- pytest-2.5.2
    plugins: cov
    collected 3 items

    mcfost/tests/test_paramfiles.py ...
    ====================== 3 passed in 1.23 seconds =======================


Descriptions of tests
-----------------------

Parameter files:

 * Test reading in the MCFOST reference parameter files for versions 2.15, 2.20, and 3.0. 
 * Test writing out the most recent parameter file (currently 3.0) and verifying that we can 
   read back in the output file with the expected values. 
 * Test the Parafile object attributes which are implemented via function calls as Python properties


There is not yet any automated testing for other parts of the code. 
