EPRsim documentation
====================

What is EPRsim ?
----------------

EPRsim is a Python 3 package for the interpretation of EPR spectra. It consists of
three major components:

EPRsim
    Functions for simulation of EPR spectra; currently limited to CW spectra.
EPRload
    A class for loading experimental EPR data; currently limited to Bruker
    and Specman data.
Tools
    Miscellaneous utilities for handling experimental and simulated EPR data.

These components follow the conventions of Stefan Stoll's `EasySpin`_ where possible.
All other compoents of this package and the documentation are private modules utilized
by the EPRsim module, and don't apply to the end user.

EPRsim was originally developed by Stephan Rein in the lab of Prof. Dr. Stefan Weber.
Although it was slightly intimidating to be the first non Ste(f|ph)an to work on
these projects, EPRload was adapted from the EasySpin implementation by Liam Twomey,
and the package is currently maintained by Liam Twomey.

.. _Easyspin: https://github.com/StollLab/EasySpin

This package is developed at the `git repository`_, and is not yet available on PyPI.
 
.. _git repository: https://github.com/Liam-Twomey/EPRsim

.. toctree::
    :maxdepth: 2
    :caption: Contents:

    eprsim
    eprload
    tools
    og-coverpage
