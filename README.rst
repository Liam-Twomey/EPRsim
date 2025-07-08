EPRsim
======

From the original PyPI: 

    EPRsim is an Open-source simulation package for cw-EPR spectra.
    EPRsim has been developed in the group of Prof. Dr. Stefan Weber
    at the University of Freiburg, Freiburg im Breisgau, Germany,
    during the last couple of years. EPRsim was developed by Stephan Rein.
    The program uses several concepts described in [#f1]_. EPRsim is open-source
    and available free of charge.

The goal of this fork is to (slowly, but surely) extend the original EPRsim package with
additional features, largely by reimplementing components of Stefan Stoll's
`Easyspin <Easyspin.org>`__ package for MATLAB.

Current capabilities
--------------------

-  Simulation of EPR spectra in fast motion, slow, motion, and solid-state
   regimes (by Stephan Rein).
-  Loading of Bruker BES3T (.DSC/.DTA) files and beta support for loading Specman (.exp/d01) files
   (adapted from the EasySpin implementation by Liam Twomey).

Documentation
-------------

Documentation for the package is hosted on `readthedocs <https://eprsim.readthedocs.io/en/latest/>`__ . Please see the development section to see how to build the docs yourself.

Installation
------------

Install from source
~~~~~~~~~~~~~~~~~~~

This project uses the ``setuptools`` build system, via the ``build``
module.

Requires Python 3; python 3.10+ is strongly recommended.

It is highly recommended (but not required) to setup a virtual
environment for EPRsim with:

.. code:: shell

   python -m venv <venvName>
   source <venvName>/bin/activate

Each time you open a terminal, you will need to run
``source <venvName>/bin/activate`` to activate the virtual environment
Then, to install the program:

.. code:: shell

   pip3 install build
   git clone https:/github.com/LiamTwomey/EPRsim
   cd EPRsim
   python -m build
   pip install .

Installing from PyPI
~~~~~~~~~~~~~~~~~~~~

This package is not yet on PyPI (the latest version on PyPI is 0.0.4,
Stefan Rein’s original version), but it’s planned to be eventually. If
it ever is, it will be installable via ``pip3 install eprsim``.

Compatibility
~~~~~~~~~~~~~

This fork has been tested with Python 3.12, Numba 0.60.0, scipy 1.14.1,
numpy 2.0.2, and matplotlib 3.9.2

Development
-----------

Installing for development
~~~~~~~~~~~~~~~~~~~~~~~~~~

To install for development, instead use:

.. code:: shell

   git clone https:/github.com/LiamTwomey/EPRsim
   cd EPRsim
   python -m build
   pip install -e '.[all]'.

This does three things:

#. Install the development dependencies ``pytest`` and ``sphinx`` for
   testing and documentation, respectively, as well as their dependencies.
#. Installs the package in editable format, so changes to the ``src``
   directory are reflected on next import, without needing to rebuild.
#. Installs the optional dependency Numba to accelerate the code using
   just-in-time compilation.

The ``all`` optional dependency group installs the ``dev``, ``doc``, and ``fast``
packages, which install dependencies for the 

Testing
~~~~~~~

This package uses ``pytest`` for testing. After installing pytest
(``pip3 install pytest``), running the terminal command ``pytest`` from
anywhere in the repository will execute the tests in ``./tests/``.

Contributing
------------
I would welcome other contributors who would like to help
make a working EPR simulation program which does not depend on
proprietary software. Please feel free to get in touch.

References
----------


References
==========

.. [#f1]
   S. Stoll, A. Schweiger, J. Magn. Reson., 2006, 178, 42-55
