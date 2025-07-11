**********************************
Packaging Documentation for EPRsim
**********************************

Packaging method
################
This package uses the following frameworks:

setuptools
    Package build system
pytest
    Testing framework
sphinx
    Documentation framework

Packaging Info
==============

Building
--------
The build system used is setuptools, which is configured by the following
parts of ``pyproject.toml``:

.. code:: toml

    [build-system]
    requires = ["setuptools"]
    build-backend = "setuptools.build_meta"

    [tool.setuptools.package-dir]
    eprsim = "src"

The package-dir block is particularly important, as it tells setuptools
that the package eprsim is located at the path ``src`` (relative to
pyproject.toml); by default, it would look for the eprsim package in the 
``eprsim`` subdirectory. We can alternately set it up to look for packages
*in* the src directory (for instance if we wanted separate ``eprsim`` and
``eprload`` packages) by instead setting:

.. code:: toml

    [tool.setuptools]
    package-dir = {"" = "src"}
    # to look for src/eprsim and src/eprload

For the end user, the package can be built with ``python -m build``, then
``pip install .``. For development, install in editable mode (``-e``) so
the files are referenced in their original locations rather than copied
to the venv, making any changes to the source files reflected in the
package when it's imported. Also install the "devel" set of optional
packages (defined under ``[project.optional-dependencies]`` in 
``pyproject.toml``), so the finall install command should be:
``pip install -e '.[devel]'``.

Editable mode does *not* extend to the package settings, so if any changes
are made in ``/pyproject.toml`` or ``/src/__init__.py``, the package will
need to be rebuilt and reinstalled.

Versioning
----------
The package and its documentation read the package version from
:code:`__version__` in :code:`src/__init__.py`. This is achieved by
the following in the ``pyproject.toml``:

.. code:: toml

    [project]
    name = "eprsim"
    dynamic = ["version"]

    [tool.setuptools.dynamic]
    version = {attr = "eprsim.__version__"}

Troubleshooting update to .toml
===============================

The package originally used the `setup.py` method to build the package. For
original setyup config, see :ref:`ogSetupCode`. To be in line with what
setuptools now expects, I needed to convert this to a `pyproject.toml` file.
This wasn't too hard (see file in package root directory), but it led to
unexpected behavior in packaging. Whereas previously, each file was treated
as subpackage of EPRsim, and so could be imported with, for instance `from
eprsim import Tools`, the `__init__.py` now appeared to treat everything as a
single module with no components.

Evidently, the `default behavior`_ for python packages is now that each
subpackage must have its own directory and `__init__.py` file. The code is in
the structure below: 

.. _default behavior: https://packaging.python.org/en/latest/guides/packaging-namespace-packages/

.. code::

    src/
    ├── __init__.py
    ├── EPRsim.py
    ├── EPRload.py
    ├── Convolutions.py
    ├── Direct_conversion_to_Field.py
    ├── FastMotion.py
    ├── Hamiltonian_Eig.py
    ├── Hamiltonian_Point_Group.py
    ├── Interpolation_lib.py
    ├── Nucdic.py
    ├── Pauli_generators.py
    ├── Presettings.py
    ├── resfield_full.py
    ├── SolidState.py
    ├── spectral_processing.py
    ├── Tools.py
    └── Validate_input_parameter.py

The outcome I *want* is to be able to `import eprsim`, or `from eprsim import
tools`. What I *get* is that the package only imports `EPRsim.py`. After a *lot*
of faffing about, I found the issue was in the `pyproject.toml`, where if no
package path is defined, it looks for a subfolder with the package name. Since
my package folder is called src, it finds nothing. I thought I dealt with it by
setting the following:

.. code :: toml

   [tool.setuptools]
   package-dir = {"" = "src"}

However, as it turns out that forces the automatic discovery system to use the
`src-layout` option. To allow using `flat-layout` without naming the source folder
`eprsim`, I had to use:

.. code :: toml

   [tool.setuptools.package_dir]
   eprsim="src"

Info Sources
~~~~~~~~~~~~

* Setuptools `package discovery documentation`_.
* A `wierd corner`_ of the setuptools documentation which I only saw on the source GitHub.
* The Python Packaging `user guide`_.

.. _package discovery documentation: https://setuptools.pypa.io/en/latest/userguide/package_discovery.html
.. _wierd corner: https://github.com/pypa/setuptools/blob/main/docs/userguide/package_discovery.rst
.. _user guide: https://packaging.python.org/en/latest/guides/writing-pyproject-toml/

.. _ogSetupCode:

Original Setup Code
-------------------

setup.py
~~~~~~~~

.. code :: python

    #! /usr/bin/env python
    # --- import -------------------------------------------------------------------------------------
    import os
    from setuptools import setup, find_packages
    # --- define -------------------------------------------------------------------------------------
    here = os.path.abspath(os.path.dirname(__file__))
    extra_files = []
    extra_files.append(os.path.join(here, "CONTRIBUTORS.txt"))
    extra_files.append(os.path.join(here, "LICENSE.txt"))
    extra_files.append(os.path.join(here, "README.md"))
    extra_files.append(os.path.join(here, "EPRsim", "VERSION"))
    # --- setup --------------------------------------------------------------------------------------
    with open(os.path.join(here, "requirements.txt")) as f:
        required = f.read().splitlines()
    with open(os.path.join(here, "EPRsim", "VERSION")) as version_file:
        version = version_file.read().strip()
    setup(
        name="EPRsim",
        version=version,
        packages=find_packages(),
        package_data={"": extra_files},
        install_requires=required,
        author="Darien Morrow",
        author_email="darienmorrow@gmail.com",
        license="GPLv3",
        url="https://github.com/darienmorrow/EPRsim",
        keywords="photophysics spectroscopy science paramagnetic resonance",
        entry_points={"console_scripts": []},
        classifiers=[
            "Development Status :: 1 - Planning",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering",
            "Natural Language :: English",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.3",
            "Programming Language :: Python :: 3.4",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
        ],
    )

requirements.txt
~~~~~~~~~~~~~~~~

.. code::  

    numpy
    scipy
    numba
    matplotlib

src/__init__.py
~~~~~~~~~~~~~~~

.. code:: python

    from .__version__ import *
    from . import Convolutions
    from . import Direct_conversion_to_Field
    from . import EPRsim
    from . import FastMotion
    from . import Hamiltonian_Eig
    from . import Hamiltonian_Point_Group
    from . import Interpolation_lib
    from . import Nucdic
    from . import Pauli_generators
    from . import Presettings
    from . import resfield_full
    from . import SolidState
    from . import spectral_processing
    from . import Tools
    from . import Validate_input_parameter

src/__version__.py
~~~~~~~~~~~~~~~~~~

.. code:: python

    """Define EPRsim version."""
    # --- import --------------------------------------------------------------------------------------
    import os
    # ---- define -------------------------------------------------------------------------------------
    here = os.path.abspath(os.path.dirname(__file__))
    __all__ = ["__version__", "__branch__"]
    # --- version -------------------------------------------------------------------------------------
    # read from VERSION file
    with open(os.path.join(here, "VERSION")) as f:
        __version__ = f.read().strip()
    # add git branch, if appropriate
    p = os.path.join(os.path.dirname(here), ".git", "HEAD")
    if os.path.isfile(p):
        with open(p) as f:
            __branch__ = f.readline().rstrip().split(r"/")[-1]
        if __branch__ != "master":
            __version__ += "-" + __branch__
    else:
        __branch__ = None

