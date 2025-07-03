.. EPRsim documentation master file, created by
   sphinx-quickstart on Sat Jun 28 13:44:13 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

EPRsim documentation
====================
.. toctree::
    :maxdepth: 1
    :hidden:

    User Guide <user/index.rst>
    Development Docs <devel/index.rst>


What is EPRsim ?
----------------

EPRsim is a Python3 package for the interpretation of EPR spectra.

EPRsim was originally developed by Stephan Rein in the lab of Prof. Dr. Stefan Weber,
following the conventions of Stefan Stoll's `Easyspin`_ .

Although it was slightly intimidating to be the first non-``Ste(f|ph)an`` to work on
these projects, EPRload was adapted from the EasySpin implementation by Liam Twomey,
and the package is currently maintained by Liam Twomey.

.. _Easyspin: https://github.com/StollLab/EasySpin

This package is developed at the `git repository`_, and is not yet available on PyPI.
 
.. _git repository: https://github.com/Liam-Twomey/EPRsim

.. grid:: 1 1 2 2 
    :gutter: 2 3 4 4

    .. grid-item-card::
        :img-top: ../source/_static/user_guide.svg
        :class-img-top: sd-card
        :text-align: center

        **User Guide**
        ^^^

        Guide to use of EPRsim and documentation of user-facing functions.

        +++

        .. button-ref:: user/index
            :expand:
            :color: secondary
            :click-parent:

            Visit User Guide

    .. grid-item-card::
        :img-top: ../source/_static/devel.svg
        :class-img-top: sd-card
        :text-align: center

        **Development Documentation**
        ^^^

        Documentation on internal functions and development processes for EPRsim 

        +++

        .. button-ref:: devel/index
            :expand:
            :color: secondary
            :click-parent:

            Visit Dev Docs 

License and Permissions
-----------------------

This project is licensed under GPLv3, meaning that it is free to use,
modify, make derivative versions, and distribute. If any modified or derivative
verison of this software is made, its source code must also be made freely
available under the GPL license.

The testing EPR data (in :code:`test/eprfiles/`) was forked from the EasySpin
project, and is licensed under the MIT license, which is GPL-compatible.
