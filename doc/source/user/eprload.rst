#######
EPRload
#######

Basic use
=========

To use the eprload class:

>>> from eprsim.EPRload import *
>>> spectrum = eprload("./testfile.dsc")

This will return an ``eprload`` object called spectrum. Following similar syntax
to EasySpin, the magnetic field abscissa can be accessed as ``spectrum.B`` or
``spectrum.Absc``, the signal as ``self.S`` or ``self.Spec``, and the experimental
parameters as ``spectrum.P`` or ``spectrum.Param``.

* If ``debug=True`` is passed when instantiating the class, it will print out debug
  info as the files are read. 
* If ``keepTmp=True`` is passed when instantiating, it will save class attributes like 
  axis information and file format, at the cost of more memory usage.

User Documentation
==================

.. module:: src.EPRload
.. autoclass:: eprload
    :members:
