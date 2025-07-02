#######
EPRload
#######

Basic use
=========

To use the eprload class:

>>> from eprsim.EPRload import *
>>> spectrum = eprload("./")

User Documentation
==================

.. module:: src.EPRload
.. autoclass:: eprload
    :members:

Internal Documentation
======================

The logic flow of this program (in pseudocode) is:

.. code:: 
    call :class:`eprload`
        if file is .dsc or .dta :method:`loadBES3T` to intiate bes3t loading
            load parameters from .dsc file
            parse parameters to determine how to load data 
            for each axis:
                if axis type is IGD, abscissa = :method:`_readNonLinearAbscissa`
                if axis type is IDX, assign abscissa as linspace defined by params
                if axis type is NTUP, error, because this is not implemented.
            assign :attr:`self.Spec` = :method:`_readBinaryDataMatrix`.
        if file is .d01 or .exp :method:`loadSpecMan`
            call :method:`_loadEXP`
                read parameters from .exp file using pyyaml
                parse parameters to correct type
                assign parameters to :attr:`self.Param`
            call :method:`_loadD01`
                test that :attr:`self.Param` is assigned.
                read initial bytes of file (# spec, # axes, # data points)
                read data into appropriately sized array


