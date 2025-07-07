#######
EPRload
#######

Module structure
================

The EPRload module consists of just the eprload class.
When instantiated on a valid filename, ``init`` calls the appropriate class method
for the file. In turn, it calls associated dundered subfunctions (``_readXXX``)
to read the files in. Note that *only the ``loadXXX`` functions mutate
the eprload object*. This is for modularity and numba compatibility. The impact of
any other method on the class object should be defined entirely by how its returns are
assigned when it is called.

The logic flow of this program (in pseudocode) is:

.. code:: 

    Instantiate :class:`eprload`, calling ``__init__``, then ``checkFileType``.
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
                read parameters from .exp file
                parse parameters to correct type
                assign parameters to :attr:`self.Param`
            call :method:`_loadD01`
                test that :attr:`self.Param` is assigned.
                read initial bytes of file (# spec, # axes, # data points)
                read data into appropriately sized array



Backend Documentation
=====================

.. autoclass:: src.EPRload.eprload
    :no-index:
    :members:
    :private-members:
    :member-order: bysource

