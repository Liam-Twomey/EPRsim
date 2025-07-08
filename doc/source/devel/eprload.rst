#######
EPRload
#######

Dev Notes
=========

Module structure
----------------

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

Notes on BES3T format
---------------------
This information is pulled from the ``BES3T-2.0.PDF`` document on our E500 PC.

Result sets (III.1.2.ii) are what causes there to be multiple data sets per DTA,
file, and are used in the *rare* case that simultaneous EPR and other data is
recorded, i.e. concurrent optical and EPR measurements. IKKF defines the number
of elements in a result set, which is stored :math:`[A_1 B_1 C_1 A_2 B_2 C_2]`
for every data point.

Basically there are three levels of data hierarchy in the ``.DTA`` file, which
when reassembled generate one :math:`N_x*N_y*N_z` array per signal.

* The ``signal``. Number of signals is determined by the length of the ``IKKF``
  property; if there is more than one value for IKKF, then some other data was 
  recorded alongside the EPR spectrum with the same number of points. The first
  ``signal`` is EPR data, all following are auxilliary data.
* The ``trace``. Each signal is composed of a Real trace or one Real and one
  Imaginary trace.
* The ``value``. Each trace is composed of values for the trace, for a specific
  data point.

Data is stored point by point, in the format :math:`S_1\{T_1[V_1,V_2],
T_2[V_1,V_2]\}, S_2\{...\}...S_N`, where only Values are actually stored.


Backend Documentation
=====================

.. autoclass:: src.EPRload.eprload
    :no-index:
    :members:
    :private-members:
    :member-order: bysource

