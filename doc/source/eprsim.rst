######
EPRsim
######

.. module:: src.EPRsim

Basic Use
=========


EPRsim User Guide
=================

The `Simulate` function
-----------------------

.. autofunction:: simulate

The `Parameters` object
-----------------------

.. autoclass:: src.EPRsim.Parameters
   :members:

EPRsim Backend Documentation
============================

Parameter Validation
--------------------

.. autoclass:: Validate_Parameters
   :members:

.. autoclass:: Simulation_Params
   :members:

Internal Simulation Utilities
------------------------------

.. autofunction:: remove_permutations
.. autofunction:: read_single_isotope_comb
.. autofunction:: remove_sub_threshold
.. autofunction:: get_indexvector_isotopes
.. autofunction:: create_nucvec
.. autofunction:: coupled_isotopes_indexvector
.. autofunction:: new_Nucsvec_and_tensors
.. autofunction:: check_eq_in_fast_motion
.. autofunction:: get_isotope_combinations
.. autofunction:: redefine_nuclear_coupling

Output utilities
----------------

.. autofunction:: check_if_instance
.. autofunction:: get_weighting_factor
.. autofunction:: warningflag
