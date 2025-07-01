######
EPRsim
######


Basic Use
=========


EPRsim User Guide
=================

The `Simulate` function
-----------------------

.. module:: src.EPRsim
.. autofunction:: simulate

..
    The `Parameters` object
    -----------------------
    
    .. module:: src.EPRsim
    
    .. autoclass:: src.EPRsim.Parameters
       :members:
    
    EPRsim Backend Documentation
    ============================
    
    Fast-Motion Simulation
    ----------------------
    
    .. automodule:: FastMotion
        :members:
    
    Solid-State simulation
    ----------------------
    
    .. automodule:: SolidState
       :members:
    
    Parameter Validation
    --------------------
    
    .. autoclass:: Validate_Parameters
       :members:
    
    .. autoclass:: Simulation_Params
       :members:
    
    Internal Simulation Utilities
    ------------------------------
    
    .. autofunction:: EPRsim.remove_permutations
    .. autofunction:: EPRsim.read_single_isotope_comb
    .. autofunction:: EPRsim.remove_sub_threshold
    .. autofunction:: EPRsim.get_indexvector_isotopes
    .. autofunction:: EPRsim.create_nucvec
    .. autofunction:: EPRsim.coupled_isotopes_indexvector
    .. autofunction:: EPRsim.new_Nucsvec_and_tensors
    .. autofunction:: EPRsim.check_eq_in_fast_motion
    .. autofunction:: EPRsim.get_isotope_combinations
    .. autofunction:: EPRsim.redefine_nuclear_coupling
    
    Output utilities
    ----------------
    
    .. autofunction:: EPRsim.check_if_instance
    .. autofunction:: EPRsim.get_weighting_factor
    .. autofunction:: EPRsim.warningflag
    
    The Tools module
    ================
    #.. module:: Tools
    #.. autoclass:: physical_constants
    #.. autofunction:: pseudo_field_modulation
    #.. autofunction:: degree_in_rad
    
