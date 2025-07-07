Writing EPRsim modules
======================


Numba support
-------------

EPRsim uses the `Numba`_ JIT compiler to accelerate operations by compiling python (specifically raw python code and numpy) to assembly via the llvm c compiler. JIT compilation will accelerate functions which are called many times the most, and loops. Read the Numba primer linked above before starting to work with Numba.

.. _Numba: https://numba.readthedocs.io/en/stable/user/5minguide.html
