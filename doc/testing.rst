.. contents::
   :depth: 3
..

Testing Information
===================

This project uses `PyTest <pytest.org>`__ for testing. Tests are
organized into files (all suffixed with ``_test.py``) by the functions
which they test.

Running the test suite
----------------------

-  Install PyTest in your environment with (``pip3 install pytest``).
-  Once PyTest is installed, ``cd`` to the ``testing`` directory in your
   terminal (not in IDLE), and run ``pytest``.
-  This will run PyTest on all unit tests in files ending with
   ``_test.py``.

Writing new unit tests
----------------------

-  Each unit test is a function whose name begins with ``test_``.
-  The remainder of the name should make it clear what the function is
   testing.
-  A docstring (``'''some description'''``) can be used to add further
   clarification on the role of the function.
-  Any files which the test functions test against should be in a
   clearly-labelled subfolder of the ``eprfiles`` folder if they are EPR
   data, or in a subdirectory of ``testing`` if they are not.
