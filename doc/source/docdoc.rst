################
Sphinx Reference
################

About the documentation
=======================

EPRsim documentation is written in the `reStructuredText`_ format, using `Sphinx`_
as the documentation builder. Categories of documentation are written in :code:`.rst`
files in :code:`doc/`, then their name is added to :code:`index.rst` to ensure Sphinx looks
for them.

The  documentation uses the `numpydoc`_ extension for formatting, and the
`furo`_ theme.

.. _numpydoc: https://numpydoc.readthedocs.io/en/latest/install.html
.. _furo: https://github.com/pradyunsg/furo

To build documentation, :code:`pip install sphinx numpydoc furo` to your project venv
(or install eprsim with the :code:`[devel]` tag), then run :code:`make html` in the
:code:`doc` directory. The html documentation will be built into
:code:`doc/_build/html`; to access them conveniently open :code:`index.html`
from that directory.

.. _reStructuredText: https://docutils.sourceforge.io/rst.html
.. _Sphinx: https://www.sphinx-doc.org/en/master/index.html

Writing EPRsim documentation
============================
The majority of the EPRsim documentation is automated via the Sphinx `autodoc`_
extension, which can parse docstrings into module/class/function/etc. documentation.
Since we also use :code:`numpydoc`, all syntax follows their `style guide`_.

.. _autodoc: https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html
.. _style guide: https://numpydoc.readthedocs.io/en/latest/format.html

To add documentation to a function, class, etc, it should have a docstring of the
following form, where :code:`<name>` denotes a place you would insert a relevant
name.

.. code:: rst 

    """
    <short description of function>

    Parameters
    ----------
    <argument 1>: <arg 1 type>
        <argument 1 description>
    <argument 2>: <arg 2 type>
        <argument 2 description>

    Returns
    -------
    <return value 1>: <return value 1 type>
        <return value 1 description>
    <return value 2>: <return value 2 type>
        <return value 2 description>

    See Also
    --------
    <Other relevant functions, same format as arguments and returns>

    Notes
    -----
    <Any notes or information about how the function works>

    """

***Important note:*** For the documentation to build correctly, every indent
in front of the docstring must be a tab character! Linters do not check for
correct rST syntax in a python file. The documentation formatting will be
all sorts of messed up if the docstring indentation is inconsistent.
The whole package uses tabs anyway; if you see a file not doing so, please fix
it! Vim :code:`:set noexpandtab` then :code:`:retab!` makes it easy.

Documentation in the :code:`doc` subdirectory. :code:`index.rst` indicates which
:code:`.rst` files to include in the documentation. There is one documentation
file per module.

To make a new documentation page, make a copy of :code:`_templates/module-template.rst`
and add it to the index file. Then, fill out each section.

#. The Basic Use section should contain a couple of `doctest`_ blocks showing how to
   use the components. 

#. The documentation block can be largely automated with :code:`..auto<type>:: <name>`
   directives. The template has an example for doing a whole module in one go, but it
   often makes sense to do non-class functions individually in their own sections.
   If you're starting at a lower level than the whole module, you'll need to specify
   the upper levels via one of the following syntaxes:

    .. code:: rst

        # define for whole file (until changed).
        .. module:: EPRsim
        .. class:: eprload
        .. autofunction:: loadBES3T

        # define for one auto statement
        .. autofunction:: EPRsim.eprload.loadBES3T
    
     Autodoc knows we're working within the eprsim package, so we don't need to define
     that. 

.. _doctest: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#doctest-blocks


Specific style considerations
-----------------------------

#. Always use tabs for indentation.
#. Never use docstring sections outside the numpydoc specification.
#. Never submit PRs with undocumented code.
#. Make sure your documentation builds without errors before submitting a PR.

Useful Links
------------

* The Sphinx `getting started`_. 
* reStructuredText `primer`_.
* Cross referencing `guide`_. 
* Python RST `documentation`_.

.. _getting started: https://www.sphinx-doc.org/en/master/usage/quickstart.html
.. _primer: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html
.. _guide: https://www.sphinx-doc.org/en/master/usage/referencing.html
.. _documentation: https://devguide.python.org/documentation/markup/

Troubleshooting
===============

Sphinx is wildly finicky about how the documents are laid out. It specifically
needs to import all the libraries in the package in order for autodoc to work,
but it won't do it! Importing the files in python works just fine, but Sphinx
is throwing a fit!

I tried reinitializing sphinx, and it's fine until the rst files with autodoc
commands in them are added to the source directory, even if they're not in the index!

* A test file with just :code:`.. module:: EPRload` is fine.
* Adding :code:`.. class:: eprload` is also fine!
* adding :code:`.. autofunction:: checkFileType` throws the error! 

Made a minimal example where the only items in the src directory are __init__.py and
EPRload.py. Introducing a file test.rst to the index:
* :code:`.. automodule:: EPRload; :members:` works
* :code:`.. module:: EPRload; .. autoclass:: eprload; :members:` it works
* :code:`.. module:: EPRload; .. class:: eprload; .. automethod:: getSpec` fails to get function
* :code:`.. module:: EPRload; .. automethod:: eprload.getSpec` fails to get function
* :code:`.. automethod:: EPRload.eprload.getSpec` fails to get function
* The patching of the sphinx `path`_ was also removed, as it is not necessary if the
package is :code:`pip install -e .`
* After copying Tools.py over :code:`.. automodule:: EPRload; :members:` no longer works
after changing files in src/!! This is fixed by re-enabling path patching to `../../src`
* :code:`.. module:: EPRload; .. automethod:: eprload.getSpec` now works?

.. _path: https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#ensuring-the-code-can-be-imported 

Bringing eprload back in, was able to get it to run. Using nepoleon:
* Sets of braces not in `:math:` sections cause issues.
* spacing on lists is very very specific: It's displayed exactly how it looks
in the docstring, so lists with multiple lines per list item *must* be spaced
in by the same number of spaces as the first character. This is in contrast to normal rst,
where lists shouldn't be indented at all!

    .. code:: rst

        #. item 1
        #. item 2
           has multiple lines of text
        #. item 4
            will look screwy because I used a tab
            * this is a nested list, unlike normal RST, do *not* add spaces before and after.
        #. and here we return to the parent list.


Also note that references in docstrings must have a space before them!ImportError: attempted relative import with no known parent package

The biggest problem is files which import other parts of the package. As it turns 
out, that issue is due to a mix of configuration issues. Specifically, for sphinx
to recognize the source as a package, it has to be pointed to the root directory
of the project, then paths to modules have to be defined from that directory. So:

.. code:: python

    # conf.py
    import sys
    from pathlib import Path
    pkg_path  = str(Path('../..').resolve())
    print('package location:', pkg_path)
    sys.path.insert(0,pkg_path) 

.. code:: rst

    In documentation rst files:

    .. module:: src.eprload

 There's also still something funky with EPRsim.py docstrings specifically.
 Even though the FastMotion can be autodocumented, EPRsim just won't. Also
 note that if there is a .rst file referencing a module/function, the build
 will throw errors, even if the file is not in the index.

 Tried again with SolidState ... it initially had no issues, but there were
 no docstrings for the two functions in it; when I added docstrings, it now
 fails.

 I have no idea what I changed (other than making sure there were no empty
 docstrings), and now it works?

 Oh my god it's the tabbing again. Somehow EPRsim.py went back to tabbing by
 spaces, and vim :code:`:retab!` fixed everything.
