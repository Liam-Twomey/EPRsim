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
needs to import all the libraries in the package in order for aut
