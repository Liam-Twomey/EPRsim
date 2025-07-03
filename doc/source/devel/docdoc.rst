######################
Documentation Overview
######################

EPRsim documentation is written in the `reStructuredText`_ format, using
`Sphinx`_ as the documentation builder. The documentation is in the ``doc/``
subdirectory of the git repo. 

.. _reStructuredText: https://docutils.sourceforge.io/rst.html
.. _Sphinx: https://www.sphinx-doc.org/en/master/index.html

Infrastructure
==============
Categories of documentation are written in ``.rst`` files in ``doc/source``,
divided up into ``user/`` for user-facing functions, and ``devel/`` for internal
documentation. Files are included in the HTML output if paths to them are added
to the ``index.rst`` file of their directory.

The  documentation uses the `napoleon`_ extension for auto-importing docstrings,
which are formatted to the `numpydoc`_ style guide, and uses the `PyData`_
sphinx theme for HTML output.

.. _napoleon: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
.. _numpydoc: https://numpydoc.readthedocs.io/en/latest/format.html
.. _PyData: https://pydata-sphinx-theme.readthedocs.io/en/stable/index.html

To work with documentation, install the package in development mode with
``pip install -e '.[devel]'``. This will install all dependencies needed for
development, as well as those needed for documentation: ``sphinx``,
``pydata_sphinx_theme``, ``sphinx-copybutton``, ``sphinx-design``, and ``mock``.

Writing EPRsim documentation
============================
The majority of the EPRsim documentation is automated via the Sphinx `autodoc`_
and ``napoleon`` extensions, which parse `numpydoc`_ -formatted docstrings in the
source code to allow us to use ReST directives to auto-index and import docstrings
into the documentation format, as described by the autodoc spec.

.. _autodoc: https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html
.. _style guide: https://numpydoc.readthedocs.io/en/latest/format.html

Docstring Syntax
----------------

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

***Important note:*** Numpydoc formatting has very strict whitespace
requirements, and they are different than those of ReST as a whole. Some
less-than-intuitive points:

#. Every indent in the docstring block must be a tab character! This is a hard
   error to catch, as the compiler will throw all sorts of wierd format errors 
   if there's a wrong space somewhere. The whole package uses tabs anyway;
   if you see a file not doing so, please fix it! Vim ``:set noexpandtab`` then
   ``:retab!`` makes it easy. This will mess up table tabbing however, so go
   back with ``v``, select, ``r<Space>`` to fix. 
#. If there are list items which take up multiple lines, the following lines
   must be indented exactly in line with the first line, or else it will throw
   an error about unexpected sections. In python docstrings, this will often 
   mean **manually using spaces to indent the following list items**.
#. Unlike normal ReST, nested lists must *not* have blank lines before and after.
#. The only allowed section headings in the docstring are the ones defined by the 
   Numpydoc spec, any others will throw an error.
#. Block directives like ``.. code::`` and ``.. math::`` *absolutely must* have 
   blank lines before and after the directive *and* the following block.

Documentation structure
-----------------------

The documentation structure is:

.. code::

  doc
  ├── make.bat
  ├── Makefile
  ├── build
  │   └─> location for built doc files
  └── source
      ├── conf.py
      ├── index.rst 
      ├──_static - 
      │   └─> location for aux files for build (css, svg, etc.)
      ├──_templates
      │   └─> to store templates, not indexed by sphinx
      ├── devel
      │   ├── index.rst
      │   └─> ReST files for dev doc pages
      ├── _templates
      │   └── module-template.rst
      └── user
          ├── index.rst
          └─> ReST files for other user pages

The ``source/index.rst`` ``.. toctree::`` directive determines which other files
will be read, in this case ``source/user/index.rst`` and ``source/devel/index.rst``.
Paths are determined from the file containing the toctree directive.

Writing New Documentation
-------------------------

To make a new documentation page, make a copy of :code:`_templates/module-template.rst`
in the user or devel directory, and and add it to the respective index file. Then,
fill out each section with info as needed.

#. The Basic Use section should contain a couple of `doctest`_ blocks showing how to
   use the components for user documentation; for devel docs, this is often omitted.
#. The documentation block can be largely automated with :code:`..auto<type>:: <name>`
   directives. The template has an example for doing a whole module in one go, but it
   often makes sense to do non-class functions individually in their own sections.
   If you're starting at a lower level than the whole module, you'll need to specify
   the upper levels via one of the following syntaxes:

   .. code:: rst

     Define for whole file (until changed).
     .. module:: EPRload
     .. autofunction:: eprload.loadBES3T

     Define full path for one auto statement
     .. autofunction:: EPRsim.eprload.loadBES3T

   Autodoc knows we're working within the eprsim package, so we don't need to define
   that.
#. For internal/development docs, make sure to set the ``:private-members:`` and
   ``:member-order: bysource`` to layout the functions in the order which they appear in
   the source file.

.. _doctest: https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html#doctest-blocks


.. important:: Never submit PRs with undocumented code, and make sure updated
   documentation builds without errors before submitting a PR.

Configuration
-------------

Everything works now, but man it was a pain to set up, so this documentation serves to
record the full setup process. In general, the `Numpy docs source`_ is a good reference for
configuration, as they use the numpydoc format (obviously), sphinx-design, and 
the same theme. For instance, the parent index file is an obvious copycat of Numpy's.

.. _Numpy docs source: https://github.com/numpy/numpy/tree/main/doc/source

#. Configuration of sphinx is done by ``/doc/source/conf.py``.

#. The PATH in ``conf.py`` must be patched to contain the *parent directory of the 
   package source files*. In our case that is ``/src``, so relative to ``doc/source/
   conf.py``, its parent directory is ``../..``.

   .. code:: python

       import sys
       from pathlib import Path
       pkg_path  = str(Path('../..').resolve())
       sys.path.insert(0,pkg_path)

   Now that the path has been patched, any files we want to reference in autodoc
   are defined *relative to that patched path*, so if we want to reference
   :class:`Parameters` from the file ``/src/EPRsim.py``, we would use
   ``.. autoclass:: src.EPRsim.Parameters``.

#. We read the version automatically from ``/source/__init__.py``, which defines
   the ``__version__`` variable. This is read by ``pyproject.toml``, and by
   ``conf.py`` via:

   .. code:: python

      from importlib.metadata import version as vn
      release = '.'.join(vn('eprsim').split('.'))

#. We depend on loading the following Sphinx extensions:
    * ``sphinx.ext.napoleon`` for numpy in autodoc
    * ``sphinx.ext.autodoc`` for automatic documentation of docstrings
    * ``sphinx.ext.mathjax`` for latex support in html 
    * ``sphinx_copybutton`` allow code copy button in examples
    * ``sphinx_design`` to enable design elements like grids

#. Finally, we set a couple of variables to get the settings how we want:

   .. code:: python

      # ignore files in this dir for indexing 
      templates_path = ['_templates']
      # set the theme for html export
      html_theme = 'pydata_sphinx_theme'
      # set directory to look for auxilliary files
      html_static_path = ['_static'] 
      # css classes for the docs, to change the theme
      html_css_files = ["eprsim.css"] 
      # Disable showing the "show source" link in the right sidebar.
      html_show_sourcelink = False 
      # Path to logo to show at top left
      html_logo = _static/EPRsim_logo.svg
      # Path to favicon for tab display
      html_favicon = _static/logo.svg
      # set html parameters for live linking to repo
      html_context = {
      	"github_user": "Liam-Twomey",
      	"github_repo" : "EPRsim",
      	"github_version" : "main",
      	"doc_path":"doc",

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
