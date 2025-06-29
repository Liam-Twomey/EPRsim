################
Sphinx Reference
################

Writing EPRsim documentation
============================

EPRsim documentation is written in the `reStructuredText`_ format, using `Sphinx`_ as the documentation builder. Categories of documentation are written in `.rst` files in `doc/`, then their name is added to `index.rst` to ensure Sphinx looks for them.

To build documentation, `pip install sphinx` to your project venv (or install eprsim with the `[devel]` tag), then run `make html` in the `doc` directory. The html documentation will be built into `doc/_build/html`; to access them conveniently open `index.html` from that directory.

.. _reStructuredText: https://docutils.sourceforge.io/rst.html
.. _Sphinx: https://www.sphinx-doc.org/en/master/index.html

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
