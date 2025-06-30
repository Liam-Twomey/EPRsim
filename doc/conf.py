# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Imports and Setup -------------------------------------------------------
import sys # for autodoc
from pathlib import Path # for autodoc
sys.path.insert(0, str(Path('../src').resolve()))

from importlib.metadata import version
release = version('eprsim')
# for example take major/minor
version = '.'.join(release.split('.'))

import mock # to allow building with import statements
# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'EPRsim'
copyright = '2025, Stefan Rein, Darien Morrow, Liam Twomey'
author = 'Stefan Rein, Darien Morrow, Liam Twomey'
release = version#__version__ 

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
	'sphinx.ext.napoleon',
]

numpydoc_class_members_toctree=False
#mock_modules = ['numpy', 'scipy', 'matplotlib', 'matplotlib.pyplot', 'Tools',
#	'Convolutions','Direct_conversion_to_field',]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
toc_object_entries  = False
html_theme = 'furo'
html_theme_options = {}
html_static_path = ['_static']
