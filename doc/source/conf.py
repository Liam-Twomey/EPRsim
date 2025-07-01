# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- imports for autofunctions
import sys
from pathlib import Path
from importlib.metadata import version as vn
#import mock

pkg_path  = str(Path('../..').resolve())
print('package location:', pkg_path)
sys.path.insert(0,pkg_path) 
ver = vn('eprsim')

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'eprsim'
copyright = '2025, Stephan Rein, Darien Morrow, Liam Twomey'
author = 'Stephan Rein, Darien Morrow, Liam Twomey'
release = '.'.join(ver.split('.'))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.napoleon']

templates_path = ['_templates']
exclude_patterns = []
numpydoc_show_class_members= False

#toMock = ['Validate_input_parameter','Tools']
#for modName in toMock:
#	sys.modules[modName] = mock.Mock()

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
