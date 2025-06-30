# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- imports for autofunctions
import sys
from pathlib import Path
sys.path.insert(0, str(Path('../..').resolve()))
from importlib.metadata import version as vn
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


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
