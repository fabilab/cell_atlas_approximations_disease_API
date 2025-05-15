# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import os
import sys
sys.path.insert(0, os.path.abspath('../../Python'))

# Hide the credits warning for atlasapprox-disease
os.environ["ATLASAPPROX_DISEASE_HIDECREDITS"] = "yes"

# -- Project information -----------------------------------------------------

project = 'Cell Atlas Approximations(Disease) API'
copyright = '2025, Fabio Zanini'
author = 'Ying Xu, Fabio Zanini'


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx_tabs.tabs",
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_gallery.gen_gallery',
]
sphinx_tabs_disable_tab_closing = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# Orders functions as they appear in the source code
autodoc_member_order = 'bysource'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for sphinx-gallery ----------------------------------------------
sphinx_gallery_conf = {
    'examples_dirs': 'python',  # Path to example scripts
    'gallery_dirs': 'auto_examples',    # Output directory for generated gallery
    'filename_pattern': '/plot_',       # Only execute files starting with 'plot_'
    'ignore_pattern': r'__init__\.py',  # Ignore specific files
    'download_all_examples': True,      # Allow downloading examples
}
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ['css/custom.css']

# -- Enable TODOs in documentation -------------------------------------------
todo_include_todos = True