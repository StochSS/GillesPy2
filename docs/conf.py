# -*- coding: utf-8 -*-
# =============================================================================
# @file  conf.py
# @brief Configuration file for the Sphinx documentation builder.
# @license Please see the file named LICENSE in the project directory
# @website https://github.com/GillesPy2/GillesPy2
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation: http://www.sphinx-doc.org/en/master/config
# =============================================================================

import os
from   os import path
import sys

# -- Path setup --------------------------------------------------------------

ROOT_DIR = path.join(path.dirname(path.abspath(__file__)), '..')

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
#   sys.path.insert(0, '/mnt/c/Users/seanm/OneDrive/Documents/Research/documentation/GillesPy2')
sys.path.append('..')


# -- Project information -----------------------------------------------------

# The following reads the variables without doing an "import gillespy2",
# because the latter will cause the python execution environment to fail if
# any dependencies are not already installed -- negating most of the reason
# we're using setup() in the first place.  This code avoids eval, for security.

gillespy2 = {}
with open(path.join(ROOT_DIR, 'gillespy2/__version__.py')) as f:
    text = f.read().rstrip().splitlines()
    vars = [line for line in text if line.startswith('__') and '=' in line]
    for v in vars:
        setting = v.split('=')
        gillespy2[setting[0].strip()] = setting[1].strip().replace("'", '')


project   = gillespy2['__title__']
copyright = gillespy2['__copyright__']
author    = gillespy2['__author__']

# The short X.Y version
version   = gillespy2['__version__']
version   = version[:version.rfind('.')]

# The full version, including alpha/beta/rc tags
release   = gillespy2['__version__']


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
]

# Disable full module names in instances when the scope is already evident.
add_module_names = False

# Use short typehint format instead of long. This omits module names from types.
autodoc_typehints_format = "short"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = None


# -- Options for HTML output -------------------------------------------------

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'GillesPy2doc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'GillesPy2.tex', 'GillesPy2 Documentation',
     'Author', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, 'gillespy2', 'GillesPy2 Documentation',
     [author], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (master_doc, 'GillesPy2', 'GillesPy2 Documentation',
     author, 'GillesPy2', 'One line description of project.',
     'Miscellaneous'),
]


# -- Options for Epub output -------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = project

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#
# epub_identifier = ''

# A unique identification for the text.
#
# epub_uid = ''

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']


# -- Extension configuration -------------------------------------------------

import alabaster

html_theme_path = [alabaster.get_path()]
extensions.append('alabaster')
html_theme = 'alabaster'

html_theme_options = {
    'github_user': 'GillesPy2',
    'github_repo': 'GillesPy2',
    'fixed_sidebar': 'true',
    'github_banner': 'true',
    'logo': 'img/gillespy2-logo.svg',
    'touch_icon': 'img/gillespy2-logo.svg',
    'show_relbar_bottom': 'true',
    'link_hover': 'purple',
    'caption_font_size': '150%',
}

html_css_files = [
    'css/gillespy2_alabaster_customizations.css',
]

# This script inserts breaks after each parameter in a function definition.
html_js_files = [
    'js/gillespy2_alabaster_customizations.js',
]

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True


# -- Additional stuff for GillesPy2 

import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
