# -*- coding: utf-8 -*-
#
# hydromt_delwaq documentation build configuration file, created by
# sphinx-quickstart on Wed Jul 24 15:19:00 2019.
#
# This file is execfile()d with the current directory set to its
# containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import shutil
import hydromt
import hydromt_delwaq
from distutils.dir_util import copy_tree

here = os.path.dirname(__file__)
sys.path.insert(0, os.path.abspath(os.path.join(here, "..")))


def remove_dir_content(path: str) -> None:
    for root, dirs, files in os.walk(path):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))
    if os.path.isdir(path):
        shutil.rmtree(path)


# -- Project information -----------------------------------------------------

project = "HydroMT delwaq"
copyright = "Deltares"
author = "Hélène Boisgontier"

# The short version which is displayed
version = hydromt_delwaq.__version__

# # -- Copy notebooks to include in docs -------
if os.path.isdir("_examples"):
    remove_dir_content("_examples")
os.makedirs("_examples")
copy_tree("../examples", "_examples")


# -- General configuration ------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_design",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.todo",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.githubpages",
    "IPython.sphinxext.ipython_directive",
    "IPython.sphinxext.ipython_console_highlighting",
    "nbsphinx",
]

autosummary_generate = True
# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]
# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"
# The master toctree document.
master_doc = "index"

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = "en"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This patterns also effect to html_static_path and html_extra_path
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False


# -- Options for HTML output ----------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"
autodoc_member_order = "bysource"  # overwrite default alphabetical sort
autoclass_content = "both"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["theme-deltares.css"]
html_theme_options = {
    "show_nav_level": 2,
    "navbar_align": "content",
    "use_edit_page_button": True,
    "icon_links": [
        {
            "name": "Deltares",
            "url": "https://deltares.nl/en/",
            "icon": "_static/deltares-white.svg",
            "type": "local",
        },
    ],
    "external_links": [
        {
            "name": "HydroMT core",
            "url": "https://deltares.github.io/hydromt/latest/index.html",
        },
    ],
}

html_context = {
    "github_url": "https://github.com",  # or your GitHub Enterprise interprise
    "github_user": "Deltares",
    "github_repo": "hydromt_delwaq",
    "github_version": "new_docs",  # FIXME
    "doc_path": "docs",
}

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# This is required for the alabaster theme
# refs: http://alabaster.readthedocs.io/en/latest/installation.html#sidebars
# html_sidebars = {
#     "**": [
#         "relations.html",  # needs 'show_related': True theme option to display
#         "searchbox.html",
#     ]
# }


# -- Options for HTMLHelp output ------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = "hydromt_delwaq_doc"


# -- Options for LaTeX output ---------------------------------------------

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
    (
        master_doc,
        "hydromt_delwaq.tex",
        "HydroMT DELWAQ plugin Documentation",
        [author],
        "manual",
    ),
]


# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (master_doc, "hydromt_delwaq", "HydroMT DELWAQ plugin Documentation", [author], 1)
]


# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        master_doc,
        "hydromt_delwaq",
        "HydroMT DELWAQ plugin Documentation",
        author,
        "HydroMT DELWAQ plugin",
        "Build and analyze DELWAQ models like a data-wizard.",
        "Miscellaneous",
    ),
]
