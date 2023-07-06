# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
import shutil
import glob

#We need this later when autobuilding API Docs
#sys.path.insert(0, os.path.abspath('../../pyiron_atomistics/'))

#Of course we should change this
project = 'pyscal_rdf'
copyright = '2023, Sarath Menon, Abril Guzman'
author = 'Sarath Menon, Abril Guzman'

def skip(app, what, name, obj, would_skip, options):
    if name in ( '__init__',):
        return False
    return would_skip

def setup(app):
    app.connect('autodoc-skip-member', skip)

#You can setup ipynb notebooks in the examples folder and they
#would be copied when the docs are built
if os.path.exists("examples"):
    shutil.rmtree("examples")
shutil.copytree("../../examples", "examples")

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    "sphinx.ext.extlinks", 
    "sphinx.ext.intersphinx", 
    "sphinx.ext.todo", 
    "sphinx.ext.viewcode",
    "myst_parser", 
    "sphinx_copybutton", 
    "sphinx_design", 
    "sphinx_inline_tabs",
    "nbsphinx",
]

myst_enable_extensions = ["dollarmath", "amsmath"]
#html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'
#html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_logo = "_static/logo.png"
html_theme_options = {
    'logo_only' : True,
    #'canonical_url' : 'https://calphy.readthedocs.io/',
}

html_extra_path = ['_static' ]

source_suffix = ['.rst', '.md']

exclude_patterns = []

