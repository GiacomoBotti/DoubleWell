# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

project = 'DoubleWell'
copyright = '2026, Giacomo Botti, Sophya Garashchuk'
author = 'Giacomo Botti, Sophya Garashchuk'
release = '1'
version = '1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinxfortran_ng.fortran_domain",
    "sphinxfortran_ng.fortran_autodoc",
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
#    "matplotlib.sphinxext.plot_directive",
]

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]

fortran_src = [os.path.abspath("../../*.f90")]

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']
