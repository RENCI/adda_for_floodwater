# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "ADCIRC Data Assimilator"
copyright = "2023, RENCI"
author = "Brian Blanton"
release = "0.2"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.githubpages",
    "sphinx.ext.viewcode",
    "sphinx.ext.todo",
    "sphinx.ext.imgmath",
    "sphinx.ext.mathjax",
    "sphinx_rtd_theme",
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "classic" # "sphinx_rtd_theme" # "classic"
html_static_path = []

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True
