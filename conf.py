# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Technical Notes'
copyright = '2025, Quemix inc'
author = 'Jun-Ichi Iwata'
version = '2025.2.14'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.mathjax', 'sphinx.ext.githubpages']
# extensions = ['sphinx.ext.mathjax', 'sphinx_rtd_theme']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'ja'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
# html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_static_path = ['_static']
html_style = 'css/my_theme.css'
#html_logo = '_static/Quloud_A_color.png'
html_title = ''
html_short_title = 'aa'
html_theme_options = {
    # 'logo_only': True,
    # 'sticky_navigation': True,
    # 'titles_only': True,
    # 'display_version': True,
    # 'flyout_display': 'attached',
    # 'version_selector': True,
    # 'language_selector': True,
    # 'collapse_navigation': True,
}
