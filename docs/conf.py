# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'eBiota'
copyright = '2024, Zhulab'
author = 'Zhulab'

release = '1.0'
version = '1.0'

# -- General configuration

extensions = [
        "sphinx_rtd_theme",
        'recommonmark',
        'sphinx_markdown_tables',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# needed by readthedocs
master_doc = 'index'