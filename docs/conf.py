# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Lible'
copyright = '2024, MihkuU'
author = 'MihkuU'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['breathe']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

breathe_default_project = "Lible"

source_suffic = '.rst'

master_doc = 'index'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']

# -- Read the Docs Shenanigans -----------------------------------------------
import subprocess, os

def configureDoxyfile(input_dir, output_dir):
	with open('Doxyfile.in', 'r') as file:
		filedata = file.read()

	filedata = filedata.replace('@doxygen_input_dir@', input_dir)
	filedata = filedata.replace('@doxygen_output_dir@', output_dir)

	print("configureDoxyfile() called:")
	print("os.listdir(): ", os.listdir())
	print("input_dir: ", input_dir,)

	with open('Doxyfile', 'w') as file:
		file.write(filedata)

# Check if we're running on 'Read the Docs' servers
read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

breathe_projects = {}

if read_the_docs_build:
	input_dir = '../src/lible/ints/'
	output_dir = 'build'
	configureDoxyfile(input_dir, output_dir)
	subprocess.call('doxygen', shell=True)
	breathe_projects['Lible'] = output_dir + '/xml'
