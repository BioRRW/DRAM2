"""Setup file for package"""
from setuptools import setup, find_namespace_packages
from dram2.annotate import __version__
from os import path

__author__ = 'rmflynn'

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="dram2",
    version=__version__,
    # scripts=['scripts/DRAM.py', 'scripts/DRAM-v.py', 'scripts/DRAM-setup.py'],
    packages=find_namespace_packages(),
    description="Distilled and Refined Annotation of Metabolism: A tool for the annotation and curation of function for"
                " microbial and viral genomes",
    long_description=long_description,
    long_description_content_type='text/markdown',  # Optional (see note above)
    package_data={'dram2': ['CONFIG', "rule_adjectives/rules.tsv"]},
    python_requires='>=3.10',
    install_requires=['scikit-bio', 'pandas', 'altair', 'sqlalchemy', 'networkx', 'openpyxl', 'numpy', 'click', 'pytest', 'biopython'],
    entry_points={
        'console_scripts': [
            'dram2 = dram2.utils.command_line:dram2',
            # 'adj = dram2.rule_adjectives:evaluate',
            # 'tree = dram2.tree_kit.dram_phylo_pipe:tree_kit',

        ],
    },
    author="Rory Flynn",
    author_email='Rory.Flynn@colostate.edu',
    url="",  # this will change
    download_url="",
    include_package_data=True  # include all files in MANIFEST.in
)
