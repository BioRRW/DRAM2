from setuptools import setup, find_packages
from dram2 import __version__ as version
from os import path

__author__ = 'shafferm, rmflynn'
__version__ = version

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="DRAM2-bio",
    version=__version__,
    # scripts=['scripts/DRAM.py', 'scripts/DRAM-v.py', 'scripts/DRAM-setup.py'],
    packages=find_packages(),
    description="Distilled and Refined Annotation of Metabolism: A tool for the annotation and curation of function for"
                " microbial and viral genomes",
    long_description=long_description,
    long_description_content_type='text/markdown',  # Optional (see note above)
    package_data={'dram2': ['CONFIG']},
    python_requires='>=3',
    install_requires=['scikit-bio', 'pandas', 'altair', 'sqlalchemy', 'networkx', 'openpyxl', 'numpy', 'click'],
    entry_points={
        'console_scripts': [
            'dram2 = dram2.annotate_bins:dram2',
            'dram2-v = dram2.annotate_vgfs:dram2_v',
            'dram2-setup = dram2.database_setup.py:dram2_setup',
        ],
    },
    author="Michael Shaffer",
    author_email='michael.t.shaffer@colostate.edu',
    url="https://github.com/shafferm/DRAM/",
    download_url="https://github.com/shafferm/DRAM/tarball/%s" % __version__
)
