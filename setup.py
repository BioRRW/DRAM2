"""Setup file for package"""
from setuptools import setup, find_namespace_packages
from os import path

__version__ = "1.0.b1"

__author__ = "rmflynn"

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="dram2",
    version=__version__,
    packages=find_namespace_packages(
        include=["dram2.*", "dram2.db_kits.*"],
    ),
    description=(
        "Distilled and Refined Annotation of Metabolism: A tool for the"
        " annotation and curation of function for"
        " microbial and viral genomes"
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",  # Optional (see note above)
    package_data={
        "dram2.rule_adjectives": ["rules.tsv"],
        "dram2.tree_kit": [
            "data/nxr_nar/color_map.tsv",
            "data/nxr_nar/nxr-nar-tree-mapping.tsv",
            "data/pmoa_amoa/pmoa_amoa-tree-mapping.tsv",
        ],
        "dram2.distill": [
            "data/amg_database.tsv",
            "data/etc_module_database.tsv",
            "data/genome_summary_form.tsv",
            "data/module_step_form.tsv",
            "data/function_heatmap_form.tsv",
        ],
        "dram2.db_kits.methyl_kit": [
            "methylotrophy.faa",
            "methylotrophy_distillate.tsv",
        ],
        "dram2.db_kits.cant_hyd_kit": [
            "data/BacMet_ExpVerified_BiocideRes_genes_SHORT.faa",
            "data/CANT_HYD_BLAST_scores.csv",
            "data/CANT_HYD_HMM_scores.csv",
            "data/engineeredsys_dram_module.tsv",
        ],
    },
    # package_dir={'': ''},
    python_requires=">=3.10",
    install_requires=[
        "scikit-bio",
        "pandas",
        "altair",
        "sqlalchemy",
        "networkx",
        "openpyxl",
        "numpy",
        "click",
        "pytest",
        "biopython",
    ],
    entry_points={
        "console_scripts": [
            "dram2 = dram2.cli:dram2",
            # 'adj = dram2.rule_adjectives:evaluate',
            # 'tree = dram2.tree_kit.dram_phylo_pipe:tree_kit',
        ],
    },
    author="Rory Flynn",
    author_email="Rory.Flynn@colostate.edu",
    url="",  # this will change
    download_url="",
    include_package_data=True,  # include all files in MANIFEST.in
)
