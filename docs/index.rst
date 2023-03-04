.. DRAM2 documentation master file, created by
   sphinx-quickstart on Wed Nov 30 14:52:41 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========================================
Welcome to the DRAM2 Docs!
=========================================

You are reading the official documentation for DRAM22 (Distilled and Refined Annotation of Metabolism, Version 2) a tool for annotating metagenomic assembled genomes. DRAM2 annotates MAGs  using a set of :ref:`databases<database_info>` including custom user databases. DRAM2 is run in two stages. First an annotation step to assign database identifiers to gene, and then a distill step to curate these annotations into useful functional categories. Additionally, viral contigs are further analyzed during to identify potential AMGs. This is done via assigning an auxiliary score and flags representing the confidence that a gene is both metabolic and viral.

For more detail on how the Science of DRAM2 works please see our `paper <>`_ or
read .

For information on how DRAM is changing, please read the most recent `release notes <https://github.com/WrightonLabCSU/DRAM/releases/latest>`_.

DRAM2 Development Note
----------------------

At this time, DRAM2 is only available internally
The DRAM development team is actively working on DRAM2. We do not anticipate adding any additional functionality to DRAM, i.e. DRAM1. Features requested for DRAM1 will be added to DRAM2, to the best of our ability and as appropriate.




.. toctree::
  :caption: Getting Started
  :name: getting_started
  :hidden:
  :maxdepth: 1

  getting_started/setup_configure_dram2
  getting_started/basic_use
  getting_started/dram1_to_dram2
  getting_started/commands
  getting_started/example_use_case

.. toctree::
  :caption: Science Guide
  :name: commands
  :hidden:
  :maxdepth: 1

  science/scientific_overview
  science/database

.. toctree::
  :caption: Command Guide
  :name: commands
  :hidden:
  :maxdepth: 1

  commands/dram2

.. toctree::
  :caption: Developer Documentaition
  :name: dev
  :hidden:
  :maxdepth: 0

  dev/docs
  dev/call_genes
  dev/annotate
  dev/db_builder
  dev/db_kits
  dev/distill
  dev/genbank
  dev/merger
  dev/rna
  dev/rule_adjectives
  dev/strainer
  dev/trash
  dev/tree_kit
  dev/utils
