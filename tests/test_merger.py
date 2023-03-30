"""
====================================
Test The Merging of DRAM Directories
====================================

Using


We need to test:

"""

import pytest
import logging
from pathlib import Path
from glob import glob

from Bio import SeqIO
from click.testing import CliRunner

from dram2.cli.context import DramContext
from dram2.annotate import annotate_cmd
from dram2.call_genes import call_genes_cmd


import pytest


def test_merger():
    pass


def test_list_dir():
    pass


def test_make_example_dram_run_to_merge(tmp_path: Path):
    context = DramContext(Path("tests/data/all_config.yaml"), tmp_path)
    infasta = Path("tests/data/NC_001422.fasta")
    out_fasta = Path(tmp_path, "split_fastas")
    out_fasta.mkdir(exist_ok=True)
    for i, line in enumerate(SeqIO.parse(open(infasta), "fasta")):
        with (out_fasta / f"NC_001422_{i}.fasta") as fasta_out:
            _ = SeqIO.write(line, fasta_out, "fasta")
        call_genes_cmd(context, glob(f"{out_fasta.as_posix}/*.fasta"))
        annotate_cmd(context)
