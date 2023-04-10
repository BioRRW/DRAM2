"""
====================================
Test The Merging of DRAM Directories
====================================

use::

    conda activate ./dram2_env
    pip install -e ./
    pip uninstall dram2
    pip uninstall dram2.db_kits
    pytest tests/test_merger.py
    exit()
    continue
    pytest tests/test_merger.py::test_merger_pipe
    dram_dirs

We need to test:

    genes += [j for i in genes_dir for j in get_genes_from_dir(i)]
"""

import pytest
from shutil import rmtree, copy
import logging
from pathlib import Path
from glob import glob

from Bio import SeqIO
from click.testing import CliRunner

from dram2.cli.context import DramContext
from dram2.annotate import annotate_pipe
from dram2.call_genes import call_genes_pipe
from dram2.utils import DramUsageError
from dram2.merger import merger_pipe
from dram2.db_builder import db_builder_pipe

import pytest


def test_merger():
    pass


def test_list_dir():
    pass


@pytest.fixture()
def cant_hyd_db(tmp_path_factory):
    out = tmp_path_factory.mktemp("ch")
    output_yaml = out / "conf.yaml"
    context = DramContext(out, None)
    db_builder_pipe(context, ["cant_hyd"], True, (), output_yaml)
    return out


@pytest.fixture()
def dram_run_to_merge(tmp_path: Path, cant_hyd_db: Path):
    infasta = Path("tests/data/Cytophaga_hutchinsonii_ATCC_33406.fasta")
    fasta_dir = tmp_path / "fastas"
    fasta_dir.mkdir(exist_ok=True)
    dram_out = tmp_path / "dram"
    dram_out.mkdir(exist_ok=True)
    for i in range(1):
        outfasta = fasta_dir / f"test_{i}.fasta"
        copy(infasta, outfasta)
        dram_dir = dram_out / f"dram_dir_{i}"
        try:
            context = DramContext(dram_dir, cant_hyd_db / "conf.yaml")
            call_genes_pipe(context, [outfasta])
            annotate_pipe(context, [], use_db=["cant_hyd"])
        except DramUsageError:
            rmtree(dram_dir)
            continue
    return dram_out


def test_merger_pipe(tmp_path: Path, dram_run_to_merge: Path, cant_hyd_db: Path):
    context = DramContext(tmp_path, cant_hyd_db)
    dram_dirs = [Path(i) for i in glob(dram_run_to_merge.as_posix()+ "/*")]
    merger_pipe(context, dram_dirs)
    breakpoint()
