"""
Test the calling of genes by dram2

TODO:

- test the coverage

"""

import pytest

import os
from io import StringIO
from functools import partial
from filecmp import cmp
from click.testing import CliRunner
from shutil import copy
from pathlib import Path

import pandas as pd
import numpy as np
from skbio.io import read as read_sequence

from dram2.call_genes import run_prodigal, filter_fasta
from dram2.utils.command_line import dram2
from dram2.utils.context import PROJECT_CONFIG_YAML_NAME

# from dram2.annotate import (
#     filter_fasta,
#     get_gene_data,
#     get_unannotated,
#     assign_grades,
#     generate_annotated_fasta,
#     create_annotated_fasta,
#     generate_renamed_fasta,
#     rename_fasta,
#     count_motifs,
#     strip_endings,
#     get_dups,
#     annotate_wraper,
#     annotate
# )


@pytest.fixture()
def logger(tmpdir):
    logger = logging.getLogger("test_log")
    setup_logger(logger)
    return logger


@pytest.fixture()
def prodigal_dir(fasta_loc, tmpdir, logger):
    prodigal_output = tmpdir.mkdir("prodigal_output")
    gff, fna, faa = run_prodigal(Path(fasta_loc), Path(prodigal_output), logger)
    return gff, fna, faa


@pytest.fixture()
def fasta_loc():
    return os.path.join("tests", "data", "NC_001422.fasta")


@pytest.fixture()
def faa_loc():
    return os.path.join("tests", "data", "NC_001422.faa")


@pytest.fixture()
def runner():
    return CliRunner()


def test_no_output(runner):
    result = runner.invoke(dram2, "call")
    assert result.exit_code != 0
    with pytest.raises(ValueError, match=r"You need to set an output directory*"):
        raise result.exception


def test_filter_fasta(fasta_loc, tmpdir):
    filt_fasta = Path(tmpdir) / "filtered_fasta.fasta"
    filter_fasta(Path(fasta_loc), output_loc=filt_fasta, min_len=5000)
    with open(filt_fasta) as f:
        next(f)
        filtered_seq = ""
        for line in f:
            filtered_seq += line.strip()
    assert len(filtered_seq) == 5386
    non_fasta = Path(tmpdir) / "filtered_non.fasta"
    filter_fasta(Path(fasta_loc), output_loc=non_fasta, min_len=6000)
    assert non_fasta.stat().st_size == 0
    filt_fasta = tmpdir.mkdir("test_filt").join("filtered_fasta.fasta")
    filter_fasta(Path(fasta_loc), min_len=5000, output_loc=Path(filt_fasta))
    assert os.path.isfile(filt_fasta)


def test_call_genes_cmd(runner, fasta_loc, tmp_path):
    result_first = runner.invoke(dram2, ["-o", tmp_path, "call", fasta_loc])
    result_colide = runner.invoke(dram2, ["-o", tmp_path, "call", fasta_loc])
    result_force = runner.invoke(dram2, ["-o", tmp_path, "call", "-f", fasta_loc])
    assert result_first.exit_code == 0
    assert result_colide.exit_code != 0
    assert result_force.exit_code == 0
    # assert result.output == 'Hello Peter!\n'
    out_ex = Path(tmp_path, PROJECT_CONFIG_YAML_NAME)
    assert out_ex.exists()
