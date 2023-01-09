"""



"""

import pytest

import os
from io import StringIO
from functools import partial
from filecmp import cmp
from click.testing import CliRunner
from shutil import copy
from pathlib import Path
from dram2.utils.tests.test_utils import logger

import pandas as pd
import numpy as np
from skbio.io import read as read_sequence

from dram2.annotate.annotate import (
    filter_fasta,
    run_prodigal,
    get_gene_data,
    get_unannotated,
    assign_grades,
    generate_annotated_fasta,
    create_annotated_fasta,
    generate_renamed_fasta,
    rename_fasta,
    count_motifs,
    strip_endings,
    get_dups,
    annotate_wraper,
    annotate
)


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

def test_no_output(runner, fasta_loc):
    result = runner.invoke(annotate_wraper, fasta_loc)
    assert result.exit_code != 0
    assert "Error: Missing option \'-o\'" in result.output


def test_just_prodigal(runner, fasta_loc, tmp_path, faa_loc, logger):
    annotate([Path(fasta_loc)], logger, output_dir=tmp_path, force=True, keep_tmp= True)
    out_ex = Path(tmp_path, "working_dir", "NC_001422", "genes.faa")
    assert out_ex.exists()
    cmp(out_ex, faa_loc)


def test_just_prodigal_no_temp(runner, logger, fasta_loc, tmp_path):
    annotate([Path(fasta_loc)], logger, output_dir=tmp_path, force=True)
    out_ex = Path(tmp_path, "working_dir")
    assert not out_ex.exists()
 

def test_just_prodigal_cmd(runner, fasta_loc, tmp_path):
    result = runner.invoke(annotate_wraper, [fasta_loc, '-o', tmp_path, '-f'])
    assert result.exit_code == 0
    # assert result.output == 'Hello Peter!\n'
    out_ex = Path(tmp_path, "working_dir")
    assert not out_ex.exists()

