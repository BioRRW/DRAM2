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

import pandas as pd
import numpy as np
from skbio.io import read as read_sequence
from dram2.db_kits.utils import (
    process_reciprocal_best_hits,
    get_sig_row,
    parse_hmmsearch_domtblout,
    do_blast_style_search,
    get_basic_description,
    make_mmseqs_db
)

from dram2.annotate.annotate import annotate_wraper, annotate

from dram2.utils.tests.test_utils import logger


@pytest.fixture()
def fasta_loc():
    return os.path.join("tests", "data", "NC_001422.fasta")


@pytest.fixture()
def faa_loc():
    return os.path.join("tests", "data", "NC_001422.faa")


@pytest.fixture()
def runner():
    return CliRunner()


def test_parse_hmmsearch_domtblout():
    parsed_hit = parse_hmmsearch_domtblout(
        Path("tests", "data", "hmmsearch_hit.txt")
    )
    assert parsed_hit.shape == (1, 23)
    assert parsed_hit.loc[0, "query_id"] == "NP_040710.1"
    assert parsed_hit.loc[0, "query_length"] == 38
    assert parsed_hit.loc[0, "target_id"] == "Microvir_J"
    assert parsed_hit.loc[0, "target_ascession"] == "PF04726.13"
    assert parsed_hit.loc[0, "full_evalue"] == 6.900000e-31


@pytest.fixture()
def processed_hits():
    forward = os.path.join(
        "tests", "data", "query_target_hits.b6"
    )
    reverse = os.path.join(
        "tests", "data", "target_query_hits.b6"
    )
    processed_hits = process_reciprocal_best_hits(forward, reverse)
    return processed_hits


def test_process_reciprocal_best_hits(processed_hits):
    assert processed_hits.shape == (7, 5)
    assert set(processed_hits.loc[processed_hits.target_RBH].index) == {
        "NC_001422.1_5",
        "NC_001422.1_4",
        "NC_001422.1_7",
        "NC_001422.1_6",
    }


def test_get_sig_row():
    names = ["target_start", "target_end", "target_length", "full_evalue"]
    assert not get_sig_row(pd.Series([1, 85, 100, 1], index=names))
    assert get_sig_row(pd.Series([1, 86, 100, 1e-16], index=names))
    assert not get_sig_row(pd.Series([1, 29, 100, 1e-20], index=names))


@pytest.fixture()
def phix_proteins():
    return os.path.join("tests", "data", "NC_001422.faa")


def test_get_basic_description():
    header_dict = {
        "YP_009015653.1": "YP_009015653.1 gp350 [Bacillus virus G]",
        "NP_077550.1": "NP_077550.1 EsV-1-65 [Ectocarpus siliculosus virus 1]",
    }
    viral_hits_data = [["NP_077550.1"], ["YP_009015653.1"]]
    viral_hits = pd.DataFrame(
        viral_hits_data, index=["gene1", "gene2"], columns=["viral_hit"]
    )
    viral_hits_add_description = get_basic_description(viral_hits, header_dict)
    assert viral_hits_add_description.shape == (2, 2)
    assert viral_hits_add_description.loc["gene1", "viral_id"] == "NP_077550.1"

class FakeDatabaseHandler:
    @staticmethod
    def get_database_names():
        return []


@pytest.fixture()
def mmseqs_db_dir(tmpdir):
    output_loc = tmpdir.mkdir('make_mmseqs_db_test')
    return output_loc



@pytest.fixture()
def phix_proteins():
    return os.path.join('tests', 'data', 'NC_001422.faa')




@pytest.fixture()
def faa_loc():
    return os.path.join("tests","data", "NC_001422.faa")

def test_just_custom_fasta(fasta_loc, tmp_path, faa_loc, logger):
    annotate(
        fasta_paths=[Path(fasta_loc)],
        output_dir=Path(tmp_path),
        force=True,
        custom_fasta_db_loc=[Path(faa_loc)],
        custom_fasta_db_name=["test"],
        logger=logger,
    )
    assert (Path(tmp_path) / "annotations.tsv").exists()


def test_just_custom_fasta_cmd(runner, fasta_loc, tmp_path, faa_loc):
    result = runner.invoke(
        annotate_wraper,
        [
            fasta_loc,
            "-o",
            tmp_path,
            "-f",
            "--custom_fasta_db_loc",
            faa_loc,
            "--custom_fasta_db_name",
            "test",
        ],
    )
    assert (Path(tmp_path) / "annotations.tsv").exists()
    out = pd.read_csv((Path(tmp_path) / "annotations.tsv"), sep="\t", index_col=0)
    expected= pd.read_csv(Path("tests", "data", "annotations_custom_fasta_NC_001422_with_self.tsv"), sep="\t", index_col=0)
    assert out.sort_index().equals(expected.sort_index())
    assert result.exit_code == 0
    assert "Getting forward best hits from test" in result.output
    assert "Getting reverse best hits from test" in result.output
    assert "Getting descriptions of hits from test" in result.output
    assert not Path(tmp_path, "working_dir").exists()



def test_custom_fasta_fail(runner, fasta_loc, tmp_path):
    result = runner.invoke(
        annotate_wraper,
        [fasta_loc, "-o", tmp_path, "-f", "--custom_fasta_db_name", "test"],
    )
    assert result.exit_code != 0

def test_custom_hmm_cmd(runner, fasta_loc, tmp_path):
    result = runner.invoke(
        annotate_wraper,
        [
            fasta_loc,
            "-o",
            tmp_path,
            "-f",
            "--custom_hmm_db_loc",
            Path("tests", "data", "camper_small.hmm"),
            "--custom_hmm_db_name",
            "test",
        ],
    )
    assert (Path(tmp_path) / "annotations.tsv").exists()
    out = pd.read_csv((Path(tmp_path) / "annotations.tsv"), sep="\t", index_col=0)
    # expected= pd.read_csv(Path("tests", "data", "annotations_custom_fasta_NC_001422_with_self.tsv"), sep="\t", index_col=0)
    #assert out.equals(expected)
    assert result.exit_code == 0
    assert "Pre-processing custom hmm database test" in result.output
    assert "Annotating custom hmm database test" in result.output
    assert not Path(tmp_path, "working_dir").exists()

