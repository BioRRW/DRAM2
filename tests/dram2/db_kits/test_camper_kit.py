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
)
from dram2.db_kits.camper_kit import (
    get_sig_row,
    hmmscan_formater,
    blast_search_formater,
)
from datetime import datetime

from dram2.annotate.annotate import annotate_wraper, annotate

from dram2.utils.tests.test_utils import logger
import json


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
    parsed_hit = parse_hmmsearch_domtblout(Path("tests", "data", "hmmsearch_hit.txt"))
    assert parsed_hit.shape == (1, 23)
    assert parsed_hit.loc[0, "query_id"] == "NP_040710.1"
    assert parsed_hit.loc[0, "query_length"] == 38
    assert parsed_hit.loc[0, "target_id"] == "Microvir_J"
    assert parsed_hit.loc[0, "target_ascession"] == "PF04726.13"
    assert parsed_hit.loc[0, "full_evalue"] == 6.900000e-31


@pytest.fixture()
def processed_hits():
    forward = os.path.join("tests", "data", "query_target_hits.b6")
    reverse = os.path.join("tests", "data", "target_query_hits.b6")
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


def test_annotate(tmpdir, runner, logger):
    tmp_out = tmpdir.mkdir("annotate_genes")
    input_faa = Path("tests", "data", "camper_test_genes.faa")
    output_dir = Path(tmp_out)

    config_loc = Path("tests/data/test_config_camper")
    config = json.loads(open(config_loc).read())
    annotate(
        [input_faa],
        logger=logger,
        output_dir=output_dir,
        config=config,
        # config={},
        # db_path="temp_db",
        use_db=["camper"],
        genes_called=True,
        force=True,
    )
    assert (output_dir / "annotations.tsv").exists()
    out = pd.read_csv((output_dir / "annotations.tsv"), sep="\t", index_col=0)
    expected = pd.read_csv(
        Path("tests", "data", "camper_test_annotations_one.tsv"), sep="\t", index_col=0
    )
    assert out.sort_index().equals(expected.sort_index())


@pytest.fixture()
def runner():
    return CliRunner()


def test_annotate_cmd(runner, tmpdir):
    tmp_out = tmpdir.mkdir("annotate_genes")
    input_faa = Path("tests", "data", "camper_test_genes.faa")
    output_dir = Path(tmp_out, "annotations")
    config = Path("tests/data/test_config_camper")
    result = runner.invoke(
        annotate_wraper,
        [
            input_faa.as_posix(),
            "-o",
            output_dir.as_posix(),
            "-f",
            "--config_loc",
            config.as_posix(),
            # "--db_path",
            # "temp_db",
            "--use_db",
            "camper",
            "--genes_called",
        ],
    )
    # uncoment soon
    # assert (output_dir / "annotations.tsv").exists()
    # out = pd.read_csv((output_dir / "annotations.tsv"), sep="\t", index_col=0)
    # expected= pd.read_csv(Path("tests", "data", "camper_test_annotations_one.tsv"), sep="\t", index_col=0)
    # assert out.sort_index().equals(expected.sort_index())
    # assert result.exit_code == 0
    # result.output


def test_hmmscan_formater():
    hits = pd.read_csv("tests/data/test_camper_hmmscan_formater_hits.csv", index_col=0)
    received = hmmscan_formater(
        hits=hits,
        db_name="CAMPER",
        hmm_info_path="tests/data/test_camper_hmmscan_formater_scores.tsv",
        top_hit=False,
    )
    expected = pd.read_csv(
        "tests/data/test_camper_hmmscan_formater_out_notop_hit.csv", index_col=0
    )
    assert received.equals(expected), "Something wrong with camper hmm"
    received = hmmscan_formater(
        hits=hits,
        db_name="CAMPER",
        hmm_info_path="tests/data/test_camper_hmmscan_formater_scores.tsv",
        top_hit=True,
    )
    expected = pd.read_csv(
        "tests/data/test_camper_hmmscan_formater_out.csv", index_col=0
    )
    assert received.equals(expected), "Something wrong with camper hmm"


def test_blast_search_formater(logger):
    received = blast_search_formater(
        hits_path="tests/data/test_camper_blast_search_formater_hits.tsv",
        db_name="CAMPER",
        logger=logger,
        info_db=pd.read_csv(
            "./tests/data/test_camper_blast_search_formater_scores.tsv",
            sep="\t",
            index_col=0,
        ),
    )
    expected = pd.read_csv(
        "./tests/data/test_camper_blast_search_formater_out.csv", index_col=0
    )
    assert received.equals(expected), "Something wrong with camper blast search"
