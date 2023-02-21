import pytest
import logging
from dram2.utils.context import DramContext
from dram2.utils.utils import setup_logger
from dram2.utils.command_line import dram2

from dram2.db_kits.camper_kit import CamperKit, hmmscan_formater
from pathlib import Path

from dram2.db_kits.utils import Fasta
from dram2.annotate import annotate, make_mmseqs_db_for_fasta
from click.testing import CliRunner
from dram2.db_kits.camper_kit import blast_search_formater
import pandas as pd


@pytest.fixture()
def dram_context(tmpdir):
    dram_context = DramContext(
        cores=10,
        db_path=tmpdir,
        config_file=Path("tests/data/all_config.yaml"),
        log_file_path=tmpdir / "log",
        output_dir=tmpdir / "out1",
        keep_tmp=False,
        verbose=0,
    )
    return dram_context

@pytest.fixture()
def dram_context_bad(tmpdir):
    dram_context = DramContext(
        cores=10,
        db_path=tmpdir,
        config_file=Path("tests/data/no_path_config.yaml"),
        log_file_path=tmpdir / "log",
        output_dir=tmpdir / "out1",
        keep_tmp=False,
        verbose=0,
    )
    return dram_context

@pytest.fixture()
def logger(tmpdir):
    logger = logging.getLogger("test_log")
    setup_logger(logger)
    return logger


@pytest.fixture()
def runner():
    return CliRunner()


@pytest.fixture()
def db_args_dict(logger, tmpdir):
    working_dir = Path(tmpdir) / "working"
    working_dir.mkdir()
    db_args = {
        "logger": logger,
        "output_dir": tmpdir,
        "working_dir": working_dir,
        "bit_score_threshold": 1,
        "rbh_bit_score_threshold": 1,
        "kofam_use_dbcan2_thresholds": False,
        "threads": 2,
        "force": False,
        "extra": False,
        "db_path": tmpdir / "data",
    }
    return db_args

@pytest.fixture()
def tmp_camper_fasta(logger, tmp_path) -> Fasta:
    input_faa = Path("tests", "data", "camper_test_genes.faa")
    fasta = Fasta("test", None, Path(tmp_path), input_faa, None, None, None)
    return make_mmseqs_db_for_fasta(fasta, logger, 1)


def test_check_setup_fail(dram_context_bad, db_args_dict):
    with pytest.raises(
        FileNotFoundError,
        match=(
            r"The file .* is not at the path"
            r" .*. Most likely you moved the DRAM"
            r" data but forgot to update the config file to"
            r" point to it. The easy fix is to set the*"
        ),
    ):
        kit = CamperKit(dram_context_bad.get_dram_config(), db_args_dict)


def test_camper_blast_search_formater(logger):
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


def test_check_setup_pass(dram_context, db_args_dict):
    kit = CamperKit(dram_context.get_dram_config(), db_args_dict)


def test_camper_hmmscan_formater():
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


def test_search(dram_context, tmp_camper_fasta, db_args_dict):

    kit = CamperKit(dram_context.get_dram_config(), db_args_dict)
    out_received = kit.search(tmp_camper_fasta)
    output_dir = dram_context.get_output_dir()
    # assert pd.read_csv(
    #     Path(output_dir, "annotations.tsv"), sep="\t", index_col=0
    # ).equals(
    #     pd.read_csv(
    #         Path("tests", "data", "camper_test_annotations.tsv"),
    #         sep="\t",
    #         index_col=0,
    #     )
    # )


def test_search_cmd(tmpdir, runner):
    input_faa = Path("tests", "data", "camper_test_genes.faa")
    input_faa_empty = Path("tests", "data", "camper_test_genes_empty.faa")
    output_dir = Path(tmpdir, "annotations")
    output_dir_empty = Path(tmpdir, "annotations_empty")
    result1 = runner.invoke(
        dram2, ["-o", output_dir, "annotate", "--use_db", "camper", input_faa]
    )
    # assert result1.exit_code == 0
    result = runner.invoke(
        dram2, ["-o", output_dir_empty, "annotate", "--use_db", "camper", input_faa]
    )
    # breakpoint()
    # assert pd.read_csv(
    #     os.path.join(output_dir, "annotations.tsv"), sep="\t", index_col=0
    # ).equals(
    #     pd.read_csv(
    #         os.path.join("tests", "data", "camper_test_annotations.tsv"),
    #         sep="\t",
    #         index_col=0,
    #     )
    # )
