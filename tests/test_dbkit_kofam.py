"""
=================
Test DBKit: KOfam
=================

DBKits have a range of functions that can be unique as you like. All of them
need test files like this one

Here we test:
    - the stand alone annotation command

"""
import pytest
import logging
from pathlib import Path
from click.testing import CliRunner

from dram2.cli.context import DramContext
from dram2.cli import dram2
from dram2.db_kits.utils import Fasta
from dram2.db_kits.kofam_kit import KOfamKit


@pytest.fixture()
def dram_context(tmpdir):
    dram_context = DramContext(
        cores=10,
        db_path=tmpdir,
        config_file=Path("tests/data/all_config.yaml"),
        log_file_path=tmpdir / "log",
        output_dir=tmpdir / "out1",
        keep_tmp=False,
        verbose=int(3),
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
        verbose=3,
    )
    return dram_context


@pytest.fixture()
def logger(tmpdir):
    logger = logging.getLogger("test_log")
    return logger


@pytest.fixture()
def runner():
    return CliRunner()


@pytest.fixture()
def test_fasta(logger, tmp_path) -> Fasta:
    input_faa = Path("tests", "data", "camper_test_genes.faa")
    fasta = Fasta("test", None, Path(tmp_path), input_faa, None, None, None)
    return make_mmseqs_db_for_fasta(fasta, logger, 1)


@pytest.fixture()
def db_args_dict(logger, tmpdir):
    working_dir = Path(tmpdir) / "working"
    working_dir.mkdir()
    db_args = {
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


def test_check_setup_fail(dram_context_bad, logger):
    with pytest.raises(
        FileNotFoundError,
        match=(
            r"The file .* is not at the path"
            r" .*. Most likely you moved the DRAM"
            r" data but forgot to update the config file to"
            r" point to it. The easy fix is to set the*"
        ),
    ):
        kit = KOfamKit(dram_context_bad.get_dram_config(logger), logger)


def test_check_setup_pass(dram_context, db_args_dict):
    kit = KOfamKit(dram_context.get_dram_config(), db_args_dict)


def test_search(dram_context, test_fasta, db_args_dict):

    kit = KOfamKit(dram_context.get_dram_config(), db_args_dict)
    out_received = kit.search(test_fasta)
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
        dram2, ["-o", output_dir, "annotate", "--use_db", "kofam", input_faa]
    )
    # assert result1.exit_code == 0
    result = runner.invoke(
        dram2, ["-o", output_dir_empty, "annotate",
                "--use_db", "kofam", input_faa]
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
