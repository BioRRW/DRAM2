"""

"""
import pytest
import logging
from dram2.utils.context import DramContext
from dram2.utils.utils import setup_logger

from dram2.db_kits.utils import (
    DRAM_DATAFOLDER_TAG,
)
from pathlib import Path

from click.testing import CliRunner


@pytest.fixture()
def dram_context(tmpdir):
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
def faa_loc():
    return Path("tests", "data", "NC_001422.faa")


@pytest.fixture()
def logger(tmpdir):
    logger = logging.getLogger("test_log")
    setup_logger(logger)
    return logger


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
def runner():
    return CliRunner()


def test_context(tmpdir):
    dram_context = DramContext(
        cores=10,
        db_path=tmpdir,
        config_file=Path("tests/data/no_path_config.yaml"),
        log_file_path=tmpdir / "log",
        output_dir=tmpdir / "out1",
        keep_tmp=False,
        verbose=0,
    )
    dram_context_absolute = DramContext(
        cores=10,
        db_path=tmpdir,
        config_file=Path("tests/data/test_dram_config_absolute.yaml"),
        log_file_path=tmpdir / "log",
        output_dir=tmpdir / "out1",
        keep_tmp=False,
        verbose=0,
    )
    config = dram_context.get_dram_config()
    config_abs = dram_context_absolute.get_dram_config()
    assert config[DRAM_DATAFOLDER_TAG] == Path("tests/data/").absolute()
    assert config_abs[DRAM_DATAFOLDER_TAG] == Path("/look/at/me/I/am/absolute")
