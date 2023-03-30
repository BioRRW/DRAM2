"""
Test the CLI and DRAM Context
=============================

conda activate dram2_env/

pytest tests/test_cli.py
pytest tests/test_cli.py::test_log_error_wraper
pytest tests/test_call_genes.py

exit

"""

import logging
from pathlib import Path

import pytest
from click.testing import CliRunner

from dram2.cli import dram2
from dram2.cli.context import (
    get_config_path,
    get_time_stamp_id,
    get_new_config_path,
    log_error_wraper,
    DramContext,
    USER_CONFIG,
    GLOBAL_CONFIG,
    LOG_FILE_NAME,
)


@pytest.fixture()
def runner():
    return CliRunner()


def test_log_error_wraper(tmp_path):
    def dev_zero(y):
        x = y / 0
        return x

    with pytest.raises(ZeroDivisionError):
        context = DramContext(tmp_path, None)
        log_error_wraper(dev_zero, context)(1)


def test_help(runner):
    result = runner.invoke(dram2, "--help")
    assert result.exit_code == 0


def test_get_config_path():
    assert Path("tests/data/all_config.yaml") == get_config_path(
        logging.getLogger("test"), Path("tests/data/all_config.yaml")
    )


def test_get_new_config_path():
    get_new_config_path()


def test_context(tmp_path):
    context = DramContext(tmp_path, None)
    logger = context.get_logger()
    assert isinstance(logger, logging.Logger)
    log_file = tmp_path / LOG_FILE_NAME
    assert log_file.exists
    assert tmp_path == context.get_dram_dir()
    assert {} == context.get_project_meta()
    config = {"test": "test"}
    context.set_project_meta(config)
    assert config == context.get_project_meta()
