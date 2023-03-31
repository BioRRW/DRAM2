"""
Test the calling of genes by dram2

Use::

    conda activate ./dram2_env
    pytest tests/test_call_genes.py
    exit
    pytest tests/test_call_genes.py::test_call_genes_pipe
    pytest tests/test_call_genes.py::test_call_genes_cmd

TODO:

- test the coverage

"""

import pytest

from io import StringIO
from click.testing import CliRunner
from pathlib import Path
import logging


from dram2.cli import dram2
from dram2.cli.context import PROJECT_META_YAML_NAME, DramContext
from dram2.call_genes import run_prodigal, filter_fasta, call_genes_pipe


TEST_FASTA = Path("tests", "data", "NC_001422.fasta")
TEST_FAA = Path("tests", "data", "NC_001422.faa")
LOGGER = logging.getLogger("test_log")


@pytest.fixture()
def prodigal_dir(TEST_FASTA, tmpdir):
    prodigal_output = tmpdir.mkdir("prodigal_output")
    gff, fna, faa = run_prodigal(Path(TEST_FASTA), Path(prodigal_output), LOGGER)
    return gff, fna, faa


@pytest.fixture()
def runner():
    return CliRunner()


def test_no_output(runner):
    result = runner.invoke(dram2, "call")
    assert result.exit_code != 0
    with pytest.raises(ValueError, match=r"You need to set an output directory*"):
        raise result.exception


def test_filter_fasta(tmpdir):
    filt_fasta = Path(tmpdir) / "filtered_fasta.fasta"
    filter_fasta(Path(TEST_FASTA), output_loc=filt_fasta, min_len=5000)
    with open(filt_fasta) as f:
        next(f)
        filtered_seq = ""
        for line in f:
            filtered_seq += line.strip()
    assert len(filtered_seq) == 5386
    non_fasta = Path(tmpdir) / "filtered_non.fasta"
    filter_fasta(Path(TEST_FASTA), output_loc=non_fasta, min_len=6000)
    assert non_fasta.stat().st_size == 0
    filt_fasta = tmpdir.mkdir("test_filt").join("filtered_fasta.fasta")
    filter_fasta(Path(TEST_FASTA), min_len=5000, output_loc=Path(filt_fasta))
    assert filt_fasta.isfile()


def test_call_genes_pipe(tmp_path):
    out_meta = call_genes_pipe(
        fasta_paths=[TEST_FASTA],
        context=DramContext(tmp_path, Path("tests/data/all_config.yaml")),
    )
    out_received = list(out_meta["genes_called"].values())[0]
    out_expected = {
        "min_contig_size": 2500,
        "prodigal_mode": "meta",
        "prodigal_trans_tables": "11",
        "working_dir": "genes",
        "fastas": ["NC_001422"],
    }
    assert out_expected == out_received
    out_fasta = out_meta["fastas"][0][1:-1]
    for i in out_fasta:
        assert (tmp_path / i).exists()


def test_call_genes_cmd(runner, tmp_path):
    result_first = runner.invoke(
        dram2, ["-d", tmp_path.as_posix(), "call", TEST_FASTA.as_posix()]
    )
    result_colide = runner.invoke(
        dram2, ["-d", tmp_path.as_posix(), "call", TEST_FASTA.as_posix()]
    )
    result_force = runner.invoke(
        dram2, ["-d", tmp_path.as_posix(), "call", "-f", TEST_FASTA.as_posix()]
    )
    assert result_first.exit_code == 0, result_first
    assert result_colide.exit_code != 0, result_colide
    assert result_force.exit_code == 0, result_force
    # assert result.output == 'Hello Peter!\n'
    out_ex = Path(tmp_path, PROJECT_META_YAML_NAME)
    assert out_ex.exists()
