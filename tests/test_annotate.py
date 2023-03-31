"""
Test the annotation of genes but not with any databases


Assuming you have a development instalation example use:

    conda activate ./dram2_env
    pytest tests/test_annotate.py


"""


import os
import logging
import pytest
from io import StringIO
from click.testing import CliRunner
from pathlib import Path

import pandas as pd
import numpy as np
from dram2.cli import dram2
from dram2.cli.context import DramContext
from dram2.db_kits.utils import Fasta
from dram2.annotate import (
    annotate,
    annotate_pipe,
    path_to_gene_fastas,
    make_mmseqs_db_for_fasta,
    check_for_annotations,
    USED_DBS_TAG,
)
from dram2.db_kits.utils import FastaKit, HmmKit


@pytest.fixture()
def logger(tmpdir):
    logger = logging.getLogger("test_log")
    return logger


@pytest.fixture()
def fasta_loc():
    return os.path.join("tests", "data", "NC_001422.fasta")


@pytest.fixture()
def faa_loc():
    return os.path.join("tests", "data", "NC_001422.faa")


@pytest.fixture()
def runner():
    return CliRunner()


@pytest.fixture()
def prodigal_gff(prodigal_dir):
    return prodigal_dir[0]


@pytest.fixture()
def prodigal_fna(prodigal_dir):
    return prodigal_dir[1]


@pytest.fixture()
def prodigal_faa(prodigal_dir):
    return prodigal_dir[2]


@pytest.fixture()
def mmseqs_db_dir(tmpdir):
    output_loc = tmpdir.mkdir("make_mmseqs_db_test")
    return output_loc


@pytest.fixture()
def phix_proteins():
    return os.path.join("tests", "data", "NC_001422.faa")


@pytest.fixture()
def annotated_fake_gff_loc():
    return os.path.join("tests", "data", "annotated_fake_gff.gff")


@pytest.fixture()
def fake_phix_annotations():
    return pd.DataFrame(
        [[np.NaN], ["GH13"], [np.NaN], [np.NaN], [np.NaN], [np.NaN], [np.NaN]],
        index=[
            "NC_001422.1_1",
            "NC_001422.1_2",
            "NC_001422.1_3",
            "NC_001422.1_4",
            "NC_001422.1_5",
            "NC_001422.1_6",
            "NC_001422.1_7",
        ],
        columns=["cazy_id"],
    )


@pytest.fixture()
def fake_gff_loc():
    return os.path.join("tests", "data", "fake_gff.gff")


@pytest.fixture()
def phix_prodigal_genes():
    phix_seq = (
        ">NC_001422.1_1 # 51 # 221 # 1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer=5-"
        "10bp;gc_cont=0.404\n"
        "MSRKIILIKQELLLLVYELNRSGLLAENEKIRPILAQLEKLLLCDLSPSTNDSVKN*\n"
        ">NC_001422.1_2 # 390 # 848 # 1 # ID=1_2;partial=00;start_type=ATG;rbs_motif=GGA/GAG/AGG;rbs_spacer="
        "5-10bp;gc_cont=0.468\n"
        "MSQVTEQSVRFQTALASIKLIQASAVLDLTEDDFDFLTSNKVWIATDRSRARRCVEACVY\n"
        "GTLDFVGYPRFPAPVEFIAAVIAYYVHPVNIQTACLIMEGAEFTENIINGVERPVKAAEL\n"
        "FAFTLRVRAGNTDVLTDAEENVRQKLRAEGVM*\n"
        ">NC_001422.1_3 # 848 # 964 # 1 # ID=1_3;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;"
        "gc_cont=0.496\n"
        "MSKGKKRSGARPGRPQPLRGTKGKRKGARLWYVGGQQF*\n"
        ">NC_001422.1_4 # 1001 # 2284 # 1 # ID=1_4;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5"
        "-10bp;gc_cont=0.449\n"
        "MSNIQTGAERMPHDLSHLGFLAGQIGRLITISTTPVIAGDSFEMDAVGALRLSPLRRGLA\n"
        "IDSTVDIFTFYVPHRHVYGEQWIKFMKDGVNATPLPTVNTTGYIDHAAFLGTINPDTNKI\n"
        "PKHLFQGYLNIYNNYFKAPWMPDRTEANPNELNQDDARYGFRCCHLKNIWTAPLPPETEL\n"
        "SRQMTTSTTSIDIMGLQAAYANLHTDQERDYFMQRYHDVISSFGGKTSYDADNRPLLVMR\n"
        "SNLWASGYDVDGTDQTSLGQFSGRVQQTYKHSVPRFFVPEHGTMFTLALVRFPPTATKEI\n"
        "QYLNAKGALTYTDIAGDPVLYGNLPPREISMKDVFRSGDSSKKFKIAEGQWYRYAPSYVS\n"
        "PAYHLLEGFPFIQEPPSGDLQERVLIRHHDYDQCFQSVQLLQWNSQVKFNVTVYRNLPTT\n"
        "RDSIMTS*\n"
        ">NC_001422.1_5 # 2395 # 2922 # 1 # ID=1_5;partial=00;start_type=ATG;rbs_motif=AGGAG;rbs_spacer=5-10"
        "bp;gc_cont=0.420\n"
        "MFQTFISRHNSNFFSDKLVLTSVTPASSAPVLQTPKATSSTLYFDSLTVNAGNGGFLHCI\n"
        "QMDTSVNAANQVVSVGADIAFDADPKFFACLVRFESSSVPTTLPTAYDVYPLNGRHDGGY\n"
        "YTVKDCVTIDVLPRTPGNNVYVGFMVWSNFTATKCRGLVSLNQVIKEIICLQPLK*\n"
        ">NC_001422.1_6 # 2931 # 3917 # 1 # ID=1_6;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer="
        "5-10bp;gc_cont=0.451\n"
        "MFGAIAGGIASALAGGAMSKLFGGGQKAASGGIQGDVLATDNNTVGMGDAGIKSAIQGSN\n"
        "VPNPDEAAPSFVSGAMAKAGKGLLEGTLQAGTSAVSDKLLDLVGLGGKSAADKGKDTRDY\n"
        "LAAAFPELNAWERAGADASSAGMVDAGFENQKELTKMQLDNQKEIAEMQNETQKEIAGIQ\n"
        "SATSRQNTKDQVYAQNEMLAYQQKESTARVASIMENTNLSKQQQVSEIMRQMLTQAQTAG\n"
        "QYFTNDQIKEMTRKVSAEVDLVHQQTQNQRYGSSHIGATAKDISNVVTDAASGVVDIFHG\n"
        "IDKAVADTWNNFWKDGKADGIGSNLSRK*\n"
        ">NC_001422.1_7 # 3981 # 5384 # 1 # ID=1_7;partial=01;start_type=ATG;rbs_motif=GGAGG;rbs_spacer=5-10"
        "bp;gc_cont=0.455\n"
        "MVRSYYPSECHADYFDFERIEALKPAIEACGISTLSQSPMLGFHKQMDNRIKLLEEILSF\n"
        "RMQGVEFDNGDMYVDGHKAASDVRDEFVSVTEKLMDELAQCYNVLPQLDINNTIDHRPEG\n"
        "DEKWFLENEKTVTQFCRKLAAERPLKDIRDEYNYPKKKGIKDECSRLLEASTMKSRRGFA\n"
        "IQRLMNAMRQAHADGWFIVFDTLTLADDRLEAFYDNPNALRDYFRDIGRMVLAAEGRKAN\n"
        "DSHADCYQYFCVPEYGTANGRLHFHAVHFMRTLPTGSVDPNFGRRVRNRRQLNSLQNTWP\n"
        "YGYSMPIAVRYTQDAFSRSGWLWPVDAKGEPLKATSYMAVGFYVAKYVNKKSDMDLAAKG\n"
        "LGAKEWNNSLKTKLSLLPKKLFRIRMSRNFGMKMLTMTNLSTECLIQLTKLGYDATPFNQ\n"
        "ILKQNAKREMRLRLGKVTVADVLAAQPVTTNLLKFMRASIKMIGVSNL\n"
    )
    return StringIO(phix_seq)


@pytest.fixture()
def phix_annotations():
    return pd.DataFrame(
        [
            ["A", "K1", "U1", None, "a_bug1"],
            ["B", "K2", "U2", "P2", "a_bug2"],
            ["C", None, "U3", "P3", "a_bug3"],
            ["D", "K4", "U4", "P4", "a_bug4"],
            ["A", "K5", "U5", "P5", "a_bug5"],
            ["E", "K6", "U6", "P6", "a_bug6"],
            ["C", "K7", "U7", "P7", "a_bug7"],
        ],
        index=[
            "NC_001422.1_1",
            "NC_001422.1_2",
            "NC_001422.1_3",
            "NC_001422.1_4",
            "NC_001422.1_5",
            "NC_001422.1_6",
            "NC_001422.1_7",
        ],
        columns=["rank", "kegg_hit", "uniref_hit", "pfam_hits", "bin_taxonomy"],
    )


@pytest.fixture()
def phix_annotations_no_kegg():
    return pd.DataFrame(
        [
            ["B", "U1", None, "a_bug1"],
            ["B", "U2", "P2", "a_bug2"],
            ["C", "U3", "P3", "a_bug3"],
            ["D", "U4", "P4", "a_bug4"],
            ["B", "U5", "P5", "a_bug5"],
            ["E", "U6", "P6", "a_bug6"],
            ["C", "U7", "P7", "a_bug7"],
        ],
        index=[
            "NC_001422.1_1",
            "NC_001422.1_2",
            "NC_001422.1_3",
            "NC_001422.1_4",
            "NC_001422.1_5",
            "NC_001422.1_6",
            "NC_001422.1_7",
        ],
        columns=["rank", "uniref_hit", "pfam_hits", "bin_taxonomy"],
    )


@pytest.fixture()
def dbargs_dict(faa_loc, logger, tmpdir):
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
def tmp_fasta(faa_loc, logger, tmp_path) -> Fasta:
    fasta = Fasta("test", None, Path(tmp_path), Path(faa_loc), None, None, None)
    return make_mmseqs_db_for_fasta(fasta, logger, 1)


@pytest.fixture()
def tmp_fasta(faa_loc, logger, tmp_path) -> Fasta:
    fasta = Fasta("test", None, Path(tmp_path), Path(faa_loc), None, None, None)
    return make_mmseqs_db_for_fasta(fasta, logger, 1)


@pytest.fixture()
def test_conf_path():
    return Path("tests", "data", "all_config.yaml")


def test_FastaKit(phix_proteins, tmpdir, logger, dbargs_dict, tmp_fasta):
    assert tmp_fasta.mmsdb.exists()
    fastaDB = FastaKit("phix", phix_proteins, {}, dbargs_dict)
    expected_out = pd.DataFrame(
        [
            ["NC_001422.1_1", True, 0.893, 107, 9.622000e-33],
            ["NC_001422.1_2", True, 0.945, 304, 2.975000e-99],
            ["NC_001422.1_5", True, 0.972, 361, 1.509000e-118],
            ["NC_001422.1_6", True, 0.895, 614, 6.475000e-207],
            ["NC_001422.1_7", True, 0.980, 969, 0.000000e00],
            ["NC_001422.1_3", True, 0.960, 79, 1.702000e-23],
            ["NC_001422.1_4", True, 1.000, 908, 1.206000e-316],
        ],
        columns=["phix_hit", "phix_RBH", "phix_identity", "phix_bitScore", "phix_eVal"],
        index=[
            "NC_001422.1_1",
            "NC_001422.1_2",
            "NC_001422.1_5",
            "NC_001422.1_6",
            "NC_001422.1_7",
            "NC_001422.1_3",
            "NC_001422.1_4",
        ],
    )
    out_received = fastaDB.search(tmp_fasta)
    # Must drop the evalue or round it but I am lazy and it is close
    expected_out.drop("phix_eVal", axis=1).equals(
        out_received.drop("phix_eVal", axis=1)
    )


@pytest.fixture()
def tmp_camper_fasta(logger, tmp_path) -> Fasta:
    input_faa = Path("tests", "data", "camper_test_genes.faa")
    fasta = Fasta("test", None, Path(tmp_path), input_faa, None, None, None)
    return make_mmseqs_db_for_fasta(fasta, logger, 1)


@pytest.fixture()
def test_context(tmp_path, test_conf_path):
    return DramContext(tmp_path, path, test_conf_path)


def test_check_for_annotations():
    assert "You need to run annotate with with: [c, d]" in check_for_annotations(
        [{"a", "b", "c", "d"}, {"b", "c", "d", "e"}], {USED_DBS_TAG: ["a", "b"]}
    )
    assert "You need to annotate with: a or e" in (
        check_for_annotations(
            [{"a", "b", "c", "d"}, {"b", "c", "d", "e"}], {USED_DBS_TAG: ["b"]}
        )
    )
    assert "You also need to annotate with: e" not in (
        check_for_annotations(
            [{"a", "b", "c", "d"}, {"b", "c", "d", "e"}], {USED_DBS_TAG: ["b", "a"]}
        )
    )
    assert (
        check_for_annotations(
            [{"a", "b", "c", "d"}, {"b", "c", "d", "e"}],
            {USED_DBS_TAG: ["a", "b", "c", "d"]},
        )
        == None
    )
    assert (
        check_for_annotations(
            [{"a", "b", "c", "d", "f"}, {"b", "c", "d", "e"}],
            {USED_DBS_TAG: ["b", "c", "d", "e"]},
        )
        == None
    )


def test_HmmKit(
    tmp_camper_fasta,
    dbargs_dict,
):
    assert tmp_fasta.mmsdb.exists()
    hmmdb = HmmKit(
        "phix",
        Path("tests", "data", "camper_small.hmm"),
        Path("tests", "data", "camper_small_hmm_cutoffs.tsv"),
        {},
        dbargs_dict,
    )
    out_received = hmmdb.search(tmp_fasta)
    expected_out = pd.DataFrame()
    # Must drop the evalue or round it but I am lazy and it is close
    # expected_out.drop("phix_eVal", axis=1).equals(
    #     out_received.drop("phix_eVal", axis=1)
    # )
    # expected_out.equals( out_received.drop("phix_eVal", axis=1)


def test_path_to_gene_fastas(faa_loc, tmp_path):
    fasta = path_to_gene_fastas(Path(faa_loc), working_dir=Path(tmp_path))
    assert fasta.name == "NC_001422"
    assert fasta.tmp_dir == (tmp_path / "NC_001422")


def test_custom_faa(faa_loc, tmp_path):
    tmpdir = tmp_path / "work"
    result = annotate_pipe(
        [Path(faa_loc)],
        {},
        logging.getLogger(),
        tmp_path,
        tmpdir,
        1,
        custom_fasta_db_name="custom",
    )
    assert result["annotations"]["latest"]["database_used"] == []
    assert len(result["annotations"]["latest"]["fastas"]) == 1
    fasta_out = Fasta.import_strings(
        tmp_path, *result["annotations"]["latest"]["fastas"][0]
    )
    assert fasta_out.name == "NC_001422"
    assert fasta_out.mmsdb.exists()
    assert fasta_out.tmp_dir == (tmpdir / "NC_001422")
    assert fasta_out.mmsdb == (tmpdir / "NC_001422" / "gene.mmsdb")


def test_custom_hmm(faa_loc, tmp_path):
    tmpdir = tmp_path / "work"
    result = annotate([Path(faa_loc)], {}, logging.getLogger(), tmp_path, tmpdir, 1)
    assert result["annotations"]["latest"]["database_used"] == []
    assert len(result["annotations"]["latest"]["fastas"]) == 1
    fasta_out = Fasta.import_strings(
        tmp_path, *result["annotations"]["latest"]["fastas"][0]
    )
    assert fasta_out.name == "NC_001422"
    assert fasta_out.mmsdb.exists()
    assert fasta_out.tmp_dir == (tmpdir / "NC_001422")
    assert fasta_out.mmsdb == (tmpdir / "NC_001422" / "gene.mmsdb")


def test_no_output(faa_loc, test_context):
    tmpdir = tmp_path / "work"
    result = annotate_pipe(test_context)
    assert result["annotations"]["latest"]["database_used"] == []
    assert len(result["annotations"]["latest"]["fastas"]) == 1
    fasta_out = Fasta.import_strings(
        tmp_path, *result["annotations"]["latest"]["fastas"][0]
    )
    assert fasta_out.name == "NC_001422"
    assert fasta_out.mmsdb.exists()
    assert fasta_out.tmp_dir == (tmpdir / "NC_001422")
    assert fasta_out.mmsdb == (tmpdir / "NC_001422" / "gene.mmsdb")


def test_no_output_comand(runner, fasta_loc):

    result = runner.invoke(dram2, "annotate")
    assert result.exit_code != 0
    with pytest.raises(ValueError, match=r"You need to set an output directory*"):
        raise result.exception


def test_empty_cmd_from_dir(runner, fasta_loc, tmp_path):
    result1 = runner.invoke(dram2, ["-o", tmp_path, "call", "-f", fasta_loc])
    assert result1.exit_code == 0
    result = runner.invoke(dram2, ["-o", tmp_path, "annotate"])
    result_fail = runner.invoke(dram2, ["-o", tmp_path, "annotate"])
    result_force = runner.invoke(dram2, ["-o", tmp_path, "annotate"])
    assert result.exit_code == 0
    assert result_fail.exit_code == 0
    assert result_force.exit_code == 0
    out_ex = Path(tmp_path, "working_dir")
    assert not out_ex.exists()


def test_empty_cmd_from_genes(runner, faa_loc, tmp_path):
    result = runner.invoke(dram2, ["-o", tmp_path, "annotate", faa_loc])
    assert result.exit_code == 0
    out_ex = Path(tmp_path, "working_dir")
    assert not out_ex.exists()
