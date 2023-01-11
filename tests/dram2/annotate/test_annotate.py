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
from dram2.utils.command_line import dram2
from dram2.annotate import annotate


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
def mmseqs_db(prodigal_faa, mmseqs_db_dir, logger):
    output_file = str(mmseqs_db_dir.join("mmseqs_db.mmsdb"))
    make_mmseqs_db(prodigal_faa, output_file, logger, True, 1)
    return output_file


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
def target_mmseqs_db(mmseqs_db_dir, phix_proteins, logger):
    output_file = str(mmseqs_db_dir.join("target.mmsdb"))
    make_mmseqs_db(phix_proteins, output_file, logger, True, 1)
    return output_file




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


# Do this better
# def test_get_gene_data(phix_prodigal_genes):
#     scaffold_gene_df = get_gene_data(phix_prodigal_genes)
#     assert scaffold_gene_df.shape == (7, 5)


# def test_get_unannotated(phix_proteins):
#     annotated_genes = [
#         "NP_040704.1",
#         "NP_040703.1",
#         "NP_040713.1",
#         "NP_040712.1",
#         "NP_040711.1",
#         "NP_040710.1",
#         "NP_040709.1",
#         "NP_040707.1",
#     ]
#     unannotated_genes = ["NP_040705.1", "NP_040706.1", "NP_040708.1"]
#     test_unannotated_genes = get_unannotated(phix_proteins, annotated_genes)
#     assert set(unannotated_genes) == set(test_unannotated_genes)


# def test_assign_grades():
#     annotations_data = [
#         [True, "K00001", False, "KOER09234OK", ["PF00001"]],
#         [False, "K00002", True, "KLODKJFSO234KL", ["PF01234"]],
#         [False, "K00003", False, "EIORWU234KLKDS", np.NaN],
#         [False, np.NaN, False, np.NaN, np.NaN],
#         [False, np.NaN, False, np.NaN, ["PF01235"]],
#     ]
#     annotations = pd.DataFrame(
#         annotations_data,
#         index=["gene1", "gene2", "gene3", "gene4", "gene5"],
#         columns=["kegg_RBH", "kegg_hit", "uniref_RBH", "uniref_hit", "pfam_hits"],
#     )
#     test_grades = assign_grades(annotations)
#     assert test_grades.loc["gene1", "rank"] == "A"
#     assert test_grades.loc["gene2", "rank"] == "B"
#     assert test_grades.loc["gene3", "rank"] == "C"
#     assert test_grades.loc["gene4", "rank"] == "E"
#     assert test_grades.loc["gene5", "rank"] == "D"
#     # test no uniref
#     annotations_data2 = [
#         [True, "K00001", ["PF00001"]],
#         [False, "K00002", ["PF01234"]],
#         [False, "K00003", np.NaN],
#         [False, np.NaN, np.NaN],
#         [False, np.NaN, ["PF01235"]],
#     ]
#     annotations2 = pd.DataFrame(
#         annotations_data2,
#         index=["gene1", "gene2", "gene3", "gene4", "gene5"],
#         columns=["kegg_RBH", "kegg_hit", "pfam_hits"],
#     )
#     test_grades2 = assign_grades(annotations2)
#     assert test_grades2.loc["gene1", "rank"] == "A"
#     assert test_grades2.loc["gene2", "rank"] == "C"
#     assert test_grades2.loc["gene3", "rank"] == "C"
#     assert test_grades2.loc["gene4", "rank"] == "E"
#     assert test_grades2.loc["gene5", "rank"] == "D"



def test_no_output(runner, fasta_loc):
    result = runner.invoke(dram2, 'annotate')
    assert result.exit_code != 0
    assert "Error: Missing option \'-o\'" in result.output


# def test_just_prodigal(ctx, runner, faa_loc, tmp_path, faa_loc, logger):
#     annotate(gene_fasta_paths=[Path(faa_loc)], logger=logger, output_dir=tmp_path, force=True, keep_tmp= True)
#     out_ex = Path(tmp_path, "working_dir", "NC_001422", "genes.faa")
#     assert out_ex.exists()
#     cmp(out_ex, faa_loc)

def test_empty_cmd_from_dir(runner, fasta_loc, tmp_path):
    result1 = runner.invoke(dram2, ['-o', tmp_path, 'call', '-f', fasta_loc])
    assert result1.exit_code == 0
    result = runner.invoke(dram2, ['-o', tmp_path, 'annotate'])
    assert result.exit_code == 0
    out_ex = Path(tmp_path, "working_dir")
    assert not out_ex.exists()

def empty_cmd_from_genes(runner, faa_loc, tmp_path):
    result = runner.invoke(dram2, ['-o', tmp_path, 'annotate', faa_loc])
    assert result.exit_code == 0
    out_ex = Path(tmp_path, "working_dir")
    assert not out_ex.exists()


# def test_generate_annotated_fasta_short(phix_prodigal_genes, phix_annotations):
#     fasta_generator_short = generate_annotated_fasta(
#         phix_prodigal_genes, phix_annotations, name="phiX"
#     )
#     short_fasta_header_dict = {
#         seq.metadata["id"]: seq.metadata["description"] for seq in fasta_generator_short
#     }
#
#     assert short_fasta_header_dict["phiX_NC_001422.1_1"] == "rank: A; K1 (db=kegg)"
#     assert short_fasta_header_dict["phiX_NC_001422.1_2"] == "rank: B; U2 (db=uniref)"
#     assert short_fasta_header_dict["phiX_NC_001422.1_4"] == "rank: D; P4 (db=pfam)"


# def test_generate_annotated_fasta_short_no_kegg(
#     phix_prodigal_genes, phix_annotations_no_kegg
# ):
#     # drop kegg hit column
#     fasta_generator_short_no_kegg = generate_annotated_fasta(
#         phix_prodigal_genes, phix_annotations_no_kegg, name="phiX"
#     )
#     short_fasta_header_dict_no_kegg = {
#         seq.metadata["id"]: seq.metadata["description"]
#         for seq in fasta_generator_short_no_kegg
#     }
#     print(short_fasta_header_dict_no_kegg)
#
#     assert (
#         short_fasta_header_dict_no_kegg["phiX_NC_001422.1_1"]
#         == "rank: B; U1 (db=uniref)"
#     )
#     assert (
#         short_fasta_header_dict_no_kegg["phiX_NC_001422.1_2"]
#         == "rank: B; U2 (db=uniref)"
#     )
#     assert (
#         short_fasta_header_dict_no_kegg["phiX_NC_001422.1_4"] == "rank: D; P4 (db=pfam)"
#     )


# def test_generate_annotated_fasta_long(phix_prodigal_genes, phix_annotations):
#     fasta_generator_long = generate_annotated_fasta(
#         phix_prodigal_genes, phix_annotations, all_annotations=True, name="phiX"
#     )
#     long_fasta_header_dict = {
#         seq.metadata["id"]: seq.metadata["description"] for seq in fasta_generator_long
#     }
#     assert (
#         long_fasta_header_dict["phiX_NC_001422.1_1"]
#         == "rank: A; K1 (db=kegg); U1 (db=uniref); a_bug1"
#     )


# def test_create_annotated_fasta(phix_prodigal_genes, phix_annotations, tmpdir):
#     filt_fasta = tmpdir.mkdir("test_annotate_fasta").join("annotated_fasta.faa")
#     create_annotated_fasta(phix_prodigal_genes, phix_annotations, str(filt_fasta))
#     assert os.path.isfile(filt_fasta)


# def test_generate_renamed_fasta(fasta_loc):
#     genome_fasta_header_dict = {
#         seq.metadata["id"]: seq.metadata["description"]
#         for seq in generate_renamed_fasta(fasta_loc, "phiX")
#     }
#     assert len(genome_fasta_header_dict) == 1
#     assert (
#         genome_fasta_header_dict["phiX_NC_001422.1"]
#         == "Coliphage phi-X174, complete genome"
#     )


# def test_rename_fasta(fasta_loc, tmpdir):
#     filt_fasta = tmpdir.mkdir("test_annotate_scaffolds").join("annotated_fasta.fasta")
#     rename_fasta(fasta_loc, str(filt_fasta), "phiX")
#     assert os.path.isfile(filt_fasta)



# def test_count_motifs(phix_proteins):
#     motif_counts = count_motifs(phix_proteins, motif="(A.A)")
#     assert len(motif_counts) == 11
#     assert motif_counts["NP_040713.1"] == 7
#     assert motif_counts["NP_040709.1"] == 0




# def test_get_dubs():
#     w_dups = [True, False, True, True]
#     assert get_dups(w_dups) == [True, True, False, False]
#     no_dups = ["a", "b", "c"]
#     assert get_dups(no_dups) == [True, True, True]


