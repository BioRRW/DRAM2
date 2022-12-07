__version__ = "2.0.0"

from dram2.utils.globals import (
    PRODIGAL_MODE_CHOICES,
    PRODIGAL_TRANS_TABLE_CHOICES,
    PRODIGAL_MODE_DFT,
    PRODIGAL_TRANS_TABLE_DFT
)
from dram2.annotate.annotate import DB_KITS
import click


@click.command("annotate")
@click.argument(
    "fasta_paths",
    type=click.Path(exists=True),
    nargs=-1,
)
@click.option(
    "-g",
    "--genes_called",
    is_flag=True,
    help="Specify that you are using called genes and not MAGs. The gene calling step will be skipped. Note that this will try to add new genes to the databases not in fact call old genes. For that you should use the --annotations_path argument.",
)
@click.option("-o", "--output_dir", help="output directory", required=True)
@click.option(
    "-a",
    "--past_annotations_path",
    help="past_annotations to append new annotations to.",
)
@click.option(
    "-m",
    "--input_mag_path",
    multiple=True,
    help="Path to a fasta file with an individual MAG or otherwise grouped uncalled genes, or a python glob compatible wildcard to point a collection of such files",
)
@click.option(
    "-g",
    "--input_faa_path",
    multiple=True,
    help="Path to a fasta file with called genes, or a python glob compatible wildcard to point a collection of such files",
)
@click.option("-o", "--output_dir", help="output directory", required=True)
@click.option(
    "-a",
    "--past_annotations_path",
    help="past_annotations to append new annotations to.",
)
@click.option(
    "--bit_score_threshold",
    type=int,
    default=60,
    help="minimum bitScore of search to retain hits",
)
@click.option(
    "--rbh_bit_score_threshold",
    type=int,
    default=350,
    help="minimum bitScore of reverse best hits to retain hits",
)
@click.option(
    "--kofam_use_dbcan2_thresholds",
    default=False,
    help="Use dbcan2 suggested HMM cutoffs for KOfam annotation instead of KOfam "
    "recommended cutoffs. This will be ignored if annotating with KEGG Genes.",
)
@click.option(
    "--custom_db_name",
    multiple=True,
    help="Names of custom databases, can be used multiple times.",
)
@click.option(
    "--custom_fasta_loc",
    multiple=True,
    help="Location of fastas to annotate against, can be used multiple times but"
    "must match nubmer of custom_db_name's",
)
@click.option(
    "--custom_hmm_name",
    multiple=True,
    help="Names of custom hmm databases, can be used multiple times.",
)
@click.option(
    "--custom_hmm_loc",
    multiple=True,
    help="Location of hmms to annotate against, can be used multiple times but"
    "must match nubmer of custom_hmm_name's",
)
@click.option(
    "--prodigal_mode",
    type=click.Choice(["train", "meta", "single"], case_sensitive=False),
    help="Mode of prodigal to use for gene calling. NOTE: normal or single mode require genomes which are high quality with low contamination and long contigs (average length >3 Kbp).",
)
@click.option(
    "--prodigal_trans_tables",
    type=click.Choice([str(i) for i in range(1, 26)], case_sensitive=False),
    default=PRODIGAL_MODE_DFT,
    help="Mode of prodigal to use for gene calling. NOTE: normal or single mode require genomes which are high quality with low contamination and long contigs (average length >3 Kbp).",
)
@click.option(
    "--dry",
    is_flag=True,
    help="Specify that this is a dry run, lists input files but dose not run them",
)
@click.option(
    "--log_file_path", help="Optional path for the log file that will document this run"
)
@click.option(
    "--custom_hmm_cutoffs_loc",
    multiple=True,
    help="Location of file with custom HMM cutoffs and descriptions, can be used "
    "multiple times.",
)
@click.option(
    "--use_db",
    multiple=True,
    type=click.Choice([i.NAME for i in DB_KITS], case_sensitive=False),
)
@click.option("--keep_tmp_dir", default=False)
@click.option(
    "--make_new_faa",
    default=None,
    help="If true the output directory will have a new genes.faa file with the"
    " anotation information apended to the headers. If false this file will not"
    " be made and time will be saved. If not specified the value will be set"
    " based on other arguments.",
)
@click.option("--threads", type=int, default=10, help="number of processors to use")
def annotate(args):
    from dram2.annotate.annotate import annotate
    annotate(**args)


@click.command("list_databases")
def list_databases():
    print([i.NAME for i in DB_KITS])
