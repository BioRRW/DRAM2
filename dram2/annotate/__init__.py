__version__ = "2.0.0"
from dram2.utils.globals import (
    PRODIGAL_MODE_CHOICES,
    PRODIGAL_TRANS_TABLE_CHOICES,
    PRODIGAL_MODE_DFT,
    PRODIGAL_TRANS_TABLE_DFT
)
from dram2.annotate.annotate import DB_KITS
import click

class RequiredIfNot(click.Option):
    def __init__(self, *args, **kwargs):
        self.required_if_not = kwargs.pop('required_if_not')
        if not self.required_if_not: raise click.UsageError("'required_if_not' parameter required")
        kwargs['help'] = (kwargs.get('help', '') +
            ' NOTE: This argument or {self.required_if_not} is requered').strip()
        super(RequiredIfNot, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        we_are_present = self.name in opts
        other_present = self.required_if_not not in opts

        if not other_present:
            if not we_are_present:
                raise click.UsageError(
                    f"Illegal usage: {self.name} or {self.required_if_not} must be used")
            else:
                self.prompt = None

        return super(RequiredIfNot, self).handle_parse_result(
            ctx, opts, args)



@click.command("annotate_bins")
@click.argument(
    "fasta_paths",
    type=click.Path(exists=True),
    multiple=True,
    nargs=-1.
    cls=RequiredIfNot,
    required_if_not="gene_faa_paths",
    help="fasta file, optionally with wildcards to point to " "individual MAGs",
)
@click.option(
    "-g",
    "--gene_faa_paths",
    required_if_not="annotate_bins",
    help="fasta file, optionally with wildcards to point to " "individual MAGs",
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
    "--prodigal_mode",
    type=click.Choice([str(i) for i in range(1, 26)], case_sensitive=False),
    default=PRODIGAL_MODE_DFT,
    help="Mode of prodigal to use for gene calling. NOTE: normal or single mode require genomes which are high quality with low contamination and long contigs (average length >3 Kbp).",
)
#     prodigal_trans_table_choices =
#     annotate_parser.add_argument('--trans_table', type=str, default='11', choices=prodigal_trans_table_choices,
#                                  help='Translation table for prodigal to use for gene calling.')
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
# @click.option('--use_uniref', default=False,
#               help='Annotate these fastas against UniRef, drastically increases run time and '
#                    'memory requirements')
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


@dram2.command("list_databases")
def list_databases():
    print([i.NAME for i in DB_KITS])
