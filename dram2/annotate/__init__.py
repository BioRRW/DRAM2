"""
Annotate Called Genes with DRAM
==================

Main control point for the annotation process

TODO
 - make the annotated fasta its own function
 -  Fix the header verbosity to be a seperate option in the annotated fasta
 _ PRODIGAL IS MONO THREADED fix that obviously
# TODO Distillate sheets is part of the config, drop it
# TODO: add ability to take into account multiple best hits as in old_code.py
# TODO: add silent mode
# TODO: add abx resistance genes
# TODO: in annotated gene faa checkout out ko id for actual kegg gene id
# TODO: add ability to handle [] in file names
# TODO: Add ability to build dbs on first run
# TODO: Set a default ouput dir on first run
# TODO: Set the working dir seperate from output
"""
import logging
from pathlib import Path
from os import path
from shutil import rmtree
import pandas as pd
from pkg_resources import resource_filename
import importlib
import pkgutil
import dram2.db_kits as db_kits
from dram2.db_kits.utils import DBKit, FastaKit, HmmKit, Fasta, make_mmseqs_db
from dram2.utils.utils import load_config
from dram2.call_genes import get_fasta_name
from typing import Sequence, Optional
# from dram2.utils.globals import




from pkg_resources import resource_filename


DEFAULT_BIT_SCORE_THRESHOLD: float = 60
DEFAULT_RBH_BIT_SCORE_THRESHOLD: float = 350
DEFAULT_THREADS: int = 10
DEFAULT_GENES_CALLED: bool = False
DEFAULT_KEEP_TMP: bool = False
DBSETS = []
import click

def get_config_loc():
    return Path(resource_filename("dram2.utils", "CONFIG"))

for i in pkgutil.iter_modules(db_kits.__path__, db_kits.__name__ + "."):
    importlib.import_module(i.name)

DB_KITS: list = [i for i in DBKit.__subclasses__() if i.selectable]

__version__ = "2.0.0"




def get_gene_fasta_names(fasta_locs:list[Path], working_dir: Path) -> Fasta:
    fasta_names = [get_fasta_name(i) for i in fasta_locs]
    # make temporary directory
    if len(fasta_names) != len(set(fasta_names)):
        raise ValueError(
            "Genome file names must be unique. At least one name appears twice in this search."
        )
    # make tmp_dirs
    tmp_dirs = [working_dir / i for i in fasta_names]
    for i in tmp_dirs:
        i.mkdir(exist_ok=True)

    return [
        Fasta(i, j, k, j, None, None, None)
        for i, j, k in zip(fasta_names, fasta_locs, tmp_dirs)
    ]


def make_mmseqs_db_for_fasta(fasta: Fasta, logger, threads):
    mmsdb = fasta.tmp_dir / "gene.mmsdb"
    make_mmseqs_db(
        fasta.faa, mmsdb.as_posix(), logger, create_index=True, threads=threads
    )
    fasta.mmsdb = mmsdb
    return fasta

def log_file_path

def annotate(
    ctx,
    gene_faa_paths: list[Path],
    logger: logging.Logger,
    output_dir: Path,
    bit_score_threshold: int = DEFAULT_BIT_SCORE_THRESHOLD,
    rbh_bit_score_threshold: int = DEFAULT_RBH_BIT_SCORE_THRESHOLD,
    past_annotations_path: str = str(None),
    use_db: Sequence = (),
    db_path: Path = None,
    custom_fasta_db_name: Sequence = (),
    custom_fasta_db_loc: Sequence = (),
    custom_hmm_db_loc: Sequence = (),
    custom_hmm_db_name: Sequence = (),
    custom_hmm_db_cutoffs_loc: Sequence = (),
    kofam_use_dbcan2_thresholds: bool = False,
    # rename_genes: bool = True,
    keep_tmp: bool = DEFAULT_KEEP_TMP,
    threads: int = DEFAULT_THREADS,
    make_new_faa: bool = bool(None),
    dry: bool = False,
    force: bool = False,
    extra=None,
    config: Optional(dict) = None,
    write_config: bool = False,
):
    cores:int = ctx.obj.cores
    timestamp_id:str = datetime.now().strftime('%Y%m%d%H%M%S')
    project_config:dict = ctx.obj.get_project_config()

    saved_faa_paths = load_faa_paths(project_config)
    gene_faa_paths

    # get assembly locations
    fastas = []
    fastas += get_gene_fasta_names(gene_faa_paths,  working_dir)
    fastas += load_saved_fastas()

    # make mmseqs_dbs
    fastas = [make_mmseqs_db_for_fasta(fasta, logger, threads) for fasta in fastas]
    fastas += ctx.obj["fastas"]

    db_args = {
        "kofam_use_dbcan2_thresholds": kofam_use_dbcan2_thresholds,
        # "fasta_paths": fasta_paths,
        "output_dir": output_dir,
        "working_dir": working_dir,
        "bit_score_threshold": bit_score_threshold,
        "rbh_bit_score_threshold": rbh_bit_score_threshold,
        "past_annotations_path": past_annotations_path,
        "threads": threads,
        "logger": logger,
        "make_new_faa": make_new_faa,
        "dry": dry,
        "db_path": db_path,
        "force": force,
        "extra": extra,
    }

    database = [i(config, db_args) for i in DB_KITS if i.name in use_db]

    # update the config
    if write_config:
        config.update({j: k for i in database for j, k in i.config.items()})

    # Add all the database that you are going to
    database += [
        FastaKit(i, j, config, db_args)
        for i, j in zip(custom_fasta_db_name, custom_fasta_db_loc)
    ]
    db_len_dif = len(custom_hmm_db_name) - len(custom_hmm_db_cutoffs_loc)
    if db_len_dif < 0:
        raise ValueError(
            f"There are more hmm cutoff files provided then custom hmm database provided"
        )

    custom_hmm_db_cutoffs_loc = list(custom_hmm_db_cutoffs_loc) + ([None] * db_len_dif)
    database += [
        HmmKit(i, j, k, config, db_args)
        for i, j, k in zip(
            custom_hmm_db_name, custom_hmm_db_loc, custom_hmm_db_cutoffs_loc
        )
    ]

    # combine those annotations
    annotations = pd.concat([j.search(i) for i in fastas for j in database])
    annotations.to_csv(output_dir / "annotations.tsv", sep="\t")
    # clean up
    if not keep_tmp:
        logger.info(f"Cleaning up temporary directory: {working_dir}")
        rmtree(working_dir)

    # # logger.info('Annotation started')
    # fasta_names = get_fasta_names()

    # if len(use_db + custom_fasta_loc + custom_hmm_loc) < 1:
    #     raise ValueError(
    #         "For some reason there are no database selected to"
    #         " annotate against. This is most likely a result of"
    #         " bad arguments."
    #     )

    # # get database locations
    # db_handler = DatabaseHandler(logger)

    # tmp_dir = path.join(output_dir, "working_dir")
    # mkdir(tmp_dir)

    # # annotate
    # annotation_locs = list()
    # faa_locs = list()
    # # TODO IO is so slow!!! We need to use more memory and make less files

    # # check if we are making a new faa
    # if make_new_faa is None:
    #     # We only make a new faa if the user asks or
    #     # we are not apending to a past annotatons
    #     make_new_faa = past_annotations_path is None

    # annotation_locs: list = [
    #     partial_annotate_one_fasta(fasta_loc) for fasta_loc in fasta_locs
    # ]
    # # merge
    # annotation_frames: list = [
    #     pd.read_csv(i, sep="\t", index_col=0) for i in annotation_locs
    # ]
    # all_annotations = pd.concat(
    #     annotation_frames,
    #     sort=False,
    # )
    # all_annotations = all_annotations.sort_values("fasta")
    # if past_annotations_path is not None:
    #     all_annotations.drop("fasta", axis=1, inplace=True)
    #     past_annotations = pd.read_csv(past_annotations_path, sep="\t", index_col=0)
    #     if set(past_annotations.columns).issubset(set(all_annotations.columns)):
    #         logger.warning(
    #             "You have passed an annotations file that containes"
    #             " columns matching the new annotations. You most"
    #             " likely are annotating with the same database again."
    #             " The falowing columns will be replaced in the new"
    #             " annotations file:\n%s"
    #             % set(past_annotations.columns).intersection(
    #                 set(all_annotations.columns)
    #             )
    #         )
    #     past_annotations = past_annotations[
    #         list(set(past_annotations.columns) - set(all_annotations.columns))
    #     ]
    #     all_annotations.index = all_annotations.index.str.removeprefix("genes_")
    #     # Note we use the old annotations fasts
    #     all_annotations = pd.merge(
    #         all_annotations,
    #         past_annotations,
    #         how="outer",
    #         left_index=True,
    #         right_index=True,
    #     )
    #     if len(all_annotations) != len(past_annotations):
    #         logger.critical(
    #             "The old and new annotations filles did not merge correctly! Check the new"
    #             " annotations file for errors. Did you use the corect genes.faa for your annotations?"
    #         )
    # all_annotations.to_csv(path.join(output_dir, "annotations.tsv"), sep="\t")
    # if len(faa_locs) > 0:
    #     merge_files(faa_locs, path.join(output_dir, "genes.faa"))




@click.command("annotate")
@click.argument(
    "fasta_paths",
    type=click.Path(exists=True, path_type=Path),
    nargs=-1,
)
@click.option(
    "-s",
    "--study_set",
    multiple=True,
    type=click.Choice(DBSETS, case_sensitive=False),
)
@click.option(
    "--use_db",
    multiple=True,
    default=[],
    type=click.Choice([i.name for i in DB_KITS], case_sensitive=False),
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
    "--custom_fasta_db_name",
    type=str,
    multiple=True,
    help="Names of custom databases, can be used multiple times.",
)
@click.option(
    "--custom_fasta_db_loc",
    multiple=True,
    type=click.Path(exists=True, path_type=Path),
    help="Location of fastas to annotate against, can be used multiple times but"
    "must match nubmer of custom_db_name's",
)
@click.option(
    "--custom_hmm_db_name",
    multiple=True,
    help="Names of custom hmm databases, can be used multiple times.",
)
@click.option(
    "--custom_hmm_db_loc",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="Location of hmms to annotate against, can be used multiple times but"
    "must match nubmer of custom_hmm_name's",
)
@click.option(
    "--custom_hmm_db_cutoffs_loc",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="Location of file with custom HMM cutoffs and descriptions, can be used "
    "multiple times.",
)
@click.option(
    "--make_new_faa",
    default=DEFAULT_KEEP_TMP,
    help="If true the output directory will have a new genes.faa file with the"
    " anotation information apended to the headers. If false this file will not"
    " be made and time will be saved. If not specified the value will be set"
    " based on other arguments.",
)
@click.pass_context
#= "~/.dram/config.yaml"
def annotate_wraper(
    ctx: list,
    gene_faa_path: list[Path],
    bit_score_threshold: int = DEFAULT_BIT_SCORE_THRESHOLD,
    rbh_bit_score_threshold: int = DEFAULT_RBH_BIT_SCORE_THRESHOLD,
    # log_file_path: str = str(None),
    past_annotations_path: str = str(None),
    use_db: Sequence = (),
    custom_fasta_db_name: Sequence = (),
    custom_fasta_db_loc: Sequence = (),
    custom_hmm_db_loc: Sequence = (),
    custom_hmm_db_name: Sequence = (),
    custom_hmm_db_cutoffs_loc: Sequence = (),
    kofam_use_dbcan2_thresholds: bool = False,
    rename_genes: bool = True,
    threads: int = DEFAULT_THREADS,
    make_new_faa: bool = bool(None),
    dry: bool = False,
    force: bool = False,
    # db_path: Path = None,
    extra=None,
    study_set: Sequence= (),
):
    try:
        logger:logging.Logger = ctx.obj.get_logger()
        output_dir:Path = ctx.obj.get_output_dir()
        keep_tmp:bool = ctx.obj.keep_tmp
        annotate(
            ctx,
            gene_faa_path,
            logger=logger,
            output_dir=output_dir,
            bit_score_threshold=bit_score_threshold,
            rbh_bit_score_threshold=rbh_bit_score_threshold,
            past_annotations_path=past_annotations_path,
            use_db=use_db,
            custom_fasta_db_name=custom_fasta_db_name,
            custom_fasta_db_loc=custom_fasta_db_loc,
            custom_hmm_db_loc=custom_hmm_db_loc,
            custom_hmm_db_name=custom_hmm_db_name,
            custom_hmm_db_cutoffs_loc=custom_hmm_db_cutoffs_loc,
            keep_tmp=keep_tmp,
            kofam_use_dbcan2_thresholds=kofam_use_dbcan2_thresholds,
            threads=threads,
            make_new_faa=make_new_faa,
            dry=dry,
            config=config,
            force=force,
            # db_path=db_path,
            extra=extra,
        )

    except Exception as e:
        logger.error(e)
        logger.exception("Fatal error in annotation")
        raise (e)


@click.command("list_databases")
def list_databases():
    print([i.name for i in DB_KITS])
