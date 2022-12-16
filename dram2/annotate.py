"""
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
# TODO: Set the woring dir seperate from output
"""
import re
import io
import time
import logging
from pathlib import Path
from functools import partial, reduce
from datetime import datetime
from skbio.io import read as read_sequence, write as write_sequence
from skbio import Sequence
from skbio.metadata import IntervalMetadata
from os import path, mkdir, stat
from shutil import rmtree, copy2
import pandas as pd
import numpy as np
from pkg_resources import resource_filename
import importlib
import pkgutil
import dram2.db_kits as db_kits
from dram2.db_kits.utils import DBKit, FastaKit, HmmKit, Fasta, make_mmseqs_db
from dram2.utils.utils import load_config


def get_config_loc():
    return Path(resource_filename("dram2.utils", "CONFIG"))


from pkg_resources import resource_filename


DEFAULT_MIN_CONTIG_SIZE: int = 2500
DEFAULT_PRODIGAL_MODE: str = "meta"
DEFAULT_TRANS_TABLE: str = "11"
DEFAULT_OUTPUT_DIR: Path = Path(".")
DEFAULT_BIT_SCORE_THRESHOLD: float = 60
DEFAULT_RBH_BIT_SCORE_THRESHOLD: float = 350
DEFAULT_THREADS: int = 10
DEFAULT_GENES_CALLED: bool = False
DEFAULT_KEEP_TMP: bool = False
DBSETS = []
import click

for i in pkgutil.iter_modules(db_kits.__path__, db_kits.__name__ + "."):
    importlib.import_module(i.name)

DB_KITS: list = [i for i in DBKit.__subclasses__() if i.selectable]

__version__ = "2.0.0"

from dram2.utils.utils import (
    run_process,
    remove_suffix,
    setup_logger,
)





def make_mmseqs_db_for_fasta(fasta: Fasta, logger, threads):
    mmsdb = fasta.tmp_dir / "gene.mmsdb"
    make_mmseqs_db(
        fasta.faa, mmsdb.as_posix(), logger, create_index=True, threads=threads
    )
    fasta.mmsdb = mmsdb
    return fasta


def annotate(
    fasta_paths: list[Path],
    logger: logging.Logger,
    output_dir: Path,
    genes_called: bool = DEFAULT_GENES_CALLED,
    bit_score_threshold: int = DEFAULT_BIT_SCORE_THRESHOLD,
    rbh_bit_score_threshold: int = DEFAULT_RBH_BIT_SCORE_THRESHOLD,
    past_annotations_path: str = str(None),
    use_db: list = (),
    db_path: Path = None,
    custom_fasta_db_name: list = (),
    custom_fasta_db_loc: list = (),
    custom_hmm_db_loc: list = (),
    custom_hmm_db_name: list = (),
    custom_hmm_db_cutoffs_loc: list = (),
    kofam_use_dbcan2_thresholds: bool = False,
    # rename_genes: bool = True,
    keep_tmp: bool = DEFAULT_KEEP_TMP,
    threads: int = DEFAULT_THREADS,
    make_new_faa: bool = bool(None),
    dry: bool = False,
    force: bool = False,
    extra=None,
    config: dict = None,
    write_config: bool = False,
):

    # if dry:
    #     logger.info(
    #         "If this was not a dry run, DRAM would have annotated {len(fastas)} assembly and {len(gene_faa_path)} gene files. all files are listed bellow"
    #     )
    #     logger(f"Assemblies: {fasta_paths[]}")
    #     return

    # Make a working_dir that may be deleted
    working_dir: Path = output_dir / "working_dir"
    working_dir.mkdir(exist_ok=True)
    # get assembly locations
    fastas = get_fasta_names(fasta_paths, genes_called, working_dir)

    # Filter and Call the genes
    if not genes_called:

    # make mmseqs_dbs
    fastas = [make_mmseqs_db_for_fasta(fasta, logger, threads) for fasta in fastas]

    # This is were you would combine the fasta data if you were so inclined
    # if len(use_db) < 1:
    #     logger.warning(
    #         f"No databases where selected to annotate with, genes will be call if not called already but nothing else will be done"
    #     )
    # else:
    #     logger.info(
    #         f"Starting the Annotation of Genes with database configuration: \n {use_db}"
    #     )

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




@click.command("annotate", chain=True)
@click.argument(
    "fasta_paths",
    type=click.Path(exists=True, path_type=Path),
    nargs=-1,
)
# @click.option(
#     "--prodigal_mode",
#     default=DEFAULT_PRODIGAL_MODE,
#     type=click.Choice(["train", "meta", "single"], case_sensitive=False),
#     help="Mode of prodigal to use for gene calling. NOTE: normal or single mode require genomes which are high quality with low contamination and long contigs (average length >3 Kbp).",
# )
# @click.option(
#     "--prodigal_trans_tables",
#     type=click.Choice([str(i) for i in range(1, 26)], case_sensitive=False),
#     default=DEFAULT_TRANS_TABLE,
#     help="Mode of prodigal to use for gene calling. NOTE: normal or single mode require genomes which are high quality with low contamination and long contigs (average length >3 Kbp).",
# )
@click.option(
    "-o", "--output_dir", help="output directory", required=True, type=click.Path()
)
# @click.option(
#     "-a",
#     "--past_annotations_path",
#     help="past_annotations to append new annotations to.",
# )
# @click.option(
#     "-g",
#     "--genes_called",
#     is_flag=True,
#     help="Specify that you are using called genes and not MAGs. The gene calling step will be skipped. Note that this will try to add new genes to the databases not in fact call old genes. For that you should use the --annotations_path argument.",
# )
# @click.option(
#     "-a",
#     "--past_annotations_path",
#     help="past_annotations to append new annotations to.",
# )
# @click.option(
#     "--kofam_use_dbcan2_thresholds",
#     default=False,
#     help="Use dbcan2 suggested HMM cutoffs for KOfam annotation instead of KOfam "
#     "recommended cutoffs. This will be ignored if annotating with KEGG Genes.",
# )
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

@click.option("-o", "--output_dir", help="output directory", required=True)
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
#= "~/.dram/config.yaml"
def annotate_wraper(
    fasta_paths: list[Path] ,
    genes_called: bool = DEFAULT_GENES_CALLED,
    min_contig_size=DEFAULT_MIN_CONTIG_SIZE,
    prodigal_mode=DEFAULT_PRODIGAL_MODE,
    prodigal_trans_tables: str = DEFAULT_TRANS_TABLE,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    bit_score_threshold: float = DEFAULT_BIT_SCORE_THRESHOLD,
    rbh_bit_score_threshold: float = DEFAULT_RBH_BIT_SCORE_THRESHOLD,
    log_file_path: str = str(None),
    past_annotations_path: str = str(None),
    use_db: list = (),
    custom_fasta_db_name: list = (),
    custom_fasta_db_loc: list = (),
    custom_hmm_db_loc: list = (),
    custom_hmm_db_name: list = (),
    custom_hmm_db_cutoffs_loc: list = (),
    kofam_use_dbcan2_thresholds: bool = False,
    rename_genes: bool = True,
    keep_tmp: bool = True,
    threads: int = DEFAULT_THREADS,
    make_new_faa: bool = bool(None),
    dry: bool = False,
    force: bool = False,
    config_path: Path = None,
    db_path: Path = None,
    extra=None,
    study_set: list = (),
):
    # Get a logger
    if log_file_path is None:
        log_file_path = path.join(output_dir, "Annotation.log")
    logger = logging.getLogger("annotation_log")

    try:
        if not output_dir.exists():
            output_dir.mkdir()
        elif not force:
            raise ValueError(
                "The output_dir already exists! try using the -f flag to overwrite"
            )
        setup_logger(logger, log_file_path)
        logger.info(f"The log file is created at {log_file_path}")
        config = load_config(config_path, logger)
        annotate(
            logger=logger,
            fasta_paths=fasta_paths,
            genes_called=genes_called,
            min_contig_size=min_contig_size,
            prodigal_mode=prodigal_mode,
            trans_table=prodigal_trans_tables,
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
            rename_genes=rename_genes,
            keep_tmp=keep_tmp,
            kofam_use_dbcan2_thresholds=kofam_use_dbcan2_thresholds,
            threads=threads,
            make_new_faa=make_new_faa,
            dry=dry,
            config=config,
            force=force,
            db_path=db_path,
            extra=extra,
        )

    except Exception as e:
        logger.error(e)
        logger.exception("Fatal error in annotation")
        raise (e)


@click.command("list_databases")
def list_databases():
    print([i.name for i in DB_KITS])
