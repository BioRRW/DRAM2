"""
Annotate Called Genes with DRAM
==================

Main control point for the annotation process, and the main way to access the annotations.tsv or what ever it becomes. When we say that is the main control point of the annotation process we mean that it calls all the available objects in the dram2.db_kits name space. And orchestrates them to make the annotations file currently a tsv. as the origin of the annotations tsv it is also in the best position to parse the annotations.tsv in order to get gene_ids used in other processes, and to check if the annotations are available 

This is the second or first step of any pipeline and I would argue it is the heart of DRAM.  

TODO
 - make the annotated FASTA its own function
 -  Fix the header verbosity to be a separate option in the annotated fasta
 - Distillate sheets is part of the config, drop it. The sheets should be updated with the dram version so they do not get out of sink with the code.
 - add ability to take into account multiple best hits as in old_code.py
 - add silent mode
 - add abx resistance genes
 - in annotated gene faa checkout out ko id for actual kegg gene id
 - Add ability to build dbs on first run
 - Set the working dir separate from output
 - replace the tsv with something faster, maybe a parquet


"""
from multiprocessing import Pool
from functools import partial

# from itertools import chain
from dram2.utils import __version__
from dram2.cli.context import DramContext, DEFAULT_KEEP_TMP, get_time_stamp_id

import logging
from pathlib import Path

# from shutil import rmtree
import pandas as pd
from pkg_resources import resource_filename
import importlib
import pkgutil
import dram2.db_kits as db_kits
from dram2.db_kits.utils import DBKit, FastaKit, HmmKit, make_mmseqs_db
from dram2.utils.utils import Fasta, DramUsageError
from dram2.utils.globals import FASTAS_CONF_TAG
from itertools import chain
from collections import Counter

# from dram2.call_genes import clean_called_genes
from typing import Sequence, Optional
from datetime import datetime

# from dram2.utils.globals import


from pkg_resources import resource_filename


ANNOTATION_RUN_TAG = "annotation"
ANNOTATION_FILE_TAG: str = "annotation_file"
DEFAULT_BIT_SCORE_THRESHOLD: float = 60
DEFAULT_RBH_BIT_SCORE_THRESHOLD: float = 350
DEFAULT_THREADS: int = 10
DEFAULT_GENES_CALLED: bool = False
DEFAULT_KEEP_TMP: bool = False
FORCE_DEFAULT: bool = False
USED_DBS_TAG: str = "used_dbs"
ANNOTATIONS_TAG = "annotations"

DISTILLATION_MIN_SET = {"kegg|kofam", "dbcan", "pfam", "motif_count", "peptidase"}
DISTILLATION_MIN_SET = {"kegg|kofam", "dbcan", "pfam", "motif_count", "peptidase"}
DBSETS = {
    "mini": {},
    "mini_kegg": {},
    "adjectives_full": {},
    "adjectives_full_kegg": {},
}
import click

# def check_columns(data, db_kits, logger):
# functions = {i:j for i,j in ID_FUNCTION_DICT.items() if i in data.columns}
# missing = [i for i in ID_FUNCTION_DICT if i not in data.columns]
# logger.info("Note: the following id fields "
#             f"were not in the annotations file and are not being used: {missing},"
#             f" but these are {list(functions.keys())}")


def get_annotation_ids_by_row(data:pd.DataFrame, db_kits:list, groupby_column:Optional[str] = None):
    if groupby_column is not None:
        data.set_index(groupby_column, inplace=True)
    out = (data.apply(lambda x: {i for i in j.get_ids(x) for j in db_kits}))
    return out


def get_all_annotation_ids(ids_by_row) -> dict:
    out = Counter(chain(*ids_by_row.values))
    return out


def get_config_loc():
    return Path(resource_filename("dram2.utils", "CONFIG"))


for i in pkgutil.iter_modules(db_kits.__path__, db_kits.__name__ + "."):
    importlib.import_module(i.name)

DB_KITS: list = [i for i in DBKit.__subclasses__() if i.selectable]

__version__ = "2.0.0"


def check_for_annotations(annotation_set: set, annotation_run: dict):
    dbs = set(annotation_run[USED_DBS_TAG])

    def give_if_mising(annot: str):
        reqs = annot.split("|")
        for i in reqs:
            if i in dbs:
                return None
        if len(reqs) > 1:
            return f"({' or '.join(reqs)})"
        return reqs

    you_need = {i for i in annotation_set if give_if_mising(i) is not None}
    if len(you_need) > 0:
        raise DramUsageError(
            f"You are trying to use a DRAM2 function that requires specific annotations which this DRAM project dose not have yet. You need to run annotation with: {','.join(you_need)}\n\n"
            f"The command to do that is like `dram2 -o this_dram_dir annotate --use_db {'--use_db'.join(you_need)}`\n"
            "You should still review the docs to make sure you are running the program correctly to get results you want."
        )


def check_fasta_nams(fastas: list[Fasta]):
    fasta_names = [i.name for i in fastas]
    if len(fasta_names) != len(set(fasta_names)):
        raise ValueError(
            "Genome file names must be unique. At least one name appears twice "
            "in this search."
        )


def path_to_gene_fastas(fasta_loc: Path, working_dir: Path) -> Fasta:
    # make temporary directory
    # make tmp_dirs
    fasta_name = fasta_loc.stem
    fasta_working_dir = working_dir / fasta_name
    fasta_working_dir.mkdir(parents=True, exist_ok=True)
    return Fasta(fasta_name, fasta_loc, fasta_working_dir, fasta_loc, None, None, None)


def make_mmseqs_db_for_fasta(fasta: Fasta, logger, threads) -> Fasta:
    if fasta.mmsdb is not None and fasta.mmsdb.exists():
        return fasta
    if fasta.tmp_dir is None:
        raise ValueError(
            "Some how a fasta was passed to the function that makes mmseqs "
            "databases which did not have an associated temporary directory in "
            "which to put that mmseqs-db. Pleas kindly file a bug report on "
            "GitHub. This indicates that the developer probably made a mistake"
        )
    if fasta.faa is None:
        raise ValueError(
            "Some how a fasta was passed to the function that makes mmseqs "
            "databases which did not have an associated faa directory in which "
            "to put that mmseqs-db. Pleas kindly file a bug report on GitHub. "
            "This indicates that the developer probably made a mistake"
        )
    mmsdb = fasta.tmp_dir / "gene.mmsdb"
    make_mmseqs_db(
        fasta.faa.absolute().as_posix(),
        mmsdb.as_posix(),
        logger,
        create_index=True,
        threads=threads,
    )
    fasta.mmsdb = mmsdb
    return fasta


def log_file_path():
    pass


def get_all_fastas(
    gene_fasta_paths: list[Path],
    # genes_runs: Optional[dict], # may need this some day
    called_fastas: Optional[dict],
    annotation_run: Optional[dict],
    working_dir: Path,
    output_dir: Path,
    use_db: Sequence[str],
    logger: logging.Logger,
    force: bool,
):
    fastas: list[Fasta] = []
    fastas += [path_to_gene_fastas(i, working_dir) for i in gene_fasta_paths]
    if called_fastas is not None:
        fastas += [Fasta.import_strings(output_dir, *j) for j in called_fastas]
    if annotation_run is not None:
        db_inter = set(annotation_run["database_used"]).intersection(set(use_db))
        fasta_inter = set(annotation_run[FASTAS_CONF_TAG]).intersection(
            {i.name for i in fastas}
        )
        if len(db_inter) > 0:
            if force:
                logger.warning(
                    f"You are re-annotating {len(fasta_inter)} of {len(fastas)} FASTAs with databases they were already annotated with: {db_inter}. The past annotations with these databases be replaced."
                )
            else:
                raise ValueError(
                    f"You are trying re-annotating genes with database they were already annotated with: {db_inter}. You need to use the force flag '-f' in order to do this."
                )
    return fastas


def get_past_annotation_run(project_config: dict) -> Optional[dict]:
    annotation_run: Optional[dict] = (
        None
        if project_config.get("annotations") is None
        else project_config["annotations"][project_config["annotations"]["latest"]]
    )
    return annotation_run


def annotate(
    gene_fasta_paths: list[Path],
    dram_config: dict,
    logger: logging.Logger,
    output_dir: Path,
    working_dir: Path,
    cores: int,
    run_id: str,
    keep_tmp: bool = DEFAULT_KEEP_TMP,
    project_config: Optional[dict] = None,
    bit_score_threshold: int = DEFAULT_BIT_SCORE_THRESHOLD,
    rbh_bit_score_threshold: int = DEFAULT_RBH_BIT_SCORE_THRESHOLD,
    past_annotations_path: str = str(None),
    use_db: Sequence[str] = (),
    db_path: Optional[Path] = None,
    custom_fasta_db_name: Sequence = (),
    custom_fasta_db_loc: Sequence = (),
    custom_hmm_db_loc: Sequence = (),
    custom_hmm_db_name: Sequence = (),
    custom_hmm_db_cutoffs_loc: Sequence = (),
    kofam_use_dbcan2_thresholds: bool = False,
    # rename_genes: bool = True,
    threads: int = DEFAULT_THREADS,
    make_new_faa: bool = bool(None),
    force: bool = False,
    extra=None,
    write_config: bool = False,
) -> dict:
    # get assembly locations
    # make mmseqs_dbs
    working_dir.mkdir(exist_ok=True)
    if project_config is None:
        project_config = {}
    # make a seperate testable function for these two
    genes_runs: Optional[dict] = project_config.get("genes_called")
    called_fastas: Optional[dict] = project_config.get(FASTAS_CONF_TAG)
    annotation_run = get_past_annotation_run(project_config)
    fastas = get_all_fastas(
        gene_fasta_paths,
        # genes_runs,
        called_fastas,
        annotation_run,
        working_dir,
        output_dir,
        use_db,
        logger,
        force,
    )
    if len(fastas) < 1:
        raise ValueError(
            "No FASTAs were passed to the annotator DRAM has nothing to do."
        )
    if cores < len(fastas):
        cores_for_sub_process = cores
        cores_for_maping = 1
    else:
        cores_for_sub_process = 1
        cores_for_maping = cores
    with Pool(cores_for_maping) as p:
        fastas = p.map(
            partial(
                make_mmseqs_db_for_fasta, logger=logger, threads=cores_for_sub_process
            ),
            fastas,
        )
    db_args = {
        "logger": logger,
        # "fasta_paths": gene_fasta_paths,
        "output_dir": output_dir,
        "working_dir": working_dir,
        "bit_score_threshold": bit_score_threshold,
        "rbh_bit_score_threshold": rbh_bit_score_threshold,
        "kofam_use_dbcan2_thresholds": kofam_use_dbcan2_thresholds,
        "threads": threads,
        "force": force,
        "extra": extra,
        "db_path": db_path,  # where to store dbs on the fly
        "keep_tmp": keep_tmp,
    }

    # initsalize all used databases
    databases = [i(dram_config) for i in DB_KITS if i.name in set(use_db)]
    # add argument for annotations
    for i in databases:
        i.set_args(**db_args)

    # update the config
    if write_config:
        dram_config.update({j: k for i in databases for j, k in i.config.items()})

    # Add all the databases that you are going to
    databases += [
        FastaKit(i, j, dram_config, db_args)
        for i, j in zip(custom_fasta_db_name, custom_fasta_db_loc)
    ]
    db_len_dif = len(custom_hmm_db_name) - len(custom_hmm_db_cutoffs_loc)
    if db_len_dif < 0:
        raise ValueError(
            f"There are more hmm cutoff files provided then custom hmm databases provided"
        )

    custom_hmm_db_cutoffs_loc = list(custom_hmm_db_cutoffs_loc) + ([None] * db_len_dif)
    databases += [
        HmmKit(i, j, k, dram_config, db_args)
        for i, j, k in zip(
            custom_hmm_db_name, custom_hmm_db_loc, custom_hmm_db_cutoffs_loc
        )
    ]
    if len(databases) < 1:
        logger.warning(
            "No databases were selected. There is nothing for DRAM to do but save progress and exit."
        )
        new_annotations = pd.DataFrame(index=[fa.name for fa in fastas])
    else:
        # combine those annotations
        number_of_fastas = len(fastas)
        _ = [j.start_counter(number_of_fastas) for j in databases]
        new_annotations = pd.concat([j.search(i) for i in fastas for j in databases])

    if annotation_run is not None:
        past_annotations = pd.read_csv(
            output_dir / annotation_run["annotation_file"],
            sep="\t",
            index_col=0,
        )
        annotations = merge_past_annotations(
            new_annotations, past_annotations, logger, force
        )
    else:
        annotations = new_annotations

    annotation_tsv = output_dir / "annotations.tsv"
    annotations.to_csv(annotation_tsv, sep="\t")

    new_project_config = {
        ANNOTATIONS_TAG: {
            "latest": run_id,
            run_id: {
                "version": __version__,
                USED_DBS_TAG: [db.name for db in databases],
                "working_dir": working_dir.relative_to(output_dir).as_posix(),
                "database_used": [i.name for i in databases],
                ANNOTATION_FILE_TAG: annotation_tsv.relative_to(output_dir).as_posix(),
                FASTAS_CONF_TAG: [i.name for i in fastas],
            },
        },
        FASTAS_CONF_TAG: [i.export(output_dir) for i in fastas],
    }
    if genes_runs is not None:
        for i in genes_runs.values():
            if "fastas" in i:
                del i["fastas"]
            i["annotated"] = True
    return new_project_config


def merge_past_annotations(
    new_annotations: pd.DataFrame,
    past_annotations: pd.DataFrame,
    logger: logging.Logger,
    force: bool,
) -> pd.DataFrame:
    colliding_columns = set(past_annotations.columns).intersection(
        set(new_annotations.columns)
    )
    if len(colliding_columns) > 0:
        if not force:
            raise ValueError(
                "There is a column name collisions in the old and new annotations file, you need to use the force flag to overwrite this data"
            )
    colliding_genes = set(new_annotations.index).intersection(
        set(past_annotations.index)
    )
    past_annotations_merge = past_annotations.loc[list(colliding_genes)].drop(
        list(colliding_columns), axis=1
    )
    past_annotations_appened = past_annotations.drop(list(colliding_genes))

    all_annotations = pd.merge(
        new_annotations,
        past_annotations_merge,
        how="outer",
        left_index=True,
        right_index=True,
    )
    if len(all_annotations) != max(len(past_annotations_merge), len(new_annotations)):
        logger.critical(
            "The old and new annotations files may not have merge correctly! Check the new"
            " annotations file for errors. Did you use the correct genes.faa for your annotations?"
        )
    all_annotations = pd.concat([all_annotations, past_annotations_appened])
    return all_annotations


@click.command(
    "annotate",
    help="Get gene identifiers from a set of databases and format them for "
    "other DRAM2 analysis tools. This would be the second step of any DRAM2 pipeline "
    "unless your genes are already called, in witch case it is the first step.",
)
@click.argument(
    "gene_fasta_paths",
    type=click.Path(exists=True, path_type=Path),
    nargs=-1,
)
@click.option(
    "-s",
    "--study_set",
    multiple=True,
    type=click.Choice(list(DBSETS.keys()), case_sensitive=False),
)
@click.option(
    "--use_db",
    multiple=True,
    default=[],
    type=click.Choice([i.name for i in DB_KITS], case_sensitive=False),
    help="Specifiy exactly which dbs to use. This argument can be used multiple times, so for example if you want to annotate with FeGenie and Camper you would have a command like `dram2 -o dram_dir annotate --use_db fegenie --use_db camper`, the options avalible are in this help.",
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
    "--tempory_dir",
    type=click.Path(path_type=Path),
    help="Location of the tempory file where the annotations will be stored, this file will still be deleaded at the end of the annotation process if the keep_tmp flag is not set.",
)
@click.option(
    "--make_new_faa",
    default=DEFAULT_KEEP_TMP,
    help="If true the output directory will have a new genes.faa file with the"
    " anotation information apended to the headers. If false this file will not"
    " be made and time will be saved. If not specified the value will be set"
    " based on other arguments.",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Remove all past annotations and annotate again.",
)
@click.pass_context
def annotate_cmd(
    ctx: click.Context,
    gene_fasta_paths: list[Path],
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
    threads: int = DEFAULT_THREADS,
    # make_new_faa: Optional[bool] = None,
    # tempory_dir: Optional[Path] = None,
    force: bool = False,
    # db_path: Path = None,
    extra=None,
    # study_set: Sequence = (),
):
    context: DramContext = ctx.obj
    run_id: str = get_time_stamp_id(ANNOTATION_RUN_TAG)
    working_dir: Path = context.get_output_dir() / run_id
    logger = logging.getLogger("dram2_log")
    context.setup_logger(logger)
    output_dir: Path = context.get_output_dir()
    keep_tmp: bool = context.keep_tmp
    cores: int = context.cores
    project_config: dict = context.get_project_config()
    dram_config = context.get_dram_config(logger)  # FIX
    try:
        new_config = annotate(
            gene_fasta_paths=gene_fasta_paths,
            dram_config=dram_config,
            project_config=project_config,
            logger=logger,
            output_dir=output_dir,
            working_dir=working_dir,
            keep_tmp=keep_tmp,
            cores=cores,
            run_id=run_id,
            bit_score_threshold=bit_score_threshold,
            rbh_bit_score_threshold=rbh_bit_score_threshold,
            past_annotations_path=past_annotations_path,
            use_db=use_db,
            custom_fasta_db_name=custom_fasta_db_name,
            custom_fasta_db_loc=custom_fasta_db_loc,
            custom_hmm_db_loc=custom_hmm_db_loc,
            custom_hmm_db_name=custom_hmm_db_name,
            custom_hmm_db_cutoffs_loc=custom_hmm_db_cutoffs_loc,
            kofam_use_dbcan2_thresholds=kofam_use_dbcan2_thresholds,
            threads=threads,
            force=force,
            # db_path=db_path,
            extra=extra,
        )
        project_config.update(new_config)
        context.set_project_config(project_config)

    except Exception as e:
        logger.error(e)
        logger.exception("Fatal error in annotation")
        raise (e)


@click.command("list_databases")
def list_databases():
    print([i.name for i in DB_KITS])
