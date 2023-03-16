"""
===============================
Annotate Called Genes with DRAM
===============================

Main control point for the annotation process, and the main way to access the
annotations.tsv or whatever it becomes. When we say that is the main control
point of the annotation process, we mean that it calls all the available objects
in the dram2.db_kits name space. And orchestrates them to make the annotations
file, currently a tsv. As the origin of the annotations TSV it is also in the
best position to parse the annotations.tsv in order to get gene IDs used in
other processes, and to check if the annotations are available.

This is the second or first step of any pipeline and I would argue it is the
heart of DRAM.

TODO
 - make the annotated FASTA its own function.
    -  Fix the header verbosity to be a separate option in the annotated FASTA.
 - Distillate sheets are part of the config, drop it. The sheets should be
   updated with the dram version, so they do not get out of sync with the code.
 - add ability to take into account multiple best hits as in old_code.py
 - add silent mode
 - add abx resistance genes
 - in annotated gene faa, checkout out ko_id for actual KEGG gene id
 - Add ability to download DBss on first run
 - Set the working-dir separate from output
 - replace the tsv with something faster, maybe a parquet
"""

import importlib
import pkgutil
import logging
from multiprocessing import Pool
from functools import partial, reduce
from pathlib import Path
from typing import Sequence, Optional
from shutil import rmtree
from itertools import chain
from collections import Counter
from pkg_resources import resource_filename

import click
import pandas as pd

from dram2 import db_kits as db_kits
from dram2.call_genes import DEFAULT_GENES_FILE
from dram2.db_kits.utils import DBKit, FastaKit, HmmKit, make_mmseqs_db
from dram2.utils import DramUsageError, Fasta
from dram2.utils.globals import FASTAS_CONF_TAG
from dram2.cli.context import (
    DramContext,
    DEFAULT_KEEP_TMP,
    get_time_stamp_id,
    __version__,
)


for i in pkgutil.iter_modules(db_kits.__path__, db_kits.__name__ + "."):
    importlib.import_module(i.name)

DB_KITS: list = [i for i in DBKit.__subclasses__() if i.selectable]


ANNOTATION_FILE_TAG: str = "annotation_file"
DEFAULT_BIT_SCORE_THRESHOLD: float = 60
DEFAULT_RBH_BIT_SCORE_THRESHOLD: float = 350
DEFAULT_THREADS: int = 10
DEFAULT_GENES_CALLED: bool = False
DEFAULT_KEEP_TMP: bool = False
FORCE_DEFAULT: bool = False
USED_DBS_TAG: str = "used_dbs"
ANNOTATIONS_TAG = "annotations"
FASTA_COL = "fasta"
GENE_ID_COL = "gene_ids"
DISTILLATION_MIN_SET_KEGG = {"kegg", "dbcan", "pfam", "heme", "peptidase"}
DISTILLATION_MIN_SET = {"kofam", "dbcan", "pfam", "heme", "peptidase"}
ADJECTIVES_SET = {
    "kofam",
    "dbcan",
    "pfam",
    "heme",
    "peptidase",
    "sulfur",
    "camper",
    "methyl",
    "fegenie",
}
ADJECTIVES_SET_KEGG = {
    "kegg",
    "dbcan",
    "pfam",
    "heme",
    "peptidase",
    "sulfur",
    "camper",
    "methyl",
    "fegenie",
}
DBSETS_COL = "db_id_sets"
DBSETS = {
    "mini": DISTILLATION_MIN_SET,
    "mini_kegg": DISTILLATION_MIN_SET,
    "adjectives": ADJECTIVES_SET,
    "adjectives_kegg": ADJECTIVES_SET_KEGG,
    "wrighton_full": {i.name for i in DB_KITS},
    # "adjectives_full": {},
    # "adjectives_full_kegg": {},
}


def get_annotation_ids_by_row(data: pd.DataFrame, db_kits: list) -> pd.DataFrame:
    # if groupby_column is not None:
    #     data.set_index(groupby_column, inplace=True)
    return data.assign(
        **{
            DBSETS_COL: lambda x: [
                {
                    i
                    for j in (k for k in db_kits if k.can_get_ids)
                    for i in j.get_ids(x)
                    if not pd.isna(i)
                }
                for _, x in data.iterrows()
            ]
        }
    )


def get_all_annotation_ids(ids_by_row) -> dict:
    out = Counter(chain(*ids_by_row[DBSETS_COL].values))
    return out


def get_config_loc():
    return Path(resource_filename("dram2.utils", "CONFIG"))


def check_for_annotations(
    annotation_sets: list[set[str]], annotation_run: dict
) -> Optional[str]:
    """
    Check for required sets of annotations intelligently
    ----------------------------------------------------
    """
    dbs_we_have = set(annotation_run[USED_DBS_TAG])
    you_need = [{j for j in i if j not in dbs_we_have}
                for i in annotation_sets]
    for i in you_need:
        if len(i) < 1:
            return None
    you_need_and = reduce(lambda x, y: x.intersection(y), you_need)
    you_need_or = [
        ", ".join(i - you_need_and) for i in you_need if len(i - you_need_and) > 0
    ]
    if len(you_need_or) < len(you_need_and):
        you_need_or = []
    error_message = "You are trying to use a DRAM2 function that requires specific annotations which this DRAM project does not have yet.\n"
    if len(you_need_and) > 0:
        error_message += (
            f"You need to run annotate with with: [{', '.join(you_need_and)}]"
            "\n"
            f"The command to do that is like `dram2 -o this_output_dir"
            f" annotate --use_db {' --use_db '.join(you_need_and)}`"
            "\n"
        )
        if len(you_need_or) > 0:
            error_message += "Also!\n"
    if len(you_need_or) > 0:
        error_message += f"You need to annotate with: {' or '.join(you_need_or)}\n\n"
    error_message += "You should still review the docs to make sure you are running the program correctly to get results you want."
    return error_message


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


def make_mmseqs_db_for_fasta(
    fasta: Fasta, logger: logging.Logger, threads: int
) -> Fasta:
    if fasta.mmsdb is not None and fasta.mmsdb.exists():
        return fasta
    if fasta.tmp_dir is None:
        raise ValueError(
            "Some how a fasta was passed to the function that makes mmseqs "
            "databases which did not have an associated temporary directory in "
            "which to put that mmseqs-db. Please kindly file a bug report on "
            "GitHub. This indicates that the developer probably made a mistake"
        )
    if fasta.faa is None:
        raise ValueError(
            "Some how a fasta was passed to the function that makes mmseqs "
            "databases which did not have an associated faa directory in which "
            "to put that mmseqs-db. Please kindly file a bug report on GitHub. "
            "This indicates that the developer probably made a mistake"
        )
    fasta.tmp_dir.mkdir(exist_ok=True, parents=True)
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


def has_dup_fasta_name(fastas: list[Fasta], logger: logging.Logger) -> int:
    duplicated = [
        item for item, count in Counter([i.name for i in fastas]).items() if count > 1
    ]
    if len(duplicated) > 0:
        logger.debug(f"duplicated names: {', '.join(duplicated)}")
        return len(duplicated)
    return 0


def get_all_fastas(
    gene_fasta_paths: list[Path],
    # genes_runs: Optional[dict], # may need this someday
    called_fastas: Optional[dict],
    annotation_run: Optional[dict],
    working_dir: Path,
    output_dir: Path,
    use_db: Sequence[str],
    logger: logging.Logger,
    force: bool,
):
    gene_fasta_obs: list[Fasta] = [
        path_to_gene_fastas(i, output_dir / DEFAULT_GENES_FILE)
        for i in gene_fasta_paths
    ]
    # get the precalled genes
    if called_fastas is not None:
        called_fastas_obs: list[Fasta] = [
            Fasta.import_strings(output_dir, *j) for j in called_fastas
        ]
        # now we check these for dups
        if (dup_count := has_dup_fasta_name(called_fastas_obs, logger)) > 0:
            raise DramUsageError(
                f"Genome file names must be unique. There is/are {dup_count} name/s that appear twice in called genes."
            )
    else:
        called_fastas_obs: list[Fasta] = []

    # Combine
    fastas = gene_fasta_obs + called_fastas_obs
    # Stop the user duping the fastas that were called
    if (dup_count := has_dup_fasta_name(fastas, logger)) > 0:
        raise DramUsageError(
            f"Genome file names must be unique. There is/are {dup_count} name/s that appear in both the called genes passed, and the called genes already in the output_dir."
        )

    if annotation_run is not None:
        db_inter = set(
            annotation_run[USED_DBS_TAG]
            if annotation_run[USED_DBS_TAG] is not None
            else []
        ).intersection(set(use_db))
        fasta_inter = set(annotation_run[FASTAS_CONF_TAG]).intersection(
            {i.name for i in fastas}
        )
        if len(db_inter) > 0:
            if force:
                logger.warning(
                    f"You are re-annotating {len(fasta_inter)} of {len(fastas)} FASTAs with databases they were already annotated with: {db_inter}. The past annotations with these databases will be replaced."
                )
            else:
                raise DramUsageError(
                    f"You are trying to re-annotate {len(fasta_inter)} of {len(fastas)} FASTAs with databases they were already annotated with: {db_inter}. You need to use the force flag '-f' in order to do this. If you use that flag the past annotations with these databases will be replaced."
                )
    return fastas


class annotation_meta(dataclass):

    def run_dict():
        pass


def get_last_annotation_run(project_config: dict) -> Optional[dict]:
    annotation_run: Optional[dict] = (
        None
        if project_config.get("annotations") is None
        else project_config["annotations"][project_config["annotations"]["latest"]]
    )
    return annotation_run


LATEST_ANNOTATION_RUN_TAG: str = "latest"


def make_new_project_config(
    old_project_config: dict,
    run_id: str,
    databases: list[DBKit],
    fastas: list[Fasta],
    working_dir: Path,
    annotation_tsv: Path,
    bit_score_threshold: int,
    rbh_bit_score_threshold: int,
    kofam_use_dbcan2_thresholds: bool,
    threads: int,
    force: bool,
    output_dir: Path,
    db_path: Path,
    extra: dict,
    keep_tmp: bool,
) -> dict:
    """
    ----------------------
    TODO:
        Remove db path
        The extra argument is for people's custom DBs it should be updated at some point
    Args:
        old_project_config:
        run_id:
        databases:
        fastas:
        working_dir:
        annotation_tsv:
        bit_score_threshold:
        rbh_bit_score_threshold:
        kofam_use_dbcan2_thresholds:
        threads:
        force:
        output_dir:
        db_path:

    Returns: A new config directory


    """
    new_project_config = old_project_config.copy()
    past_annotation_run = get_past_annotation_run(old_project_config)
    database_names: set = {db.name for db in databases}
    if past_annotation_run is not None and USED_DBS_TAG in past_annotation_run:
        database_names.update(past_annotation_run[USED_DBS_TAG])
    else:
        new_project_config[ANNOTATIONS_TAG] = {}
    new_project_config[ANNOTATIONS_TAG][run_id] = {
        "version": __version__,
        "working_dir": working_dir.relative_to(output_dir).as_posix(),
        FASTAS_CONF_TAG: [i.name for i in fastas],
        ANNOTATION_FILE_TAG: annotation_tsv.relative_to(output_dir).as_posix(),
        "bit_score_threshold": bit_score_threshold,
        "rbh_bit_score_threshold": rbh_bit_score_threshold,
        "kofam_use_dbcan2_thresholds": kofam_use_dbcan2_thresholds,
        "threads": threads,
        "force": force,
        "keep_tmp": keep_tmp,
        USED_DBS_TAG: list(database_names),
    }
    new_project_config[ANNOTATIONS_TAG][LATEST_ANNOTATION_RUN_TAG] = run_id
    new_project_config[FASTAS_CONF_TAG] = [
        i.export(output_dir) for i in fastas]
    return new_project_config


def search_fasta_with_database(databases: list[DBKit], fasta: Fasta) -> pd.DataFrame:
    data = (
        reduce(
            partial(pd.merge, left_index=True, right_index=True, how="outer"),
            [j.search(fasta) for j in databases],
        )
        .assign(**{FASTA_COL: fasta.name})
        .assign(**{GENE_ID_COL: lambda x: x.index})
    )
    data.index = [f"{fasta.name}_{j}" for j in data.index.values]
    return data


def annotate(
    gene_fasta_paths: list[Path],
    dram_config: dict,
    logger: logging.Logger,
    output_dir: Path,
    working_dir: Path,
    cores: int,
    run_id: str,
    project_config: dict,
    keep_tmp: bool = DEFAULT_KEEP_TMP,
    bit_score_threshold: int = DEFAULT_BIT_SCORE_THRESHOLD,
    rbh_bit_score_threshold: int = DEFAULT_RBH_BIT_SCORE_THRESHOLD,
    # past_annotations_path: str = str(None),
    use_db: Sequence[str] = (),
    db_path: Optional[Path] = None,
    custom_fasta_db_name: Sequence = (),
    custom_fasta_db_loc: Sequence = (),
    custom_hmm_db_loc: Sequence = (),
    custom_hmm_db_name: Sequence = (),
    custom_hmm_db_cutoffs_loc: Sequence = (),
    kofam_use_dbcan2_thresholds: bool = False,
    # rename_genes: bool = True,
    # make_new_faa: bool = bool(None),
    force: bool = False,
    extra=None,
    write_config: bool = False,
) -> dict:
    # get assembly locations
    # make mmseqs_dbs
    working_dir.mkdir(exist_ok=True)
    # make a separate testable function for these two
    genes_runs: Optional[dict] = project_config.get("genes_called")
    called_fastas: Optional[dict] = project_config.get(FASTAS_CONF_TAG)
    past_annotation_run = get_past_annotation_run(project_config)
    fastas = get_all_fastas(
        gene_fasta_paths,
        # genes_runs,
        called_fastas,
        past_annotation_run,
        working_dir,
        output_dir,
        use_db,
        logger,
        force,
    )
    logger.info(f"Started annotation with databases: {','.join(use_db)}")
    # initialize all used databases
    databases = [i(dram_config, logger)
                 for i in DB_KITS if i.name in set(use_db)]
    if len(fastas) < 1:
        raise DramUsageError(
            "No FASTAs were passed to the annotator DRAM has nothing to do."
        )

    # cores_for_sub_process:int = cores
    # fasta = [make_mmseqs_db_for_fasta(i, logger=logger, threads=cores_for_sub_process) for i in fastas]

    if cores > len(fastas):
        cores_for_sub_process: int = cores
        cores_for_maping: int = 1
    else:
        cores_for_sub_process: int = 1
        cores_for_maping: int = cores

    with Pool(cores_for_maping) as p:
        fastas = p.map(
            partial(
                make_mmseqs_db_for_fasta, logger=logger, threads=cores_for_sub_process
            ),
            fastas,
        )

    db_args = {
        # "fasta_paths": gene_fasta_paths,
    }

    # add argument for annotations
    for i in databases:
        i.load_dram_config()
        i.set_args(
            output_dir=output_dir,
            working_dir=working_dir,
            bit_score_threshold=bit_score_threshold,
            rbh_bit_score_threshold=rbh_bit_score_threshold,
            kofam_use_dbcan2_thresholds=kofam_use_dbcan2_thresholds,
            threads=cores,
            force=force,
            extra=extra,
            db_path=db_path,  # where to store dbs on the fly
            keep_tmp=keep_tmp,
        )

    # update the config
    if write_config:
        dram_config.update(
            {j: k for i in databases for j, k in i.config.items()})

    # Add all the databases that you are going to
    databases += [
        FastaKit(i, j, dram_config, db_args)
        for i, j in zip(custom_fasta_db_name, custom_fasta_db_loc)
    ]
    db_len_dif = len(custom_hmm_db_name) - len(custom_hmm_db_cutoffs_loc)
    if db_len_dif < 0:
        raise DramUsageError(
            f"There are more hmm cutoff files provided then custom hmm databases provided"
        )

    custom_hmm_db_cutoffs_loc = list(
        custom_hmm_db_cutoffs_loc) + ([None] * db_len_dif)
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
        new_annotations = pd.concat(
            [search_fasta_with_database(databases, fasta=i) for i in fastas]
        )

    # ADD DESCRIPTIONS
    # with Pool(cores) as p:
    #     descriptions: list[pd.DataFrame] = p.map(get_descriptions_for_annotations, databases)
    descriptions: list[pd.DataFrame] = [
        i.get_descriptions(new_annotations) for i in databases
    ]
    new_annotations = reduce(
        partial(pd.merge, left_index=True, right_index=True, how="outer"),
        descriptions + [new_annotations],
    )

    # make a tsv even through that is stupid
    annotation_tsv = output_dir / "annotations.tsv"
    # could be a match statement
    if past_annotation_run is not None and ANNOTATION_FILE_TAG in past_annotation_run:
        logger.info(
            "Found past annotations in project config, DRAM will attempt to merge new annotations"
        )
        past_annotations = pd.read_csv(
            output_dir / past_annotation_run[ANNOTATION_FILE_TAG],
            sep="\t",
            index_col=0,
        )
        annotations = merge_past_annotations(
            new_annotations, past_annotations, logger, force
        )
        # The only case we update past dbs
    elif force and annotation_tsv.exists():
        logger.info(
            "Found past annotations in the output path, DRAM will attempt to force-fully merge new annotations"
        )
        past_annotations = pd.read_csv(
            annotation_tsv,
            sep="\t",
            index_col=0,
        )
        annotations = merge_past_annotations(
            new_annotations, past_annotations, logger, force
        )
    else:
        annotations = new_annotations
    annotations.to_csv(annotation_tsv, sep="\t")
    if genes_runs is not None:
        for i in genes_runs.values():
            i["annotated"] = True
    if not keep_tmp:
        logger.info(f"Removing the temporary directory: {working_dir}.")
        rmtree(working_dir)
    return make_new_project_config(
        project_config,
        run_id,
        databases,
        fastas,
        **db_args,
        annotation_tsv=annotation_tsv,
    )


def merge_past_annotations(
    new_annotations: pd.DataFrame,
    past_annotations: pd.DataFrame,
    logger: logging.Logger,
    force: bool,
) -> pd.DataFrame:
    known_colliders = {FASTA_COL, GENE_ID_COL}
    colliding_columns = set(past_annotations.columns).intersection(
        set(new_annotations.columns)
    )
    problem_colliders = colliding_columns - known_colliders
    if len(problem_colliders) > 0:
        if not force:
            raise DramUsageError(
                f"There is a name collision s for the column/columns:"
                f" ({', '.join(problem_colliders)}). \n\n You need to"
                f" use the force flag to overwrite this data in the"
                f" old annotations file, "
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
            "The old and new annotation files may not have merged correctly!"
            " Check the new annotations file for errors. Did you use the"
            " correct genes.faa for your annotations?"
        )
    all_annotations = pd.concat([all_annotations, past_annotations_appened])
    return all_annotations


@click.command(
    "annotate",
)
@click.argument(
    "gene_fasta_paths",
    type=click.Path(exists=True, path_type=Path),
    nargs=-1,
)
@click.option(
    "-s",
    "--use_dbset",
    multiple=True,
    type=click.Choice(list(DBSETS.keys()), case_sensitive=False),
)
@click.option(
    "--use_db",
    multiple=True,
    default=[],
    type=click.Choice([i.name for i in DB_KITS], case_sensitive=False),
    help=("Specify exactly which DBs to use. This argument can be used"
          "  multiple times, so for example if you want to annotate with"
          "  FeGenie and Camper you would have a command like `dram2 - o"
          "  output/dir annotate --use_db fegenie --use_db camper`,"
          "  the options available are in this help."),
)
@ click.option(
    "--bit_score_threshold",
    type=int,
    default=60,
    help="minimum bit score of search to retain hits",
)
@ click.option(
    "--rbh_bit_score_threshold",
    type=int,
    default=350,
    help="minimum bit score of reverse best hits to retain hits",
)
@ click.option(
    "--custom_fasta_db_name",
    type=str,
    multiple=True,
    help="Names of custom databases, can be used multiple times.",
)
@ click.option(
    "--custom_fasta_db_loc",
    multiple=True,
    type=click.Path(exists=True, path_type=Path),
    help="""
    Location of fastas to annotate against, can be used multiple times but
    "must match number of custom_db_name's
    """,
)
@ click.option(
    "--custom_hmm_db_name",
    multiple=True,
    help="Names of custom hmm databases, can be used multiple times.",
)
@ click.option(
    "--custom_hmm_db_loc",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="Location of HMMs to annotate against, can be used multiple times but"
    "must match number of custom_hmm_name's",
)
@ click.option(
    "--custom_hmm_db_cutoffs_loc",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="""
    Location of file with custom HMM cutoffs and descriptions, can be used
    multiple times.
    """,
)
@ click.option(
    "--tempory_dir",
    type=click.Path(path_type=Path),
    help="""
    Location of the temporary file where the annotations will be stored, this
    file will still be defeated at the end of the annotation process if the
    keep tmp flag is not set.
    """,
)
@ click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Remove all past annotations and annotate again.",
)
@ click.pass_context
def annotate_cmd(
    ctx: click.Context,
    gene_fasta_paths: list[Path],
    use_db: list[str],
    bit_score_threshold: int = DEFAULT_BIT_SCORE_THRESHOLD,
    rbh_bit_score_threshold: int = DEFAULT_RBH_BIT_SCORE_THRESHOLD,
    # log_file_path: str = str(None),
    # past_annotations_path: str = str(None),
    use_dbset: Sequence = (),
    custom_fasta_db_name: Sequence = (),
    custom_fasta_db_loc: Sequence = (),
    custom_hmm_db_loc: Sequence = (),
    custom_hmm_db_name: Sequence = (),
    custom_hmm_db_cutoffs_loc: Sequence = (),
    kofam_use_dbcan2_thresholds: bool = False,
    # make_new_faa: Optional[bool] = None,
    tempory_dir: Optional[Path] = None,
    force: bool = False,
    # db_path: Path = None,
    extra=None,
    # study_set: Sequence = (),
):
    """
    Annotate Genes with Gene Database
    - --

    Get gene identifiers from a set of databases and format them for other
    DRAM2 analysis tools. To use this tool, your genes should already be
    called.


    The annotation process depends on the user's selection. You can use the
    --use_db argument to select a set of databases, or use the use_dbset
    argument to use a pre-configured set of databases.

    """
    context: DramContext = ctx.obj
    run_id: str = get_time_stamp_id(ANNOTATIONS_TAG)
    if tempory_dir is None:
        working_dir: Path = context.get_output_dir() / run_id
    else:
        working_dir: Path = tempory_dir
    # logger = logging.getLogger("dram2_log")
    logger = context.get_logger()
    output_dir: Path = context.get_output_dir()
    keep_tmp: bool = context.keep_tmp
    cores: int = context.cores
    project_config: dict = context.get_project_config()
    dram_config = context.get_dram_config(logger)  # FIX
    if len(use_dbset) > 0:
        use_db = list(use_db) + [i for j in use_dbset for i in DBSETS[j]]

    use_db = list(set(use_db))
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
            # past_annotations_path=past_annotations_path,
            use_db=use_db,
            custom_fasta_db_name=custom_fasta_db_name,
            custom_fasta_db_loc=custom_fasta_db_loc,
            custom_hmm_db_loc=custom_hmm_db_loc,
            custom_hmm_db_name=custom_hmm_db_name,
            custom_hmm_db_cutoffs_loc=custom_hmm_db_cutoffs_loc,
            kofam_use_dbcan2_thresholds=kofam_use_dbcan2_thresholds,
            force=force,
            extra=extra,
        )
        project_config.update(new_config)
        context.set_project_config(project_config)

    except Exception as e:
        logger.error(e)
        logger.exception("Fatal error in annotation")
        raise (e)


@ click.command("list_databases")
def list_databases():
    """
    List available database
    - --


    This is a simple tool to list the databases that are available to use in dram commands, mostly annotations. This includes both annotations' formal name, the key that identifies it to dram, always lowercase and all one word, and the citation if it exists.
    """
    for i in DB_KITS:
        print(
            f'{i.formal_name}:\n   Use the key name "{i.name}" to select\n    Citation: {i.citation}\n\n'
        )
    """

    Args:
        gene_fasta_paths:
        dram_config:
        logger:
        output_dir:
        working_dir:
        cores:
        run_id:
        project_config:
        keep_tmp:
        bit_score_threshold:
        rbh_bit_score_threshold:
        use_db:
        db_path:
        custom_fasta_db_name:
        custom_fasta_db_loc:
        custom_hmm_db_loc:
        custom_hmm_db_name:
        custom_hmm_db_cutoffs_loc:
        kofam_use_dbcan2_thresholds:
        threads:
        force:
        extra:
        write_config:
        new_annotations:
        past_annotations:
        logger:
        force:
        ctx:
        gene_fasta_paths:
        use_db:
        bit_score_threshold:
        rbh_bit_score_threshold:
        use_dbset:
        custom_fasta_db_name:
        custom_fasta_db_loc:
        custom_hmm_db_loc:
        custom_hmm_db_name:
        custom_hmm_db_cutoffs_loc:
        kofam_use_dbcan2_thresholds:
        threads:
        tempory_dir:
        force:
        extra:

    Raises:
        ValueError:

    """
    """

    Args:
        gene_fasta_paths:
        dram_config:
        logger:
        output_dir:
        working_dir:
        cores:
        run_id:
        project_config:
        keep_tmp:
        bit_score_threshold:
        rbh_bit_score_threshold:
        use_db:
        db_path:
        custom_fasta_db_name:
        custom_fasta_db_loc:
        custom_hmm_db_loc:
        custom_hmm_db_name:
        custom_hmm_db_cutoffs_loc:
        kofam_use_dbcan2_thresholds:
        threads:
        force:
        extra:
        write_config:
        new_annotations:
        past_annotations:
        logger:
        force:
        ctx:
        gene_fasta_paths:
        use_db:
        bit_score_threshold:
        rbh_bit_score_threshold:
        use_dbset:
        custom_fasta_db_name:
        custom_fasta_db_loc:
        custom_hmm_db_loc:
        custom_hmm_db_name:
        custom_hmm_db_cutoffs_loc:
        kofam_use_dbcan2_thresholds:
        threads:
        tempory_dir:
        force:
        extra:

    Raises:
        ValueError:

    """
