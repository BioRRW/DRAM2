"""
===============================
Annotate Called Genes with DRAM
===============================

Main control point for the annotation process, and the main way to access the
annotations.tsv or whatever it becomes. When we say that is the main control
point of the annotation process, we mean that it calls all the available
objects in the dram2.db_kits name space. And orchestrates them to make the
annotations file, currently a tsv. As the origin of the annotations TSV it is
also in the best position to parse the annotations.tsv in order to get gene IDs
used in other processes, and to check if the annotations are available.

This is the second or first step of any pipeline and I would argue it is the
heart of DRAM.

example use::

    dram2 annotate --help

Todo:
----
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
from dataclasses import dataclass
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
FORCE_DEFAULT: bool = False
USED_DBS_TAG: str = "used_dbs"
ANNOTATIONS_TAG = "annotations"
LATEST_TAG: str = "latest"
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
VERSION_TAG = "version"
WORKING_DIR_TAG = "working_dir"
BIT_SCORE_THRESHOLD_TAG = "bit_score_threshold"
RBH_BIT_SCORE_THRESHOLD_TAG = "rbh_bit_score_threshold"
FORCE_TAG = "force"
KEEP_TMP_TAG = "keep_tmp"
KOFAM_USE_DBCAN2_THRESHOLDS_TAG = "kofam_use_dbcan2_thresholds"


@dataclass
class AnnotationMeta:

    output_dir: Path
    used_dbs: set[str]
    fasta_names: set[str]
    working_dir: Path
    annotation_tsv: Path
    bit_score_threshold: int
    rbh_bit_score_threshold: int
    kofam_use_dbcan2_thresholds: bool
    force: bool
    keep_tmp: bool
    version: str

    def get_dict(self) -> dict:
        conf_dict = {
            VERSION_TAG: self.version,
            WORKING_DIR_TAG: (self.working_dir.relative_to(self.output_dir).as_posix()),
            FASTAS_CONF_TAG: list(self.fasta_names),
            ANNOTATION_FILE_TAG: (
                self.annotation_tsv.relative_to(self.output_dir).as_posix()
            ),
            BIT_SCORE_THRESHOLD_TAG: self.bit_score_threshold,
            RBH_BIT_SCORE_THRESHOLD_TAG: self.rbh_bit_score_threshold,
            FORCE_TAG: self.force,
            KEEP_TMP_TAG: self.keep_tmp,
            KOFAM_USE_DBCAN2_THRESHOLDS_TAG: self.kofam_use_dbcan2_thresholds,
            USED_DBS_TAG: list(self.used_dbs),
        }
        return conf_dict

    pass


def get_annotation_ids_by_row(data: pd.DataFrame, db_kits: list) -> pd.DataFrame:
    """
    Extract the annotaion IDs from each row.

    Extract the annotaion IDs from each row and return a data frame with a new column
    with name DBSETS_COL containing these sets.

    :param data:
    :param db_kits:
    :returns:
    """
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
    """
    Extract all annotation ids/counts.

    Take the output from get annotation_ids_by_row and combine all the sets there in
    and return the ids and the counts. this will be the count of ids found in the full
    set. This is typicaly used with a groupby.

    :param ids_by_row: the output from get_annotation_ids_by_row
    :returns:
    """
    out = Counter(chain(*ids_by_row[DBSETS_COL].values))
    return out


def check_for_annotations(
    annotation_sets: list[set[str]], annotation_meta: AnnotationMeta
) -> Optional[str]:
    """
    Check for Required Sets of Annotations Intelligently
    ----------------------------------------------------

    This method is not intended for the annotations themselves but for
    downstream processes that depend on these annotations. It is used to check if the list of annotations passed with the
    annotation_sets argument are present.

    :param annotation_sets:
    :param annotation_meta:
    :returns:
    """
    dbs_we_have: set = annotation_meta.used_dbs
    you_need = [{j for j in i if j not in dbs_we_have} for i in annotation_sets]
    for i in you_need:
        if len(i) < 1:
            return None
    you_need_and = reduce(lambda x, y: x.intersection(y), you_need)
    you_need_or = [
        ", ".join(i - you_need_and) for i in you_need if len(i - you_need_and) > 0
    ]
    if len(you_need_or) < len(you_need_and):
        you_need_or = []
    error_message = (
        "You are trying to use a DRAM2 function that requires"
        "specific annotations which this DRAM project does not"
        "have yet.\n"
    )
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
        error_message += (
            f"You need to annotate with:" f" {' or '.join(you_need_or)}\n\n"
        )
    error_message += (
        "You should still review the docs to make sure you are"
        " running the program correctly to get results you want."
    )
    return error_message


def check_fasta_names(fastas: list[Fasta]):
    """
    Are the fasta names unique?

    :param fastas:
    :raises ValueError:
    """
    fasta_names = [i.name for i in fastas]
    if len(fasta_names) != len(set(fasta_names)):
        raise ValueError(
            "Genome file names must be unique. At least one name appears twice"
            " in this search."
        )


def path_to_gene_fastas(fasta_loc: Path, working_dir: Path) -> Fasta:
    """
    Take a path and make a genes fasta object.

    :param fasta_loc:
    :param working_dir:
    :returns:
    """
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
            "databases which did not have an associated temporary directory in"
            " which to put that mmseqs-db. Please kindly file a bug report on "
            "GitHub. This indicates that the developer probably made a mistake"
        )
    if fasta.faa is None:
        raise ValueError(
            "Some how a fasta was passed to the function that makes mmseqs"
            " databases which did not have an associated faa directory in"
            " which to put that mmseqs-db. Please kindly file a bug report on"
            " GitHub. This indicates that the developer probably made a"
            " mistake"
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
    called_fastas: Optional[dict],
    annotation_meta: AnnotationMeta | None,
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
                f"Genome file names must be unique. There is/are {dup_count}"
                f" name/s that appear twice in called genes."
            )
    else:
        called_fastas_obs: list[Fasta] = []

    # Combine
    fastas = gene_fasta_obs + called_fastas_obs
    # Stop the user duping the fastas that were called
    if (dup_count := has_dup_fasta_name(fastas, logger)) > 0:
        raise DramUsageError(
            f"Genome file names must be unique. There is/are {dup_count}"
            " name/s that appear in both the called genes passed, and the"
            "called genes already in the output_dir."
        )

    if annotation_meta is not None:
        db_inter = annotation_meta.used_dbs.intersection(set(use_db))
        fasta_inter = annotation_meta.fasta_names.intersection({i.name for i in fastas})
        if len(db_inter) > 0:
            if force:
                logger.warning(
                    f"You are re-annotating {len(fasta_inter)} of"
                    f" {len(fastas)} FASTAs with databases they were already"
                    f" annotated with: {db_inter}. The past annotations with"
                    " these databases will be replaced."
                )
            else:
                raise DramUsageError(
                    f"You are trying to re-annotate {len(fasta_inter)} of"
                    f" {len(fastas)} FASTAs with databases they were already"
                    f" annotated with: {db_inter}. You need to use the force"
                    f" flag '-f' in order to do this. If you use that flag"
                    f" the past annotations with these databases will be"
                    f" replaced."
                )
    return fastas


def dict_to_annotation_meta(annotation_dict: dict, output_dir: Path) -> AnnotationMeta:
    used_dbs: set[str] = set(annotation_dict[USED_DBS_TAG])
    fasta_names: set[str] = set(annotation_dict[FASTAS_CONF_TAG])
    working_dir: Path = Path(annotation_dict[WORKING_DIR_TAG])
    annotation_tsv: Path = Path(annotation_dict[ANNOTATION_FILE_TAG])
    if ANNOTATION_FILE_TAG not in annotation_dict:
        raise DramUsageError(
            """
            There is no annotations.tsv recorded in the project_config
            provided.\n\n It must be the case that the DRAM directory does not
            contain the result of a successful annotation call.\n Run `dram2
            get_status` to see if annotations have been run on this dram
            directory  or if it is valid at all. If you have called genes but
            not run annotations then run `dram2 - o this_output_dir annotate
            db_set 'distill'` in order to get the minimal annotation set for
            distillation.\n Review the documentation to learn more about the
            required pipeline needed  to run dram distill.
                """
        )
    annotations_path = output_dir / Path(annotation_dict[ANNOTATION_FILE_TAG])
    if not annotations_path.exists():
        raise DramUsageError(
            f"""
            The path to annotations exists but it does not point to a
            annotations file that exists in the dram_directory make sure the
            path to your annotations is at the relive path
            {relative_annotation_path} with respect to the dram_directory:
            {output_dir}.
            """
        )
    if not annotation_tsv.exists():
        raise DramUsageError(
            """
            There an annotations.tsv recorded in the project_config provided,
            but it dose not point to an existing files in the project/output
            directory.\n\n
            If you have called genes and need to re-run annotations then run
            `dram2 - o this_output_dir annotate db_set mini` in order to get
            the minimal annotation set for distillation. \n
            Review the documentation to learn more about the required pipeline
            needed to run dram distill.
            """
        )
    bit_score_threshold: int = annotation_dict[BIT_SCORE_THRESHOLD_TAG]
    rbh_bit_score_threshold: int = annotation_dict[RBH_BIT_SCORE_THRESHOLD_TAG]
    version: str = annotation_dict[VERSION_TAG]
    force: bool = annotation_dict[FORCE_TAG]
    keep_tmp: bool = annotation_dict[KEEP_TMP_TAG]
    kofam_use_dbcan2_thresholds: bool = annotation_dict[KOFAM_USE_DBCAN2_THRESHOLDS_TAG]

    return AnnotationMeta(
        output_dir=output_dir,
        kofam_use_dbcan2_thresholds=kofam_use_dbcan2_thresholds,
        used_dbs=used_dbs,
        fasta_names=fasta_names,
        working_dir=working_dir,
        annotation_tsv=annotation_tsv,
        bit_score_threshold=bit_score_threshold,
        rbh_bit_score_threshold=rbh_bit_score_threshold,
        version=version,
        force=force,
        keep_tmp=keep_tmp,
    )


def get_last_annotation_meta(
    project_meta: dict, output_dir: Path
) -> AnnotationMeta | None:
    annotaions_dicts = project_meta.get(ANNOTATIONS_TAG)
    if annotaions_dicts is None:
        return None
    run_id = annotaions_dicts[LATEST_TAG]
    latest_dict = annotaions_dicts[run_id]
    annotation_meta = dict_to_annotation_meta(latest_dict, output_dir)
    return annotation_meta


def make_new_meta_data(
    project_meta: dict,
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
    - ---------------------
    TODO:
        Remove db path
        The extra argument is for people's custom DBs it should be updated at
        some point.

    Returns: A new config directory
    """
    new_annotation_meta = AnnotationMeta(
        output_dir=output_dir,
        version=__version__,
        used_dbs={i.name for i in databases},
        fasta_names={i.name for i in fastas},
        working_dir=working_dir,
        annotation_tsv=annotation_tsv,
        bit_score_threshold=bit_score_threshold,
        rbh_bit_score_threshold=rbh_bit_score_threshold,
        kofam_use_dbcan2_thresholds=kofam_use_dbcan2_thresholds,
        force=force,
        keep_tmp=keep_tmp,
    )
    past_annotation_meta = get_last_annotation_meta(project_meta, output_dir)
    if past_annotation_meta is not None:
        new_annotation_meta.used_dbs.update(past_annotation_meta.used_dbs)
    if ANNOTATIONS_TAG in project_meta:
        project_meta[ANNOTATIONS_TAG].update(new_annotation_meta.get_dict())
    else:
        project_meta[ANNOTATIONS_TAG] = new_annotation_meta.get_dict()
    project_meta[ANNOTATIONS_TAG][LATEST_TAG] = run_id
    project_meta[FASTAS_CONF_TAG] = [i.export(output_dir) for i in fastas]
    return project_meta


def search_fasta_with_database(databases: list[DBKit], fasta: Fasta) -> pd.DataFrame:
    search_results: list[pd.DataFrame] = [j.search(fasta) for j in databases]
    data = (
        reduce(
            partial(pd.merge, left_index=True, right_index=True, how="outer"),
            search_results,
        )
        .assign(**{FASTA_COL: fasta.name})
        .assign(**{GENE_ID_COL: lambda x: x.index})
    )
    data.set_index([f"{fasta.name}_{j}" for j in data.index.values])
    return data


def annotate(
    gene_fasta_paths: list[Path],
    dram_config: dict,
    logger: logging.Logger,
    output_dir: Path,
    working_dir: Path,
    cores: int,
    run_id: str,
    project_meta: dict,
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
    called_fastas: Optional[dict] = project_meta.get(FASTAS_CONF_TAG)
    past_annotation_meta = get_last_annotation_meta(project_meta, output_dir)
    fastas = get_all_fastas(
        gene_fasta_paths,
        called_fastas,
        past_annotation_meta,
        output_dir,
        use_db,
        logger,
        force,
    )
    logger.info(f"Started annotation with databases: {','.join(use_db)}")
    # initialize all used databases
    databases = [i(dram_config, logger) for i in DB_KITS if i.name in set(use_db)]
    if len(fastas) < 1:
        raise DramUsageError(
            "No FASTAs were passed to the annotator DRAM has nothing to do."
        )

    # cores_for_sub_process:int = cores
    # fasta = [make_mmseqs_db_for_fasta(i, logger=logger,
    # threads=cores_for_sub_process) for i
    # in fastas]

    if cores > len(fastas):
        cores_for_sub_process: int = cores
        cores_for_maping: int = 1
    else:
        cores_for_sub_process: int = 1
        cores_for_maping: int = cores

    with Pool(cores_for_maping) as p:
        fastas = p.map(
            partial(
                make_mmseqs_db_for_fasta,
                logger=logger,
                threads=cores_for_sub_process,
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
        dram_config.update({j: k for i in databases for j, k in i.config.items()})

    # Add all the databases that you are going to
    databases += [
        FastaKit(i, j, dram_config, db_args)
        for i, j in zip(custom_fasta_db_name, custom_fasta_db_loc)
    ]
    db_len_dif = len(custom_hmm_db_name) - len(custom_hmm_db_cutoffs_loc)
    if db_len_dif < 0:
        raise DramUsageError(
            """
            There are more hmm cutoff files provided then custom hmm
            databases provided.
            """
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
            """
            No databases were selected. There is nothing for DRAM to do but
            save progress and exit."""
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
    #     descriptions: list[pd.DataFrame] =
    #     p.map(get_descriptions_for_annotations, databases)
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
    if past_annotation_meta is not None:
        logger.info(
            """
                    Found past annotations in project config, DRAM will attempt
                    to merge new annotations.
                    """
        )
        past_annotations = pd.read_csv(
            output_dir / past_annotation_meta.annotation_tsv,
            sep="\t",
            index_col=0,
        )
        annotations = merge_past_annotations(
            new_annotations, past_annotations, logger, force
        )
        # The only case we update past dbs
    elif force and annotation_tsv.exists():
        logger.info(
            """
                    Found past annotations in the output path, DRAM will
                    attempt to force-fully merge new annotations.
                    """
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
    # TODO: make this a genes modual fuction to do this
    genes_runs: Optional[dict] = project_meta.get("genes_called")
    if genes_runs is not None:
        for i in genes_runs.values():
            i["annotated"] = True
    if not keep_tmp:
        logger.info(f"Removing the temporary directory: {working_dir}.")
        rmtree(working_dir)
    return make_new_meta_data(
        project_meta,
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
                f"""
                There is a name collision s for the column/columns:
                ({', '.join(problem_colliders)}). \n\n You need to
                use the force flag to overwrite this data in the
                old annotations file,
                """
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
            """
            The old and new annotation files may not have merged correctly!
            Check the new annotations file for errors. Did you use the correct
            genes.faa for your annotations?
            """
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
    help=(
        """
        Specify exactly which DBs to use. This argument can be used multiple
        times, so for example if you want to annotate with FeGenie and Camper
        you would have a command like `dram2 - o output/dir annotate --use_db
        fegenie --use_db camper`, the options available are in this help.
        """
    ),
)
@click.option(
    "--bit_score_threshold",
    type=int,
    default=60,
    help="minimum bit score of search to retain hits",
)
@click.option(
    "--rbh_bit_score_threshold",
    type=int,
    default=350,
    help="minimum bit score of reverse best hits to retain hits",
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
    help="""
    Location of fastas to annotate against, can be used multiple times but
    "must match number of custom_db_name's
    """,
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
    help="Location of HMMs to annotate against, can be used multiple times but"
    "must match number of custom_hmm_name's",
)
@click.option(
    "--custom_hmm_db_cutoffs_loc",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="""
    Location of file with custom HMM cutoffs and descriptions, can be used
    multiple times.
    """,
)
@click.option(
    "--tempory_dir",
    type=click.Path(path_type=Path),
    help="""
    Location of the temporary file where the annotations will be stored, this
    file will still be defeated at the end of the annotation process if the
    keep tmp flag is not set.
    """,
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
    - -use_db argument to select a set of databases, or use the use_dbset
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
    project_meta: dict = context.get_project_meta()
    dram_config = context.get_dram_config(logger)  # FIX
    if len(use_dbset) > 0:
        use_db = list(use_db) + [i for j in use_dbset for i in DBSETS[j]]

    use_db = list(set(use_db))
    try:
        new_meta = annotate(
            gene_fasta_paths=gene_fasta_paths,
            dram_config=dram_config,
            project_meta=project_meta,
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
        project_meta.update(new_meta)
        context.set_project_meta(project_meta)

    except Exception as e:
        logger.error(e)
        logger.exception("Fatal error in annotation")
        raise (e)


@click.command("list_databases")
def list_databases():
    """
    List available database
    - --


    This is a simple tool to list the databases that are available to use in dram commands, mostly annotations. This includes both annotations' formal name, the key that identifies it to dram, always lowercase and all one word, and the citation if it exists.
    """
    for i in DB_KITS:
        print(
            f"""
            {i.formal_name}:\n   Use the key name \"{i.name}\" to select\n
            Citation: {i.citation}\n\n
            """
        )
