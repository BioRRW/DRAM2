"""
merge DRAM Directories


example use::
    conda activate ./dram2_env
    dram2 merge --help
    dram2 -h

Todo
----
    - Add more checks for annotations like bit-score excetera
    - Maybe bring back the force flag.

Test::
    exit()
    continue
    pytest tests/test_merger.py::test_merger_pipe
"""

import logging
from pathlib import Path
from shutil import copytree

import click
import pandas as pd

from dram2.cli.context import (
    DramContext,
    __version__,
    get_time_stamp_id,
    log_error_wraper,
)
from dram2.utils import DramUsageError, Fasta
from dram2.call_genes import DEFAULT_GENES_FILE
from dram2.utils.globals import FASTAS_CONF_TAG
from dram2.annotate import (
    AnnotationMeta,
    merge_past_annotations,
    ANNOTATIONS_TAG,
    get_last_annotation_meta,
)
from functools import partial, reduce
from typing import Iterable


# def merge_gtdb_taxonomy(
#     annotations: pd.DataFrame, gtdb_taxonomy_paths: list[Path], logger: logging.Logger
# ):
#     """
#     This needs to be updated
#
#     :param annotations:
#     :param gtdb_taxonomy_paths:
#     :param logger:
#     """
#     gtdb_taxonomy = pd.concat(
#         [pd.read_csv(i, sep="\t", index_col=0) for i in gtdb_taxonomy_paths]
#     )
#     taxonomy = list()
#     taxonomy_missing_bins = list()
#     for i in annotations.fasta:
#         # add taxonomy
#         if i in gtdb_taxonomy.index:
#             taxonomy.append(gtdb_taxonomy.loc[i, "classification"])
#         else:
#             taxonomy.append(i)
#             taxonomy_missing_bins.append(i)
#     for i in set(taxonomy_missing_bins):
#         logger.warning(
#             f"Bin {i} was not found in taxonomy file, replaced with bin name."
#         )
#     annotations["bin_taxonomy"] = taxonomy


# use merge_files to merge these data:
# genes_fna faa scaffolds gff trna rrnas
# make output gbk dir
# gbk_dir = path.join(output_dir, "genbank")
# mkdir(gbk_dir)
#  for anno in annotations_list:
#       if anno is not None:
#            # TODO: make annotate_fasta generate a genbank dir and
#            #        then copy it's contents, get rid of Annotation.name
#            if path.isfile(anno.gbk_loc):
#                 copy2(anno.gbk_loc, path.join(gbk_dir, "%s.gbk" % anno.name))
#             else:
#                 for gbk_loc in glob(path.join(anno.gbk_loc, "*.gbk")):
#                     copy2(gbk_loc, gbk_dir)


# def merge_checkm_quality(
#     annotations: pd.DataFrame, checkm_quality_path: Path, logger: logging.Logger
# ):
#     checkm_quality = pd.concat(
#         [pd.read_csv(i, sep="\t", index_col=0) for i in checkm_quality_path]
#     )
#     checkm_quality.index = [
#         (str(i).removesuffix(".fa").removesuffix(".fasta").removesuffix(".fna"))
#         for i in checkm_quality.index
#     ]
#
#     completeness = list()
#     contamination = list()
#     quality_missing_bins = list()
#     for i in annotations.fasta:
#         # add completeness and contamination
#         if i in checkm_quality.index:
#             completeness.append(checkm_quality.loc[i, "Completeness"])
#             contamination.append(checkm_quality.loc[i, "Contamination"])
#         else:
#             completeness.append(0)
#             contamination.append(100)
#             quality_missing_bins.append(i)
#         for j in set(quality_missing_bins):
#             logger.warning(
#                 f"Bin {j} was not found in quality file, "
#                 "replaced with completeness 0 and contamination 100."
#             )
#     annotations["bin_completeness"] = completeness
#     annotations["bin_contamination"] = contamination


@click.group(
    "merge",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.argument(
    "dram_dirs",
    type=click.Path(exists=True, path_type=Path, file_okay=False),
    nargs=-1,
)
@click.option(
    "-g",
    "--gene_dir",
    multiple=True,
    type=click.Path(exists=True, path_type=Path, file_okay=False),
    help="""
    Use this flag to point a directory of prodigal called genes. You can
    use this flag as many times as you want to point to as many
    directoryies as you want.
    These derectories will not be checked like files pulled from dram_dir metadata.
    You do not need this if you point to a dram project directory that has metadata
    intact.
        """,
)
@click.option(
    "-a",
    "--annotations_tsv",
    multiple=True,
    type=click.Path(exists=True, path_type=Path, dir_okay=False),
    help="""
    Use this flag to point an annotaitons tsv. You can use this flag as many
    times as you want to point to as many directoryies as you want.
    These derectories will not be checked like files pulled from dram_dir metadata.
    You do not need this if you point to a dram project directory that has metadata
    intact.
    """,
)
@click.option(
    "-r",
    "--rrna",
    multiple=True,
    type=click.Path(exists=True, path_type=Path, file_okay=False),
    help="""
        """,
)
@click.option(
    "-t",
    "--trna",
    multiple=True,
    type=click.Path(exists=True, path_type=Path, file_okay=False),
    help="""
        """,
)
@click.pass_context
def merger_cmd(
    ctx: click.Context,
    dram_dirs: list[Path],
    annotations_tsv: list[Path],
    genes_dir: list[Path],
    rrna: list[Path],
    trna: list[Path],
):
    """
    Merge DRAM Projects
    ----

    You may have separate dram projects that you need to merge. This command
    will be able to do so with variable levels of safety. Merging dram runs is
    tricky because there is no guaranty that the annotations will be the same in all
    the directoryies you want to merge.

    Example of use:"
        conda activate ./dram2_env
        dram2 merge -h
        dram2 -h
    pytest tests/test_merger.py
    """
    context: DramContext = ctx.obj
    log_error_wraper(merger_pipe, context, "merging annotations and genes")(
        context=context,
        dram_dirs=dram_dirs,
        genes_dir=genes_dir,
        rrna=rrna,
        trna=trna,
        annotations_tsv=annotations_tsv,
    )


def merger_pipe(
    context: DramContext,
    dram_dirs: list[Path],
    annotations_tsv: list[Path] | None = None,
    genes_dir: list[Path] | None = None,
    rrna: list[Path] | None = None,
    trna: list[Path] | None = None,
):
    """
    The pipe line for merging.


    test::

            pytest tests/test_merger.py::test_merger_pipe

    :param context:
    :param dram_dirs:
    :param annotations_tsv:
    :param genes_dir:
    :param rrna:
    :param trna:
    :raises DramUsageError:
    """
    if rrna is None:
        rrna: list[Path] = []
    if trna is None:
        trna: list[Path] = []
    if genes_dir is None:
        genes_dir: list[Path] = []
    if annotations_tsv is None:
        annotations_tsv: list[Path] = []
    if len(rrna) + len(trna) > 0:
        raise DramUsageError("Sorry rRNA and tRNA ")
    genes: list[Fasta] = []
    annotations: list[AnnotationMeta] = []
    if len(dram_dirs) + len(annotations_tsv) + len(genes_dir) < 1:
        raise DramUsageError(
            "No DRAM directories were specified, no annotaions were specified, "
            "and no genes directorys were specified. You must use use one of"
            " these aguments successfully."
        )
    genes += [j for i in dram_dirs for j in get_genes_from_meta(i)]
    genes += [j for i in genes_dir for j in get_genes_from_dir(i)]
    logger = context.get_logger()
    annotations += [i for i in get_annotation_from_meta(dram_dirs, logger)]
    annotations += [get_annotation_from_file(i) for i in annotations_tsv]
    project_meta: dict = context.get_project_meta()
    output_dir: Path = context.get_dram_dir()
    make_new_meta(context, project_meta, output_dir, genes, annotations)
    merger(genes, annotations, output_dir, logger)
    breakpoint()


def make_new_meta(
    context: DramContext,
    project_meta: dict,
    output_dir: Path,
    fastas: list[Fasta],
    annotation_meta: list[AnnotationMeta],
):
    run_id: str = get_time_stamp_id(ANNOTATIONS_TAG)
    # called_fastas: list | None = project_meta.get(FASTAS_CONF_TAG)

    # out_annotation_meta = get_last_annotation_meta(project_meta, output_dir)
    # if out_annotation_meta is not None:
    #     for i in annotation_meta:
    #         out_annotation_meta.used_dbs.update(i.used_dbs)
    # if ANNOTATIONS_TAG in project_meta:
    #     project_meta[ANNOTATIONS_TAG].update(new_annotation_meta.get_dict())
    # else:
    #     project_meta[ANNOTATIONS_TAG] = new_annotation_meta.get_dict()
    # project_meta[ANNOTATIONS_TAG][LATEST_TAG] = run_id
    # project_meta[FASTAS_CONF_TAG] = [i.export(output_dir) for i in fastas]
    # return project_meta
    # project_meta[FASTAS_CONF_TAG] = [i.export(output_dir) for i in genes]
    # context.set_dram_config(project_meta)


def get_meta(dir: Path) -> dict:
    """
    Get metadata

    :param dir: Dram dir
    :returns: the meta data dict
    """
    return DramContext(dir, None).get_project_meta()


def get_annotation_from_meta(
    dirs: list[Path], logger: logging.Logger
) -> Iterable[AnnotationMeta]:
    """
    Get all the annotations that we are merging, as specified by meta.

    :param project_meta:
    """
    for dir in dirs:
        project_meta: dict = get_meta(dir)
        past_annotation_meta = get_last_annotation_meta(project_meta, dir)
        if past_annotation_meta is None:
            logger.warning(f"There is no annotations in the DRAM-dir: {dir.as_posix()}")
            continue
        yield past_annotation_meta


def get_annotation_from_file(annotations_tsv: Path) -> AnnotationMeta:
    """
    Get all the annotations that we are merging, as specified by path.

    :param annotations_tsv:
    """
    return AnnotationMeta(
        output_dir=annotations_tsv.parent,
        version=__version__,
        used_dbs=set(),
        fasta_names=set(),
        working_dir=annotations_tsv.parent,
        annotation_tsv=annotations_tsv.relative_to(annotations_tsv.parent),
        bit_score_threshold=-1,
        rbh_bit_score_threshold=-1,
        kofam_use_dbcan2_thresholds=False,
        force=False,
        keep_tmp=False,
    )


def path_to_gene_fastas(fasta_dir: Path) -> Fasta:
    """
    Take a path and make a genes fasta object.

    Todo:
    ----

        - Make this into a universal for merger and this
    :param fasta_loc:
    :param working_dir:
    :returns:
    """
    fasta_name = fasta_dir.stem
    fasta_working_dir = fasta_dir.relative_to(fasta_dir.parent)
    return Fasta(
        fasta_name,
        fasta_dir / "gene.faa",
        fasta_working_dir,
        fasta_working_dir / "gene.faa",
        fasta_working_dir / "gene.fna",
        fasta_working_dir / "gene.gff",
        fasta_working_dir / "gene.msdb",
    )


def get_genes_from_dir(
    dir: Path,
) -> list[Fasta]:
    """
    Get all the genes that we are merging as specified by path

    :param genes_dir:
    :param project_meta:
    """
    gene_fasta_obs: list[Fasta] = [path_to_gene_fastas(i) for i in dir.iterdir()]
    return gene_fasta_obs


def get_genes_from_meta(dir: Path) -> list[Fasta]:
    """
    Get all the genes that we are merging as specified by meta

    :param genes_dir:
    :param project_meta:
    """
    project_meta: dict = get_meta(dir)
    called_fastas: list | None = project_meta.get(FASTAS_CONF_TAG)
    if called_fastas is None:
        raise DramUsageError(
            f"No genes data found, for DRAM directory {dir}. Try mergeing another way"
        )
    called_fastas_obs: list[Fasta] = [
        Fasta.import_strings(dir, *j) for j in called_fastas
    ]
    return called_fastas_obs


def get_rrna_to_merge() -> list[Path]:
    """
    Get all the rRNA that we are merging

    :param genes_dir:
    :param project_meta:
    """
    return []


def get_trna_to_merge() -> list[Path]:
    """
    Get all the tRNA that we are merging

    :param genes_dir:
    :param project_meta:
    """
    return []


def merger(
    fastas: list[Fasta],
    annotations_meta: list[AnnotationMeta],
    output_dir: Path,
    logger: logging.Logger,
):
    merge_annotations([i.annotation_tsv for i in annotations_meta], logger)
    genes_dir = output_dir / DEFAULT_GENES_FILE
    merge_genes_dir([i.tmp_dir for i in fastas], genes_dir)


def merge_annotations(annotations_paths: list[Path], logger):
    merge_two = partial(
        merge_past_annotations,
        logger=logger,
        force=True,
    )
    return reduce(merge_two, annotations_paths)


def merge_genes_dir(annotations_path: list[Path], output):
    for i in annotations_path:
        copytree(i, output)
