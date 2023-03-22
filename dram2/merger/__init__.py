"""
merge DRAM Directories


example use::
    conda activate ./dram2_env
    dram2 merge --help
"""

import logging
from pathlib import Path

import click
import pandas as pd

from dram2.cli.context import (
    DramContext, DEFAULT_KEEP_TMP, __version__, get_time_stamp_id)
from dram2.utils import merge_files


def merge_gtdb_taxonomy(annotations: pd.DataFrame,
                        gtdb_taxonomy_paths: list[Path],
                        logger: logging.Logger):
    gtdb_taxonomy = pd.concat(
        [pd.read_csv(i, sep="\t", index_col=0)
         for i in gtdb_taxonomy_paths]
    )
    taxonomy = list()
    taxonomy_missing_bins = list()
    for i in annotations.fasta:
        # add taxonomy
        if i in gtdb_taxonomy.index:
            taxonomy.append(gtdb_taxonomy.loc[i, "classification"])
        else:
            taxonomy.append(i)
            taxonomy_missing_bins.append(i)
    for i in set(taxonomy_missing_bins):
        logger.warning(
            f"Bin {i} was not found in taxonomy file, replaced with bin name."
        )
    annotations["bin_taxonomy"] = taxonomy


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


def merge_checkm_quality(annotations: pd.DataFrame,
                         checkm_quality_path: Path,
                         logger: logging.Logger
                         ):
    checkm_quality = pd.concat(
        [pd.read_csv(i, sep="\t", index_col=0) for i in checkm_quality_path]
    )
    checkm_quality.index = [
        (str(i)
         .removesuffix(".fa")
         .removesuffix(".fasta")
         .removesuffix(".fna")
         ) for i in checkm_quality.index]

    completeness = list()
    contamination = list()
    quality_missing_bins = list()
    for i in annotations.fasta:
        # add completeness and contamination
        if i in checkm_quality.index:
            completeness.append(checkm_quality.loc[i, "Completeness"])
            contamination.append(checkm_quality.loc[i, "Contamination"])
        else:
            completeness.append(0)
            contamination.append(100)
            quality_missing_bins.append(i)
        for j in set(quality_missing_bins):
            logger.warning(
                f"Bin {j} was not found in quality file, "
                "replaced with completeness 0 and contamination 100."
            )
    annotations["bin_completeness"] = completeness
    annotations["bin_contamination"] = contamination


@click.group('merge',
             )
@click.option('--safe/--un_safe', is_flag=True)
@click.argument(
    "gene_fasta_paths",
    type=click.Path(exists=True, path_type=Path),
    nargs=-1,
)
@click.option(
    "--gene_dir",
    multiple=True,
    type=click.Path(exists=True, path_type=Path, file_okay=False),
    help="""
        Use this flag to point a directory of prodigal called genes. You can
        use this flag as many times as you want to point to as many
        directoryies as you want.
        """,
)
@click.option(
    "--annotations_tsv",
    multiple=True,
    type=click.Path(exists=True, path_type=Path, file_okay=False),
    help="""
        Use this flag to point an annotaitons tsv. You can
        use this flag as many times as you want to point to as many
        directoryies as you want.
        """,
)
@ click.pass_context
def merger_cmd(
    ctx: click.Context,
):
    """
    Merge DRAM Projects
    ----

    You may have separate dram projects that you need to merge. This command
    will be able to do so with variable levels of safety. Merging dram runs is
    tricky because there is no guaranty that the annotations will the same.
    """
    context: DramContext = ctx.obj
    run_id: str = get_time_stamp_id(ANNOTATIONS_TAG)
    logger = context.get_logger()
    output_dir: Path = context.get_output_dir()
    project_config: dict = context.get_project_config()
    try:
        merger()
    except Exception as e:
        logger.error(e)
        logger.exception("Fatal error in merging")
        raise (e)


def get_merge_files_from_meta(project_meta: {}):
    pass


def merger():
    pass
