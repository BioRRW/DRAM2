from dram2.tree_kit.dram_phylo_pipe import (
    phylo_tree,
    MIN_DIF_LEN_RATIO_DFLT,
    MAX_LEN_TO_LABEL_DFLT,
)

import click
import logging
from pathlib import Path
from typing import Optional
import os

from dram2.cli.context import (
    DramContext,
    __version__,
    get_time_stamp_id,
)
from dram2.annotate import ANNOTATIONS_TAG
from dram2.utils.utils import Fasta, DramUsageError
from dram2.annotate import (
    get_past_annotation_run,
    ANNOTATION_FILE_TAG,
    check_for_annotations,
    DISTILLATION_MIN_SET,
    DISTILLATION_MIN_SET_KEGG,
)
from dram2.utils.globals import FASTAS_CONF_TAG


@click.command("phylotree")
# @click.version_option(__version__)
@click.option(
    "-a",
    "annotations",
    type=click.Path(exists=True),
    required=False,
    help="The DRAM annotations file, not necessary"
    "  if you use the dram_directory option",
)
@click.option(
    "-g",
    "genes",
    type=click.Path(exists=True, type=Optional[Path]),
    default=None,
    required=False,
    help="The gene fasta file, genes.faa file from dram combine genes.",
)
@click.option(
    "-d",
    "dram_directory",
    type=click.Path(exists=True),
    required=False,
    help="The dram input file, with no names changed so it contains annotations.txt and genes.faa, genes.fna",
)
# @click.option(
#     "-o",
#     "--output_dir",
#     default="./",
#     type=click.Path(exists=False),
#     required=False,
#     help="The output directory, includes new annotations, phylo_xml files, and log",
# )
# @click.option(
#     "-o",
#     "output_dir",
#     type=click.Path(exists=False),
#     required=False,
#     help="The output annotations file",
# )
# @click.option(
#     "-c",
#     "--cores",
#     type=int,
#     required=False,
#     help="The number of cores to use",
#     default=10,
# )
@click.option(
    "--min_dif_len_ratio",
    type=int,
    required=False,
    help="The minimum ratio, in distances, between the nearest labeled gene and the nearest differently labeled gene for a placed gene to adopt that label. This does not apply to genes labeled based on placement in a labeled clad. So if the first distance is d1 and the second longer distance is d2 if d2/d1 - 1 < min_dif_len_ratio then the paced gene will fail to be labeled.",
    default=MIN_DIF_LEN_RATIO_DFLT,
)
@click.option(
    "--max_len_to_label",
    type=int,
    required=False,
    help="The maximum distance that a placed gene can be from the nearest labeled node and adopt that label. This does not apply to genes labeled based on placement in a labeled clad.",
    default=MAX_LEN_TO_LABEL_DFLT,
)
# @click.option(
#     "--annotate_all",
#     is_flag=True,
#     show_default=True,
#     default=False,
#     help="Don't place just the ambiguous genes, place all of them",
# )
@click.option(
    "--keep_temp",
    is_flag=True,
    show_default=True,
    default=False,
    help="Keep the tempory directory where these trees are kept",
)
# @click.option(
#     "-f",
#     "--force",
#     is_flag=True,
#     show_default=True,
#     default=False,
# )
# @click.option(
#     "--output_dir",
#     default="./",
#     help="Don't place uncertain genes, place all of them",
# )
def get_annotations_and_genes_path(
    project_config: dict,
    annotations: Optional[Path],
    genes: Optional[Path],
    force: bool,
    output_dir: Path,
) -> tuple[Path, Path | list]:
    if annotations is not None:
        if not force:
            raise DramUsageError(
                "You must use the --force flag if you pass an annotations and don't use an output directory with a project config file in it. This is to be shure you know what you are doing, and know that normal checks will be skiped."
            )
        if genes is None:
            raise DramUsageError(
                "You must use the --genes option if you pass an annotations file. The genes dir is used for phylo_genetic placement."
            )
        return annotations, genes
    if genes is not None:
        raise DramUsageError(
            "You must use the --annotations option to point to an annotations file if you pass an genes instead of use an output directory with a project config file in it."
        )
    if (annotation_run := get_past_annotation_run(project_config)) is None:
        raise DramUsageError(
            "The project_config dose not exist or you have not annotated. Have you called genes and annotations? You must at least annotate with KOfam or KEGG in order to run this script."
        )

    if (
        db_error := check_for_annotations(
            [set(["kegg"]), set(["kofam"])], annotation_run
        )
    ) is not None:
        raise DramUsageError(db_error)
    relative_annotation_path: Optional[str] = annotation_run.get(ANNOTATION_FILE_TAG)
    relative_genes_path: Optional[str] = annotation_run.get(ANNOTATION_FILE_TAG)
    if relative_annotation_path is None:
        raise DramUsageError(
            "There is no annotations.tsv recorded in the project_config provided.\n\n"
            "It must be the case that the DRAM directory dose not contain the result of "
            "a successful annotation.\n"
            "Run `dram2 get_status` to see if annotations have been run on this dram directory "
            "or if it is valid at all. "
        )
    if relative_genes_path is None:
        raise DramUsageError(
            "There is no annotations.tsv recorded in the project_config provided.\n\n"
            "It must be the case that the DRAM directory dose not contain the result of "
            "a successful annotation or DRAM2 called genes.\n"
        )
    annotations_path = output_dir / relative_annotation_path
    genes_list: Optional[dict] = project_config.get(FASTAS_CONF_TAG)
    if genes_list is None:
        raise DramUsageError("There are no called genes in the cofig file")

    if not annotations_path.exists():
        raise DramUsageError(
            f"The path to annotations exists but it dose not point to a annotations file that exists in the dram_directory make sure the path to your annotations is at the relive path {relative_annotation_path} with respect to the dram_directory: {output_dir}."
        )
    if force:
        logging.warning(
            "Skipping the normal checks for needed annotations "
            "because the force flag was passed."
        )
    else:
        if (
            db_error := check_for_annotations(
                [DISTILLATION_MIN_SET_KEGG, DISTILLATION_MIN_SET], annotation_run
            )
        ) is not None:
            raise DramUsageError(db_error)

    return annotations_path, genes_list


DEFAULT_COMBINED_GENES_NAME = "genes.faa"
from Bio import SeqIO
from collections.abc import Iterator


def rename_fastas(fastas) -> Iterator:
    for fasta in fastas:
        with open(fasta.faa, "w") as faa_feed:
            for record in SeqIO.parse(faa_feed, "fasta"):
                record.id = f"{fasta.name}_{record.id}"
                yield record


def combine_genes(
    genes_list: list,
    output_dir: Path,
    genes_out_path_name: str | Path = DEFAULT_COMBINED_GENES_NAME,
) -> Path:
    fastas = [Fasta.import_strings(output_dir, *j) for j in genes_list]
    if isinstance(genes_out_path_name, str):
        genes_out: Path = output_dir / genes_out_path_name
    else:
        genes_out: Path = output_dir / genes_out_path_name
    seqs = [i for i in rename_fastas(fastas)]
    with open(genes_out, "w") as faa_feed:
        SeqIO.write(seqs, faa_feed, "fasta")
    return genes_out


TREE_TAG: str = "phylo_trees"


@click.pass_context
def phylo_tree_cmd(
    ctx: click.Context,
    annotations: Optional[Path],
    genes: Optional[Path],
    # annotate_all: bool = False,
    keep_temp: bool,
    force: bool,
    max_len_to_label,
    min_dif_len_ratio,
):
    """
    Philogenetic trees to DRAM
    ___

    Some genes can perform more than one function, phyloginy offers a system to identifie what the metibolic function an abiguase gene is perfoming.

    In order to use this command you must first complete the falowing for all genes in your project:
    - Call the Genes with `dram2 call`
    - annotate the genes with`dram2 annotate using `dram2 annotate` be shure you use the fallowing database:
        - KEGG or KOfam

    Note that at the moment NXR/NAR is the only tree active but PMOA/AMOA is in the works.

    """
    context: DramContext = ctx.obj
    logger = context.get_logger()
    output_dir: Path = context.get_output_dir()
    cores: int = context.cores
    project_config: dict = context.get_project_config()
    try:
        annotations, genes_path_list = get_annotations_and_genes_path(
            project_config, annotations, genes, force, output_dir
        )
        if isinstance(genes_path_list, list):
            genes_fasta: Path = combine_genes(genes_path_list, output_dir)
        else:
            genes_fasta: Path = genes_path_list
        run_id: str = get_time_stamp_id(TREE_TAG)
        tree_names = phylo_tree(
            annotations=annotations,
            gene_fasta=genes_fasta,
            output_dir=output_dir,
            logger=logger,
            # annotate_all =  annotate_all,
            keep_temp=keep_temp,
            cores=cores,
            # force=force,
            max_len_to_label=max_len_to_label,
            min_dif_len_ratio=min_dif_len_ratio,
        )
        tree_names
        project_config.update({TREE_TAG: {run_id: {"names": tree_names}}})

    except Exception as e:
        logger.error(e)
        logger.exception("Fatal error in DRAM2 Tree_Kit")
        raise (e)
