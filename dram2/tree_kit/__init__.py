"""The main controle point where trees are loaded and functions are defined."""
import logging
from pathlib import Path
from typing import Optional
import os

import click

from dram2.cli.context import (
    DramContext,
    __version__,
    get_time_stamp_id,
)
from dram2.utils import get_package_path, DramUsageError

from dram2.annotate import (
    get_last_annotation_meta,
    AnnotationMeta,
    DB_KITS,
    check_for_annotations,
)
from dram2.tree_kit.dram_phylo_pipe import (
    phylo_tree,
    MIN_DIF_LEN_RATIO_DFLT,
    MAX_LEN_TO_LABEL_DFLT,
)
from dram2.utils.globals import FASTAS_CONF_TAG
from dram2.tree_kit.pplacer import DramTree

TREE_TAG: str = "phylo_trees"

DATA_PATH: Path = get_package_path(Path("tree_kit", "data"))
NXR_NAR_TREE = DramTree(
    name="nxr_nar",
    pplacer_profile=DATA_PATH / "nxr_nar" / "nxr_nar.refpkg",
    target_ids=["K11180", "dsrA", "dsrB", "K11181"],
    target_dbs=[{"kegg", "sulfur"}, {"kofam", "sulfur"}],
    reference_seq=os.path.join(
        DATA_PATH, "nxr_nar", "nxr-nar_seqs_for_tree_aligned.faa"
    ),
    gene_mapping_path=DATA_PATH / "nxr_nar" / "nxr-nar-tree-mapping.tsv",
    color_mapping_path=DATA_PATH / "nxr_nar" / "color_map.tsv",
)
AMOA_PMOA_TREE = DramTree(
    name="amoa_pmoa",
    pplacer_profile=DATA_PATH / "nxr_nar" / "nxr_nar.refpkg",
    target_ids=["K11180", "dsrA", "dsrB", "K11181"],
    target_dbs=[{"kegg"}, {"kofam"}],
    reference_seq=os.path.join(
        DATA_PATH, "nxr_nar", "nxr-nar_seqs_for_tree_aligned.faa"
    ),
    gene_mapping_path=DATA_PATH / "nxr_nar" / "nxr-nar-tree-mapping.tsv",
    color_mapping_path=DATA_PATH / "nxr_nar" / "color_map.tsv",
)
TREES = [NXR_NAR_TREE]
LATEST_TREE_RUN_TAG: str = "latest"


def get_annotations_and_genes_path(
    annotation_meta: Optional[AnnotationMeta],
    genes_list: Optional[list],
    annotations: Optional[Path],
    genes: Optional[Path],
    force: bool,
    logger: logging.Logger,
) -> tuple[Path, Path | list]:
    """
    Select the most apropriate annotations and genes files to use based on user input.

    Todo:
    ----
      - Make this a match statment

    :param annotation_meta:
    :param genes_list:
    :param annotations:
    :param genes:
    :param force:
    :param logger:
    :returns:
    :raises DramUsageError:

    """
    if force:
        logger.warning(
            "Skipping the normal checks for needed annotations "
            "because the force flag was passed."
        )
    elif annotation_meta is not None:
        if (
            db_error := check_for_annotations(NXR_NAR_TREE.target_dbs, annotation_meta)
        ) is not None:
            raise DramUsageError(db_error)
    else:
        raise DramUsageError(
            """
            The project_config does not exist, or you have not annotated. \n Have you
            called genes and annotations?\n\n You must at least annotate with KOfam or
            KEGG in order to run this script.\n\nIf you have done the steps but don't
            have a project_config for whatever reason then you can use the force option
            to skip this check.
            """
        )
    if isinstance(annotations, Path):
        annotation_to_use: Path = annotations
    elif annotation_meta is None:
        raise DramUsageError(
            """
            There is no project_config in the output directory and no annotation.tsv
            was provided in its place. \nIf your output directory has no meta data, you
            must use the - -annotations option to point to an annotations.tsv.
            """
        )
    else:
        annotation_to_use: Path = annotation_meta.annotation_tsv
    if isinstance(genes, Path):
        return annotation_to_use, genes
    elif isinstance(genes_list, list):
        return annotation_to_use, genes_list
    else:
        raise DramUsageError(
            """
            There is no genes directory recorded in the project_config file provided,
            and no genes.faa was provided in its place. \nIf your output directory has
            no project config file, you must use the --genes option to point to a
            genes.faa.
            """
        )


@click.command("phylotree")
# @click.version_option(__version__)
@click.option(
    "-a",
    "annotations",
    type=click.Path(exists=True),
    required=False,
    help="""
    The DRAM annotations.tsv file, it must meet all the requirement for this tool.
    pointing to this not necessary if you use the output directory option. If there is
    no project config in the output directory you will need to use force.
    """,
)
@click.option(
    "-g",
    "genes",
    type=click.Path(exists=True, path_type=Path),
    context_settings=dict(help_option_names=["-h", "--help"]),
    default=None,
    required=False,
    help="The gene fasta file, genes.faa file from dram combine genes.",
)
@click.option(
    "--min_dif_len_ratio",
    type=int,
    required=False,
    help="""
    The minimum ratio, in distances, between the nearest labeled gene and the nearest
    differently labeled gene for a placed gene to adopt that label. This does not apply
    to genes labeled based on placement in a labeled clad. So if the first distance is
    d1 and the second longer distance is d2 if d2/d1 - 1 < min_dif_len_ratio then the
    placed gene will fail to be labeled.
    """,
    default=MIN_DIF_LEN_RATIO_DFLT,
)
@click.option(
    "--max_len_to_label",
    type=int,
    required=False,
    help="""
    The maximum distance that a placed gene can be from the nearest labeled node and
    adopt that label. This does not apply to genes labeled based on placement in a
    labeled clad.
    """,
    default=MAX_LEN_TO_LABEL_DFLT,
)
@click.option(
    "--keep_temp",
    is_flag=True,
    show_default=True,
    default=False,
    help="""
    Keep the temporary directory where pplacer files and other intermediaries are.
    """,
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    show_default=True,
    default=False,
    help="""
    Skip the normal checks on the dram project_config and just try to use the
    annotations and genes provided.
    """,
)
@click.option(
    "--use_tree",
    multiple=True,
    default=[t.name for t in TREES],
    type=click.Choice([t.name for t in TREES], case_sensitive=False),
    help="""
    Specify exactly which trees to use. This argument can be used multiple times, but
    for now there is only one option.""",
)
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
    use_tree: list[str],
):
    """
    Phylogenetic trees to DRAM
    ---

    Some genes can perform more than one function, phylogeny offers a system
    to identify what the metabolic function, an ambiguous gene is performing.

    In order to use this command, you must first complete the following for
    all genes in your project:

    - Call the Genes with `dram2 call` (and don’t remove the genes directory!!)
    - annotate the genes with`dram2 annotate using `dram2 annotate` be
    shure you use the following database:
        - Sulfur and KEGG or KOfam

    Note that at the moment NXR/NAR is the only tree active but PMOA/AMOA is
    in the works.

    """
    context: DramContext = ctx.obj
    logger = context.get_logger()
    output_dir: Path = context.get_output_dir()
    cores: int = context.cores
    project_config: dict = context.get_project_meta()
    dram_config: dict = context.get_dram_config(logger)  # FIX
    try:
        genes_list: Optional[list] = project_config.get(FASTAS_CONF_TAG)
        annotation_meta: Optional[AnnotationMeta] = get_last_annotation_meta(
            project_config, output_dir
        )
        annotations_path, genes_path_list = get_annotations_and_genes_path(
            annotation_meta, genes_list, annotations, genes, force, logger
        )
        run_id: str = get_time_stamp_id(TREE_TAG)
        db_kits_with_ids = [i for i in DB_KITS if i.selectable and i.can_get_ids]
        if annotation_meta is None and not force:
            raise DramUsageError(
                "There is no project_config or it is missing the data from a"
                " annotation run, but you have not used the force flag to skip"
                " checks. We can't check the annotations without the"
                " annotation run data. If you know what you are doing, you can"
                " use the force flag to continue without checks"
            )
        elif annotation_meta is not None and not force:
            dbs_we_have_ano = set(annotation_meta.used_dbs)
            db_kits_with_ids = [
                i for i in db_kits_with_ids if i.name in dbs_we_have_ano
            ]
        else:
            logger.warning("Skipping the normal checks because of the force flag")
        db_kits = [i(dram_config, logger) for i in db_kits_with_ids]
        trees = [t for t in TREES if t.name in use_tree]
        tree_names, tree_paths = phylo_tree(
            annotations_path=annotations_path,
            gene_fasta_list=genes_path_list,
            output_dir=output_dir,
            logger=logger,
            # annotate_all =  annotate_all,
            keep_temp=keep_temp,
            cores=cores,
            # force=force,
            max_len_to_label=max_len_to_label,
            min_dif_len_ratio=min_dif_len_ratio,
            db_kits=db_kits,
            trees=trees,
        )
        if TREE_TAG in project_config:
            project_config[TREE_TAG].update(
                {
                    LATEST_TREE_RUN_TAG: run_id,
                    run_id: {i: j for i, j in zip(tree_names, tree_paths)},
                }
            )
        else:
            project_config.update(
                {
                    TREE_TAG: {
                        LATEST_TREE_RUN_TAG: run_id,
                        run_id: {i: j for i, j in zip(tree_names, tree_paths)},
                    }
                }
            )
        context.set_project_meta(project_config)

    except Exception as e:
        logger.error(e)
        logger.exception("Fatal error in DRAM2 Tree_Kit")
        raise (e)
