__version__ = "0.0.1"

from dram2.tree_kit.dram_phylo_pipe import (
    tree_kit,
    MIN_DIF_LEN_RATIO_DFLT,
    MAX_LEN_TO_LABEL_DFLT,
)

import click


@click.command("phylo_tree")
@click.version_option(__version__)
@click.option(
    "-a",
    "dram_annotations",
    type=click.Path(exists=True),
    required=False,
    help="The DRAM annotations file, not necessary"
    "  if you use the dram_directory option",
)
@click.option(
    "-g",
    "gene_fasta",
    type=click.Path(exists=True),
    required=False,
    help="The gene fasta file, genes.faa file from dram output.",
)
@click.option(
    "-d",
    "dram_directory",
    type=click.Path(exists=True),
    required=False,
    help="The dram input file, with no names changed so it contains annotations.txt and genes.faa, genes.fna",
)
@click.option(
    "-o",
    "--output_dir",
    default="./",
    type=click.Path(exists=False),
    required=False,
    help="The output directory, includes new annotations, phylo_xml files, and log",
)
@click.option(
    "-o",
    "output_dir",
    type=click.Path(exists=False),
    required=False,
    help="The output annotations file",
)
@click.option(
    "-c",
    "--cores",
    type=int,
    required=False,
    help="The number of cores to use",
    default=10,
)
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
@click.option(
    "--annotate_all",
    is_flag=True,
    show_default=True,
    default=False,
    help="Don't place just the ambiguous genes, place all of them",
)
@click.option(
    "--keep_temp",
    is_flag=True,
    show_default=True,
    default=False,
    help="",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    show_default=True,
    default=False,
    help="Don't place just the ambiguous genes, place all of them",
)
@click.option(
    "--output_dir",
    default="./",
    help="Don't place uncertain genes, place all of them",
)
def phylo_tree(**args):
    tree_kit(**args)
