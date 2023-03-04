"""
This program adds profiles based on phylogenetic trees into dram. The key to this process is pplacer which places leaves into pre-existing trees

NOTE pplacer uses about 1/4 of the memory when placing on a FastTree tree as compared to a RAxML tree inferred with GTRGAMMA. If your reads are short and in a fixed region, the memory used by pplacer v1.1 alpha08 (or later) scales with respect to the total number of non-gap columns in your query alignment. You can also make it use less memory (and run faster) by cutting down the size of your reference tree.

TODO Allow this to use the genes directory instead of this
TODO Update the clade data type to be more representative
TODO Maybe switch to FastTree -- or don't
TODO Make database check tree specific
"""

# import tempfile
# import matplotlib.pyplot as plt
# from io import StringIO
# from mag_annotator.utils import setup_logger
# from mag_annotator.pull_sequences import pull_sequences
# from pyvis.network import Network
# import networkx as nx
import os
from collections.abc import Iterator
from collections import namedtuple
from functools import partial

import logging
import pandas as pd
from tempfile import TemporaryDirectory
from Bio import Phylo as phy
from Bio.Phylo.BaseTree import Clade
from Bio.Phylo.Newick import Tree
from Bio import SeqIO
from pathlib import Path

from dram2.annotate import get_annotation_ids_by_row, DBSETS_COL
from dram2.tree_kit.pplacer import DramTree
from dram2.utils.utils import run_process, Fasta
from dram2.db_kits.utils import DBKit

DEFAULT_COMBINED_GENES_NAME = "genes.faa"


def read_rename_fasta(fasta: Fasta) -> Iterator:
    if fasta.faa is None or not fasta.faa.exists():
        raise ValueError("Missing faa error")
    with open(fasta.faa, "r") as faa_feed:
        for record in SeqIO.parse(faa_feed, "fasta"):
            record.id = f"{fasta.name}_{record.id}"
            yield record


def combine_genes(
    genes_list: list,
    output_dir: Path,
    work_dir: Path,
    genes_out_path_name: str | Path = DEFAULT_COMBINED_GENES_NAME,
) -> Path:
    fastas = [Fasta.import_strings(output_dir, *j) for j in genes_list]
    genes_out: Path = work_dir / genes_out_path_name
    with open(genes_out, "w") as faa_feed:
        for fasta in fastas:
            SeqIO.write([i for i in read_rename_fasta(fasta)], faa_feed, "fasta")
    return genes_out


UNPLACE_LABEL = "UNPLACEABLE"

# pplacer/guppy columns in one place for reference
PLACEMAT_NAME_1_COL = "name1"
PLACEMAT_NAME_2_COL = "name2"
PLACEMAT_DISTANCE_COL = "distance"
MIN_DIF_LEN_RATIO_DFLT: float = 0.10
MAX_LEN_TO_LABEL_DFLT: float = 8
# prefixes for the notes section of the output
PROXIMITY_INFO_PREFIX: str = "Placed based on proximity to labeled genes:"
CLADE_INFO_PREFIX: str = "Placed based destination clade:"
UNPLACE_PREFIX: str = "Can't be placed:"


def phylo_tree(
    annotations_path: Path,
    gene_fasta_list: Path | list,
    output_dir: Path,
    logger: logging.Logger,
    # annotate_all: bool,
    keep_temp: bool,
    cores: int,
    # force: bool,
    max_len_to_label: float,
    min_dif_len_ratio: float,
    db_kits: list[DBKit],
    trees: list[DramTree],
) -> tuple[list[str], list[str]]:
    logger.info("Start phylogenetic tree disambiguation")

    logger.info("Processing annotations")
    annotations = pd.read_csv(annotations_path, sep="\t", index_col=0)
    annotation_ids: pd.Series = get_annotation_ids_by_row(annotations, db_kits)[
        DBSETS_COL
    ]
    tree_paths = []
    with TemporaryDirectory(dir=output_dir, prefix="tmp_") as work_dir_str:
        work_dir: Path = Path(work_dir_str)
        if isinstance(gene_fasta_list, list):
            gene_fasta: Path = combine_genes(gene_fasta_list, output_dir, work_dir)
        else:
            gene_fasta: Path = gene_fasta_list
        for tree in trees:
            # For now we need to combine genes if they are separate. It would be super cool if we did not.
            logger.info(f"Performing phylogenetic disambiguation with tree {tree.name}")
            tree.set_logger(logger)
            logger.info(f"Made temporary files in {work_dir}")
            ids_keep = set(
                annotation_ids[
                    annotation_ids.apply(
                        lambda x: len(x.intersection(tree.target_ids)) > 0
                    )
                ].index
            )
            if len(ids_keep) < 1:
                logger.info(
                    f"No enigmatic genes were found relating to {tree.name}. Skipping this tree."
                )
                tree_paths += [None]
                continue
            trimed_fa = extract_enigmatic_genes(ids_keep, gene_fasta, work_dir, logger)
            logger.info("Placing enigmatic genes")
            jplace_file = tree.pplacer_place_sequences(
                trimed_fa.as_posix(), work_dir_str, threads=cores
            )
            treeph = read_phtree(jplace_file, work_dir_str, logger)
            edpl = read_edpl(jplace_file, work_dir_str, logger)
            known_terminals, placed_terminals = color_known_termininals(
                treeph, annotations, tree.mapping
            )

            logger.info("Applying labels based on location")
            # treeph, path_df = color_tree_paths(treeph, tree)
            get_clade_info(treeph.root, tree)
            find_all_nearest(treeph, known_terminals, placed_terminals)
            logger.info("Applying labels based on proximity")
            tree_df = make_df_of_tree(
                placed_terminals, tree.name, edpl, max_len_to_label, min_dif_len_ratio
            )
            logger.info("Writing output, to {output_dir.as_posix()}")
            write_files(
                tree.name,
                tree_df,
                treeph,
                jplace_file,
                output_dir,
                work_dir_str,
                keep_temp,
            )
            tree_paths += [f"{tree.name}_tree_data.tsv"]
            logger.info(end_message(tree_df, tree.name))
    tree_names = [tree.name for tree in trees]
    return tree_names, tree_paths


def extract_enigmatic_genes(
    ids_keep: set[str],
    gene_fasta: Path,
    work_dir: Path,
    logger: logging.Logger,
) -> Path:
    """
    :param ids_keep: Enigmatic ids from annotations from dram run
    :param gene_fasta: faa from dram run
    :param work_dir: Temp files here
    :param logger: Standard DRAM logger
    :returns: The path to the ambiguous genes in a fasta file

    Takes in a fasta file of genes and a list of ids in the dram annotation, and returns a filtered fasta to match.
    """
    logger.info("Finding enigmatic genes")
    output_fasta = work_dir / "trim.faa"
    with open(gene_fasta, "r") as faa_feed, open(output_fasta, "w") as out_faa:
        output_fasta_generator = (
            i for i in SeqIO.parse(faa_feed, "fasta") if i.id in ids_keep
        )
        SeqIO.write(output_fasta_generator, out_faa, "fasta")
    return output_fasta


def read_phtree(
    jplace_file: str,
    work_dir: str,
    logger: logging.Logger,
):
    placed_tree_file = os.path.join(work_dir, "placed_tree.nh")
    _ = run_process(
        [
            "guppy",
            "tog",
            "-o",
            placed_tree_file,
            jplace_file,
        ],
        logger,
        capture_stdout=True,
    )
    treeph = phy.read(placed_tree_file, format="newick")
    # Combine gene mapping and color mapping into an omni mapping file
    return treeph


def find_all_nearest(
    treeph: Tree, known_terminals: list[Clade], placed_terminals: list[Clade]
):
    apply_names(treeph.root)  # giving names helps debug the tree
    parents = all_parents(treeph)
    for clade in placed_terminals:
        clade.nearest = []
        nearest = find_distances(clade, known_terminals, parents)
        clade.nearest.append(nearest)
        if clade.label is None:
            clade.nearest.append(
                find_distances(clade, known_terminals, parents, nearest.end.label)
            )


def color_known_termininals(treeph, annotations: pd.DataFrame, mapping: pd.DataFrame):
    terminals = treeph.get_terminals()
    known_terminals = [i for i in terminals if (i.name in mapping.index)]
    placed_terminals = [i for i in terminals if (i.name in annotations.index)]
    for cl in known_terminals:
        cl.color = (
            None if cl.name not in mapping.index else mapping.loc[cl.name, "color"]
        )
    return known_terminals, placed_terminals


def read_edpl(
    jplace_file: str,
    work_dir: str,
    logger: logging.Logger,
):
    """
    There are not as many options for certainty as first thought, but in any case I will for now go with the EDPL.

    So the deal with uncertainty is that it may only mean something if the distance to the contrary node is smaller than some ratio of the EDPL. Of cores, it is not clear if the distances would even be distributed equally along all paths.

    In any case, there is no reason to think that these values are not a measure of uncertainty, and the user will probably have use for them but it may mean that in the future we provide additional pplacer output.

    The best solution and one that we could certainly achieve, is to use pplacer to place and then


    I think [the manual](http://matsen.github.io/pplacer/generated_rst/guppy_edpl.html#guppy-edpl) describes it best, and it says this with EDPL:

        The expected distance between placement locations (EDPL) is a means of understanding the uncertainty of a placement using placer. The motivation for using such a metric comes from when there are a number of closely-related sequences present in the reference alignment. In this case, there may be considerable uncertainty about which edge is best as measured by posterior probability or likelihood weight ratio. However, the actual uncertainty as to the best region of the tree for that query sequence may be quite small. For instance, we may have a number of very similar subspecies of a given species in the alignment, and although it may not be possible to be sure to match a given query to a subspecies, one might be quite sure that it is one of them.

        The EDPL metric is one way of resolving this problem by considering the distances between the possible placements for a given query. It works as follows. Say the query bounces around to the different placement positions according to their posterior probability; i.e. the query lands with location one with probability p_1, location two with probability p_2, and so on. Then the EDPL value is simply the expected distance it will travel in one of those bounces (if you don’t like probabilistic language, it’s simply the average distance it will travel per bounce when allowed to bounce between the placements for a long time with their assigned probabilities). Here’s an example, with three hypothetical locations for a given query sequence:

    The [pplacer paper]() has even more to say. It will discuss the use of place vis which is a tool that will most likely make its way into our work also.

        Quantifying uncertainty in placement location

        Pplacer calculates edge uncertainty via posterior probability and the likelihood weight ratio. These methods quantify uncertainty on an edge-by-edge basis by comparing the best placement locations on each edge. Such quantities form the basis of an understanding of placement uncertainty.

        The Expected Distance between Placement Locations (EDPL) is used to overcome difficulties in distinguishing between local and global uncertainty, which is a complication of relying on confidence scores determined on an edge-by-edge basis. This quantity is computed as follows for a given query sequence. Pplacer first determines the top-scoring collection of edges; the optimal placement on each edge is assigned a probability defining confidence, which is the likelihood weight ratio (in ML mode) or the posterior probability (in Bayesian mode). The EDPL uncertainty is the weighted-average distance between those placements (Figure 4), i.e. the sum of the distances between the optimal placements weighted by their probability (4). The EDPL thus uses distances on the tree to distinguish between cases where nearby edges appear equally good, versus cases when a given query sequence does not have a clear position in the tree. These measures of uncertainty can then be viewed with placeviz as described below.
    """
    placed_edpl_file = os.path.join(work_dir, "edpl.txt")
    _ = run_process(
        [
            "guppy",
            "edpl",
            "--csv",
            "-o",
            placed_edpl_file,
            jplace_file,
        ],
        logger,
        capture_stdout=True,
    )
    edpl = pd.read_csv(placed_edpl_file, names=["gene", "edpl"])
    return edpl
    # = phy.read(placed_edpl_file, format="newick")
    # Combine gene mapping and color mapping into an omni mapping file
    # return t


"""
from this, I want to get:
    names for each of the clades in the figure
    labels for each known clade, unknown clade
    distance to the next node
"""


def apply_labels(clade: Clade, label: str):
    clade.label = label
    for i in clade:
        apply_labels(i, label)


def breath_search(clade: Clade):
    clades = [i for i in clade]
    while len(clades) > 0:
        for i in clades:
            yield i
        clades = [j for i in clades for j in i]


def apply_names(clade: Clade, name: str = "multiple"):
    clade.name = name
    names_num = {}
    for i, c in enumerate(breath_search(clade)):
        if c.name is not None or c.label is None:
            continue
        if c.label.startswith(name):
            c.name = f"{name}-{i}"
        else:
            name = c.label
            if name not in names_num:
                names_num[name] = 0
            else:
                names_num[name] += 1
            apply_names(c, f"{name}-{names_num[name]}")


def all_parents(tree: Tree):
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            parents[child] = clade
    return parents


PathNode = namedtuple("PathNode", ["end", "len"])


def find_distances(
    clade: Clade, known_terminals: set[str], parents: dict, skip_label: str = str(None)
) -> PathNode:
    past = clade
    present = PathNode(parents[clade], clade.branch_length)
    min_dist = float("inf")
    nearest = PathNode(None, None)
    while present.len < min_dist:
        paths = [
            PathNode(i, i.branch_length + present.len) for i in present.end if i != past
        ]
        while len(paths) > 0:
            for i in paths:
                if i.end in known_terminals and i.len < min_dist:
                    if skip_label is not None and i.end.label == skip_label:
                        continue
                    min_dist = i.len
                    nearest = i
            paths = [
                PathNode(j, i.len + j.branch_length)
                for i in paths
                if i.len <= min_dist
                for j in i.end
            ]
        past = present.end
        if past not in parents:  # the only clade that has no parents is the root
            break  # if this is the root we are done
        present = PathNode(parents[past], present.len + past.branch_length)
    return nearest


def get_clade_info(clade, tree: DramTree):
    """
    label the clades for the default tree

    Uses breadth first search to find all monolithic clades, and multiple clades.
    After this all clades will have a color and all clades will have a label.

    :param clade:
    :param tree:
    :raises ValueError:
    """
    if len(childs := list(clade)) == 0:  # The clade is a terminal node
        if (
            clade.name is None
        ):  # Needs to raise error, unnamed terminals should be impossible
            raise ValueError("Terminal clade without name")
        if (
            clade.name not in tree.mapping.index
        ):  # If the clade is not in the mapping it must be one that was added
            clade.label = None
        else:
            clade.label = tree.mapping.loc[clade.name, "call"]
    else:
        for i in childs:
            get_clade_info(i, tree)
        labels: set[str] = {
            j
            for i in childs
            if i.label is not None
            for j in (
                [i.label]
                if not i.label.startswith("multiple: ")
                else i.label.strip("multiple: ").split(", ")
            )
        }
        if len(labels) == 1:
            clade.label = "".join(labels)
            for i in childs:
                if i.color is not None and i.color != (0, 0, 0):
                    clade.color = i.color
                    break
        elif len(labels) < 1:
            clade.label = None
        else:
            clade.label = f"multiple: {', '.join(labels)}"
            for i in childs:
                if i.label is not None and not i.label.startswith("multiple"):
                    apply_labels(i, i.label)


def color_paths_by_location(treeph):
    for cl in treeph.get_terminals():
        add_color = None
        for i in treeph.get_path(cl):
            if i.color is not None and i.color != (0, 0, 0):
                add_color = i.color
            i.color = add_color
        cl = add_color
    return treeph


def pull_labes_from_tree():
    """this may not be useful"""
    pass


"""
Notes:
   the to_networkx command works, but it can get complicated this may be needed later
       for example, phy.to_networkx(phy.read("color_tree_branch.xml", "newick"))
   This makes a tree that is totally unreadable
       phy.draw_ascii(treeph)
    Need to consider offering a way to output the visualizations from pplacer

# make a test set soon
tree = NXR_NAR_TREE
logger = logging.getLogger('dram_tree_log')
annotations = pd.read_csv('example_one/all_bins_combined_3217db_ACTIVE_GENES_annotations.txt', sep='\t', index_col=0)
jplace_file = 'example_place_output.jplace'
"""


def clade_info_to_series(
    clade: Clade, tree_name: str, max_len_to_label: float, min_dif_len_ratio: float
) -> pd.DataFrame:
    """
    Note that we use labeled nodes for distance, so the distance is not to the root of a clade but to its nearest endpoint in such a clade this could have unexpected consequences, but it means that we ground our choices in known genes and not in emergent behavior of clades and the labeling algorithm.

    """
    delta: float = None
    if clade.label is not None:
        label = clade.label
        place_info = f"{CLADE_INFO_PREFIX} Nearest gene is {clade.nearest[0].end.name}"
        dist = clade.nearest[0].len
    else:
        if (dist := clade.nearest[0].len) > max_len_to_label:
            label = UNPLACE_LABEL
            place_info = f"{UNPLACE_PREFIX} The distance to the nearest labeled node is {dist}, which is more than {max_len_to_label} (the max_len_to_label filter)."
        elif (
            delta := clade.nearest[0].len / clade.nearest[1].len - 1
        ) > min_dif_len_ratio:
            label = UNPLACE_LABEL
            place_info = f"{UNPLACE_PREFIX} The difference between the nearest labeled and alternatively labeled nodes, is {delta}, which is less than {max_len_to_label} (the max_len_to_label filter)."

        else:
            first = clade.nearest[0]
            second = clade.nearest[1]
            clade.label = label = first.end.label
            place_info = f"{PROXIMITY_INFO_PREFIX} the closes labeled gene was {first.end.name}, at distance {first.len}, the nearest counter label was {second.end.label} on gene {second.end.name} at distance {second.len}"
    return pd.DataFrame(
        {
            f"{tree_name}_labels": label,
            f"distance_to_nearest_label": dist,
            f"difference_to_nearest_alt_label": delta,
            f"{tree_name}_placement_info": place_info,
        },
        index=[clade.name],
    )


def make_df_of_tree(
    placed_terminals: list[Clade],
    tree_name: str,
    edpl: pd.DataFrame,
    max_len_to_label: float,
    min_dif_len_ratio: float,
) -> pd.DataFrame:
    to_searies = partial(
        clade_info_to_series,
        tree_name=tree_name,
        max_len_to_label=max_len_to_label,
        min_dif_len_ratio=min_dif_len_ratio,
    )
    data = pd.concat([to_searies(i) for i in placed_terminals], axis=0).merge(
        edpl, left_index=True, right_index=True, how="left", copy=False
    )
    return data


def write_files(
    tree_name: str,
    tree_df: pd.DataFrame,
    treeph: Tree | None,
    jplace_file: str | None,
    output_dir: Path,
    work_dir: str,
    keep_temp: bool,
):
    """
    Write all the Tree files
    ---------------_


    FIX THIS
    """
    tree_df.to_csv(os.path.join(output_dir, f"{tree_name}_tree_data.tsv"), sep="\t")
    if treeph is not None:
        phy.write(
            treeph,
            os.path.join(output_dir, f"{tree_name}_labeled_tree.xml"),
            "phyloxml",
        )
    if jplace_file is not None:
        os.rename(
            jplace_file,
            os.path.join(output_dir, f"{tree_name}.jplace"),
        )
    if keep_temp:
        os.rename(
            work_dir,
            os.path.join(output_dir, f"{tree_name}_working_dir"),
        )


def end_message(tree_df: pd.DataFrame, tree_name: str) -> str:
    info_col = f"{tree_name}_placement_info"
    full_len = len(tree_df)
    clade_len = sum(tree_df[info_col].str.startswith(CLADE_INFO_PREFIX))
    prox_len = sum(tree_df[info_col].str.startswith(PROXIMITY_INFO_PREFIX))
    fail_len = sum(tree_df[info_col].str.startswith(UNPLACE_PREFIX))
    return (
        "Run phylogenetic disambiguation complete."
        f"\nOf the {full_len} that request phylogenetic"
        f" placement, {clade_len} genes where placed"
        f" based on the clade they fell into, {prox_len}"
        f" were classified based on the relative distance"
        f" to labeled nodes. There were {fail_len} genes"
        f" that could not be placed and so remain ambiguous"
    )


# if __name__ == "__main__":
#     dram_tree_kit()
