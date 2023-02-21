"""
This program add profiles based on phylogentic trees into dram. The key to this process is pplacer which places leaves into pre-exitsting trees

NOTE pplacer uses about 1/4 of the memory when placing on a FastTree tree as compared to a RAxML tree inferred with GTRGAMMA. If your reads are short and in a fixed region, the memory used by pplacer v1.1 alpha08 (or later) scales with respect to the total number of non-gap columns in your query alignment. You can also make it use less memory (and run faster) by cutting down the size of your reference tree.
TODO switch to FastTree
"""

# import tempfile
# import matplotlib.pyplot as plt
# from io import StringIO
# from mag_annotator.utils import setup_logger
# from mag_annotator.pull_sequences import pull_sequences
# from pyvis.network import Network
# import networkx as nx
import os
from importlib.resources import files
from collections import namedtuple
from functools import partial
import logging
import pandas as pd
from tempfile import TemporaryDirectory
from dram2.tree_kit import __version__
from dram2.tree_kit.pplacer import DramTree
from dram2.utils.utils import run_process
from dram2.annotate import get_annotation_ids_by_row
from skbio import write as write_sq
from skbio import read as read_sq
from Bio import Phylo as phy
from Bio.Phylo.BaseTree import Clade
from Bio.Phylo.Newick import Tree


"""
import os
os.system('dram2 annotate --help')
os.system('ls /home/projects-wrighton-2/DRAM/development_flynn/dram2_dev/nov_30_22_trees_and_adjectives/DRAM2/dram2/tree_kit/data/')
/home/projects-wrighton-2/DRAM/development_flynn/dram2_dev/nov_30_22_trees_and_adjectives/DRAM2/dram2/tree_kit/data/nxr_nar/nxr-nar-tree-mapping.tsv
"""

data_path = files('dram2.tree_kit').joinpath('data')
NXR_NAR_TREE = DramTree(
    name="nxr_nar",
    pplacer_profile=os.path.join(data_path, "nxr_nar", "nxr_nar.refpkg"),
    target_ids=["K11180", "dsrA", "dsrB", "K11181"],
    reference_seq=os.path.join(data_path, "nxr_nar", "nxr-nar_seqs_for_tree_aligned.faa"),
    gene_mapping_path=os.path.join(data_path, "nxr_nar", "nxr-nar-tree-mapping.tsv"),
    color_mapping_path=os.path.join(data_path, "nxr_nar", "color_map.tsv")
)
TREES = [NXR_NAR_TREE]
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


def tree_kit(
    dram_annotations,
    gene_fasta: str,
    dram_directory: str,
    logger: logging.Logger,
    output_dir: str = "./",
    annotate_all: bool = False,
    keep_temp: bool = False,
    cores: int = 10,
    # logg_path: str = "phylo_tree.log",
    force: bool = False,
    max_len_to_label: float = MAX_LEN_TO_LABEL_DFLT,
    min_dif_len_ratio: float = MIN_DIF_LEN_RATIO_DFLT,
):
    # logger = logging.getLogger("dram_tree_log")
    # eetup_logger(logger, logg_path)
    tree = NXR_NAR_TREE
    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    elif not force:
        raise ValueError(
            "The output_dir already exists! try using the -f flag to overwrite"
        )
    logger.info("Start phylogenetic tree disambiguation")

    logger.info("Processing annotations")
    annotations = pd.read_csv(dram_annotations, sep="\t", index_col=0)
    annotation_ids = get_annotation_ids_by_row(annotations, logger)
    with TemporaryDirectory(dir=output_dir, prefix="tmp_") as work_dir:
        for tree in TREES:
            logger.info(f"Performing phylogenetic disambiguation with tree {tree.name}")
            tree.set_logger(logger)
            logger.info(f"Made temporary files in {work_dir}")
            trimed_fa = extract_enigmatic_genes(
                annotation_ids, gene_fasta, work_dir, tree.target_ids, logger
            )
            logger.info("Placing enigmatic genes")
            jplace_file = tree.pplacer_place_sequences(
                trimed_fa, work_dir, threads=cores
            )
            treeph = read_phtree(jplace_file, work_dir, logger)
            edpl = read_edpl(jplace_file, work_dir, logger)
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
            logger.info("Writing output, to {output_dir}")
            write_files(
                tree.name, tree_df, treeph, jplace_file, output_dir, work_dir, keep_temp
            )
            logger.info(end_message(tree_df, tree.name))


"""
import os
os.system('test_env/bin/python3 dram_tree_kit/dram_phylo_tree.py  -f -a tests/data/mini_nxr_nar_annotations.tsv -g tests/data/mini_nxr_nar_genes.faa -c 30')
dram_tree_kit()
"""


def extract_enigmatic_genes(
    annotation_ids: pd.Series,
    gene_fasta: str,
    work_dir: str,
    target_ids: set,
    logger: logging.Logger,
) -> str:
    """
    :param annotation_ids: annotations from dram run
    :param gene_fasta: faa from dram run
    :param work_dir: Temp files here
    :param target_ids: ID set needing phylo info, used to filter genes
    :param logger: Standard DRAM logger
    :returns: The path to the ambiguous genes in a fasta file

    Takes in a fasta file of genes and a list of ids in the dram annotation, and returns a filtered fasta to match.
    """
    logger.info("Finding enigmatic genes")
    output_fasta = os.path.join(work_dir, "trim.faa")
    ids_keep = annotation_ids[
        annotation_ids.apply(lambda x: len(x.intersection(target_ids)) > 0)
    ]
    output_fasta_generator = (
        i
        for i in read_sq(gene_fasta, format="fasta")
        if i.metadata["id"] in ids_keep.index
    )
    write_sq(output_fasta_generator, format="fasta", into=output_fasta)
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
    # Combine gene maping and color maping into an omni maping file
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


def color_known_termininals(treeph, annotations, mapping):
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
    There are not as many options for certainty as first thought but is in any case I will for now go with the EDPL.

    So the deal with uncertainty is that it may only mean something if the distance to the contrary node is smaller than some ratio of the edpl. Of cores, it is not clear if the distance is would even be distributed equally along all paths.

    In any case, there is no reason to think that these values are not a measure of uncertainty, and the user will probability have use for them but it may mean that in the future we provide additional pplacer output.

    The best solutution and one that we could cirtanly acheave, is to use pplacer to place and then


    I think [the manual](http://matsen.github.io/pplacer/generated_rst/guppy_edpl.html#guppy-edpl) describes it best, and it says this with EDPL:

        The expected distance between placement locations (EDPL) is a means of understanding the uncertainty of a placement using placer. The motivation for using such a metric comes from when there are a number of closely-related sequences present in the reference alignment. In this case, there may be considerable uncertainty about which edge is best as measured by posterior probability or likelihood weight ratio. However, the actual uncertainty as to the best region of the tree for that query sequence may be quite small. For instance, we may have a number of very similar subspecies of a given species in the alignment, and although it may not be possible to be sure to match a given query to a subspecies, one might be quite sure that it is one of them.

        The EDPL metric is one way of resolving this problem by considering the distances between the possible placements for a given query. It works as follows. Say the query bounces around to the different placement positions according to their posterior probability; i.e. the query lands with location one with probability p_1, location two with probability p_2, and so on. Then the EDPL value is simply the expected distance it will travel in one of those bounces (if you don’t like probabilistic language, it’s simply the average distance it will travel per bounce when allowed to bounce between the placements for a long time with their assigned probabilities). Here’s an example, with three hypothetical locations for a given query sequence:

    The [pplacer paper]() has even more to say. It will discus the use of place vis which is a tool that will most likely make its whay into our work also.

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
    # Combine gene maping and color maping into an omni maping file
    # return t


"""
from this, I want to get:
    names for each of the clades in the figure
    labels for each known clade, unkown clade
    distace to the next node
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

    Uses breath first search to find all monolithic clades, and multiple clade.
    After this all clades will have a color and all clade will have a label.

    :param clade:
    :param tree:
    :raises ValueError:
    """
    if len(childs := list(clade)) == 0:  # The clade is a terminal node
        if (
            clade.name is None
        ):  # Needs to raise error, unamed terminals should be imposible
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
    """this may not be usefull"""
    pass


"""
Notes:
   the to_networkx comand works but it can get complicated this may be needed later
       for example phy.to_networkx(phy.read("color_tree_branch.xml", "newick"))
   This makes a tree that is totaly unreadable
       phy.draw_ascii(treeph)
    Need to condider offering a way to output the visualizations from pplacer

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
    Note that we use labeled nodes for distace, so the distace is not to the root of a clade but to its nearest endpoint in such a clade this could have un-expected conciquences but it means that we ground our choices in known genes and not in emergent behavure of clades and the lableing algorithm.

    """
    delta: float = None
    if clade.label is not None:
        label = clade.label
        place_info = f"{CLADE_INFO_PREFIX} Nearest gene is {clade.nearest[0].end.name}"
        dist = clade.nearest[0].len
    else:
        if (dist := clade.nearest[0].len) > max_len_to_label:
            label = UNPLACE_LABEL
            place_info = f"{UNPLACE_PREFIX} The distance to the nearest labeled node, is {dist}, which is more than {max_len_to_label} (the max_len_to_label filter)."
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
    treeph: Tree,
    jplace_file: str,
    output_dir: str,
    work_dir: str,
    keep_temp: bool,
):
    tree_df.to_csv(os.path.join(output_dir, f"{tree_name}_tree_data.tsv"), sep="\t")
    phy.write(
        treeph,
        os.path.join(output_dir, f"{tree_name}_labeled_tree.xml"),
        "phyloxml",
    )
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


if __name__ == "__main__":
    dram_tree_kit()
