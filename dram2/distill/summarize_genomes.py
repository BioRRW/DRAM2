"""
DRAM Distillate
---------------

This is the script that distills the genomes from the annotations step.
It requires annotations and the dram distillation step

"""
import logging
from collections import Counter
from itertools import chain
from typing import Optional

import pandas as pd
from collections import Counter, defaultdict
from pathlib import Path
import altair as alt
import networkx as nx
from itertools import tee
import re
import numpy as np
from dram2.utils.globals import FASTAS_CONF_TAG

from dram2.db_kits.utils import DBKit, FastaKit, HmmKit, DRAM_DATAFOLDER_TAG
from dram2.utils.utils import (
    get_ordered_uniques,
    DramUsageError,
    Fasta,
    get_package_path,
)

# from dram2.annotate import  get_ids_from_annotations_all
from dram2.annotate import (
    get_annotation_ids_by_row,
    get_all_annotation_ids,
    DB_KITS,
    get_past_annotation_run,
    ANNOTATION_FILE_TAG,
    check_for_annotations,
    USED_DBS_TAG,
    DISTILLATION_MIN_SET,
    DISTILLATION_MIN_SET_KEGG,
)

__version__ = "2.0.0"
DEFAULT_GROUPBY_COLUMN = "fasta"
DEFAULT_SHOW_DISTILLATE_GENE_NAMES = False
DEFAULT_GENOMES_PER_PRODUCT = 1000
GENOMES_PRODUCT_LIMIT = 2000
DEFAULT_MAKE_BIG_HTML = False
SUMMARIZE_METABOLISM_TAG: str = "summarize_metabolism"
MAKE_GENOME_STATS_TAG: str = "make_genome_stats"
MAKE_PRODUCT_TAG: str = "make_product"
PRODUCT_TAG = "product"
DISTILLATE_MODULES: tuple[str, str, str] = (
    SUMMARIZE_METABOLISM_TAG,
    MAKE_GENOME_STATS_TAG,
    MAKE_PRODUCT_TAG,
)
DRAM_SHEET_TAG = "dram_sheets"
GENOME_SUMMARY_FORM_TAG = "genome_summary_form"
MODULE_STEPS_FORM_TAG = "module_step_form"
FUNCTION_HEATMAP_FORM_TAG = "function_heatmap_form"
ETC_MODULE_DF_TAG = "etc_module_database"
FORM_TAGS: set = {
    GENOME_SUMMARY_FORM_TAG,
    MODULE_STEPS_FORM_TAG,
    FUNCTION_HEATMAP_FORM_TAG,
    ETC_MODULE_DF_TAG,
}

FILES_NAMES: dict[str, Path] = {
    GENOME_SUMMARY_FORM_TAG: Path("distill", "data", "genome_summary_form.tsv"),
    MODULE_STEPS_FORM_TAG: Path("distill", "data", "module_step_form.tsv"),
    FUNCTION_HEATMAP_FORM_TAG: Path("distill", "data", "function_heatmap_form.tsv"),
    ETC_MODULE_DF_TAG: Path("distill", "data", "etc_module_database.tsv"),
}

# DB_KITS = [i for i in DBKit.__subclasses__() if i.selectable and i.has_genome_summary]


# TODO: add RBH information to output
# TODO: add flag to output table and not xlsx

FRAME_COLUMNS = [
    "gene_id",
    "gene_description",
    "module",
    "sheet",
    "header",
    "subheader",
]
RRNA_TYPES = ["5S rRNA", "16S rRNA", "23S rRNA"]
HEATMAP_MODULES = [
    "M00001",
    "M00004",
    "M00008",
    "M00009",
    "M00012",
    "M00165",
    "M00173",
    "M00374",
    "M00375",
    "M00376",
    "M00377",
    "M00422",
    "M00567",
]
HEATMAP_CELL_HEIGHT = 10
HEATMAP_CELL_WIDTH = 10
KO_REGEX = r"^K\d\d\d\d\d$"
ETC_COVERAGE_COLUMNS = [
    "module_id",
    "module_name",
    "complex",
    "genome",
    "path_length",
    "path_length_coverage",
    "percent_coverage",
    "genes",
    "missing_genes",
    "complex_module_name",
]
TAXONOMY_LEVELS = ["d", "p", "c", "o", "f", "g", "s"]
COL_HEADER, COL_SUBHEADER, COL_MODULE, COL_GENE_ID, COL_GENE_DESCRIPTION = (
    "header",
    "subheader",
    "module",
    "gene_id",
    "gene_description",
)
CONSTANT_DISTILLATE_COLUMNS = [
    COL_GENE_ID,
    COL_GENE_DESCRIPTION,
    COL_MODULE,
    COL_HEADER,
    COL_SUBHEADER,
]
DISTILATE_SORT_ORDER_COLUMNS = [COL_HEADER, COL_SUBHEADER, COL_MODULE, COL_GENE_ID]
EXCEL_MAX_CELL_SIZE = 32767
LOCATION_TAG = "location"


def get_distillate_sheet(form_tag: str, dram_config: dict, logger: logging.Logger):
    """
    Paths in the config can be complicated. Here is a function that will get
    you the absolute path, the relative path, or whatever. Specifically for
    distillate sheets This should be more
    formalized and the config

    should actually be managed in its own structure. With a data file class
    that can use this function.
    """
    if (
        (dram_sheets := dram_config.get(DRAM_SHEET_TAG)) is None
        or dram_sheets.get(form_tag) is None
        or (sheet_path_str := dram_sheets[form_tag].get(LOCATION_TAG)) is None
    ):
        sheet_path: Path = get_package_path(Path(FILES_NAMES[form_tag]))
        logger.debug(
            f"Using the default distillation sheet for {form_tag} with its location at {sheet_path}. This information is only important if you intended to use a custom distillation sheet"
        )
    else:
        sheet_path = Path(sheet_path_str)
        if not sheet_path.is_absolute():
            dram_data_folder: Optional[str] = config.get(DRAM_DATAFOLDER_TAG)
            if dram_data_folder is None:
                raise DramUsageError(
                    f"In the DRAM2 config File, the path {form_tag} is a relative path and the {DRAM_DATAFOLDER_TAG} is not set!"
                )
            sheet_path = Path(dram_data_folder) / sheet_path
        if not sheet_path.exists():
            raise DramUsageError(
                f"The file {form_tag} is not at the path"
                f" {sheet_path}. Most likely you moved the DRAM"
                f" data but forgot to update the config file to"
                f" point to it. The easy fix is to set the"
                f" {DRAM_DATAFOLDER_TAG} variable in the config"
                f" like:\n"
                f" {DRAM_DATAFOLDER_TAG}: the/path/to/my/file"
                f" If you are using full paths and not the"
                f" {DRAM_DATAFOLDER_TAG} you may want to revue the"
                f" Configure Dram section of the documentation to"
                f" make sure your config will work with dram."
                f" rememberer that the config must be a valid yaml"
                f" file to work. Also you can always use"
                f" db_builder to remake your databases and the"
                f" config file if you don't feel up to editing it"
                f" yourself."
            )

    return pd.read_csv(sheet_path, sep="\t")


def fill_genome_summary_frame(
    annotation_ids_by_row,
    genome_summary_frame,
    groupby_column,
    logger: logging.Logger,
):
    genome_summary_id_sets = [
        set([k.strip() for k in j.split(",")]) for j in genome_summary_frame["gene_id"]
    ]
    for genome, frame in annotation_ids_by_row.groupby(groupby_column, sort=False):
        id_dict = get_all_annotation_ids(frame)
        counts = list()
        for i in genome_summary_id_sets:
            identifier_count = 0
            for j in i:
                if j in id_dict:
                    identifier_count += id_dict[j]
            counts.append(identifier_count)
        genome_summary_frame[genome] = counts
    return genome_summary_frame


def fill_genome_summary_frame_gene_names(
    annotation_ids_by_row, genome_summary_frame, groupby_column, logger
):
    genome_summary_id_sets = [
        set([k.strip() for k in j.split(",")]) for j in genome_summary_frame["gene_id"]
    ]
    for genome, frame in annotation_ids_by_row.groupby(groupby_column, sort=False):
        # make dict of identifiers to gene names
        id_gene_dict = defaultdict(list)
        for gene, ids in frame.items():
            for id_ in ids:
                id_gene_dict[id_].append(gene)
        # fill in genome summary_frame
        values = list()
        for id_set in genome_summary_id_sets:
            this_value = list()
            for id_ in id_set:
                this_value += id_gene_dict[id_]
            values.append(",".join(this_value))
        genome_summary_frame[genome] = values
    return genome_summary_frame


def summarize_rrnas(rrnas_df, groupby_column=DEFAULT_GROUPBY_COLUMN):
    genome_rrna_dict = dict()
    for genome, frame in rrnas_df.groupby(groupby_column):
        genome_rrna_dict[genome] = Counter(frame["type"])
    row_list = list()
    for rna_type in RRNA_TYPES:
        row = [
            rna_type,
            "%s ribosomal RNA gene" % rna_type.split()[0],
            "rRNA",
            "rRNA",
        ]
        for genome, rrna_dict in genome_rrna_dict.items():
            row.append(genome_rrna_dict[genome].get(rna_type, 0))
        row_list.append(row)
    rrna_frame = pd.DataFrame(
        row_list, columns=FRAME_COLUMNS + list(genome_rrna_dict.keys())
    )
    return rrna_frame


def summarize_trnas(trnas_df, groupby_column=DEFAULT_GROUPBY_COLUMN):
    # first build the frame
    combos = {(line.Type, line.Codon, line.Note) for _, line in trnas_df.iterrows()}
    frame_rows = list()
    for combo in combos:
        if combo[2] == "pseudo":
            gene_id = "%s, pseudo (%s)"
            gene_description = "%s pseudo tRNA with %s Codon"
        else:
            gene_id = "%s (%s)"
            gene_description = "%s tRNA with %s Codon"
        gene_id = gene_id % (combo[0], combo[1])
        gene_description = gene_description % (combo[0], combo[1])
        module_description = "%s tRNA" % combo[0]
        frame_rows.append(
            [gene_id, gene_description, module_description, "tRNA", "tRNA", ""]
        )
    trna_frame = pd.DataFrame(frame_rows, columns=FRAME_COLUMNS)
    trna_frame = trna_frame.sort_values("gene_id")
    # then fill it in
    trna_frame = trna_frame.set_index("gene_id")
    for group, frame in trnas_df.groupby(groupby_column):
        gene_ids = list()
        for _, line in frame.iterrows():
            if line.Note == "pseudo":
                gene_id = "%s, pseudo (%s)"
            else:
                gene_id = "%s (%s)"
            gene_ids.append(gene_id % (line.Type, line.Codon))
        trna_frame[group] = pd.Series(Counter(gene_ids))
    trna_frame = trna_frame.reset_index()
    trna_frame = trna_frame.fillna(0)
    return trna_frame


def make_genome_summary(
    genome_summary_frame: pd.DataFrame,
    annotation_ids_by_row: pd.DataFrame,
    logger: logging.Logger,
    trna_frame,
    rrna_frame,
    groupby_column: str,
):
    summary_frames = list()
    # get ko summaries
    summary_frames.append(
        fill_genome_summary_frame(
            annotation_ids_by_row,
            genome_summary_frame.copy(),
            groupby_column,
            logger,
        )
    )

    # add rRNAs
    if rrna_frame is not None:
        summary_frames.append(summarize_rrnas(rrna_frame, groupby_column))

    # add tRNAs
    if trna_frame is not None:
        summary_frames.append(summarize_trnas(trna_frame, groupby_column))

    # merge summary frames
    summarized_genomes = pd.concat(summary_frames, sort=False)
    return summarized_genomes


def split_column_str(names):
    if len(names) < EXCEL_MAX_CELL_SIZE:
        return [names]
    out = [""]
    name_list = names.split(",")
    j = 0
    for i in name_list:
        if len(out[j]) + len(i) + 1 < EXCEL_MAX_CELL_SIZE:
            out[j] = ",".join([out[j], i])
        else:
            j += 1
            out += [""]
    return out


def split_names_to_long(col: pd.Series):
    dex = col.index
    splits = [split_column_str(i) for i in col.values]
    ncols = max([len(i) for i in splits])
    col_names = [col.name if i == 0 else f"{col.name}[{i + 1}]" for i in range(ncols)]
    return pd.DataFrame(splits, columns=col_names, index=dex).fillna("")


def write_summarized_genomes_to_xlsx(summarized_genomes, output_file):
    # turn all this into an xlsx
    with pd.ExcelWriter(output_file) as writer:
        for sheet, frame in summarized_genomes.groupby("sheet", sort=False):
            frame = frame.sort_values(DISTILATE_SORT_ORDER_COLUMNS)
            frame = frame.drop(["sheet"], axis=1)
            gene_columns = list(set(frame.columns) - set(CONSTANT_DISTILLATE_COLUMNS))
            split_genes = pd.concat(
                [split_names_to_long(frame[i].astype(str)) for i in gene_columns],
                axis=1,
            )
            frame = pd.concat([frame[CONSTANT_DISTILLATE_COLUMNS], split_genes], axis=1)
            frame.to_excel(writer, sheet_name=sheet, index=False)


# TODO: add assembly stats like N50, longest contig, total assembled length etc
# TODO: add in assembly stats
def make_genome_stats(
    annotations: pd.DataFrame,
    rrna_frame: Optional[pd.DataFrame] = None,
    trna_frame: Optional[pd.DataFrame] = None,
    groupby_column: str = DEFAULT_GROUPBY_COLUMN,
):
    rows = list()
    columns = ["genome"]
    if "scaffold" in annotations.columns:
        columns.append("number of scaffolds")
    if "bin_taxonomy" in annotations.columns:
        columns.append("taxonomy")
    if "bin_completeness" in annotations.columns:
        columns.append("completeness score")
    if "bin_contamination" in annotations.columns:
        columns.append("contamination score")
    if rrna_frame is not None:
        columns += RRNA_TYPES
    if trna_frame is not None:
        columns.append("tRNA count")
    if (
        "bin_completeness" in annotations.columns
        and "bin_contamination" in annotations.columns
        and rrna_frame is not None
        and trna_frame is not None
    ):
        columns.append("assembly quality")
    for genome, frame in annotations.groupby(groupby_column, sort=False):
        row = [genome]
        if "scaffold" in frame.columns:
            row.append(len(set(frame["scaffold"])))
        if "bin_taxonomy" in frame.columns:
            row.append(frame["bin_taxonomy"][0])
        if "bin_completeness" in frame.columns:
            row.append(frame["bin_completeness"][0])
        if "bin_contamination" in frame.columns:
            row.append(frame["bin_contamination"][0])
        has_rrna = list()
        if rrna_frame is not None:
            genome_rrnas = rrna_frame.loc[rrna_frame.fasta == genome]
            for rrna in RRNA_TYPES:
                sixteens = genome_rrnas.loc[genome_rrnas.type == rrna]
                if sixteens.shape[0] == 0:
                    row.append("")
                    has_rrna.append(False)
                elif sixteens.shape[0] == 1:
                    row.append(
                        "%s (%s, %s)"
                        % (
                            sixteens["scaffold"].iloc[0],
                            sixteens.begin.iloc[0],
                            sixteens.end.iloc[0],
                        )
                    )
                    has_rrna.append(True)
                else:
                    row.append("%s present" % sixteens.shape[0])
                    has_rrna.append(False)
        if trna_frame is not None:
            # TODO: remove pseudo from count?
            row.append(trna_frame.loc[trna_frame[groupby_column] == genome].shape[0])
        if "assembly quality" in columns and trna_frame is not None:
            if (
                frame["bin_completeness"][0] > 90
                and frame["bin_contamination"][0] < 5
                and np.all(has_rrna)
                and len(set(trna_frame.loc[trna_frame[groupby_column] == genome].Type))
                >= 18
            ):
                row.append("high")
            elif (
                frame["bin_completeness"][0] >= 50
                and frame["bin_contamination"][0] < 10
            ):
                row.append("med")
            else:
                row.append("low")
        rows.append(row)
    genome_stats = pd.DataFrame(rows, columns=columns)
    return genome_stats


def build_module_net(module_df):
    """Starts with a data from including a single module"""
    # build net from a set of module paths
    num_steps = max([int(i.split(",")[0]) for i in set(module_df["path"])])
    module_net = nx.DiGraph(
        num_steps=num_steps,
        module_id=list(module_df["module"])[0],
        module_name=list(module_df["module_name"])[0],
    )
    # go through all path/step combinations
    for module_path, frame in module_df.groupby("path"):
        split_path = [int(i) for i in module_path.split(",")]
        step = split_path[0]
        module_net.add_node(module_path, kos=set(frame["ko"]))
        # add incoming edge
        if step != 0:
            module_net.add_edge("end_step_%s" % (step - 1), module_path)
        # add outgoing edge
        module_net.add_edge(module_path, "end_step_%s" % step)
    return module_net


def get_module_step_coverage(kos, module_net):
    # prune network based on what kos were observed
    pruned_module_net = module_net.copy()
    module_kos_present = set()
    for node, data in module_net.nodes.items():
        if "kos" in data:
            ko_overlap = data["kos"] & kos
            if len(ko_overlap) == 0:
                pruned_module_net.remove_node(node)
            else:
                module_kos_present = module_kos_present | ko_overlap
    # count number of missing steps
    missing_steps = 0
    for node, data in pruned_module_net.nodes.items():
        if ("end_step" in node) and (pruned_module_net.in_degree(node) == 0):
            missing_steps += 1
    # get statistics
    num_steps = pruned_module_net.graph["num_steps"] + 1
    num_steps_present = num_steps - missing_steps
    coverage = num_steps_present / num_steps
    return num_steps, num_steps_present, coverage, sorted(module_kos_present)


def make_module_coverage_df(annotation_df, module_nets):
    kos_to_genes = defaultdict(list)
    ko_id: Optional(str) = None
    ko_id_names: list[str] = ["kegg_id", "kofam_id", "ko_id"]
    for id in ko_id_names:
        if id in annotation_df:
            ko_id = id
            break
    if ko_id is None:
        raise ValueError(
            f"No KEGG or KOfam id column could be found. These names were tried: {', '.join(ko_id_names)}"
        )
    for gene_id, ko_list in annotation_df[ko_id].items():
        if type(ko_list) is str:
            for ko in ko_list.split(","):
                kos_to_genes[ko].append(gene_id)
    coverage_dict = {}
    for i, (module, net) in enumerate(module_nets.items()):
        (
            module_steps,
            module_steps_present,
            module_coverage,
            module_kos,
        ) = get_module_step_coverage(set(kos_to_genes.keys()), net)
        module_genes = sorted([gene for ko in module_kos for gene in kos_to_genes[ko]])
        coverage_dict[module] = [
            net.graph["module_name"],
            module_steps,
            module_steps_present,
            module_coverage,
            len(module_kos),
            ",".join(module_kos),
            ",".join(module_genes),
        ]
    coverage_df = pd.DataFrame.from_dict(
        coverage_dict,
        orient="index",
        columns=[
            "module_name",
            "steps",
            "steps_present",
            "step_coverage",
            "ko_count",
            "kos_present",
            "genes_present",
        ],
    )
    return coverage_df


def make_module_coverage_frame(
    annotations, module_nets, groupby_column=DEFAULT_GROUPBY_COLUMN
):
    # go through each scaffold to check for modules
    module_coverage_dict = dict()
    for group, frame in annotations.groupby(groupby_column, sort=False):
        module_coverage_dict[group] = make_module_coverage_df(frame, module_nets)
    module_coverage = pd.concat(module_coverage_dict)
    module_coverage.index = module_coverage.index.set_names(["genome", "module"])
    return module_coverage.reset_index()


def make_module_coverage_heatmap(module_coverage, mag_order=None):
    num_mags_in_frame = len(set(module_coverage["genome"]))
    c = (
        alt.Chart(module_coverage, title="Module")
        .encode(
            x=alt.X(
                "module_name",
                title=None,
                sort=None,
                axis=alt.Axis(labelLimit=0, labelAngle=90),
            ),
            y=alt.Y("genome", title=None, sort=mag_order, axis=alt.Axis(labelLimit=0)),
            tooltip=[
                alt.Tooltip("genome", title="Genome"),
                alt.Tooltip("module_name", title="Module Name"),
                alt.Tooltip("steps", title="Module steps"),
                alt.Tooltip("steps_present", title="Steps present"),
            ],
        )
        .mark_rect()
        .encode(
            color=alt.Color(
                "step_coverage",
                legend=alt.Legend(title="% Complete"),
                scale=alt.Scale(domain=(0, 1)),
            )
        )
        .properties(
            width=HEATMAP_CELL_WIDTH * len(HEATMAP_MODULES),
            height=HEATMAP_CELL_HEIGHT * num_mags_in_frame,
        )
    )
    return c


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


def first_open_paren_is_all(str_):
    """Go through string and return true"""
    curr_level = 1
    for char in str_[1:-1]:
        if char == ")":
            curr_level -= 1
        elif char == "(":
            curr_level += 1
        if curr_level == 0:
            return False
    return True


def split_into_steps(definition, split_char=" "):
    """Very fancy split on string of chars"""
    curr_level = 0
    step_starts = [-1]
    for i, char in enumerate(definition):
        if char == "(":
            curr_level += 1
        if char == ")":
            curr_level -= 1
        if (curr_level == 0) and (char in split_char):
            step_starts.append(i)
    step_starts.append(len(definition))
    steps = list()
    for a, b in pairwise(step_starts):
        step = definition[a + 1 : b]
        if step.startswith("(") and step.endswith(")"):
            if first_open_paren_is_all(step):
                step = step[1:-1]
        steps.append(step)
    return steps


def is_ko(ko):
    return re.match(KO_REGEX, ko) is not None


def make_module_network(
    definition, network: nx.DiGraph = None, parent_nodes=("start",)
):
    # TODO: Figure out how to add 'end' step to last step at end
    if network is None:
        network = nx.DiGraph()
    last_steps = []
    for step in split_into_steps(definition, ","):
        prev_steps = parent_nodes
        for substep in split_into_steps(step, "+"):
            if is_ko(substep):
                for prev_step in prev_steps:
                    network.add_edge(prev_step, substep)
                prev_steps = [substep]
            else:
                network, prev_steps = make_module_network(substep, network, prev_steps)
        last_steps += prev_steps
    return network, last_steps


def get_module_coverage(module_net: nx.DiGraph, genes_present: set):
    max_coverage = -1
    max_coverage_genes = list()
    max_coverage_missing_genes = list()
    max_path_len = 0
    for net_path in nx.all_simple_paths(module_net, source="start", target="end"):
        net_path = set(net_path[1:-1])
        overlap = net_path & genes_present
        coverage = len(overlap) / len(net_path)
        if coverage > max_coverage:
            max_coverage = coverage
            max_coverage_genes = overlap
            max_coverage_missing_genes = net_path - genes_present
            max_path_len = len(net_path)
    return (
        max_path_len,
        len(max_coverage_genes),
        max_coverage,
        max_coverage_genes,
        max_coverage_missing_genes,
    )


def make_etc_coverage_df(
    etc_module_df,
    annotation_ids_by_row: pd.DataFrame,
    logger: logging.Logger,
    groupby_column=DEFAULT_GROUPBY_COLUMN,
):
    etc_coverage_df_rows = list()
    for _, module_row in etc_module_df.iterrows():
        definition = module_row["definition"]
        # remove optional subunits
        definition = re.sub(r"-K\d\d\d\d\d", "", definition)
        module_net, _ = make_module_network(definition)
        # add end node
        no_out = [
            node for node in module_net.nodes() if module_net.out_degree(node) == 0
        ]
        for node in no_out:
            module_net.add_edge(node, "end")
        # go through each genome and check pathway coverage
        for group, frame in annotation_ids_by_row.groupby(groupby_column):
            # get annotation genes
            grouped_ids = set(get_all_annotation_ids(frame).keys())
            (
                path_len,
                path_coverage_count,
                path_coverage_percent,
                genes,
                missing_genes,
            ) = get_module_coverage(module_net, grouped_ids)
            complex_module_name = "Complex %s: %s" % (
                module_row["complex"].replace("Complex ", ""),
                module_row["module_name"],
            )
            etc_coverage_df_rows.append(
                [
                    module_row["module_id"],
                    module_row["module_name"],
                    module_row["complex"].replace("Complex ", ""),
                    group,
                    path_len,
                    path_coverage_count,
                    path_coverage_percent,
                    ",".join(sorted(genes)),
                    ",".join(sorted(missing_genes)),
                    complex_module_name,
                ]
            )
    return pd.DataFrame(etc_coverage_df_rows, columns=ETC_COVERAGE_COLUMNS)


def make_etc_coverage_heatmap(etc_coverage, mag_order=None, module_order=None):
    num_mags_in_frame = len(set(etc_coverage["genome"]))
    charts = list()
    for i, (etc_complex, frame) in enumerate(etc_coverage.groupby("complex")):
        # if this is the first chart then make y-ticks otherwise none
        c = (
            alt.Chart(frame, title=etc_complex)
            .encode(
                x=alt.X(
                    "module_name",
                    title=None,
                    axis=alt.Axis(labelLimit=0, labelAngle=90),
                    sort=module_order,
                ),
                y=alt.Y(
                    "genome",
                    axis=alt.Axis(title=None, labels=False, ticks=False),
                    sort=mag_order,
                ),
                tooltip=[
                    alt.Tooltip("genome", title="Genome"),
                    alt.Tooltip("module_name", title="Module Name"),
                    alt.Tooltip("path_length", title="Module Subunits"),
                    alt.Tooltip("path_length_coverage", title="Subunits present"),
                    alt.Tooltip("genes", title="Genes present"),
                    alt.Tooltip("missing_genes", title="Genes missing"),
                ],
            )
            .mark_rect()
            .encode(
                color=alt.Color(
                    "percent_coverage",
                    legend=alt.Legend(title="% Complete"),
                    scale=alt.Scale(domain=(0, 1)),
                )
            )
            .properties(
                width=HEATMAP_CELL_WIDTH * len(set(frame["module_name"])),
                height=HEATMAP_CELL_HEIGHT * num_mags_in_frame,
            )
        )
        charts.append(c)
    concat_title = alt.TitleParams("ETC Complexes", anchor="middle")
    return alt.hconcat(*charts, spacing=5, title=concat_title)


def make_functional_df(
    annotations,
    function_heatmap_form,
    logger,
    annotation_ids_by_row: pd.DataFrame,
    groupby_column=DEFAULT_GROUPBY_COLUMN,
):
    # clean up function heatmap form
    function_heatmap_form = function_heatmap_form.apply(
        lambda x: x.str.strip() if x.dtype == "object" else x
    )
    function_heatmap_form = function_heatmap_form.fillna("")
    # build dict of ids per genome
    genome_to_id_dict = dict()
    for genome, frame in annotation_ids_by_row.groupby(groupby_column, sort=False):
        id_list = get_all_annotation_ids(frame).keys()
        genome_to_id_dict[genome] = set(id_list)
    # build long from data frame
    rows = list()
    for function, frame in function_heatmap_form.groupby("function_name", sort=False):
        for bin_name, id_set in genome_to_id_dict.items():
            presents_in_bin = list()
            functions_present = set()
            for _, row in frame.iterrows():
                function_id_set = set(
                    [i.strip() for i in row.function_ids.strip().split(",")]
                )
                present_in_bin = id_set & function_id_set
                functions_present = functions_present | present_in_bin
                presents_in_bin.append(len(present_in_bin) > 0)
            function_in_bin = np.all(presents_in_bin)
            row = frame.iloc[0]
            rows.append(
                [
                    row.category,
                    row.subcategory,
                    row.function_name,
                    ", ".join(functions_present),
                    "; ".join(get_ordered_uniques(frame.long_function_name)),
                    "; ".join(get_ordered_uniques(frame.gene_symbol)),
                    bin_name,
                    function_in_bin,
                    "%s: %s" % (row.category, row.function_name),
                ]
            )
    return pd.DataFrame(
        rows,
        columns=list(function_heatmap_form.columns)
        + ["genome", "present", "category_function_name"],
    )


def make_functional_heatmap(functional_df, mag_order=None):
    # build heatmaps
    charts = list()
    for i, (group, frame) in enumerate(functional_df.groupby("category", sort=False)):
        # set variables for chart
        function_order = get_ordered_uniques(list(frame.function_name))
        num_mags_in_frame = len(set(frame["genome"]))
        chart_width = HEATMAP_CELL_WIDTH * len(function_order)
        chart_height = HEATMAP_CELL_HEIGHT * num_mags_in_frame
        # if this is the first chart then make y-ticks otherwise none
        if i == 0:
            y = alt.Y(
                "genome",
                title=None,
                sort=mag_order,
                axis=alt.Axis(
                    labelLimit=0, labelExpr="replace(datum.label, /_\d*$/gi, '')"
                ),
            )
        else:
            y = alt.Y(
                "genome",
                axis=alt.Axis(title=None, labels=False, ticks=False),
                sort=mag_order,
            )
        # set up colors for chart
        rect_colors = alt.Color(
            "present",
            legend=alt.Legend(
                title="Function is Present", symbolType="square", values=[True, False]
            ),
            sort=[True, False],
            scale=alt.Scale(range=["#2ca25f", "#e5f5f9"]),
        )
        # define chart
        # TODO: Figure out how to angle title to take up less space
        c = (
            alt.Chart(frame, title=alt.TitleParams(group))
            .encode(
                x=alt.X(
                    "function_name",
                    title=None,
                    axis=alt.Axis(labelLimit=0, labelAngle=90),
                    sort=function_order,
                ),
                tooltip=[
                    alt.Tooltip("genome", title="Genome"),
                    alt.Tooltip("category", title="Category"),
                    alt.Tooltip("subcategory", title="Subcategory"),
                    alt.Tooltip("function_ids", title="Function IDs"),
                    alt.Tooltip("function_name", title="Function"),
                    alt.Tooltip("long_function_name", title="Description"),
                    alt.Tooltip("gene_symbol", title="Gene Symbol"),
                ],
            )
            .mark_rect()
            .encode(y=y, color=rect_colors)
            .properties(width=chart_width, height=chart_height)
        )
        charts.append(c)
    # merge and return
    function_heatmap = alt.hconcat(*charts, spacing=5)
    return function_heatmap


# TODO: refactor this to handle splitting large numbers of genomes into multiple heatmaps here
def fill_product_dfs(
    annotations,
    module_nets,
    etc_module_df,
    function_heatmap_form,
    logger,
    annotation_ids_by_row: pd.DataFrame,
    groupby_column=DEFAULT_GROUPBY_COLUMN,
):
    module_coverage_frame = make_module_coverage_frame(
        annotations, module_nets, groupby_column
    )

    # make ETC frame
    etc_coverage_df = make_etc_coverage_df(etc_module_df, annotation_ids_by_row, groupby_column)

    # make functional frame
    function_df = make_functional_df(
        annotations,
        function_heatmap_form,
        logger,
        annotation_ids_by_row,
        groupby_column,
    )

    return module_coverage_frame, etc_coverage_df, function_df


def rename_genomes_to_taxa(function_df, labels, mag_order):
    function_df = function_df.copy()
    new_genome_column = [labels[i] for i in function_df["genome"]]
    function_df["genome"] = pd.Series(new_genome_column, index=function_df.index)
    mag_order = [labels[i] for i in mag_order]
    return function_df, mag_order


def make_product_heatmap(
    module_coverage_frame, etc_coverage_df, function_df, mag_order=None, labels=None
):
    module_coverage_heatmap = make_module_coverage_heatmap(
        module_coverage_frame, mag_order
    )
    etc_heatmap = make_etc_coverage_heatmap(etc_coverage_df, mag_order=mag_order)
    if labels is not None:
        function_df, mag_order = rename_genomes_to_taxa(function_df, labels, mag_order)
    function_heatmap = make_functional_heatmap(function_df, mag_order)

    liquor = alt.hconcat(
        alt.hconcat(module_coverage_heatmap, etc_heatmap), function_heatmap
    )
    return liquor


def make_product_df(module_coverage_frame, etc_coverage_df, function_df):
    liquor_df = pd.concat(
        [
            module_coverage_frame.pivot(
                index="genome", columns="module_name", values="step_coverage"
            ),
            etc_coverage_df.pivot(
                index="genome", columns="complex_module_name", values="percent_coverage"
            ),
            function_df.pivot(
                index="genome", columns="category_function_name", values="present"
            ),
        ],
        axis=1,
        sort=False,
    )
    return liquor_df


def get_phylum_and_most_specific(taxa_str):
    taxa_ranks = [i[3:] for i in taxa_str.split(";")]
    phylum = taxa_ranks[1]
    most_specific_rank = TAXONOMY_LEVELS[sum([len(i) > 0 for i in taxa_ranks]) - 1]
    most_specific_taxa = taxa_ranks[sum([len(i) > 0 for i in taxa_ranks]) - 1]
    if most_specific_rank == "d":
        return "d__%s;p__" % most_specific_taxa
    if most_specific_rank == "p":
        return "p__%s;c__" % most_specific_taxa
    else:
        return "p__%s;%s__%s" % (phylum, most_specific_rank, most_specific_taxa)


def make_strings_no_repeats(genome_taxa_dict):
    labels = dict()
    seen = Counter()
    for genome, taxa_string in genome_taxa_dict.items():
        final_taxa_string = "%s_%s" % (taxa_string, str(seen[taxa_string]))
        seen[taxa_string] += 1
        labels[genome] = final_taxa_string
    return labels


def get_past_annotations(
    annotation_run: Optional[dict], output_dir, force
) -> pd.DataFrame:
    if annotation_run is None:
        raise DramUsageError(
            "There is no annotations recorded in the project_config provided.\n\n"
            "It must be the case that the DRAM directory does not contain the result of "
            "a successful annotation call.\n"
            "Run `dram2 -o this_dram_dir get_status` to see if annotations have been run on this dram directory "
            "or if it is valid at all. "
            "If you have called genes but not run annotations then run "
            "`dram2 -o this_output_dir annotate db_set 'distill'` in order to get the minimal "
            "annotation set for distillation.\n"
            " Review the documentation to learn more about the required pipeline needed "
            "to run dram distill."
        )
    relative_annotation_path: Optional[str] = annotation_run.get(ANNOTATION_FILE_TAG)
    if relative_annotation_path is None or annotation_run is None:
        raise dramusageerror(
            "There is no annotations.tsv recorded in the project_config provided.\n\n"
            "It must be the case that the DRAM directory does not contain the result of "
            "a successful annotation call.\n"
            "Run `dram2 get_status` to see if annotations have been run on this dram directory "
            "or if it is valid at all. "
            "If you have called genes but not run annotations then run "
            "`dram2 -o this_output_dir annotate db_set 'distill'` in order to get the minimal "
            "annotation set for distillation.\n"
            " Review the documentation to learn more about the required pipeline needed "
            "to run dram distill."
        )
    annotations_path = output_dir / relative_annotation_path
    if not annotations_path.exists():
        raise DramUsageError(
            f"The path to annotations exists but it does not point to a annotations file that exists in the dram_directory make sure the path to your annotations is at the relive path {relative_annotation_path} with respect to the dram_directory: {output_dir}."
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

    return pd.read_csv(annotations_path, sep="\t", index_col=0)


def get_dbkit_genome_summary_forms(
    annotation_run: Optional[dict],
    dram_config: dict,
    force: bool,
    use_db_distilate: Optional[list],
) -> list[Path]:
    db_kits_with_distilate = [
        i for i in DB_KITS if i.selectable and i.has_genome_summary
    ]
    if annotation_run is None and use_db_distilate is None:
        return []
    elif annotation_run is not None and use_db_distilate is None:
        dbs = set(annotation_run[USED_DBS_TAG])
        dist_paths = [
            i(dram_config).get_genome_summary()
            for i in db_kits_with_distilate
            if i.name in dbs
        ]
        return [i.get_genome_summary() for i in dist_paths]
    elif annotation_run is not None and use_db_distilate is not None and not force:
        dbs = set(annotation_run[USED_DBS_TAG])
        missing = [i for i in use_db_distilate if i not in annotation_run]
        if len(missing) > 0:
            raise DramUsageError(
                "There are databases missing for the distillate selected, the missing databases are: {missing}"
            )
        dist_paths = [
            i(dram_config).get_genome_summary()
            for i in db_kits_with_distilate
            if i.name in use_db_distilate
        ]
        return [i.get_genome_summary() for i in dist_paths]
    elif annotation_run is None and use_db_distilate is not None and force:
        logging.debug("Skipping normal checks for annotations")
        dist_paths = [
            i(dram_config).get_genome_summary()
            for i in db_kits_with_distilate
            if i.name in use_db_distilate
        ]
        return [i.get_genome_summary() for i in dist_paths]
    else:
        raise ValueError(
            "The developer failed to support the use case for get_dbkit_genome_summary_forms, contact the developer"
        )


def test_get_dbkit_genome_summary_forms():
    class TestDBKit(DBKit):
        name = "test"
        has_genome_summary = True

        def get_genome_summary(self):
            return Path("test/path")

        def search(self):
            pass

    get_dbkit_genome_summary_forms({USED_DBS_TAG: "test"}, {}, True, set(), {"test"})


def distill(
    run_id: str,
    logger: logging.Logger,
    output_dir: Path,
    project_config: dict,
    dram_config: dict,
    # not from context
    annotations_tsv_path: Optional[Path],
    modules: tuple[str, str, str],
    force: bool,
    trna_path: Optional[Path],
    rrna_path: Optional[Path],
    custom_summary_form: Optional[Path],
    show_gene_names: bool = DEFAULT_SHOW_DISTILLATE_GENE_NAMES,
    genomes_per_product: int = DEFAULT_GENOMES_PER_PRODUCT,
    make_big_html: bool = DEFAULT_MAKE_BIG_HTML,
    groupby_column: str = DEFAULT_GROUPBY_COLUMN,
    use_db_distilate: Optional[list] = None,
) -> tuple[Optional[list[Fasta]], Optional[Path], Optional[Path], Optional[Path]]:

    # Out_paths for config
    genome_stats_path: Optional[Path] = None
    metabolism_summary_output_path: Optional[Path] = None
    product_tsv_output: Optional[Path] = None

    annotation_run: Optional[dict] = None
    called_fastas: Optional[dict] = project_config.get(FASTAS_CONF_TAG)
    fastas = (
        None
        if called_fastas is None
        else [Fasta.import_strings(output_dir, *j) for j in called_fastas]
    )
    if annotations_tsv_path is not None:
        if not force:
            raise ValueError(
                "You need to pass the -f/--force flag also, in order to bypass the checks, that would rely on the project config, in order to ensure that you "
            )
        logger.warning("You are u")
        annotations = pd.read_csv(annotations_tsv_path, sep="\t", index_col=0)
    else:
        annotation_run = get_past_annotation_run(project_config)
        annotations = get_past_annotations(annotation_run, output_dir, force)
    if "bin_taxnomy" in annotations:
        annotations = annotations.sort_values("bin_taxonomy")

    trna = None if trna_path is None else pd.read_csv(trna_path, sep="\t")
    rrna = None if rrna_path is None else pd.read_csv(rrna_path, sep="\t")
    # SUMMARIZE_METABOLISM_TAG, MAKE_GENOME_STATS_TAG, MAKE_PRODUCT_TAG

    db_kits_with_ids = [i for i in DB_KITS if i.selectable and i.can_get_ids]
    if not force:
        dbs_we_have_ano = set(annotation_run[USED_DBS_TAG])
        db_kits_with_ids = [i for i in db_kits_with_ids if i.name in dbs_we_have_ano]
    else:
        logger.warning("Skipping the normal checks because of the force flag")
    logger.debug("Loading dram_config into db_kits")
    db_kits = [i(dram_config, logger) for i in db_kits_with_ids]
    annotation_ids_by_row: pd.dataframe = get_annotation_ids_by_row(
        annotations, db_kits
    )
    module_set: set[str] = set(modules)
    if MAKE_GENOME_STATS_TAG in module_set:
        genome_stats_path = output_dir / "genome_stats.tsv"
        make_genome_stats_file(
            genome_stats_path,
            annotations,
            trna,
            rrna,
            logger,
            groupby_column,
        )
    if SUMMARIZE_METABOLISM_TAG in module_set:
        # load the base summary form from config of package data
        genome_summary_form = get_distillate_sheet(
            GENOME_SUMMARY_FORM_TAG, dram_config, logger
        )
        # load runtime custom summary form
        if custom_summary_form is not None:
            genome_summary_form = pd.concat(
                [genome_summary_form, pd.read_csv(custom_summary_form, sep="\t")]
            )
        db_kit_forms = get_dbkit_genome_summary_forms(
            annotation_run, dram_config, force, use_db_distilate
        )
        # Merge base, custom,  and db_kit summary forms into one
        genome_summary_form = pd.concat(
            [genome_summary_form] + [pd.read_csv(i, sep="\t") for i in db_kit_forms]
        )
        # TODO add more information here
        # output location is hard coded for now
        metabolism_summary_output_path = output_dir / "metabolism_summary.xlsx"

        make_metabolism_summary(
            metabolism_summary_output_path,
            annotations,
            annotation_ids_by_row,
            trna,
            rrna,
            genome_summary_form,
            logger,
            groupby_column,
            show_gene_names,
        )
    if MAKE_PRODUCT_TAG in module_set:
        module_steps_form = get_distillate_sheet(
            MODULE_STEPS_FORM_TAG, dram_config, logger
        )
        etc_module_df = get_distillate_sheet(ETC_MODULE_DF_TAG, dram_config, logger)
        function_heatmap_form = get_distillate_sheet(
            FUNCTION_HEATMAP_FORM_TAG, dram_config, logger
        )
        product_tsv_output = output_dir / "product.tsv"
        product_html_base_path = output_dir / "product"
        make_product(
            product_tsv_output,
            product_html_base_path,
            annotations,
            logger,
            module_steps_form,
            etc_module_df,
            function_heatmap_form,
            annotation_ids_by_row,
            groupby_column,
            genomes_per_product,
            make_big_html,
        )

        # NOTE: We don't document the html because it is not used as input and we can't guarantee it gets made
    return fastas, genome_stats_path, metabolism_summary_output_path, product_tsv_output


def make_genome_stats_file(
    output_path: Path,
    annotations: pd.DataFrame,
    trna: Optional[pd.DataFrame],
    rrna: Optional[pd.DataFrame],
    logger: logging.Logger,
    groupby_column: str = DEFAULT_GROUPBY_COLUMN,
):
    # make genome stats
    genome_stats = make_genome_stats(
        annotations, rrna, trna, groupby_column=groupby_column
    )
    genome_stats.to_csv(output_path, sep="\t", index=False)
    logger.info("Calculated genome statistics")


def make_metabolism_summary(
    output_path: Path,
    annotations: pd.DataFrame,
    annotation_ids_by_row: pd.DataFrame,
    trna: Optional[pd.DataFrame],
    rrna: Optional[pd.DataFrame],
    genome_summary_form: pd.DataFrame,
    logger: logging.Logger,
    groupby_column: str = DEFAULT_GROUPBY_COLUMN,
    show_distillate_gene_names: bool = DEFAULT_SHOW_DISTILLATE_GENE_NAMES,
):
    # make genome metabolism summary

    if show_distillate_gene_names:
        summarized_genomes = fill_genome_summary_frame_gene_names(
            annotations, genome_summary_form, groupby_column, logger
        )
    else:
        summarized_genomes = make_genome_summary(
            genome_summary_form,
            annotation_ids_by_row,
            logger,
            trna,
            rrna,
            groupby_column,
        )
    write_summarized_genomes_to_xlsx(summarized_genomes, output_path)
    logger.info("Generated genome metabolism summary")


def make_product(
    tsv_output: Path,
    html_base_path: Path,
    annotations: pd.DataFrame,
    logger: logging.Logger,
    module_steps_form: pd.DataFrame,
    etc_module_df: pd.DataFrame,
    function_heatmap_form: pd.DataFrame,
    annotation_ids_by_row: pd.DataFrame,
    groupby_column: str = DEFAULT_GROUPBY_COLUMN,
    genomes_per_product: int = DEFAULT_GENOMES_PER_PRODUCT,
    make_big_html: bool = DEFAULT_MAKE_BIG_HTML,
):

    # make product
    if "bin_taxonomy" in annotations:
        genome_order = get_ordered_uniques(
            annotations.sort_values("bin_taxonomy")[groupby_column]
        )
        # if gtdb format then get phylum and most specific
        if all(
            [
                i[:3] == "d__" and len(i.split(";")) == 7
                for i in annotations["bin_taxonomy"].fillna("")
            ]
        ):
            taxa_str_parser = get_phylum_and_most_specific
        # else just throw in what is there
        else:
            taxa_str_parser = lambda x: x
        labels = make_strings_no_repeats(
            {
                row[groupby_column]: taxa_str_parser(row["bin_taxonomy"])
                for _, row in annotations.iterrows()
            }
        )
    else:
        genome_order = get_ordered_uniques(
            annotations.sort_values(groupby_column)[groupby_column]
        )
        labels = None

    # make module coverage frame
    module_nets = {
        module: build_module_net(module_df)
        for module, module_df in module_steps_form.groupby("module")
        if module in HEATMAP_MODULES
    }

    if len(genome_order) > genomes_per_product and (
        make_big_html or len(genome_order) < GENOMES_PRODUCT_LIMIT
    ):

        module_coverage_dfs = list()
        etc_coverage_dfs = list()
        function_dfs = list()
        # generates slice start and slice end to grab from genomes and labels from 0 to end of genome order
        pairwise_iter = pairwise(
            list(range(0, len(genome_order), genomes_per_product)) + [len(genome_order)]
        )
        for i, (start, end) in enumerate(pairwise_iter):
            genomes = genome_order[start:end]
            annotations_subset = annotations.loc[
                [genome in genomes for genome in annotations[groupby_column]]
            ]
            dfs = fill_product_dfs(
                annotations_subset,
                module_nets,
                etc_module_df,
                function_heatmap_form,
                logger,
                annotation_ids_by_row,
            )
            module_coverage_df_subset, etc_coverage_df_subset, function_df_subset = dfs
            module_coverage_dfs.append(module_coverage_df_subset)
            etc_coverage_dfs.append(etc_coverage_df_subset)
            function_dfs.append(function_df_subset)
            product = make_product_heatmap(
                module_coverage_df_subset,
                etc_coverage_df_subset,
                function_df_subset,
                genomes,
                labels,
            )
            product.save(f"{html_base_path.as_posix()}_{i}.html")
            logger.info(f"Generated product heatmap {i}")
        product_df = make_product_df(
            pd.concat(module_coverage_dfs),
            pd.concat(etc_coverage_dfs),
            pd.concat(function_dfs),
        )
    else:
        module_coverage_df, etc_coverage_df, function_df = fill_product_dfs(
            annotations,
            module_nets,
            etc_module_df,
            function_heatmap_form,
            logger,
            annotation_ids_by_row,
        )
        product_df = make_product_df(module_coverage_df, etc_coverage_df, function_df)
        if make_big_html or len(genome_order) < GENOMES_PRODUCT_LIMIT:
            product = make_product_heatmap(
                module_coverage_df, etc_coverage_df, function_df, genome_order, None
            )
            product.save(f"{html_base_path.as_posix()}.html")
            logger.info(f"Generated product heatmap")
    product_df.to_csv(tsv_output, sep="\t")
    logger.info("Generated product table")

    logger.info("Completed distillation")


"""
import os

os.system("rm -r ./test_15soil/distillation")
os.system("DRAM.py distill -i ./test_15soil/annotations.tsv -o ./test_15soil/distillation")

"""

