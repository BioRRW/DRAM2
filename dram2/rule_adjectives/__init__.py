"""Contains main entry point to the program and helper functions"""
import os
from pathlib import Path
import logging

import pandas as pd
import click

from dram2.rule_adjectives.rule_graph import RuleParser, get_positive_genes
from dram2.rule_adjectives.annotations import Annotations

from dram2.utils.utils import get_package_path, DramUsageError
from dram2.cli.context import (
    DramContext,
    __version__,
    OrderedGroup,
    get_time_stamp_id,
    __version__,
)
from dram2.annotate import (
    ANNOTATION_FILE_TAG,
    get_past_annotation_run,
    LATEST_ANNOTATION_RUN_TAG,
    check_for_annotations,
)
from dram2.annotate import DB_KITS
from dram2.db_kits.utils import DBKit

from dram2.tree_kit import (
    TREE_TAG,
    TREES,
    LATEST_TREE_RUN_TAG,
    DramTree,
    NXR_NAR_TREE,
    AMOA_PMOA_TREE,
)
from dram2.tree_kit.pplacer import DramTree


ADJECTIVES_TSV_NAME: str = "adjectives.tsv"
RULES_TSV_PATH: Path = get_package_path(Path("rule_adjectives", "rules.tsv"))

# class PythonLiteralOption(click.Option):
#     def type_cast_value(self, ctx, value):
#         try:
#             return ast.literal_eval(value)
#         except:
#             raise click.BadParameter(value)


def list_adjectives(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    rules = RuleParser(RULES_TSV_PATH)
    print("In the current rules file, these adjectives are available:")
    for i in rules.data.index[~rules.data["name"].isna()].unique():
        print(i)


def list_adjective_name(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    rules = RuleParser(RULES_TSV_PATH)
    print("In the current rules file, these adjectives are available:")
    for i in rules.data["name"].unique():
        print(i)


def show_rules_path(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print(RULES_TSV_PATH)


@click.group(
    "adjectives",
    cls=OrderedGroup,
    help="Describe Gene Features (Adjectives)\n\n The DRAM2 Adjectives are features, defined by a set of rules which themselves are based on key genes and selective use of phylogenetic trees.",
)
@click.pass_context
def adjectives_cmd(
    ctx: click.Context,
):
    """
    Describe Gene Features (Adjectives)
    ---

    The DRAM2 Adjectives are features, defined by a set of rules which themselves are based on key genes and selective use of phylogenetic trees.

    In order to use this command, you must first complete the following for all genes in your project:
    - Call the Genes with `dram2 call`
    - annotate the genes with`dram2 annotate using `dram2 annotate` be sure you use the following database:
        - KEGG or KOfam
        - Pfam
        - CAMPER
        - FeGenie
        - Sulfur
    - Add phylogenetic tree information to the annotations with `dram2 philotrees`

    This will make the adjectives.tsv
    """
    pass


@click.command("eval")
@click.option(
    "--annotations_tsv_path",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Location of an annotations.tsv. You don't need to use this option if you are using the output_dir for dram with a project_config. If you use this option, you must also use the force flag to bypass the safeguards that prevent you from running distill with insufficient data",
)
@click.option(
    "--adjectives_tsv_path",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Location of the output adjectives.tsv. if you leave this blank the adjectives.tsv file will be put in the output directory",
)
@click.option(
    "-a",
    "--adjectives",
    multiple=True,
    default=[],
    help="A list of adjectives, by name, to evaluate. This limits the number of adjectives that are evaluated and is faster.",
)
@click.option(
    "-p",
    "--plot_adjectives",
    multiple=True,
    default=[],
    help="A list of adjectives, by name, to plot. This limits the number of adjectives that are plotted and is probably needed for speed.",
)
@click.option(
    "-g",
    "--plot_genomes",
    multiple=True,
    default=[],
)
@click.option(
    "--plot_path",
    type=click.Path(exists=False),
    default=None,
    help="will become a folder of output plots, no path no plots.",
)
@click.option(
    "--strainer_tsv",
    type=click.Path(exists=False),
    default=None,
    help="The path for a tsv that will pass to strainer to filter genes. The only option at this time is ‘pgtb’ for positive genes that are on true bugs.",
)
@click.option(
    "--strainer_type",
    type=click.Path(exists=False),
    default=None,
    help="The type of process that should make the strainer file.",
)
@click.option(
    "--debug_ids_by_fasta_to_tsv",
    type=click.Path(exists=False),
    default=None,
    help="This is a tool to debug the list of IDs found by DRAM it is mostly for experts",
)
@click.option(
    "--user_rules_tsv",
    type=click.Path(exists=True),
    default=RULES_TSV_PATH,
    help="This is an optional path to a rules file with strict formatting. It will overwrite the original rules file that is stored with the script.",
)
@click.option(
    "--show_rules_path",
    is_flag=True,
    callback=show_rules_path,
    expose_value=False,
    is_eager=True,
    help="Show the path to the default rules path.",
)
@click.option(
    "--list_name",
    is_flag=True,
    callback=list_adjective_name,
    expose_value=False,
    is_eager=True,
    help="List the names for all adjectives_tsv that are"
    " available, you can pass these names to limit the"
    " adjectives that are evaluated",
)
@click.option(
    "--list_id",
    is_flag=True,
    callback=list_adjectives,
    expose_value=False,
    is_eager=True,
    help="List the names for all adjectives_tsv that are"
    " available, you can pass these names to limit the"
    " adjectives that are evaluated",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    show_default=True,
    default=False,
    help="Skip the normal checks on the dram project_config and just try to use the annotations and trees provided if available.",
)
@click.pass_context
def evaluate_cmd(
    ctx: click.Context,
    user_rules_tsv: Path,
    adjectives: list[str],
    plot_adjectives: list[str],
    plot_genomes: list[str],
    plot_path: str,
    debug_ids_by_fasta_to_tsv: str,
    strainer_tsv: Path,
    strainer_type: str,
    annotations_tsv_path: Path | None,
    adjectives_tsv_path: Path | None,
    force: bool,
):
    """
     Evaluate Genes and Describe Their Features (Adjectives)
    ---

    The DRAM2 Adjectives are features, defined by a set of rules which themselves are based on key genes and selective use of phylogenetic trees.

    In order to use this command you must first complete the falowing for all genes in your project:
    - Call the Genes with `dram2 call` (the genes dir itself is not needed)
    - annotate the genes with`dram2 annotate using `dram2 annotate` be sure you use the following database:
        - KEGG or KOfam
        - Pfam
        - CAMPER
        - FeGenie
        - Sulfur
    - Add phylogenetic tree information to the annotations with `dram2 philotrees`

    This will make the adjectives.tsv
    """

    context: DramContext = ctx.obj
    logger = context.get_logger()
    output_dir: Path = context.get_output_dir()
    cores: int = context.cores
    project_config: dict = context.get_project_config()
    dram_config = context.get_dram_config(logger)  # FIX
    try:

        run_id: str = get_time_stamp_id(TREE_TAG)
        annotation_run: dict | None = get_past_annotation_run(project_config)
        if adjectives_tsv_path is not None:
            adjectives_tsv = adjectives_tsv_path
        else:
            adjectives_tsv = output_dir / ADJECTIVES_TSV_NAME

        annotations_tsv = get_annotations_path(
            annotation_run, annotations_tsv_path, force, output_dir
        )
        if isinstance(user_rules_tsv, Path):
            rules_tsv = user_rules_tsv
        else:
            rules_tsv: Path = RULES_TSV_PATH
        if (tree_conf := project_config.get(TREE_TAG)) is not None:
            tree_dict = tree_conf[tree_conf[LATEST_TREE_RUN_TAG]]
            tree_paths:dict[str, None | Path] = {i: (None if j is None else output_dir / j) for i, j in tree_dict.items()}
        else:
            tree_paths = {}

        database = [i(dram_config, logger) for i in DB_KITS]
        evaluate(
            annotations_tsv=annotations_tsv,
            adjectives_tsv=adjectives_tsv,
            rules_tsv=rules_tsv,
            adjectives=adjectives,
            plot_adjectives=plot_adjectives,
            plot_genomes=plot_genomes,
            plot_path=plot_path,
            strainer_tsv=strainer_tsv,
            strainer_type=strainer_type,
            debug_ids_by_fasta_to_tsv=debug_ids_by_fasta_to_tsv,
            tree_paths=tree_paths,
            logger=logger,
            cores=cores,
            database=database
        )
        if TREE_TAG in project_config:
            project_config[TREE_TAG].update({run_id: {"names": tree_names}})
        else:
            project_config.update({TREE_TAG: {run_id: {"names": tree_names}}})

    except Exception as e:
        logger.error(e)
        logger.exception("Fatal error in DRAM2 Adjectives")
        raise (e)


def get_annotations_path(
    annotation_run: dict | None,
    annotations: Path | None,
    force: bool,
    output_dir: Path,
) -> tuple[Path, Path | list]:
    if force:
        logging.warning(
            "Skipping the normal checks for needed annotations "
            "because the force flag was passed."
        )
    elif annotation_run is not None:
        if (
            db_error := check_for_annotations(
                [{"kofam"}, {"kegg"}],
                annotation_run,
            )
        ) is not None:
            raise DramUsageError(db_error)
    else:
        raise DramUsageError(
            "The project_config does not exist or you have not annotated. \n"
            "Have you called genes and annotations?\n\n You must at least annotate"
            " with KOfam or KEGG in order to run this script.\n\nIf you have done"
            " the steps but don't have a project_config for whatever reason then "
            "you can use the force option to skip this check."
        )
    if isinstance(annotations, Path):
        annotation_to_use: Path = annotations
    elif annotation_run is None:
        raise DramUsageError(
            "There is no project_config in the output directory and no"
            " annotation.tsv was provided in its place."
            "\nIf your output directory has no project config file,"
            "you must use the --annotations option to point to an annotations.tsv."
        )
    else:
        annotation_to_use_str: str | None = annotation_run.get(ANNOTATION_FILE_TAG)
        if annotation_to_use_str is None:
            raise DramUsageError(
                "There is no annotations.tsv recorded in the project_config "
                "provided.\n\nIt must be the case that the DRAM directory"
                " does not contain the result of a successful annotation, "
                " because it was not run or it was moved."
                # "\nRun `dram2 get_status` to see if annotations have been "
                # "run on this dram directory or if it is valid at all. "
            )
        else:
            annotation_to_use: Path = output_dir / annotation_to_use_str
    return annotation_to_use


def evaluate(
    annotations_tsv: Path,
    adjectives_tsv: Path,
    rules_tsv: Path,
    adjectives: list[str],
    plot_adjectives: list[str],
    plot_genomes: list[str],
    plot_path: str,
    debug_ids_by_fasta_to_tsv: str,
    strainer_tsv: Path,
    strainer_type: str,
    tree_paths: dict[str, Path],
    logger: logging.Logger,
    cores: int,
    database: list[DBKit]
):
    """
    Using a DRAM annotations file, make a table of adjectives.

    :param annotations_tsv: Path to a DRAM annotations file.
    :param adjectives_tsv: Path for the output true false table.
    :param rules_tsv: Path to a rules file with strict formatting, this is optional.
    :param adjectives: Adjectives to evaluate.
    :param plot_adjectives: Adjectives to plot
    :param plot_genomes: Genomes to plot.
    :param plot_path: The path that will become a folder of output plots, no path no plots.
    :param strainer_tsv: The path for a tsv that will pass to strainer to filter genes.
    :param strainer_type: The type of process that should make the strainer file.
        the only option at this time is ‘pgtb’ for positive genes that are on true bugs.
    """

    annotations = Annotations(annotations_tsv.absolute().as_posix(), db_kits= database)
    rules = RuleParser(rules_tsv, adjectives=set(adjectives))
    tree_data: dict[str, pd.Series] = {
        i: (pd.Series() if j is None else pd.read_csv(j, sep="\t", index_col=0)[i])
        for i, j in tree_paths.items()
    }
    adjectives = rules.check_genomes(annotations, logger, tree_data)
    if debug_ids_by_fasta_to_tsv is not None:
        annotations.ids_by_fasta.to_csv(debug_ids_by_fasta_to_tsv, sep="\t")
        exit()
    adjectives.to_csv(adjectives_tsv, sep="\t")
    # annotations.ids_by_fasta.iloc[1]['annotations']
    if plot_path is not None:
        rules.plot_cause(
            plot_path,
            adjectives=plot_adjectives,
            genomes=plot_genomes,
            show_steps=False,
        )
    if strainer_tsv is not None:
        strainer_data = get_positive_genes(rules, annotations, adjectives)
        strainer_data.to_csv(strainer_tsv, sep="\t")


@click.command("adjectives_plot")
@click.argument(
    "plot_path", type=click.Path(exists=False), default=None
)  # , help='will become a folder of output plots, no path no plots.')
@click.option(
    "-a",
    "--adjectives",
    multiple=True,
    default=[],
    help="A list of adjectives, by name, to evaluate. This limits the number of adjectives that are evaluated and is faster.",
)
@click.option(
    "--rules_tsv",
    type=click.Path(exists=True),
    default=RULES_TSV_PATH,
    help="The path that will become a folder of output plots, no path no plots.",
)  # , help='The rules file which adhere to strict formatting' )
@click.option(
    "--list_name",
    is_flag=True,
    callback=list_adjective_name,
    expose_value=False,
    is_eager=True,
    help="List the names for all adjectives_tsv that are"
    " available, you can pass these names to limit the"
    " adjectives that are evaluated",
)
@click.option(
    "--show_rules_path",
    is_flag=True,
    callback=show_rules_path,
    expose_value=False,
    is_eager=True,
    help="Show the path to the default rules path.",
)
@click.option(
    "--list_id",
    is_flag=True,
    callback=list_adjectives,
    expose_value=False,
    is_eager=True,
    help="List the names for all adjectives_tsv that are"
    " available, you can pass these names to limit the"
    " adjectives that are evaluated",
)
@click.pass_context
def plot_rules_cmd(
    ctx: click.Context,
    rules_tsv: str = RULES_TSV_PATH,
    adjectives: list = None,
    plot_path: str = None,
):
    """
    Make a Set of Plots to Show How Rules are Evaluated
    ---

    Understanding how annotations use its rules to call a given gene is non-trivial. This function provides a method to name these genes.

    """
    context: DramContext = ctx.obj


def rule_plot(
    rules_tsv: Path = RULES_TSV_PATH,
    adjectives: list = None,
    plot_path: str = None,
):
    """

    :param annotations_tsv: Path to a DRAM annotations file.
    :param adjectives_tsv: Path for the output true false table.
    :param rules_tsv: Path to a rules file with strict formatting, this is optional.
    :param adjectives: Adjectives to evaluate.
    :param plot_adjectives: Adjectives to plot
    :param plot_genomes: Genomes to plot.
    :param plot_path: The path that will become a folder of output plots, no path no plots.
    :param strainer_tsv: The path for a tsv that will pass to strainer to filter genes.
    :param strainer_type: The type of process that should make the strainer file.
        the only option at this time is `pgtb` for positive genes that are on true bugs.
    """
    rules = RuleParser(rules_tsv, adjectives=adjectives)
    rules.plot_rule(plot_path)


adjectives_cmd.add_command(evaluate_cmd)
adjectives_cmd.add_command(plot_rules_cmd)

