"""
DRAM Distillate
---------------

This is the script that distills the genomes from the annotations step.
It requires annotations and the dram distillation step

"""
import click
from pathlib import Path
from dram2.cli.context import DramContext, get_time_stamp_id, __version__


from typing import Optional
import logging


from dram2.distill.summarize_genomes import (
    distill,
    DEFAULT_GROUPBY_COLUMN,
    DEFAULT_SHOW_DISTILLATE_GENE_NAMES,
    DEFAULT_GENOMES_PER_PRODUCT,
    DEFAULT_SHOW_DISTILLATE_GENE_NAMES,
    GENOMES_PRODUCT_LIMIT,
    DISTILLATE_MODULES,
    DB_KITS,
    LOCATION_TAG,
)

from dram2.annotate import FASTAS_CONF_TAG

DISTILLATE_RUN_TAG = "distill"
GENOME_STATS_TAG = "genome_stats"
METABOLISM_SUMMARY_TAG = "metabolism_summary"
DISTILLATE_TAG = "distillation"
PRODUCT_TAG = "product"

COMMAND_NAME = "distill"


@click.command(COMMAND_NAME)
@click.option(
    "--annotations_tsv_path",
    type=click.Path(exists=True, path_type=Path),
    help="Location of an annotations.tsv. You don't need to use this option if you are using the output_dir for dram with a project config. If you use this option, you must also use the force flag to bypass the safeguards that prevent you from running distill with insufficient data",
)
@click.option(
    "-m",
    "--modules",
    help="What distillate module to run. It can be time consuming to run all the distillate module for all projects.",
    type=click.Choice(DISTILLATE_MODULES, case_sensitive=False),
    default=DISTILLATE_MODULES,
    multiple=True,
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="Remove skip the normal checks.",
)
# @click.option(
#    "--groupby_column",
#    default=DEFAULT_GROUPBY_COLUMN,
#    help="Column from the annotations tsv to group as organism units",
# )
# These arguments are only for the Metabolism Summary, and the Genome Summary
@click.option(
    "--rrna_path",
    type=click.Path(exists=True, path_type=Path),
    help="rRNA output from a dram RNA script. You don't need to explicitly give this path if you are using an output_dir from dram with a project config file. The rRNA run will be automatically detected if you have a project config file.",
)
@click.option(
    "--trna_path",
    type=click.Path(exists=True, path_type=Path),
    help="tRNA output from a dram annotation. You don't need to explicitly give this path if you are using an output_dir from dram with a project_config file. The tRNA run will be automatically detected if you have a project config file.",
)
# These arguments are only for the Metabolism Summary
@click.option(
    "--show_gene_names",
    is_flag=True,
    default=DEFAULT_SHOW_DISTILLATE_GENE_NAMES,
    help="Give names of genes instead of counts in genome metabolism summary. This tool is not fully supported, and may run into the limits of Excel. Use with caution.",
)
@click.option(
    "--use_db_distilate",
    multiple=True,
    default=None,
    type=click.Choice([i.name for i in DB_KITS], case_sensitive=False),
    help="Specify exactly which db specific distillate to use. If you know what you are doing it may be useful to force the output of the program. If you already annotated with a database that has an associated distillate file eg:methyl and still have your project config, there should be no need for this command. If you use this command, you should have a good idea what you are doing and use the force command also.",
)
@click.option(
    "--custom_summary_form",
    type=click.Path(exists=True, path_type=Path),
    help="Custom distillate form to add your own modules to the metabolism summary. You will need to read the docs to find the format that this tsv file must take.",
)
# These arguments are only for the product
@click.option(
    "--genomes_per_product",
    help=(
        "Number of genomes per product.html output. Decrease "
        "value if getting JavaScript Error: Maximum call stack "
        "size exceeded when viewing product.html in browser. "
        "Note that by default the product html will not be created "
        f"if the number of genomes is over {GENOMES_PRODUCT_LIMIT}. "
        "You must pass the make_big_html flag in order to make that html"
    ),
    default=DEFAULT_GENOMES_PER_PRODUCT,
    type=int,
)
@click.option(
    "--make_big_html",
    is_flag=True,
    help=f"It is felt that if the number of genomes is over {GENOMES_PRODUCT_LIMIT} "
    "that product may be of limited use because of the size and the number of html "
    "files that will be made. In order to avoid the "
    "large amount of time it will take to make these distillates it makes sense to just make the product html",
)
# TODO let the product html be make from the product itself
@click.pass_context
def distill_cmd(
    ctx: click.Context,
    annotations_tsv_path: Optional[Path],
    modules: tuple[str, str, str],
    force: bool,
    # groupby_column: str,
    trna_path: Optional[Path],
    rrna_path: Optional[Path],
    show_gene_names: bool,
    genomes_per_product: int,
    make_big_html: bool,
    custom_summary_form: Optional[Path],
    use_db_distilate: Optional[list],
):
    """
    DRAM Distillate
    ---

    This command can generate three outputs. The user can select the files they want to make with the -m/--module options, by default it generates all of them.

    The output files are:

    - The genome_summary.xlsx which contains a summary of metabolisms present in each genome. It gives gene by gene information across various metabolisms for every genome in your dataset.
     - The genome_statistics.tsv file contains all measures required by the MIMAG about each fasta used as input.
     - The product.html and product.tsv allow the user to understand the metabolic potential of each MAG or Collection at a glance. The product.html is an interactive html that allows users to hover over each box to see what genes prompted the box color (Example here) and was manually curated to consider alternate genes for pathways and single processes. This heat map allows the user to quickly profile ecosystem relevant processes across hundreds of genomes. The product.tsv provides the same information in a tsv format. If you have many FASTA's your product html may be split into many output files.

    """
    context: DramContext = ctx.obj
    run_id: str = get_time_stamp_id(DISTILLATE_RUN_TAG)
    logger: logging.Logger = logging.getLogger("dram2_log")
    context.get_logger()
    output_dir: Path = context.get_output_dir()
    project_config: dict = context.get_project_config()
    dram_config: dict = context.get_dram_config(logger)  # FIX
    try:
        (
            fastas_distilled,
            genome_stats_path,
            metabolism_summary_output_path,
            product_tsv_output,
        ) = distill(
            run_id=run_id,
            logger=logger,
            output_dir=output_dir,
            project_config=project_config,
            dram_config=dram_config,
            annotations_tsv_path=annotations_tsv_path,
            modules=modules,
            force=force,
            # groupby_column=groupby_column,
            trna_path=trna_path,
            rrna_path=rrna_path,
            custom_summary_form=custom_summary_form,
            show_gene_names=show_gene_names,
            genomes_per_product=genomes_per_product,
            make_big_html=make_big_html,
            use_db_distilate=use_db_distilate,
        )
        new_project_config = {
            "latest": run_id,
            run_id: {
                "version": __version__,
                FASTAS_CONF_TAG: None
                if fastas_distilled is None
                else [i.name for i in fastas_distilled],
                GENOME_STATS_TAG: {
                    LOCATION_TAG: path_str_for_project_config(
                        output_dir, genome_stats_path
                    )
                },
                METABOLISM_SUMMARY_TAG: {
                    LOCATION_TAG: path_str_for_project_config(
                        output_dir, metabolism_summary_output_path
                    )
                },
                PRODUCT_TAG: {
                    LOCATION_TAG: path_str_for_project_config(
                        output_dir, product_tsv_output
                    )
                },
            },
        }
        if DISTILLATE_TAG in project_config:
            project_config[DISTILLATE_TAG].update(new_project_config)
        else:
            project_config.update({DISTILLATE_TAG: new_project_config})
        context.set_project_config(project_config)
    except Exception as e:
        logger.error(e)
        logger.exception(f"Fatal error in {COMMAND_NAME}")
        raise (e)


def path_str_for_project_config(output_dir, out_path: Optional[Path]) -> Optional[str]:
    return out_path if out_path is None else out_path.relative_to(output_dir).as_posix()

