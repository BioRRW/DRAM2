"""
One file to centralize all DRAM commands

This is an interesting file, it brakes some python conventions in the interest of speed that is why as many imports as possible are in the functions so they load there, and don't slow things down.

TODO:

 - Wrap the full comand instead of individulas and pull the log from the click context https://click.palletsprojects.com/en/8.1.x/advanced/#global-context-access
 - working dir too
"""
#! /bin/python3
import click
import logging

from pathlib import Path

from dram2.annotate.annotate import annotate, list_databases
from dram2.rule_adjectives import  evaluate, rule_plot
from dram2.tree_kit import phylo_tree
from dram2.distill import distill
from dram2.genbank import generate_genbank
from dram2.merger import merger
from dram2.strainer import strainer
from dram2.rna import pull_trna, pull_rrna
from dram2.db_builder import db_builder
from dram2.amg_summary import amg_summary

from dram2.utils.globals import DEFAULT_FORCE, DEFAULT_OUTPUT_DIR

@click.group(chain=True)
@click.option("--verbose/--quiet", default=False)
@click.option(
    "-d" "--db_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Specify a directory Path for any databases that you need, but have not setup to live. If you don't specify this path and request a databases that has not yet been bilt this will fail. Dram will not use files that are in this database if they are not also specified in the config file you pass.",
)
@click.option(
    "--config_file",
    type=click.Path(path_type=Path),
    help="Point to a config file that you would like to use.",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    show_default=True,
    default=DEFAULT_FORCE,
    help="Don't place just the ambiguous genes, place all of them",
)
@click.option("-o", "--output_dir",
    type=click.Path(path_type=Path),
              help="output directory", required=True)
@click.option("--threads", type=int, default=10, help="number of processors to use")
@click.pass_context
def dram2(ctx, db_path:Path=None, config_file:Path = None, log_file_path: Path = None, output_dir: Path = None, force:bool= DEFAULT_FORCE, verbose=None):
    # Get a logger
    logger = logging.getLogger("dram2_log")
    if log_file_path is None:
        log_file_path = output_dir / "dram2.log"
    logger.info(f"The log file is created at {log_file_path}")
    setup_logger(logger, log_file_path)
    cdx.obj['logger']
    if not output_dir.exists():
        output_dir.mkdir()
    elif not force:
        raise ValueError(
            "The output_dir already exists! try using the -f flag to overwrite"
            )


def dram2_logged_entry():
    click.get_current_context()
    pass

dram2.add_command(annotate)
dram2.add_command(list_databases)
dram2.add_command(evaluate)
dram2.add_command(rule_plot)
dram2.add_command(phylo_tree)
dram2.add_command(distill)
dram2.add_command(amg_summary)
dram2.add_command(generate_genbank)
dram2.add_command(merger)
dram2.add_command(strainer)
dram2.add_command(pull_rrna)
dram2.add_command(pull_trna)
dram2.add_command(db_builder)


if __name__ == '__main__':
    dram2()
