"""
The Comandline interface
========================

One file to centralize all DRAM commands

This is an interesting file, it brakes some python conventions in the interest of speed that is why as many imports as possible are in the functions so they load there, and don't slow things down.

TODO:

 - Wrap the full comand instead of individulas and pull the log from the click context https://click.palletsprojects.com/en/8.1.x/advanced/#global-context-access
 - working dir too
 - consider adding the ability to use help in sub comands and overide the group, this is eassy to do. but It would make the subcomands less indepentndt and and invokes inharatences over composition without necesisty https://stackoverflow.com/questions/54106603/python-click-display-help-for-subcommand-even-if-exception-is-raised-in-group
 - Move the context definion to a different section
 - add a safe and replace force

"""
#! /bin/python3
import click
from pathlib import Path
from typing import Optional

from dram2.cli.context import DramContext, DEFAULT_KEEP_TMP
from dram2.call_genes import call_genes_cmd
from dram2.annotate import annotate_cmd, list_databases

from dram2.utils.globals import DEFAULT_FORCE, DEFAULT_OUTPUT_DIR

from datetime import datetime

DEFAULT_VERBOSE = 3 # maps to logging levels



@click.group(
    chain=True,
    help=(
        "\b\n"  # this \b tells click to respect your formating till \n\n
        "____________  ___  ___  ___\n"
        "|  _  \\ ___ \\/ _ \\ |  \\/  |\n"
        "| | | | |_/ / /_\\ \\| .  . |\n"
        "| | | |    /|  _  || |\\/| |\n"
        "| |/ /| |\\ \\| | | || |  | |\n"
        "|___/ \\_| \\_\\_| |_/\\_|  |_/\n\n"
        "Welcome to DRAM2, the premiere metabolic annotator and analyzer. "
        "Bellow you will find a set of commands that you can use to annualize "
        "your MAGS V-mags or other genome collections. If you are new why not "
        "start with the pre defined study commands which will take you from "
        "raw data to defined distillate in a mater of seconds. \n\nThere are "
        "many advanced options to choose from so pleas look over the DRAM2 "
        "documentation at: www.ADDTHISWHENREADTHEDOCSISPUBLISHED.com"
    ),
)
@click.option(
    "-d",
    "--db_path",
    type=click.Path(path_type=Path),
    default=None,
    help=(
        "Specify a directory Path for any databases that you need, but have "
        "not setup to live. If you don't specify this path and request a "
        "databases that has not yet been built this will fail. Dram will not "
        "use files that are in this database if they are not also specified in "
        "the config file you pass."
    ),
)
@click.option(
    "--config_file",
    type=click.Path(path_type=Path),
    help="Point to a config file that you would like to use.",
)
@click.option(
    "--keep_tmp",
    is_flag=True,
    show_default=True,
    default=DEFAULT_KEEP_TMP,
    help="Keep all temporary files",
)
@click.option(
    "-o", "--output_dir", type=click.Path(path_type=Path), help="output directory"
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    default=DEFAULT_VERBOSE,
    help=(
        "Verbosity of the logging output, the number of 'v's maps to the default"
        " logging level of pythons logging modual. the default is -vvv. The "
        "mapping is 0 = CRITICAL ,-v = ERROR, -vv = WARNING, -vvv = INFO, "
        "-vvvv = DEBUG, -vvvvv... = NOTSET "
    ),
)
@click.option("-c", "--cores", type=int, default=10, help="number of processors to use")
@click.pass_context
def dram2(
    ctx,
    cores: int,
    db_path: Optional[Path] = None,
    config_file: Optional[Path] = None,
    log_file_path: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    # force: bool = DEFAULT_FORCE,
    keep_tmp: bool = DEFAULT_KEEP_TMP,
    verbose: int= DEFAULT_VERBOSE,
):
    ctx.obj = DramContext(
        cores=cores,
        db_path=db_path,
        config_file=config_file,
        log_file_path=log_file_path,
        output_dir=output_dir,
        verbose=verbose,
        keep_tmp=keep_tmp,
    )


def dram2_logged_entry():
    click.get_current_context()
    pass

# These are the essential components of dram, and call may not be
dram2.add_command(annotate_cmd)
dram2.add_command(call_genes_cmd)
dram2.add_command(list_databases)

try:
    from dram2.rule_adjectives import evaluate, rule_plot
    dram2.add_command(evaluate)
    dram2.add_command(rule_plot)
except ImportError:
    pass

try:
   from dram2.distill import distill_cmd
   dram2.add_command(distill_cmd)
except ImportError:
    pass

try:
   from dram2.genbank import generate_genbank
   dram2.add_command(generate_genbank)
except ImportError:
    pass

try:
   from dram2.merger import merger
   dram2.add_command(merger)
except ImportError:
    pass

try:
   from dram2.strainer import strainer
   dram2.add_command(strainer)
except ImportError:
    pass

try:
   from dram2.rna import pull_trna, pull_rrna
   dram2.add_command(pull_rrna)
   dram2.add_command(pull_trna)
except ImportError:
    pass

try:
   from dram2.db_builder import db_builder
   dram2.add_command(db_builder)
except ImportError:
    pass

try:
   from dram2.tree_kit import phylo_tree
   dram2.add_command(phylo_tree)
except ImportError:
    pass

try:
   from dram2.amg_summary import amg_summary
   dram2.add_command(amg_summary)
except ImportError:
    pass


if __name__ == "__main__":
    dram2()
