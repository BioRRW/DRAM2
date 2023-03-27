"""
he Comandline interface
========================

One file to centralize all DRAM commands

This is an interesting file, it brakes some python conventions in the interest of speed
that is why as many imports as possible are in the functions so they load there, and
don't slow things down.
conda activate ./dram2_env
dram2 --help
dram2 -h

TODO:

 - Wrap the full comand instead of individulas and pull the log from the click context
   https://click.palletsprojects.com/en/8.1.x/advanced/#global-context-access.
 - consider adding the ability to use help in sub comands and overide the group, this
   is eassy to do. but It would make the subcomands less indepentndt and and invokes
   inharatences over composition without necesisty
   https://stackoverflow.com/questions/54106603/python-click-display-help-for-subcommand-even-if-exception-is-raised-in-group.
 - Move the context definion to a different section
 - add a safe and replace force

"""
from dram2.rule_adjectives import adjectives_cmd
from dram2.tree_kit import phylo_tree_cmd
from dram2.neighbors import neighbors_cmd
from dram2.strainer import strainer_cmd
from dram2.merger import merger_cmd
from dram2.genbank import generate_genbank_cmd
from dram2.distill import distill_cmd
from dram2.rna import pull_trna_cmd, pull_rrna_cmd
from pathlib import Path
from typing import Optional

import click

from dram2.cli.context import DramContext, DEFAULT_KEEP_TMP, __version__, OrderedGroup
from dram2.call_genes import call_genes_cmd
from dram2.annotate import annotate_cmd, list_databases


DEFAULT_VERBOSE = 4  # maps to logging levels


@click.group(
    cls=OrderedGroup,
    help=(
        "\b\n"  # this \b tells click to respect your formating till \n\n
        "____________  ___  ___  ___\n"
        "|  _  \\ ___ \\/ _ \\ |  \\/  |\n"
        "| | | | |_/ / /_\\ \\| .  . |\n"
        "| | | |    /|  _  || |\\/| |\n"
        "| |/ /| |\\ \\| | | || |  | |\n"
        "|___/ \\_| \\_\\_| |_/\\_|  |_/"
        """


        Welcome to DRAM2, the premiere metabolic annotator and analyzer.
        Bellow you will find a set of commands that you can use to annualize
        your MAGS V-mags or other genome collections.\n\n If you are new why not
        start with the pre defined protocal commands which will take you from
        raw data to refined distillate in a mater of seconds. \n\nThere are
        many advanced options to choose from so pleas look over the DRAM2
        documentation at: www.ADDTHISWHENREADTHEDOCSISPUBLISHED.com
        """),
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@ click.option(
    "-d",
    "--dram_dir",
    type=click.Path(path_type=Path),
    help="""
    This is both the location where the output of new DRAM2 actions will go, and
    also the location where the ouputs of past DRAM2 actions and metadata can be found.
    """,
)
@ click.option(
    "--config_file",
    type=click.Path(path_type=Path),
    help="""
    Point to a config file that you would like to use for this action specificaly. With
    out this argument DRAM2 looks for the config at first in
    `USER_HOME/.config/dram_config.yaml` then in `/etc/dram2_config.yaml`.
    """,
)
@ click.version_option(__version__)
@ click.option(
    "-v",
    "--verbose",
    count=True,
    type=int,
    default=DEFAULT_VERBOSE,
    help="""
    Verbosity of the logging output, the number of 'v's maps to the default
    logging level of pythons logging modual. the default is 4 aka - vvvv. The mapping is
    -v=CRITICAL, -vv=ERROR, -vvv=WARNING, -vvvv=INFO, -vvvvv=DEBUG,
    -vvvvvv...=NOTSET.
    """,
)
@ click.option(
    "-t",
    "--threads",
    type=int,
    default=10,
    help="""
    number of threads to run or the number of processors to use. This command is passed
    to external tools called by DRAM2, as well as being used in DRAM2 its self. Note
    that increasing the threads may increase the amount of memory required.
    """,
)
@ click.option(
    "--keep_tmp",
    is_flag=True,
    show_default=True,
    default=DEFAULT_KEEP_TMP,
    help="Keep all temporary files",
)
# @click.option(
#     "-d",
#     "--db_path",
#     type = click.Path(path_type=Path),
#     default = None,
#     help = (
#         "Specify a directory Path for any databases that you need, but have "
#         "not setup to live. If you don't specify this path and request a "
#         "databases that has not yet been built this will fail. Dram will not "
#         "use files that are in this database if they are not also specified in "
#         "the config file you pass."
#     ),
# )
@ click.pass_context
def dram2(
    ctx: click.Context,
    threads: int,
    db_path: Optional[Path] = None,
    config_file: Optional[Path] = None,
    log_file_path: Optional[Path] = None,
    output_dir: Optional[Path] = None,
    # force: bool = DEFAULT_FORCE,
    keep_tmp: bool = DEFAULT_KEEP_TMP,
    verbose: int = DEFAULT_VERBOSE,
):
    ctx.obj = DramContext(
        cores=threads,
        db_path=db_path,
        config_file=config_file,
        log_file_path=log_file_path,
        output_dir=output_dir,
        verbose=verbose,
        keep_tmp=keep_tmp,
    )


# These are the essential components of dram, and call may not be
dram2.add_command(call_genes_cmd)
dram2.add_command(annotate_cmd)
dram2.add_command(list_databases)
dram2.add_command(pull_rrna_cmd)
dram2.add_command(pull_trna_cmd)

# Optional arguments, you can remove any of these and they will just magically disappear

dram2.add_command(distill_cmd)

try:
    from dram2.amg_summary import amg_summary_cmd

    dram2.add_command(amg_summary_cmd)
except ImportError:
    pass


dram2.add_command(generate_genbank_cmd)

dram2.add_command(merger_cmd)

dram2.add_command(strainer_cmd)

dram2.add_command(neighbors_cmd)


dram2.add_command(phylo_tree_cmd)

dram2.add_command(adjectives_cmd)


try:
    from dram2.db_builder import db_builder_cmd

    dram2.add_command(db_builder_cmd)
except ImportError:
    pass


if __name__ == "__main__":
    dram2()
