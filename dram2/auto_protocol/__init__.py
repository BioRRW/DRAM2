"""
===============
Auto Protocalls
===============
"""
import click
from dram2.cli.context import DramContext, __version__


@click.command("protocal",
               context_settings=dict(help_option_names=["-h", "--help"]),
               )
# @click.argument(
#     "fasta_paths",
#     type=click.Path(exists=True, path_type=Path),
#     nargs=-1,
#     help="The path/paths to FASTAs representing mags or other genome collections of uncalled genes."
# )
@click.pass_context
def auto_protocol(
    ctx: click.Context,
):
    context: DramContext = ctx.obj
    """
    From FASTAs to Final Output in One Command 
    ___

    This command allows you to give a set of FASTAs and get out the result of a pre-defined protocol with no other arguments.

    Use `dram2 protocol <protocol_command> --help` to learn more about the protocols

    """
    pass
