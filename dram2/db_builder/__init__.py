import click


from dram2.cli.context import DramContext, DEFAULT_KEEP_TMP, __version__


def db_builder():
    print("This comand is comming soon")


@click.command('build_database_set',
               context_settings=dict(help_option_names=["-h", "--help"]),
               )
@click.pass_context
def db_builder_md():
    """
    Build Your Own Custom DRAM Database (Not Ready)
    ___

    A tool to build your own custom dram database. Download the raw data and format it for use in dram. This is only for advanced users but it is necessary for our kegg users.

    This command should only call process that are in the db_kits. That is where downloading should be handled.
    TODO:

    - Fully Implement the Download options in the DB_Kits
    - Fully Implement the process options in the DB_KITS

    """
    print("Sorry this command is not implemented yet")
