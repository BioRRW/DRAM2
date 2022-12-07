import click

__version__ = '2.0.0'


@click.command('build_database_set')
@click.version_option(__version__)
def db_builder():
    print("This comand is comming soon")
