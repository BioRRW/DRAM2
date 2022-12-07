import click

__version__ = '2.0.0'


@click.command('summarize_amgs')
@click.version_option(__version__)
def amg_summary():
    print("This comand is comming soon")
