import click

__version__ = '2.0.0'

@click.command('pull_rrna')
@click.version_option(__version__)
def pull_rrna():
    print("This comand is comming soon")

@click.command('pull_trna')
@click.version_option(__version__)
def pull_trna():
    print("This comand is comming soon")
