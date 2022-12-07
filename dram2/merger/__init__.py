import click

__version__ = '2.0.0'

@click.command('merge')
@click.version_option(__version__)
def merger():
    print("This comand is comming soon")
