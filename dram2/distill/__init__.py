import click

__version__ = '2.0.0'

@click.command('distill')
@click.version_option(__version__)
def distill():
    print("This comand is comming soon")
