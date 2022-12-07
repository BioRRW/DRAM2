import click

__version__ = '2.0.0'


@click.command('strainer')
@click.version_option(__version__)
def strainer():
    print("This comand is comming soon")
