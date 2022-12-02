"""
One file to centralize all DRAM commands

This is an interesting file, it brakes some python conventions in the interest of speed that is why as many imports as possible are in the functions so they load there, and don't slow things down.
"""
#! /bin/python3
import click
from dram2.annotate import annotate, list_databases


@click.group(chain=True)
@click.option("--verbose/--quiet", default=False)
def dram2(verbose):
    pass



dram2.add_command(annotate)
dram2.add_command(list_databases)

if __name__ == '__main__':
    dram2()
