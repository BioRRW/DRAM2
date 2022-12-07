"""
One file to centralize all DRAM commands

This is an interesting file, it brakes some python conventions in the interest of speed that is why as many imports as possible are in the functions so they load there, and don't slow things down.
"""
#! /bin/python3
import click
from dram2.annotate import annotate, list_databases
from dram2.rule_adjectives import  evaluate, rule_plot
from dram2.tree_kit import phylo_tree
from dram2.distill import distill
from dram2.genbank import generate_genbank
from dram2.merger import merger
from dram2.strainer import strainer
from dram2.rna import pull_trna, pull_rrna
from dram2.db_builder import db_builder
from dram2.amg_summary import amg_summary


@click.group(chain=True)
@click.option("--verbose/--quiet", default=False)
def dram2(verbose):
    pass


dram2.add_command(annotate)
dram2.add_command(list_databases)
dram2.add_command(evaluate)
dram2.add_command(rule_plot)
dram2.add_command(phylo_tree)
dram2.add_command(distill)
dram2.add_command(amg_summary)
dram2.add_command(generate_genbank)
dram2.add_command(merger)
dram2.add_command(strainer)
dram2.add_command(pull_rrna)
dram2.add_command(pull_trna)
dram2.add_command(db_builder)


if __name__ == '__main__':
    dram2()
