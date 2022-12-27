
from os import path, stat
import tarfile
from shutil import move, rmtree
from dram2.db_kits.fegenie_kit import process
from dram2.utils.utils import download_file, run_process
from dram2.db_kits.utils import (
    make_mmseqs_db,
    run_hmmscan,
    get_best_hits,
    BOUTFMT6_COLUMNS,
    DBKit,
    get_sig_row, Fasta,
)

from pathlib import Path
from functools import partial
import logging
import pandas as pd

def get_uniref_description(uniref_hits, header_dict):
    """Gets UniRef ID's, taxonomy and full string from list of UniRef IDs for output in annotations"""
    gene_description = list()
    uniref_list = list()
    gene_taxonomy = list()
    for uniref_hit in uniref_hits.uniref_hit:
        header = header_dict[uniref_hit]
        gene_description.append(header)
        uniref_list.append(header[header.find("RepID=") + 6 :])
        gene_taxonomy.append(re.search(r"Tax=(.*?) (\S*?)=", header).group(1))
    new_df = pd.DataFrame(
        [uniref_list, gene_description, gene_taxonomy],
        index=["uniref_id", "uniref_hit", "uniref_taxonomy"],
        columns=uniref_hits.index,
    )
    return pd.concat(
        [new_df.transpose(), uniref_hits.drop("uniref_hit", axis=1)], axis=1, sort=False
    )
