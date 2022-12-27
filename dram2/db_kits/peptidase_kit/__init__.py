
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

def get_peptidase_description(peptidase_hits, header_dict):
    peptidase_list: list[str] = list()
    peptidase_family: list[str] = list()
    peptidase_descirption: list[str] = list()
    for peptidase_hit in peptidase_hits.peptidase_hit:
        header = header_dict[peptidase_hit]
        peptidase_list.append(peptidase_hit)
        peptidase_family.append(re.search(r"#\w*.#", header).group()[1:-1])
        peptidase_descirption.append(header)
    new_df = pd.DataFrame(
        [peptidase_list, peptidase_family, peptidase_descirption],
        index=["peptidase_id", "peptidase_family", "peptidase_hit"],
        columns=peptidase_hits.index,
    )
    return pd.concat(
        [new_df.transpose(), peptidase_hits.drop("peptidase_hit", axis=1)],
        axis=1,
        sort=False,
    )
