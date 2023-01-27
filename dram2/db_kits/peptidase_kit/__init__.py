from os import path, stat
import re
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
    get_sig_row,
    Fasta,
)

from pathlib import Path
from functools import partial
import logging
import pandas as pd

CITATION = ("N. D. Rawlings, A. J. Barrett, P. D. Thomas, X. Huang, A"
                      ". Bateman, and R. D. Finn, \"The merops database of prot"
                      "eolytic enzymes, their substrates and inhibitors in 2017"
                      " and a comparison with peptidases in the panther databas"
                      "e,\" Nucleic acids research, vol. 46, no. D1, pp. D624–D"
                      "632, 2018."
                      )

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


class PeptidaseKit(DBKit):

    name = "peptidases"
    formal_name: str = "Peptidase"
    citation: str = CITATION

    def check_setup(self):
        pass

    def search(self):
        pass

    def get_descriptions(self):
        pass

    @classmethod
    def get_ids(cls, annotatons):
        pass
