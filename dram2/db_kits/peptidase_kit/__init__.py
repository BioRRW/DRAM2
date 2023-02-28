from os import path, stat
import re
import tarfile
from shutil import move, rmtree

from dram2.db_kits.fegenie_kit import process
from dram2.utils.utils import download_file, run_process
from dram2.db_kits.utils import (
    DBKit,
    do_blast_style_search,
    get_basic_descriptions,
)
from dram2.utils.utils import Fasta
from pathlib import Path
from functools import partial
import logging
import pandas as pd

from sqlalchemy import Column, String
from dram2.db_kits.utils.sql_descriptions import SQLDescriptions, BASE


CITATION = (
    "N. D. Rawlings, A. J. Barrett, P. D. Thomas, X. Huang, A"
    '. Bateman, and R. D. Finn, "The merops database of prot'
    "eolytic enzymes, their substrates and inhibitors in 2017"
    " and a comparison with peptidases in the panther databas"
    'e," Nucleic acids research, vol. 46, no. D1, pp. D624–D'
    "632, 2018."
)


class PeptidaseDescription(BASE):
    __tablename__ = "peptidase_description"

    id = Column(String(10), primary_key=True, nullable=False, index=True)

    description = Column(String(1000))

    @property
    def serialize(self):
        return {
            "peptidase_id": self.id,
            "peptidase_description": self.description,
        }


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

    def load_dram_config(self):
        self.mmsdb = self.get_config_path("mmsdb")
        self.description_db = SQLDescriptions(
            self.get_config_path("description_db"),
            self.logger,
            PeptidaseDescription,
            self.name,
        )
    def setup(self):
        pass


    def search(self, fasta: Fasta) -> pd.DataFrame | pd.Series:
        tmp_dir = self.working_dir / fasta.name
        tmp_dir.mkdir()
        return do_blast_style_search(
            fasta.mmsdb,
            self.mmsdb,
            tmp_dir,
            self.logger,
            self.name,
            self.bit_score_threshold,
            self.rbh_bit_score_threshold,
            self.threads,
        )

    def get_descriptions(self, hits) -> pd.DataFrame:
        header_dict = self.description_db.get_descriptions(
            hits[f"{self.name}_hit"], f"{self.name}_description"
        )
        return get_basic_descriptions(hits, header_dict, self.name)

    def get_ids(self, annotations: pd.Series) -> list:
        main_id = "peptidase_family"
        if main_id in annotations:
            return  [j for j in annotations[main_id].split(";")]
        return []
