from os import path, stat
import re
import tarfile
import logging
from shutil import move, rmtree
from pathlib import Path
from functools import partial

import pandas as pd
from sqlalchemy import Column, String

from dram2.db_kits.utils.sql_descriptions import SQLDescriptions, BASE
from dram2.utils import download_file, run_process, Fasta
from dram2.db_kits.utils import (
    DBKit,
    do_blast_style_search,
)

CITATION = (
    "N. D. Rawlings, A. J. Barrett, P. D. Thomas, X. Huang, A"
    '. Bateman, and R. D. Finn, "The merops database of prot'
    "eolytic enzymes, their substrates and inhibitors in 2017"
    " and a comparison with peptidases in the panther databas"
    'e," Nucleic acids research, vol. 46, no. D1, pp. D624–D'
    "632, 2018."
)


class MeropsDescription(BASE):
    __tablename__ = "merops_description"

    id = Column(String(10), primary_key=True, nullable=False, index=True)

    description = Column(String(1000))

    @property
    def serialize(self):
        return {
            "merops_id": self.id,
            "merops_description": self.description,
        }


# def get_merops_description(merops_hits, header_dict):
#     merops_list: list[str] = list()
#     merops_family: list[str] = list()
#     merops_descirption: list[str] = list()
#     for merops_hit in merops_hits.merops_id:
#         header = header_dict[merops_id]
#         merops_list.append(merops_id)
#         merops_family.append(re.search(r"#\w*.#", header).group()[1:-1])
#         merops_descirption.append(header)
#     new_df = pd.DataFrame(
#         [merops_list, merops_family, merops_descirption],
#         index=["merops_id", "merops_family", "merops_hit"],
#         columns=merops_hits.index,
#     )
#     return pd.concat(
#         [new_df.transpose(), merops_hits.drop("merops_hit", axis=1)],
#         axis=1,
#         sort=False,
#     )


def get_merops_descriptions(
    hits: pd.DataFrame, header_dict: dict[str, str], db_name: str
) -> pd.DataFrame:
    """
    Get viral gene full descriptions based on headers (text before first space)
    """
    descriptions: pd.DataFrame = pd.DataFrame(hits[f"{db_name}_hit"].dropna()).rename(
        columns={f"{db_name}_hit": f"{db_name}_id"}
    )
    descriptions[f"{db_name}_description"] = descriptions[f"{db_name}_id"].apply(
        lambda x: header_dict[x]
    )
    descriptions[f"{db_name}_family"] = descriptions[f"{db_name}_description"].apply(
        lambda x: re.search(r"#\w*.#", x).group()[1:-1]
    )
    return descriptions


class MeropsKit(DBKit):
    """
    The DBKit object for peptiadse

    :attribute mmsdb:
    :attribute description_db:
    """

    name = "merops"
    formal_name: str = "MEROPS Peptidases"
    citation: str = CITATION

    def load_dram_config(self):
        self.mmsdb = self.get_config_path("mmsdb")
        self.description_db = SQLDescriptions(
            self.get_config_path("description_db"),
            self.logger,
            MeropsDescription,
            self.name,
        )

    def setup(self):
        pass

    def search(self, fasta: Fasta) -> pd.DataFrame | pd.Series:
        tmp_dir = self.working_dir / fasta.name
        tmp_dir.mkdir(exist_ok=True, parents=True)
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
            hits[f"{self.name}_hit"], f"description"
        )
        return get_merops_descriptions(hits, header_dict, self.name)

    def get_ids(self, annotations: pd.Series) -> list:
        main_id = "merops_family"
        if main_id not in annotations:
            self.logger.debug(
                f"Expected {main_id} to be in annotations,  but it was not found"
            )
        elif not pd.isna(annotations[main_id]):
            return [j for j in str(annotations[main_id]).split(";")]
        return []
