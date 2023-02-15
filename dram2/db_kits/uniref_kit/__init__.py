import re

# from dram2.utils.utils import download_file, run_process
from dram2.db_kits.utils import (
    DBKit,
    do_blast_style_search,
    get_basic_descriptions,
)
from dram2.utils.utils import Fasta

from pathlib import Path
import logging
import pandas as pd

from sqlalchemy import Column, String
from dram2.db_kits.utils.sql_descriptions import SQLDescriptions, BASE

UNIREF_CITATION = (
    "Y. Wang, Q. Wang, H. Huang, W. Huang, Y. Chen, P. B. McGarv"
    'ey, C. H. Wu, C. N. Arighi, and U. Consortium, "A crowdsour'
    'cing open platform for literature curation in uniprot," PLo'
    "S Biology, vol. 19, no. 12, p. e3001464, 2021."
)


class UniRefDescription(BASE):
    __tablename__ = "uniref_description"
    id = Column(String(40), primary_key=True, nullable=False, index=True)

    description = Column(String(1000))

    @property
    def serialize(self):
        return {
            "kegg_id": self.id,
            "kegg_description": self.description,
        }


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


class UniRefKit(DBKit):
    name = "uniref"
    formal_name: str = "UniRef"
    version: str = ""
    citation: str = UNIREF_CITATION

    def setup(self):
        pass


    def load_dram_config(self):
        self.mmsdb = self.get_config_path("mmsdb")
        self.description_db = SQLDescriptions(
            self.get_config_path("description_db"),
            self.logger,
            UniRefDescription,
            self.name,
        )

    def search(self, fasta: Fasta) -> pd.DataFrame | pd.Series:
        return do_blast_style_search(
            fasta.mmsdb,
            self.mmsdb,
            self.working_dir,
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

    @classmethod
    def get_ids(cls, annotatons):
        pass
