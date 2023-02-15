from os import path, stat
import tarfile
from shutil import move, rmtree
from dram2.db_kits.fegenie_kit import process
from dram2.utils.utils import download_file, run_process
from dram2.db_kits.utils import (
    DBKit,
    run_mmseqs_profile_search,
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
    "J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G. A. Sal"
    "azar, E. L. Sonnhammer, S. C. Tosatto, L. Paladin, S. Raj, L."
    ' J. Richardson et al., "Pfam: The protein families database '
    'in 2021," Nucleic acids research, vol. 49, no. D1, pp. D412–'
    "D419, 2021."
)


class PfamDescription(BASE):
    __tablename__ = "pfam_description"

    id = Column(String(12), primary_key=True, nullable=False, index=True)

    description = Column(String(1000))

    @property
    def serialize(self):
        return {
            "pfam_id": self.id,
            "pfam_description": self.description,
        }


class PfamKit(DBKit):
    name = "pfam"
    formal_name: str = "Pfam"
    version: str = ""
    citation: str = CITATION

    def setup(self):
        pass

    def load_dram_config(self):
        self.mmsdb = self.get_config_path("mmsdb")
        self.description_db = SQLDescriptions(
            self.get_config_path("description_db"),
            self.logger,
            PfamDescription,
            self.name,
        )

    def search(self, fasta: Fasta):
        run_mmseqs_profile_search(
            fasta.mmsdb,
            self.mmsdb,
            self.working_dir,
            self.logger,
            output_prefix=self.name,
            db_handler=self.description_db,
            threads=self.threads,
        )

    def get_descriptions(self, hits) -> pd.DataFrame:
        header_dict = self.description_db.get_descriptions(
            hits[f"{self.name}_hit"], f"{self.name}_description"
        )
        return get_basic_descriptions(hits, header_dict, self.name)

    @classmethod
    def get_ids(cls, annotatons):
        pass
