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
    get_sig_row,
    Fasta,
)

from pathlib import Path
from functools import partial
import logging
import pandas as pd

from sqlalchemy import Column, String
from dram2.db_kits.utils.sql_descriptions import SQLDescriptions, BASE



CITATION = ("J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G. A. Sal"
                 "azar, E. L. Sonnhammer, S. C. Tosatto, L. Paladin, S. Raj, L."
                 " J. Richardson et al., \"Pfam: The protein families database "
                 "in 2021,\" Nucleic acids research, vol. 49, no. D1, pp. D412–"
                 "D419, 2021."
                 )

class PfamDescription(BASE):
    __tablename__ = 'pfam_description'

    id = Column(String(12), primary_key=True, nullable=False, index=True)

    description = Column(String(1000))

    @property
    def serialize(self):
        return {
            'pfam_id': self.id,
            'pfam_description': self.description,
        }

class PfamKit(DBKit):
    name = "pfam"
    formal_name: str = "Pfam"
    version: str = ""
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
