from os import path, stat
import tarfile
from shutil import move, rmtree
from dram2.db_kits.fegenie_kit import process
from dram2.utils.utils import download_file, run_process
from dram2.db_kits.utils import (
    DBKit,
    get_basic_descriptions,
    BOUTFMT6_COLUMNS
)
from dram2.utils.utils import Fasta

from pathlib import Path
from functools import partial
import logging
import pandas as pd
import re

from sqlalchemy import Column, String
from dram2.db_kits.utils.sql_descriptions import SQLDescriptions, BASE


CITATION = (
    "J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G. A. Sal"
    "azar, E. L. Sonnhammer, S. C. Tosatto, L. Paladin, S. Raj, L."
    ' J. Richardson et al., "Pfam: The protein families database '
    'in 2021," Nucleic acids research, vol. 49, no. D1, pp. D412–'
    "D419, 2021."
)

def run_mmseqs_profile_search(
    query_db,
    pfam_profile,
    output_loc,
    logger,
    output_prefix="mmpro_results",
    db_handler=None,
    threads=10,
) -> pd.DataFrame:
    """Use mmseqs to run a search against pfam, currently keeping all hits and not doing any extra filtering"""
    tmp_dir = path.join(output_loc, "tmp")
    output_db = path.join(output_loc, "%s.mmsdb" % output_prefix)
    run_process(
        [
            "mmseqs",
            "search",
            query_db,
            pfam_profile,
            output_db,
            tmp_dir,
            "-k",
            "5",
            "-s",
            "7",
            "--threads",
            str(threads),
        ],
        logger,
    )
    output_loc = path.join(output_loc, "%s_output.b6" % output_prefix)
    run_process(
        ["mmseqs", "convertalis", query_db, pfam_profile, output_db, output_loc],
        logger,
    )
    pfam_results = pd.read_csv(
        output_loc, sep="\t", header=None, names=BOUTFMT6_COLUMNS
    )
    return pd.DataFrame(columns=[f"{output_prefix}_hits"])

class PfamDescription(BASE):
    __tablename__ = "pfam_description"

    id = Column(String(12), primary_key=True, nullable=False, index=True)

    description = Column(String(1000))

    @property
    def serialize(self):
        return {
            "pfam_id": self.id,
            "description": self.description,
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
        tmp_dir = self.working_dir / fasta.name
        tmp_dir.mkdir()
        return run_mmseqs_profile_search(
            fasta.mmsdb.absolute().as_posix(),
            self.mmsdb.absolute().as_posix(),
            tmp_dir.absolute().as_posix(),
            self.logger,
            output_prefix=self.name,
            db_handler=self.description_db,
            threads=self.threads,
        )

    def get_descriptions(self, hits) -> pd.DataFrame:
        self.logger.warning(f"Descriptions are not implimented for {self.name} we may not need them")
        return pd.DataFrame()
        header_dict: dict = self.description_db.get_descriptions( hits[f"{self.name}_hit"], "description")
        descriptions = get_basic_descriptions(hits, header_dict, self.name)
        pfam_dict = dict()
        for gene, pfam_frame in pfam_results.groupby("qId"):
            if len(pfam_descriptions) < 1:
                pfam_dict[gene] = "; ".join(pfam_frame.tId)
            else:
                pfam_dict[gene] = "; ".join(
                    [
                        "%s [%s]" % (pfam_descriptions[ascession], ascession)
                        for ascession in pfam_frame.tId
                    ]
                )
        return pd.DataFrame(pfam_dict, index=[f"{output_prefix}_hits"]).T

    def get_ids(self, annotations: pd.Series) -> list:
        main_id = "pfam_hits"
        if main_id not in annotations:
            self.logger.debug(f"Expected {main_id} to be in annotations,  but it was not found")
        elif not pd.isna(annotations[main_id]):
            return [
                j[1:-1].split(".")[0]
                for j in re.findall(r"\[PF\d\d\d\d\d.\d*\]", str(annotations[main_id]))
            ]
        return []
