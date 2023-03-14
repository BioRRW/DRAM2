import re

from dram2.db_kits.utils import do_blast_style_search, DBKit
from dram2.utils import Fasta
import logging
import pandas as pd

CITATION = (
    " M. Kanehisa, M. Furumichi, Y. Sato, M. Ishiguro-Watanabe, and"
    ' M. Tanabe, "Kegg: integrating viruses and cellular organisms'
    '," Nucleic acids research, vol. 49, no. D1, pp. D545–D551, 20'
    "21."
)

from sqlalchemy import Column, String
from dram2.db_kits.utils.sql_descriptions import SQLDescriptions, BASE


class KeggDescription(BASE):
    __tablename__ = "kegg_description"

    id = Column(String(20), primary_key=True, nullable=False, index=True)

    description = Column(String(100000))

    @property
    def serialize(self):
        return {
            "id": self.id,
            "description": self.description,
        }


class KeggKit(DBKit):

    name = "kegg"
    formal_name: str = "KEGG"
    citation: str = CITATION

    def setup(self):
        pass

    def search(self, fasta: Fasta):
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

    def load_dram_config(self):
        self.mmsdb = self.get_config_path("mmsdb")
        self.description_db = SQLDescriptions(
            self.get_config_path("description_db"),
            self.logger,
            KeggDescription,
            self.name,
        )

    def get_descriptions(self, hits):
        """
        Gets the KEGG IDs, and full KEGG hits from list of KEGG IDs for output in annotations.
        """
        header_dict = self.description_db.get_descriptions(
            hits["%s_hit" % self.name], "description"
        )
        new_df = pd.DataFrame(hits.kegg_hit.dropna()).rename(
            columns={"kegg_hit": "kegg_genes_id"}
        )
        new_df["kegg_description"] = new_df["kegg_genes_id"].apply(
            lambda x: header_dict[x]
        )
        new_df["kegg_id"] = new_df["kegg_description"].apply(
            lambda x: ",".join(re.findall(r"(K\d\d\d\d\d)", x))
        )
        return new_df

    def get_ids(self, annotations: pd.Series) -> list:
        ko_id = "kegg_id"
        ec_id = "kegg_description"
        ids = []
        # OLD ID TO REMOVE
        if ko_id not in annotations:
            self.logger.debug(f"Expected {ko_id} to be in annotations but not found")
        if ko_id in annotations and not pd.isna(annotations[ko_id]):
            ids += [j for j in str(annotations[ko_id]).split(",")]
        if ec_id in annotations and not pd.isna(annotations[ec_id]):
            ids += [
                i[1:-1]
                for i in re.findall(r"\[EC:\d*.\d*.\d*.\d*\]", str(annotations[ec_id]))
            ]
        return ids


"/home/projects-wrighton-2/DRAM/development_flynn/dram2_dev/jan_26_23_main_pipeline/dram_db/kegg.sqlite"
"/home/projects-wrighton-2/DRAM/development_flynn/dram2_dev/jan_26_23_main_pipeline/split_sql_dbs/kegg.sqlite"
"""
from pathlib import Path

conf = {
    "dram_data_folder": Path("/home/Database/DRAM/sep_12_22_dram1.4.0_full_db/"),
    "db_kits":{
        "kegg": {
            "mmsdb": {'location': Path("kegg.20221012.mmsdb")},
            'description_db':{ 'location': Path("/home/projects-wrighton-2/DRAM/development_flynn/dram2_dev/jan_26_23_main_pipeline/split_sql_dbs/kegg.sqlite")}}}}


data = pd.read_csv("soil/test1/annotations.tsv", sep='\t', index_col=0)


test = KeggKit(conf, logging.getLogger())
test.load_dram_config()
test.get_descriptions(data)
test.get_descriptions(pd.DataFrame({'kegg_hit': ['chu:CHU_1316']}))

data.columns
test.description_db.get_descriptions(data["kegg_hit"], "description")
"""
