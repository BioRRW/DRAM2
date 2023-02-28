import re
from dram2.db_kits.utils import do_blast_style_search, DBKit
from dram2.utils.utils import Fasta
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
            "kegg_id": self.id,
            "kegg_description": self.description,
        }


class KeggKit(DBKit):

    name = "kegg"
    formal_name: str = "KEGG"
    citation: str = CITATION

    def setup(self):
        pass


    def search(self, fasta: Fasta):
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
            hits["%s_hit" % self.name], "%s_description" % self.name
        )
        gene_description = list()
        ko_list = list()
        for kegg_hit in hits.kegg_hit:
            header = header_dict[kegg_hit]
            gene_description.append(header)
            kos = re.findall(r"(K\d\d\d\d\d)", header)
            if len(kos) == 0:
                ko_list.append("")
            else:
                ko_list.append(",".join(kos))
        # TODO: change kegg_id to kegg_genes_id so that people get an error and not the wrong identifier
        new_df = pd.DataFrame(
            [hits["kegg_hit"].values, ko_list, gene_description],
            index=["kegg_genes_id", "ko_id", "kegg_hit"],
            columns=hits.index,
        )
        return pd.concat(
            [new_df.transpose(), hits.drop("kegg_hit", axis=1)], axis=1, sort=False
        )

    def get_ids(self, annotations: pd.Series) -> list:
        ko_id = "ko_id"
        kegg_id = "kegg_hit"
        ec_id = "kegg_id"
        ids = []
        if ko_id in annotations:
            ids += [j for j in annotations[ko_id].split(",")],
        if kegg_id in annotations:
            ids += [j for j in annotations[kegg_id].split(",")],
        if ec_id in annotations:
            ids += [i[1:-1] for i in re.findall(r"\[EC:\d*.\d*.\d*.\d*\]", annotations[ec_id])]
        return ids
