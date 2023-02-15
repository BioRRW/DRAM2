
import pandas as pd
from dram2.db_kits.utils import DBKit, get_basic_descriptions, do_blast_style_search
from dram2.utils.utils import Fasta
from sqlalchemy import Column, String
from dram2.db_kits.utils.sql_descriptions import SQLDescriptions, BASE

VIRAL_REFSEQ_CITATION = ("J. R. Brister, D. Ako-Adjei, Y. Bao, and O. Blinkova, "
                         "\"Ncbi viral genomes resource,\" Nucleic acids researc"
                         "h, vol. 43, no. D1, pp. D571–D577, 2015. [3] M. Kanehi"
                         "sa, M. Furumichi, Y. Sato, M. Ishiguro-Watanabe, and M"
                         ". Tan-abe, \"Kegg: integrating viruses and cellular or"
                         "ganisms,\" Nucleic acids research, vol. 49, no. D1, pp"
                         ". D545–D551, 2021."
                 )

class ViralDescription(BASE):
    __tablename__ = 'viral_description'
    id = Column(String(14), primary_key=True, nullable=False, index=True)

    description = Column(String(1000))

    @property
    def serialize(self):
        return {
            'viral_id': self.id,
            'viral_description': self.description,
        }
class RefSeqViralKit(DBKit):
    name = 'viral'
    formal_name: str = 'RefSeq-Viral'
    citation: str = VIRAL_REFSEQ_CITATION

    def setup(self):
        pass

    def load_dram_config(self):
        self.mmsdb = self.get_config_path("mmsdb")
        self.description_db = SQLDescriptions(
            self.get_config_path("description_db"),
            self.logger,
            ViralDescription,
            self.name,
        )

    def search(self, fasta: Fasta)-> pd.DataFrame| pd.Series:
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

    def get_descriptions(self, hits)-> pd.DataFrame:
        header_dict = self.description_db.get_descriptions(
            hits[f"{self.name}_hit"], f"{self.name}_description"
        )
        return get_basic_descriptions(hits, header_dict, self.name)

    @classmethod
    def get_ids(cls, annotatons):
        pass

