"""
Count Heme Regulatory Motifs
----------------------------

This
"""
from typing import Optional
from dram2.db_kits.utils import DBKit
from dram2.utils.utils import Fasta
import pandas as pd
from skbio import read as read_sequence
from Bio import SeqIO

MOTIF="(C..CH)"

class CountMotifsKit(DBKit):
    name = "heme"
    formal_name: str = "Heme Regulatory Motifs Counts"
    version: Optional[str] = None
    citation:str = "This database is so simple \"(C..CH)\" it dose not warrant a citation."
    date: Optional[str] = None
    max_threads = 1

    def load_dram_config(self):
        pass
    def setup(self):
        pass


    def search(self, fasta: Fasta):
        with open(fasta.faa, 'r') as faa:
            return pd.DataFrame(
                {
                    seq.metadata["id"]: len(list(seq.find_with_regex(MOTIF)))
                    for seq in read_sequence(faa, format="fasta")
                },
                index=["heme_regulatory_motif_count"],
            ).T

    def get_descriptions(self, hits) -> pd.DataFrame:
        return hits

