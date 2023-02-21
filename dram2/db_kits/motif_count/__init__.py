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

MOTIF="(C..CH)"

class CountMotifsKit(DBKit):
    name = "motif"
    formal_name: str = "Heme Regulatory Motifs Counts"
    version: Optional[str] = None
    citation: Optional[str] = None
    date: Optional[str] = None

    def load_dram_config(self):
        pass
    def setup(self):
        pass


    def search(self, fasta: Fasta):
        return pd.DataFrame(
            {
                seq.metadata["id"]: len(list(seq.find_with_regex(MOTIF)))
                for seq in read_sequence(fasta.faa, format="fasta")
            },
            index=["heme_regulatory_motif_count"],
        ).T

    def get_descriptions(self, hits) -> pd.DataFrame:
        return hits

