"""
Count Heme Regulatory Motifs
----------------------------

This
"""
from typing import Optional
from dram2.db_kits.utils import DBKit, Fasta
import pandas as pd
from skbio import read as read_sequence


def count_motifs(gene_faa, motif="(C..CH)"):
    return pd.DataFrame(
        {
            seq.metadata["id"]: len(list(seq.find_with_regex(motif)))
            for seq in read_sequence(gene_faa, format="fasta")
        },
        index=["heme_regulatory_motif_count"],
    ).T

class CountMotifs(DBKit):
    name = "motif"
    formal_name: str = "Heme Regulatory Motifs Counts"
    version: Optional[str] = None
    citation: Optional[str] = None
    date: Optional[str] = None

    def check_setup(self):
        pass

    def search(self):
        pass

    def get_descriptions(self):
        pass

    @classmethod
    def get_ids(cls, annotatons):
        pass
