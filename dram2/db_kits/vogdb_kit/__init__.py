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
)
from dram2.utils.utils import Fasta
from pathlib import Path
from functools import partial
import logging
import pandas as pd

from sqlalchemy import Column, String
from dram2.db_kits.utils.sql_descriptions import SQLDescriptions, BASE

VOGDB_CITATION = (
    "J. Thannesberger, H.-J. Hellinger, I. Klymiuk, M.-T. Kastner"
    ", F. J. Rieder, M. Schneider, S. Fister, T. Lion, K. Kosulin"
    ', J. Laengle et al., "Viruses comprise an extensive pool of'
    " mobile genetic elements in eukaryote cell cultures and huma"
    'n clinical samples," The FASEB Journal, vol. 31, no. 5, pp.'
    " 1987–2000, 2017."
)

class VOGDBDescription(BASE):
    __tablename__ = 'vogdb_description'

    id = Column(String(10), primary_key=True, nullable=False, index=True)

    description = Column(String(1000))

    @property
    def serialize(self):
        return {
            'vogdb_id': self.id,
            'vogdb_description': self.description,
        }

def vogdb_hmmscan_formater(hits: pd.DataFrame, db_name: str, db_handler=None):
    hits_sig = hits[hits.apply(get_sig_row, axis=1)]
    if len(hits_sig) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame()
    # Get the best hits
    hits_best = hits_sig.sort_values("full_evalue").drop_duplicates(subset=["query_id"])
    if db_handler is None:
        hits_df = hits_best[["target_id", "query_id"]].copy()
    else:
        # get_descriptions
        desc_col = f"{db_name}_hits"
        descriptions = pd.DataFrame(
            db_handler.get_descriptions(
                hits_best["target_id"].unique(), f"{db_name}_description"
            ),
            index=[desc_col],
        ).T
        categories = descriptions[desc_col].apply(lambda x: x.split("; ")[-1])
        descriptions[f"{db_name}_categories"] = categories.apply(
            lambda x: ";".join(set([x[i : i + 2] for i in range(0, len(x), 2)]))
        )
        descriptions["target_id"] = descriptions.index
        hits_df = pd.merge(
            hits_best[["query_id", "target_id"]], descriptions, on=f"target_id"
        )
    hits_df.set_index("query_id", inplace=True, drop=True)
    hits_df.rename_axis(None, inplace=True)
    hits_df.rename(columns={"target_id": f"{db_name}_id"}, inplace=True)
    return hits_df


class VogDBKit(DBKit):
    name = "vogdb"
    formal_name: str = "VogDB"
    citation: str = VOGDB_CITATION

    def setup(self):
        pass


    def load_dram_config(self):
        self.hmm = self.get_config_path("hmmdb")
        self.description_db = SQLDescriptions(
            self.get_config_path("description_db"),
            self.logger,
            VOGDBDescription,
            self.name,
        )

    def search(self, fasta: Fasta)-> pd.DataFrame| pd.Series:
        run_hmmscan(
            genes_faa=fasta.faa.as_posix(),
            db_loc=self.hmm.as_posix(),
            db_name=self.name,
            threads=self.threads,
            output_loc=self.working_dir.as_posix(),
            formater=partial(
                vogdb_hmmscan_formater,
                db_name=self.name,
                db_handler=self.description_db,
            ),
            logger=self.logger,
        )

    def get_descriptions(self, hits):
        "fix"
        return hits

    @classmethod
    def get_ids(cls, annotatons):
        pass
