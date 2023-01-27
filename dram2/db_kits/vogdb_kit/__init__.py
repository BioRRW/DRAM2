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

VOGDB_CITATION = (
    "J. Thannesberger, H.-J. Hellinger, I. Klymiuk, M.-T. Kastner"
    ", F. J. Rieder, M. Schneider, S. Fister, T. Lion, K. Kosulin"
    ', J. Laengle et al., "Viruses comprise an extensive pool of'
    " mobile genetic elements in eukaryote cell cultures and huma"
    'n clinical samples," The FASEB Journal, vol. 31, no. 5, pp.'
    " 1987–2000, 2017."
)


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


class VogDB(DBKit):
    name = "vogdb"
    formal_name: str = "VogDB"
    citation: str = VOGDB_CITATION

    def check_setup(self):
        pass

    def search(self):
        pass

    def get_descriptions(self):
        pass

    @classmethod
    def get_ids(cls, annotatons):
        pass
