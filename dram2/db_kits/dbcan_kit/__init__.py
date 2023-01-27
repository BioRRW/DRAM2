from os import path, stat
import re
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

VERSION = '11'
DATE = '08062022'
CITATION = ("Y. Yin, X. Mao, J. Yang, X. Chen, F. Mao, and Y. Xu, \"dbcan"
                  ": a web resource for automated carbohydrate-active enzyme an"
                  "notation,\" Nucleic acids research, vol. 40, no. W1, pp. W44"
                  "5–W451, 2012."
                  )

def find_best_dbcan_hit(genome: str, group: pd.DataFrame):
    group["perc_cov"] = group.apply(
        lambda x: (x["target_end"] - x["target_start"]) / x["target_length"], axis=1
    )
    group.sort_values("perc_cov", inplace=True)
    group.columns
    group.sort_values("full_evalue", inplace=True)
    return group.iloc[0]["target_id"]


def dbcan_hmmscan_formater(hits: pd.DataFrame, db_name: str, db_handler=None):
    """
    format the ouput of the dbcan database.

    Note these changes
    introduce new column for best hit from cazy database- this will be the best hit above the already established threshold (0.35 coverage, e-18 evalue) - then for the distillate pull info from the best hit only
    introduce new column for corresponding EC number information from sub families (EC numbers are subfamily ECs)
    Make sure that ids and descriptions are separate (descriptions these are family based)
    :param hits:
    :param db_name:
    :param db_handler:
    :returns:
    """
    hits_sig = hits[hits.apply(partial(get_sig_row, evalue_lim=1e-18), axis=1)]
    if len(hits_sig) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame()
    hit_groups = hits_sig.groupby("query_id")
    all_hits = hit_groups.apply(
        lambda x: "; ".join(x["target_id"].apply(lambda y: y[:-4]).unique())
    )
    hits_df = pd.DataFrame(all_hits)
    hits_df.columns = [f"{db_name}_ids"]

    def description_pull(x: str):
        id_list = ([re.findall("^[A-Z]*[0-9]*", str(x))[0] for x in x.split("; ")],)
        id_list = [y for x in id_list for y in x if len(x) > 0]
        raise ValueError("You need to impliment descriptions")
        # description_list = db_handler.get_descriptions(
        #     id_list, "dbcan_description"
        # ).values()
        description_str = "; ".join(description_list)
        return description_str

    if db_handler is not None:
        hits_df[f"{db_name}_hits"] = hits_df[f"{db_name}_id"].apply(description_pull)
        hits_df[f"{db_name}_subfam_ec"] = hits_df[f"{db_name}_ids"].apply(
            lambda x: "; ".join(
                db_handler.get_descriptions(
                    x.split("; "), "dbcan_description", description_name="ec"
                ).values()
            )
        )
    hits_df[f"{db_name}_best_hit"] = [find_best_dbcan_hit(*i) for i in hit_groups]
    hits_df.rename_axis(None, inplace=True)
    hits_df.columns
    return hits_df


class dbCANKit(DBKit):

    name = "dbcan"
    formal_name: str = "dbCAN"
    version: str = VERSION
    citation: str = CITATION
    date: str = DATE

    def check_setup(self):
        pass

    def search(self):
        pass

    def get_descriptions(self):
        pass

    @classmethod
    def get_ids(cls, annotatons):
        pass
