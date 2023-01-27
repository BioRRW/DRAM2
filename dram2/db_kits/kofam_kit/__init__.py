
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
    get_sig_row, Fasta,
)

from pathlib import Path
from functools import partial
import logging
import pandas as pd

KOFAM_CITATION = ("T. Aramaki, R. Blanc-Mathieu, H. Endo, K. Ohkubo, M. Kanehisa"
                  ", S. Goto, and H. Ogata, \"Kofamkoala: Kegg ortholog assignme"
                  "nt based on profile hmm and adaptive score threshold,\" Bioin"
                  "formatics, vol. 36, no. 7, pp. 2251–2252, 2020."
                  )

def kofam_hmmscan_formater(
    hits: pd.DataFrame,
    hmm_info_path: str = None,
    use_dbcan2_thresholds: bool = False,
    top_hit: bool = True,
):
    hmm_info = pd.read_csv(hmm_info_path, sep="\t", index_col=0)
    if use_dbcan2_thresholds:
        hits_sig = hits[hits.apply(get_sig_row, axis=1)]
    else:
        hits_sig = sig_scores(hits, hmm_info)
    # if there are any significant results then parse to dataframe
    if len(hits_sig) == 0:
        return pd.DataFrame()
    kegg_dict = dict()
    for gene, frame in hits_sig.groupby("query_id"):
        # TODO: take top hit for full length genes and all hits for domains?
        # TODO: if top hit then give all e-value and bitscore info
        if top_hit:
            best_hit = frame[frame.full_evalue == frame.full_evalue.min()]
            ko_id = best_hit["target_id"].iloc[0]
            kegg_dict[gene] = [ko_id, hmm_info.loc[ko_id, "definition"]]
        else:
            kegg_dict[gene] = [
                ",".join([i for i in frame.target_id]),
                "; ".join([hmm_info.loc[i, "definition"] for i in frame.target_id]),
            ]
    return pd.DataFrame(kegg_dict, index=["ko_id", "kegg_hit"]).transpose()

class KOfamKit(DBKit):
    pass
