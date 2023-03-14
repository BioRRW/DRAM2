"""
KEGG DBKit
==========

This sub packege holds the data releated to kegg. 

"""

from os import path, stat
import tarfile
from shutil import move, rmtree
from glob import glob

from pathlib import Path
from functools import partial
import logging
import pandas as pd


from dram2.utils import download_file, run_process, merge_files, Fasta
from dram2.db_kits.utils import (
    sig_scores,
    run_hmmscan,
    DBKit,
    get_sig_row,
)

CITATION = (
    "T. Aramaki, R. Blanc-Mathieu, H. Endo, K. Ohkubo, M. Kanehisa"
    ', S. Goto, and H. Ogata, "Kofamkoala: Kegg ortholog assignme'
    'nt based on profile hmm and adaptive score threshold," Bioin'
    "formatics, vol. 36, no. 7, pp. 2251–2252, 2020."
)
def process_kofam_hmm(
        kofam_profile_tar_gz:Path,
        output_dir:Path, 
        logger: logging.Logger
        ) -> dict:
    kofam_profiles = output_dir /'kofam_profiles'
    kofam_profiles.mkdir()
    run_process(['tar',
                 '-xzf',
                 kofam_profile_tar_gz.absolute().as_posix(), 
                 '-C',
                 kofam_profiles.absolute().as_posix()],
                logger)
    merged_kofam_profiles = path.join(output_dir, 'kofam_profiles.hmm')
    merge_files(glob((kofam_profiles /'profiles'/'*.hmm').as_posix()), 
                merged_kofam_profiles)
    run_process(['hmmpress', '-f', merged_kofam_profiles], logger)
    logger.info('KOfam database processed')
    return {'kofam_hmm': merged_kofam_profiles}


def download_kofam_ko_list(
        output_dir:Path, 
        logger: logging.Logger
        ) -> Path:
    kofam_ko_list_gz = output_dir / 'kofam_ko_list.tsv.gz'
    url = 'ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz'
    url_http = 'https://www.genome.jp/ftp/db/kofam/ko_list.gz'
    download_file(url, kofam_ko_list_gz.as_posix(), logger, alt_urls=[url_http])
    return kofam_ko_list_gz

def download_kofam_hmm(output_dir:Path, logger:logging.Logger):
    kofam_profile_tar_gz = path.join(output_dir, 'kofam_profiles.tar.gz')
    url = 'ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz'
    url_http = 'https://www.genome.jp/ftp/db/kofam/profiles.tar.gz'
    download_file(url, kofam_profile_tar_gz, logger, alt_urls=[url_http])
    return kofam_profile_tar_gz



def process_kofam_ko_list(
        kofam_ko_list_gz:Path, 
        output_dir:Path, 
        logger: logging.Logger
        ) -> dict:
    # TODO: fix this so that it is gunzipped to the path
    kofam_ko_list = path.join(output_dir, 'kofam_ko_list.tsv')
    run_process(['gunzip', '-c', kofam_ko_list_gz], logger, save_output=kofam_ko_list)
    logging.info('KOfam ko list processed')
    return {'kofam_ko_list': kofam_ko_list}


def kofam_hmmscan_formater(
    hits: pd.DataFrame,
    hmm_info_path: Path,
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
    return pd.DataFrame(kegg_dict, index=["kofam_id", "kfam_hit"]).transpose()



class KOfamKit(DBKit):
    name = "kofam"
    formal_name: str = "KOfam"
    max_threads:int = 2
    citation: str = CITATION

    def search(self, fasta: Fasta) -> pd.DataFrame | pd.Series:
        return run_hmmscan(
            genes_faa=fasta.faa.as_posix(),
            db_loc=self.hmm.as_posix(),
            db_name=self.name,
            output_loc=self.working_dir.as_posix(),  # check_impliments
            threads=self.threads,  # check_impliments
            logger=self.logger,
            formater=partial(
                kofam_hmmscan_formater,
                hmm_info_path=self.kofam_ko_list,
                top_hit=True,
                use_dbcan2_thresholds=self.kofam_use_dbcan2_thresholds,
            ),
        )

    def setup(self):
        pass

    def load_dram_config(self):
        self.hmm: Path = self.get_config_path("hmmdb")
        self.kofam_ko_list: Path= self.get_config_path("kofam_ko_list")


    def get_ids(self, annotations: pd.Series) -> list:
        ids = []
        # OLD ID TO REMOVE
        ko_id = f"ko_id"
        if ko_id in annotations:
            ids += [annotations[ko_id]]
        # OLD ID TO REMOVE
        main_id = f"ko_id"
        if main_id in annotations:
            ids += [annotations[main_id]]
        return ids
