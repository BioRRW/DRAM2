from os import path, stat
import re
import tarfile
from shutil import move, rmtree
from dram2.utils import download_file, run_process, make_mmseqs_db, \
    run_hmmscan, get_best_hits, BOUTFMT6_COLUMNS, do_blast_style_search
from functools import partial
import logging
import pandas as pd

VERSION = '1.0.0-beta.1'
NAME = 'kegg'
NAME_FORMAL = 'KEGG'
CITATION = "CAMPER has no citeation and is in beta so you should not be using it."
SETTINGS = { 
    'search_databases':{
    "kegg": {
      "location": None,
      "name": "KEGG db",
      "description_db_updated": "Unknown, or Never",
      "citation": "Kanehisa, M., Furumichi, M., Sato, Y., Ishiguro-Watanabe, M., and Tanabe, M.; KEGG: integrating viruses and cellular organisms. Nucleic Acids Res. 49, D545-D551 (2021)."
    }
    }
}
# the format is input file: options
DOWNLOAD_OPTIONS = {'camper_tar_gz': {'version': VERSION}}
PROCESS_OPTIONS = {'camper_tar_gz': {'version': VERSION}}


def get_kegg_description(kegg_hits, header_dict):
    """Gets the KEGG IDs, and full KEGG hits from list of KEGG IDs for output in annotations"""
    gene_description = list()
    ko_list = list()
    for kegg_hit in kegg_hits.kegg_hit:
        header = header_dict[kegg_hit]
        gene_description.append(header)
        kos = re.findall(r'(K\d\d\d\d\d)', header)
        if len(kos) == 0:
            ko_list.append('')
        else:
            ko_list.append(','.join(kos))
    # TODO: change kegg_id to kegg_genes_id so that people get an error and not the wrong identifier
    new_df = pd.DataFrame([kegg_hits['kegg_hit'].values, ko_list, gene_description],
                          index=['kegg_genes_id', 'ko_id', 'kegg_hit'], columns=kegg_hits.index)
    return pd.concat([new_df.transpose(), kegg_hits.drop('kegg_hit', axis=1)], axis=1, sort=False)


def search(query_db:str, gene_faa:str, tmp_dir:str, logger:logging.Logger, 
           threads:str, verbose:str, db_handler, bit_score_threshold, 
           rbh_bit_score_threshold, **args, ):
    logger.info(f"Annotating genes with {NAME_FORMAL}.")
    hits = do_blast_style_search(query_db, 
                          db_handler.config["search_databases"]['kegg']['location'], 
                          tmp_dir,
                          db_handler, 
                          get_kegg_description, 
                          logger,
                          'kegg', 
                          bit_score_threshold, 
                          rbh_bit_score_threshold, 
                          threads)
    return hits
