from os import path, stat
import tarfile
from shutil import move, rmtree
from mag_annotator.utils import do_blast_style_search, get_basic_description, \
    run_process, make_mmseqs_db, multigrep, process_reciprocal_best_hits
from functools import partial
import logging
import pandas as pd
# get_basic_description('/home/projects-wrighton-2/DRAM/scratch_space_flynn/aug_9_methyl/test_output/dram2/final/working_dir/genes/gene_methyl.minbitscore60.tophit.swapped.mmsdb', db_name='methyl')
VERSION = '0.1.0'
NAME = 'methyl'
NAME_FORMAL = 'Methyl'
CITATION = "Methyl is a in house db mostly make by McKayla Borton"
SETTINGS = { 
    "search_databases": {
        'methyl_fa_db': {'location': None, 'citation': CITATION, 'name': 'CAMPER FASTA db'},
    },
  "dram_sheets": {
    'methyl_distillate': {'location': None,'citation': CITATION, 'name': 'CAMPER Distillate form'}
  }
}

'''
import os
drampath  = '/home/projects-wrighton-2/DRAM/scratch_space_flynn/aug_9_methyl/snakeSpeed/results/dram'
os.system(f'dram2 annotate_genes -i {drampath}/genes.faa -a {drampath}/annotations.tsv -o ./test_output --use_db methyl')


'''
# target_db = "/home/projects-wrighton-2/DRAM/dram_data/dram1.4_final_06_07_22/methyl.mmsdb"
# multigrep(hits['%s_hit' % db_name], '%s_h' % target_db, '\x00', working_dir)
PROCESS_OPTIONS = {}

def process(methyl_fa, output_dir, logger, threads=10, verbose=False) -> dict:
    methyl_fa_db = path.join(output_dir, 'methyl.mmsdb')
    make_mmseqs_db(methyl_fa, methyl_fa_db, logger)
    return {'methyl_fa_db': methyl_fa_db}
# process('/home/projects-wrighton-2/DRAM/dram_data/methyl/sep0122_methylotrophy_all_dist.faa',
#         '/home/projects-wrighton-2/DRAM/dram_data/dram1.4_final_06_07_22/', logging.getLogger())

# in the future the database will get the same input as was given in the data
def search(query_db:str, gene_faa:str, tmp_dir:str, logger:logging.Logger, 
           threads:str, verbose:str, db_handler, bit_score_threshold, rbh_bit_score_threshold, **args):
     logger.info(f"Annotating genes with {NAME_FORMAL}.")
     get_custom_description = partial(get_basic_description, db_name=NAME)
     hits= do_blast_style_search( query_db, db_handler.config["search_databases"]['methyl_fa_db']['location'], 
                           tmp_dir, db_handler, get_custom_description, logger, NAME,
         bit_score_threshold, rbh_bit_score_threshold, threads,
         verbose)
     return hits 

