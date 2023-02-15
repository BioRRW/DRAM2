from os import path, stat
import tarfile
from shutil import move, rmtree
from dram2.db_kits.utils import (
    do_blast_style_search,
    get_basic_descriptions,
    make_mmseqs_db,
    DBKit,
    multigrep,
)
from dram2.utils.utils import Fasta
from functools import partial
import logging
import pandas as pd

VERSION = "0.1.0"
NAME = "methyl"
NAME_FORMAL = "Methyl"
CITATION = "Methyl is a in house db mostly make by McKayla Borton"
SETTINGS = {
    "search_databases": {
        "methyl_fa_db": {
            "location": None,
            "citation": CITATION,
            "name": "CAMPER FASTA db",
        },
    },
    "dram_sheets": {
        "methyl_distillate": {
            "location": None,
            "citation": CITATION,
            "name": "CAMPER Distillate form",
        }
    },
}

"""
import os
drampath  = '/home/projects-wrighton-2/DRAM/scratch_space_flynn/aug_9_methyl/snakeSpeed/results/dram'
os.system(f'dram2 annotate_genes -i {drampath}/genes.faa -a {drampath}/annotations.tsv -o ./test_output --use_db methyl')


"""
# target_db = "/home/projects-wrighton-2/DRAM/dram_data/dram1.4_final_06_07_22/methyl.mmsdb"
# multigrep(hits['%s_hit' % db_name], '%s_h' % target_db, '\x00', working_dir)
PROCESS_OPTIONS = {}


def process(methyl_fa, output_dir, logger, threads=10, verbose=False) -> dict:
    methyl_fa_db = path.join(output_dir, "methyl.mmsdb")
    make_mmseqs_db(methyl_fa, methyl_fa_db, logger)
    return {"methyl_fa_db": methyl_fa_db}


# process('/home/projects-wrighton-2/DRAM/dram_data/methyl/sep0122_methylotrophy_all_dist.faa',
#         '/home/projects-wrighton-2/DRAM/dram_data/dram1.4_final_06_07_22/', logging.getLogger())

# in the future the database will get the same input as was given in the data
def search(
    query_db: str,
    gene_faa: str,
    tmp_dir: str,
    logger: logging.Logger,
    threads: str,
    verbose: str,
    db_handler,
    bit_score_threshold,
    rbh_bit_score_threshold,
    **args,
):
    logger.info(f"Annotating genes with {NAME_FORMAL}.")


class MethylKit(DBKit):
    name = NAME
    formal_name = NAME_FORMAL
    citation: str = CITATION

    def setup(self):
        # somthing like
        # make_mmseqs_db("/home/Database/DRAM/methyl/methylotrophy.faa", "/home/Database/DRAM/methyl/methylotrophy.mmsdb", logging.getLogger())
        pass

    def load_dram_config(self):
        self.mmsdb = self.get_config_path("mmsdb")

    def search(self, fasta: Fasta):
        get_custom_description = partial(get_basic_descriptions, db_name=NAME)
        hits = do_blast_style_search(
            fasta.mmsdb.as_posix(),
            self.mmsdb.as_posix(),
            self.working_dir,
            self.logger,
            self.name,
            self.bit_score_threshold,
            self.rbh_bit_score_threshold,
            self.threads,
        )
        return hits

    def get_descriptions(self, hits):
        return hits

    @classmethod
    def get_ids(cls, annotatons):
        pass
