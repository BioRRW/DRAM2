import pytest
from dram_tree_kit

logger = logging.getLogger('dram_tree_log')
tree = DramTree(pplacer_profile='./data/dram_trees/nxr_nar/nxr_nar.refpkg', 
                reference_seq=('./data/dram_trees/nxr_nar/'
                               'nxr-nar_seqs_for_tree_aligned.faa'),
                logger=logger)
output_dir = './'

query_faa = "./../second_try/query.faa"
query_faa = '../first_try/example_one/all_bins_combined_3217db_genes.faa'

with TemporaryDirectory(dir=output_dir) as temp_dir:
    pplacerjson = tree.pplacer_place_sequences(query_faa=query_faa, working_dir=temp_dir, threads=20)
