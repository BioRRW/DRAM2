
def test_extract_enigmatic_genes(temp_dir, logger):
    temp_dir = "../temp"
    annotations = pd.DataFrame({"ko_id": ["K00001", "K00002"]}, index=["one", "two"])
    annotation_ids = get_ids_from_annotations_by_row(annotations, logger)
    start_fa = os.path.join(temp_dir, "in.faa")
    end_fa = os.path.join(temp_dir, "trim.faa")
    with open(start_fa, "w") as out:
        out.write(">one\naabb\n>two\nbbaa")
    extract_enigmatic_genes(annotations, start_fa, temp_dir, {"K00001"}, logger)
    with open(end_fa, "r") as file:
        result = file.read()
    assert result == ">one\naabb\n"




"""
import os
os.system('test_env/bin/python3 dram_tree_kit/dram_phylo_tree.py  -f -a tests/data/mini_nxr_nar_annotations.tsv -g tests/data/mini_nxr_nar_genes.faa -c 30')
dram_tree_kit()
"""
