from dram2.utils.tests.test_utils import annotations
from dram2.strainer.pull_sequences import get_genes_from_identifiers

def test_get_genes_from_identifiers(annotations):
    # filter by genes
    genes_to_keep = ["gene1", "gene4", "gene8"]
    kept_genes = get_genes_from_identifiers(annotations, genes=genes_to_keep)
    assert set(kept_genes) == set(genes_to_keep)
    # filter by scaffolds
    scaffolds_to_keep = ["scaffold_2"]
    kept_genes = get_genes_from_identifiers(annotations, scaffolds=scaffolds_to_keep)
    assert set(kept_genes) == {"gene5"}
    # filter by bin
    bins_to_keep = ["bin.2"]
    kept_genes = get_genes_from_identifiers(annotations, fastas=bins_to_keep)
    assert set(kept_genes) == {"gene6", "gene7", "gene8"}
    # filter by bin and scaffold
    kept_genes = get_genes_from_identifiers(
        annotations, fastas=["bin.2"], scaffolds=["scaffold_2"]
    )
    assert set(kept_genes) == {"gene5", "gene6", "gene7", "gene8"}
    # filter by ko
    kept_genes = get_genes_from_identifiers(annotations, identifiers=["K00001"])
    assert set(kept_genes) == {"gene2", "gene8"}
    # filter by bin and ko
    kept_genes = get_genes_from_identifiers(
        annotations, fastas=["bin.1"], identifiers=["K00001"]
    )
    assert set(kept_genes) == {"gene2"}
