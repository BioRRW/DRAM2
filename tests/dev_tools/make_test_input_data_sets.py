
from Bio import SeqIO


with open("/home/Database/KEGG/kegg-all-orgs_20220129/kegg-all-orgs_unique_reheader_20220129.pep", "r") as faa_feed, open("./tests/dev_tools/nxr_nar_tree_fasta.faa", "w") as out_faa:
    output_fasta_generator = ( i for i in SeqIO.parse(faa_feed, "fasta") if ("K11180" in i.description) or ("K11181" in i.description) )
    SeqIO.write(output_fasta_generator, out_faa, "fasta")


