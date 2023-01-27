
from dram2.db_kits.utils import DBKit

VIRAL_REFSEQ_CITATION = ("J. R. Brister, D. Ako-Adjei, Y. Bao, and O. Blinkova, "
                         "\"Ncbi viral genomes resource,\" Nucleic acids researc"
                         "h, vol. 43, no. D1, pp. D571–D577, 2015. [3] M. Kanehi"
                         "sa, M. Furumichi, Y. Sato, M. Ishiguro-Watanabe, and M"
                         ". Tan-abe, \"Kegg: integrating viruses and cellular or"
                         "ganisms,\" Nucleic acids research, vol. 49, no. D1, pp"
                         ". D545–D551, 2021."
                 )

class RefSeqViral(DBKit):
    name = 'viral'
    formal_name: str = 'RefSeq-Viral'
    citation: str = VIRAL_REFSEQ_CITATION

    def check_setup(self):
        pass

    def search(self):
        pass

    def get_descriptions(self):
        pass

    @classmethod
    def get_ids(cls, annotatons):
        pass

