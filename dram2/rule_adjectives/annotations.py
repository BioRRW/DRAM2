"""
Tool to parse the annotations file, and store it in
convenient transformations.
"""
import re
import pandas as pd
from collections import Counter
from itertools import chain
# import dask.dataframe as dd
from dram2.annotate import get_annotation_ids_by_row, DB_KITS
from dram2.db_kits.utils import DBKit

ID_SET:set[str] = {"sulfur", "kegg   ", "kofam  ", "pfam   ", "camper ", "fegenie"}
SULFUR_ID = 'sulfur_id'
FEGENIE_ID = 'fegenie_id'
COLUMN_SET= {
    'camper_id',
    FEGENIE_ID,
    SULFUR_ID,
    'kegg_genes_id',
    'nxr_nar_labels',
    'ko_id',
    'kegg_id',
    'kegg_hit',
    'peptidase_family',
    'cazy_hits',
    'cazy_subfam_ec',
    'pfam_hits',
}


def get_ids_from_annotation(frame):
    print('getting ids from annotations')
    return Counter(chain(*frame.apply(get_ids_from_row, axis=1).values))


class Annotations():

    def __init__(self, annotations_tsv:str, db_kits: list[DBKit]):
        # self.ids_by_fasta = None
        # self.ids_by_row = None
        self.data = pd.read_csv(annotations_tsv, sep='\t', index_col=0, low_memory=False)
        self.db_kits = db_kits
        self.data.set_index(['fasta', self.data.index], inplace=True)
        self.set_annotation_ids_by_row()
        self.set_annotations(annotations_tsv)

    def set_annotations(self, annotations_tsv:str):
        data = self.ids_by_row.copy()
        data['annotations'] = data['annotations'].apply(list)
        annot_fasta_ids = data.groupby('fasta')
        annot_fasta_ids = annot_fasta_ids.apply(lambda x: Counter(chain(*x['annotations'].values)))
        annot_fasta_ids = pd.DataFrame(annot_fasta_ids, columns=['annotations'])
        self.ids_by_fasta = annot_fasta_ids

    def set_annotation_ids_by_row(self):
        # self.raw_annotations = anno_data
        print("generating IDs by index, this may take some time")
        # data = self.data.copy()
        # data.set_index(['fasta', data.index], inplace=True)
        # annot_ids = get_ids_from_annotations_by_row(self.data)
        annot_ids = get_annotation_ids_by_row(self.data, db_kits=self.db_kits)
        annot_ids = pd.DataFrame(annot_ids['db_id_sets'].apply(list))
        annot_ids.columns=['annotations']
        self.ids_by_row = annot_ids

