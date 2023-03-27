"""

"""

from skbio import read as read_sequence
from skbio import write as write_sequence
from os import mkdir, path
import warnings
import logging

# from dram2.distill.summarize_vgfs import filter_to_amgs
import click
import pandas as pd

from dram2.cli.context import DramContext, DEFAULT_KEEP_TMP

# TODO: filter by taxonomic level, completeness, contamination
# TODO: filter scaffolds file, gff file
# TODO: add negate, aka pull not from given list
LOGGER = logging.getLogger("pull_sequences.log")


def get_genes_from_identifiers(
    annotations,
    genes=None,
    fastas=None,
    scaffolds=None,
    identifiers=None,
    categories=None,
    custom_distillate=None,
):
    specific_genes_to_keep = list()
    # filter fastas
    if fastas is not None:
        for fasta in fastas:
            specific_genes_to_keep += list(
                annotations.loc[annotations["fasta"] == fasta].index
            )
    # filter scaffolds
    # TODO: remove this functionality or modify since scaffolds are guaranteed unique
    if scaffolds is not None:
        for scaffold in scaffolds:
            specific_genes_to_keep += list(
                annotations.loc[annotations["scaffold"] == scaffold].index
            )
    # filter genes
    if genes is not None:
        specific_genes_to_keep += genes
    # filter down annotations based on specific genes
    if len(specific_genes_to_keep) > 0:
        annotations_to_keep = annotations.loc[specific_genes_to_keep]
    else:
        annotations_to_keep = annotations

    # filter based on annotations
    if (identifiers is not None) or (categories is not None):
        annotation_genes_to_keep = list()

        # make a dictionary of genes to
        gene_to_ids = dict()
        for i, row in annotations_to_keep.iterrows():
            row_ids = get_ids_from_row(row)
            if len(row_ids) > 0:
                gene_to_ids[i] = set(row_ids)

        # get genes with ids
        if identifiers is not None:
            identifiers = set(identifiers)
            for gene, ids in gene_to_ids.items():
                if len(set(ids) & set(identifiers)) > 0:
                    annotation_genes_to_keep.append(gene)

        # get genes from distillate categories
        # if categories is not None:
        #     database_handler = DatabaseHandler()
        #     genome_summary_form = pd.read_csv(
        #         database_handler.db_locs["genome_summary_form"], sep="\t"
        #     )
        #     if custom_distillate is not None:
        #         genome_summary_form = pd.concat(
        #             [genome_summary_form, pd.read_csv(custom_distillate, sep="\t")]
        #         )
        #     for level in ["module", "sheet", "header", "subheader"]:
        #         for category, frame in genome_summary_form.loc[
        #             ~pd.isna(genome_summary_form[level])
        #         ].groupby(level):
        #             if category in categories:
        #                 for gene, ids in gene_to_ids.items():
        #                     if len(ids & set(frame["gene_id"])) > 0:
        #                         annotation_genes_to_keep.append(gene)
    else:
        annotation_genes_to_keep = list(annotations_to_keep.index)
    return annotation_genes_to_keep


def strainer(
    input_annotations,
    input_fasta,
    output_fasta,
    fastas=None,
    scaffolds=None,
    genes=None,
    identifiers=None,
    categories=None,
    taxonomy=None,
    completeness=None,
    contamination=None,
    amg_flags=None,
    aux_scores=None,
    virsorter_category=None,
    putative_amgs=False,
    max_auxiliary_score=3,
    remove_transposons=False,
    remove_fs=False,
    custom_distillate=None,
):
    annotations = pd.read_csv(input_annotations, sep="\t", index_col=0)

    # first filter based on specific names
    annotation_genes_to_keep = get_genes_from_identifiers(
        annotations,
        genes,
        fastas,
        scaffolds,
        identifiers,
        categories,
        custom_distillate,
    )
    annotations = annotations.loc[annotation_genes_to_keep]
    if len(annotations) == 0:
        raise ValueError("Categories or identifiers provided yielded no annotations")

    # DRAM specific filtering
    if taxonomy is not None:
        taxonomy = set(taxonomy)
        annotations = annotations.loc[
            [len(set(i.split(";")) & taxonomy) > 0 for i in annotations["bin_taxonomy"]]
        ]
    if completeness is not None:
        annotations = annotations.loc[
            annotations["bin_completeness"].astype(float) > completeness
        ]
    if contamination is not None:
        annotations = annotations.loc[
            annotations["bin_contamination"].astype(float) < contamination
        ]
    if len(annotations) == 0:
        raise ValueError("DRAM filters yielded no annotations")

    # DRAM-v specific filtering
    # AMGnonsense if putative_amgs:  # get potential amgs
    # AMGnonsense     # annotations = filter_to_amgs(annotations.fillna(''), max_aux=max_auxiliary_score,
    # AMGnonsense     #                              remove_transposons=remove_transposons, remove_fs=remove_fs, remove_js=remove_js)
    # AMGnonsense     annotations = filter_to_amgs(
    # AMGnonsense         annotations.fillna(""),
    # AMGnonsense         max_aux=max_auxiliary_score,
    # AMGnonsense         remove_transposons=remove_transposons,
    # AMGnonsense         remove_fs=remove_fs,
    # AMGnonsense     )
    # AMGnonsense else:
    # AMGnonsense     # filter based on virsorter categories
    # AMGnonsense     if virsorter_category is not None:
    # AMGnonsense         annotations = annotations.loc[
    # AMGnonsense             [i in virsorter_category for i in annotations.virsorter]
    # AMGnonsense         ]
    # AMGnonsense     # filter based on aux scores
    # AMGnonsense     if aux_scores is not None:
    # AMGnonsense         annotations = annotations.loc[
    # AMGnonsense             [i in aux_scores for i in annotations.auxiliary_score]
    # AMGnonsense         ]
    # AMGnonsense     # filter based on amg flags
    # AMGnonsense     if amg_flags is not None:
    # AMGnonsense         amg_flags = set(amg_flags)
    # AMGnonsense         annotations = annotations.loc[
    # AMGnonsense             [
    # AMGnonsense                 len(set(i) & amg_flags) > 0 if not pd.isna(i) else False
    # AMGnonsense                 for i in annotations.amg_flags
    # AMGnonsense             ]
    # AMGnonsense         ]
    if len(annotations) == 0:
        raise ValueError("DRAM-v filters yielded no annotations")

    # make output
    output_fasta_generator = (
        i
        for i in read_sequence(input_fasta, format="fasta")
        if i.metadata["id"] in annotations.index
    )
    write_sequence(output_fasta_generator, format="fasta", into=output_fasta)


@click.command(
    "strainer",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.pass_context
def strainer_cmd(
    ctx: click.Context,
):
    """
    Strain to Genes of Interest (Not Ready)
    ___

    After you have completed your annotation and distillation, you may want to further analyze genes of interest by making trees or functional modeling. To pull the genes you can use DRAM.py strainer. For example, if you want to pull all pmoa/amoa genes based on KEGG annotations:

    DRAM.py strainer --identifiers K10944 -i annotations.tsv -f genes.faa -o amoa_pmoa_genes.faa

    Or you might want to blast a few specific genes:

    DRAM.py strainer --genes bin.2_scaffold_2_3 bin.4_scaffold_12_42 -i annotations.tsv -f genes.dna -on my_genes.fna

    Or maybe you only want to see genes that are involved in glycolysis or the TCA cycle that are from bins from the Roseburia genus:

    DRAM.py strainer -i hmp_bins/annotations.tsv -f hmp_bins/genes.fna -o genes.roseburia.glycoloysis_tca.fna --taxonomy g__Roseburia --categories glycolysis TCA

    DRAM strainer parameters

        Input files:
            Annotations: annotations.tsv file generated during the annotate step
            Input fasta: genes fasta file (.faa or .fna) to be filtered
            Output fasta: location to save filtered fasta file
        Default Parameters:
            Fastas: None
            Scaffolds: None
            Genes: None
            Identifiers: None

    """
    print("This command requires more work to work in dram2")

    # strainer(
    #     input_annotations,
    #     input_fasta,
    #     output_fasta,
    #     fastas=None,
    #     scaffolds=None,
    #     genes=None,
    #     identifiers=None,
    #     categories=None,
    #     taxonomy=None,
    #     completeness=None,
    #     contamination=None,
    #     amg_flags=None,
    #     aux_scores=None,
    #     virsorter_category=None,
    #     putative_amgs=False,
    #     max_auxiliary_score=3,
    #     remove_transposons=False,
    #     remove_fs=False,
    #     custom_distillate=None,
    # )
