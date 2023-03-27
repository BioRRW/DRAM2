"""
Gene Neighborhoods
___________________

"""
import logging

import click

from dram2.cli.context import DramContext, DEFAULT_KEEP_TMP, __version__


def find_neighborhoods(
    annotations, genes_from_ids, distance_bp=None, distance_genes=None
):
    # get neighborhoods as dataframes
    neighborhood_frames = list()
    for neighborhood_number, gene in enumerate(genes_from_ids):
        gene_row = annotations.loc[gene]
        scaffold_annotations = annotations.loc[
            annotations["scaffold"] == gene_row["scaffold"]
        ]
        # get neighbors based on bp
        if distance_bp is not None:
            right_dist = gene_row["end_position"] + distance_bp
            left_dist = gene_row["start_position"] - distance_bp
            neighborhood_annotations = scaffold_annotations.loc[
                (scaffold_annotations["end_position"] >= left_dist)
                & (scaffold_annotations["start_position"] <= right_dist)
            ]
        else:
            neighborhood_annotations = scaffold_annotations
        # get neighbors based on annotations
        if distance_genes is not None:
            right_genes = gene_row["gene_position"] + distance_genes
            left_genes = gene_row["gene_position"] - distance_genes
            neighborhood_annotations = scaffold_annotations.loc[
                (scaffold_annotations["gene_position"] >= left_genes)
                & (scaffold_annotations["gene_position"] <= right_genes)
            ]
        # add neighborhood number and neighborhood center as columns
        neighborhood_annotations["neighborhood_number"] = neighborhood_number
        neighborhood_annotations["neighborhood_center"] = [
            i == gene for i in neighborhood_annotations.index
        ]
        neighborhood_frames.append(neighborhood_annotations)
        if len(neighborhood_annotations) == 0:
            warnings.warn("")

    # merge data frames and write to file
    return pd.concat(neighborhood_frames)


def get_gene_neighborhoods(
    input_file,
    output_dir,
    logger: logging.Logger,
    genes=None,
    identifiers=None,
    categories=None,
    genes_loc=None,
    scaffolds_loc=None,
    distance_genes=None,
    distance_bp=None,
    custom_distillate=None,
):
    # check inputs, make output
    if distance_genes is None and distance_bp is None:
        raise ValueError("Must provide distance away in bp, genes or both.")

    # get data
    annotations = pd.read_csv(input_file, sep="\t", index_col=0)
    genes_from_ids = get_genes_from_identifiers(
        annotations,
        genes=genes,
        identifiers=identifiers,
        categories=categories,
        custom_distillate=custom_distillate,
    )
    if len(genes_from_ids) == 0:
        raise ValueError(
            "No genes were found based on your filtering parameters. No neighborhoods will be generated."
        )

    mkdir(output_dir)

    neighborhood_all_annotations = find_neighborhoods(
        annotations, genes_from_ids, distance_bp, distance_genes
    )
    neighborhood_all_annotations.to_csv(
        path.join(output_dir, "neighborhood_annotations.tsv"), sep="\t"
    )
    logging.info("Neighborhood Annotations witten to tsv")
    # filter files if given
    if genes_loc is not None:
        output_fasta_generator = (
            i
            for i in read_sequence(genes_loc, format="fasta")
            if i.metadata["id"] in neighborhood_all_annotations.index
        )
        # TODO: potentially generate one fasta file per neighborhood
        write_sequence(
            output_fasta_generator,
            format="fasta",
            into=path.join(
                output_dir, "neighborhood_genes.%s" % genes_loc.split(".")[-1]
            ),
        )
        logging.info("Gene Neighborhood fasta generated")
    if scaffolds_loc is not None:
        neighborhood_all_annotations["scaffold_mod"] = [
            "%s_%s" % (row["fasta"], row["scaffold"])
            for i, row in neighborhood_all_annotations.iterrows()
        ]
        neighborhood_scaffolds = list()
        for scaffold in read_sequence(scaffolds_loc, format="fasta"):
            if (
                scaffold.metadata["id"]
                in neighborhood_all_annotations["scaffold_mod"].values
            ):
                scaffold_frame = neighborhood_all_annotations.loc[
                    neighborhood_all_annotations["scaffold_mod"]
                    == scaffold.metadata["id"]
                ]
                for neighborhood, neighborhood_frame in scaffold_frame.groupby(
                    "neighborhood_number"
                ):
                    neighborhood_frame = neighborhood_frame.sort_values(
                        "start_position"
                    )
                    neighborhood_scaffolds.append(
                        scaffold[
                            neighborhood_frame["start_position"][
                                0
                            ]: neighborhood_frame["end_position"][-1]
                        ]
                    )
        write_sequence(
            (i for i in neighborhood_scaffolds),
            format="fasta",
            into=path.join(output_dir, "neighborhood_scaffolds.fna"),
        )
        logging.info("Scaffolds Neighborhood fasta generated")


@click.command(
    "neighbors",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.pass_context
def neighbors_cmd(
    ctx: click.Context,
):
    """
    Pull Genes based on Their Neighborhoods
    ___

    DRAM2 Can pull genes based on their proximity to other genes. I have not even written documentation for this yet.
    """
    print("This command requires more work to function in dram2")
