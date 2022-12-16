import click

__version__ = '2.0.0'

@click.command('merge')
@click.version_option(__version__)
def merger():
    print("This comand is comming soon")


def merge_gtdb_taxonomy(annotations, gtdb_taxonomy):
    gtdb_taxonomy = pd.concat(
        [pd.read_csv(i, sep="\t", index_col=0) for i in gtdb_taxonomy]
    )
    taxonomy = list()
    taxonomy_missing_bins = list()
    for i in annotations.fasta:
        # add taxonomy
        if i in gtdb_taxonomy.index:
            taxonomy.append(gtdb_taxonomy.loc[i, "classification"])
        else:
            taxonomy.append(i)
            taxonomy_missing_bins.append(i)
    for i in set(taxonomy_missing_bins):
        logger.warning(
            "Bin %s was not found in taxonomy file, replaced with bin name." % i
        )
    annotations["bin_taxonomy"] = taxonomy


def merge_annotations(annotations_list, output_dir, write_annotations=False):
    # merge annotation dicts
    all_annotations = pd.concat(
        [i.get_annotations() for i in annotations_list if i is not None], sort=False
    )
    all_annotations = all_annotations.sort_values(
        ["fasta", "scaffold", "gene_position"]
    )

    # merge gene files
    merge_files(
        [i.genes_fna_loc for i in annotations_list if i is not None],
        path.join(output_dir, "genes.fna"),
    )
    merge_files(
        [i.genes_faa_loc for i in annotations_list if i is not None],
        path.join(output_dir, "genes.faa"),
    )
    merge_files(
        [i.scaffolds_loc for i in annotations_list if i is not None],
        path.join(output_dir, "scaffolds.fna"),
    )
    merge_files(
        [i.gff_loc for i in annotations_list if i is not None],
        path.join(output_dir, "genes.gff"),
        True,
    )
    trnas_locs = [
        i.trnas_loc
        for i in annotations_list
        if i is not None
        if i.trnas_loc is not None
    ]
    if len(trnas_locs) > 0:
        merge_files(trnas_locs, path.join(output_dir, "trnas.tsv"), True)
    rrnas_locs = [
        i.rrnas_loc
        for i in annotations_list
        if i is not None
        if i.rrnas_loc is not None
    ]
    if len(rrnas_locs) > 0:
        merge_files(rrnas_locs, path.join(output_dir, "rrnas.tsv"), True)

    # make output gbk dir
    gbk_dir = path.join(output_dir, "genbank")
    mkdir(gbk_dir)
    for anno in annotations_list:
        if anno is not None:
            # TODO: make annotate_fasta generate a genbank dir and then copy it's contents, get rid of Annotation.name
            if path.isfile(anno.gbk_loc):
                copy2(anno.gbk_loc, path.join(gbk_dir, "%s.gbk" % anno.name))
            else:
                for gbk_loc in glob(path.join(anno.gbk_loc, "*.gbk")):
                    copy2(gbk_loc, gbk_dir)
    if write_annotations:
        all_annotations.to_csv(path.join(output_dir, "annotations.tsv"), sep="\t")
    else:
        return all_annotations

def merge_checkm_quality(annotations, checkm_quality):
    checkm_quality = pd.concat(
        [pd.read_csv(i, sep="\t", index_col=0) for i in checkm_quality]
    )
    checkm_quality.index = [
        strip_endings(i, [".fa", ".fasta", ".fna"]) for i in checkm_quality.index
    ]

    completeness = list()
    contamination = list()
    quality_missing_bins = list()
    for i in annotations.fasta:
        # add completeness and contamination
        if i in checkm_quality.index:
            completeness.append(checkm_quality.loc[i, "Completeness"])
            contamination.append(checkm_quality.loc[i, "Contamination"])
        else:
            completeness.append(0)
            contamination.append(100)
            quality_missing_bins.append(i)
        for j in set(quality_missing_bins):
            logger.warning(
                "Bin %s was not found in quality file, "
                "replaced with completeness 0 and contamination 100." % j
            )
    annotations["bin_completeness"] = completeness
    annotations["bin_contamination"] = contamination

