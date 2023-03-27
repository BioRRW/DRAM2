import click

from dram2.cli.context import DramContext, DEFAULT_KEEP_TMP, __version__

# TODO: make it take an input and output gff location and not overwrite
# TODO: for some reason 1 is getting added to intervals when added to gff


def add_intervals_to_gff(
    annotations_loc, gff_loc, len_dict, interval_function, groupby_column, logger
):
    # get fasta length dict so we can merge, I'd love to be able to get this from somewhere else
    annotation_frame = pd.read_csv(annotations_loc, sep="\t")
    # process trnas to intervals
    annotation_dict = dict()
    for scaffold, frame in annotation_frame.groupby(groupby_column):
        if type(scaffold) is not str:
            scaffold = str(scaffold)
        else:
            scaffold = scaffold.strip()
        im = IntervalMetadata(len_dict[scaffold])
        for i, (_, row) in enumerate(frame.iterrows()):
            i += 1
            begin, end, metadata = interval_function(scaffold, row, i)
            begin -= 1
            if begin < 0:
                begin = 0
                logger.warning(
                    "Interval added to gff started less than zero, set to zero"
                )
            im.add(bounds=[(begin, end)], metadata=metadata)
        annotation_dict[scaffold] = im
    # add trna intervals to gff
    gff = list(read_sequence(gff_loc, format="gff3"))
    with open(gff_loc, "w") as f:
        for scaffold, gff_intervals in gff:
            gff_intervals = IntervalMetadata(len_dict[scaffold], gff_intervals)
            if scaffold in annotation_dict:
                gff_intervals.merge(annotation_dict[scaffold])
                gff_intervals.sort()
            f.write(
                gff_intervals.write(
                    io.StringIO(), format="gff3", seq_id=scaffold
                ).getvalue()
            )


# TODO: add annotations that don't end with '_id'
def annotate_gff(input_gff, output_gff, annotations, prefix=None):
    """Go through a gff and add a prefix to the scaffold and gene number for all ID's"""
    f = open(input_gff)
    o = open(output_gff, "w")
    for line in f:
        line = line.strip()
        if not line.startswith("#"):
            # replace id with new name
            old_scaffold = line.split("\t")[0]
            match = re.search(r"ID=\d*_\d*;", line)
            gene_number = match.group().split("_")[-1][:-1]
            old_gene_name = "%s_%s" % (old_scaffold, gene_number)
            if prefix is not None:  # add prefix to line name and gene name if given
                line = "%s_%s" % (prefix, line)
                gene_name = "%s_%s" % (prefix, old_gene_name)
            else:
                gene_name = old_gene_name
            line = re.sub(r"ID=\d*_\d*;", "ID=%s;" % gene_name, line)
            # get annotations to add from annotations file and add to end of line
            annotations_to_add = {
                strip_endings(i, ["_id"]): annotations.loc[old_gene_name, i]
                for i in annotations.columns
                if i.endswith("_id")
            }
            database_information = [
                'Dbxref="%s:%s"' % (key.strip(), value.strip().replace(":", "_"))
                for key, values in annotations_to_add.items()
                if not pd.isna(values)
                for value in values.split("; ")
                if value != ""
            ]
            if len(database_information) > 0:
                line += "%s;" % ";".join(database_information)
        o.write("%s\n" % line)


def get_gene_data(fasta_loc):
    """Take the prodigal gene headers and get the scaffold that it came from
    Based on idba_ud 'scaffold_#' scaffold names with gene name after
    """
    df_dict = dict()
    for seq in read_sequence(fasta_loc, format="fasta"):
        split_label = seq.metadata["id"].split("_")
        scaffold = "_".join(split_label[:-1])
        gene_position = split_label[-1]
        start_position, end_position, strandedness = seq.metadata["description"].split(
            "#"
        )[1:4]
        df_dict[seq.metadata["id"]] = [
            scaffold,
            int(gene_position),
            int(start_position),
            int(end_position),
            int(strandedness),
        ]
    return pd.DataFrame.from_dict(
        df_dict,
        orient="index",
        columns=[
            "scaffold",
            "gene_position",
            "start_position",
            "end_position",
            "strandedness",
        ],
    )


def get_unannotated(fasta_loc, annotations):
    """Get the genes from the fasta which did not get any annotations"""
    return [
        seq.metadata["id"]
        for seq in read_sequence(fasta_loc, format="fasta")
        if seq.metadata["id"] not in annotations
    ]


def assign_grades(annotations):
    """Grade genes based on reverse best hits to KEGG, UniRef and Pfam"""
    grades = dict()
    for gene, row in annotations.iterrows():
        if row.get("kegg_RBH") is True:
            rank = "A"
        elif row.get("uniref_RBH") is True:
            rank = "B"
        elif not pd.isna(row.get("kegg_hit")) or not pd.isna(row.get("uniref_hit")):
            rank = "C"
        elif (
            not pd.isna(row.get("pfam_hits"))
            or not pd.isna(row.get("cazy_hits"))
            or not pd.isna(row.get("peptidase_hit"))
        ):
            rank = "D"
        else:
            rank = "E"
        grades[gene] = rank
    return pd.DataFrame(grades, index=["rank"]).T


def generate_annotated_fasta(
    input_fasta, annotations, name=None, all_annotations: bool = False
):
    """Generates fasta entries with added annotation information to the header of a fasta
    either add best annotation (based on grade) (verbosity = short) or all annotations (verbosity = long)
    """
    for seq in read_sequence(input_fasta, format="fasta"):
        if not seq.metadata["id"] in annotations.index:
            seq.metadata["description"] = ""
            yield seq
            continue
        annotation = annotations.loc[seq.metadata["id"]]
        if "rank" in annotations.columns and not all_annotations:
            annotation_str = "rank: %s" % annotation["rank"]
            if annotation["rank"] == "A":
                annotation_str += "; %s (db=%s)" % (annotation.kegg_hit, "kegg")
            elif annotation["rank"] == "B":
                annotation_str += "; %s (db=%s)" % (annotation.uniref_hit, "uniref")
            elif annotation["rank"] == "C":
                if "kegg_hit" in annotation:
                    annotation_str += "; %s (db=%s)" % (annotation.kegg_hit, "kegg")
                if "uniref_hit" in annotation:
                    annotation_str += "; %s (db=%s)" % (annotation.uniref_hit, "uniref")
            elif annotation["rank"] == "D":
                annotation_str += "; %s (db=%s)" % (annotation.pfam_hits, "pfam")
            else:
                pass
        else:
            annotation_list = []
            if "rank" in annotations.columns:
                annotation_list += ["rank: %s" % annotation["rank"]]
            if "kegg_hit" in annotations.columns:
                if not pd.isna(annotation.kegg_hit):
                    annotation_list += ["%s (db=%s)" % (annotation.kegg_hit, "kegg")]
            if "uniref_hit" in annotations.columns:
                if not pd.isna(annotation.uniref_hit):
                    annotation_list += [
                        "%s (db=%s)" % (annotation.uniref_hit, "uniref")
                    ]
            if "pfam_hits" in annotations.columns:
                if not pd.isna(annotation.pfam_hits):
                    annotation_list += ["%s (db=%s)" % (annotation.pfam_hits, "pfam")]
            if "bin_taxonomy" in annotations.columns:
                if not pd.isna(annotation.bin_taxonomy):
                    annotation_list += [annotation.bin_taxonomy]
            annotation_str = "; ".join(annotation_list)
        if name is not None:
            seq.metadata["id"] = "%s_%s" % (name, seq.metadata["id"])
        seq.metadata["description"] = annotation_str
        yield seq


def create_annotated_fasta(input_fasta, annotations, output_fasta, name=None):
    """For use with genes files, added annotations"""
    write_sequence(
        generate_annotated_fasta(input_fasta, annotations, name),
        format="fasta",
        into=output_fasta,
    )


def generate_renamed_fasta(input_fasta, prefix):
    """For use with scaffolds files, merges together bins with fasta name added as a prefix to the file"""
    for seq in read_sequence(input_fasta, format="fasta"):
        seq.metadata["id"] = "%s_%s" % (prefix, seq.metadata["id"])
        yield seq


def rename_fasta(input_fasta, output_fasta, prefix):
    """See above"""
    write_sequence(
        generate_renamed_fasta(input_fasta, prefix), format="fasta", into=output_fasta
    )


def make_gbk_from_gff_and_fasta(
    gff_loc="genes.gff", fasta_loc="scaffolds.fna", faa_loc="genes.faa", output_gbk=None
):
    # filter scaffolds out of fasta which are not in genes.gff
    gff_frame = pd.read_csv(gff_loc, sep="\t", comment="#", header=None)
    gff_scaffolds = set(gff_frame[0])
    capture_fasta = io.StringIO()
    write_sequence(
        (
            i
            for i in read_sequence(fasta_loc, format="fasta")
            if i.metadata["id"] in gff_scaffolds
        ),
        into=capture_fasta,
        format="fasta",
    )
    # scikit-bio can make a genbank for a concatenated gff and fasta file so make this first
    gff = open(gff_loc)
    aa_records = {i.metadata["id"]: i for i in read_sequence(faa_loc, format="fasta")}
    # the gff with fasta format is the gff then ##FASTA then the fasta
    concat_gff = "%s\n##FASTA\n%s" % (gff.read(), capture_fasta.getvalue())

    # now we can only ready one record from the gff/fasta hybrid at a time
    # read a record at a time till end of the fasta
    genbank_records = ""
    # +1 because there is no \n> for the first line in the file and +1 because we are indexing starting at 1
    for i in range(1, capture_fasta.getvalue().count("\n>") + 1 + 1):
        seq = read_sequence(
            io.StringIO(concat_gff), format="gff3", into=Sequence, seq_num=i
        )
        seq.metadata["LOCUS"] = {
            "locus_name": seq.metadata["id"],
            "size": len(seq),
            "unit": "bp",
            "mol_type": "DNA",
            "shape": "linear",
            "division": "ENV",
            "date": time.strftime("%d-%b-%Y").upper(),
        }
        seq_bounds = (
            seq.interval_metadata.lower_bound,
            seq.interval_metadata.upper_bound,
        )
        for interval in seq.interval_metadata.query([seq_bounds]):
            interval.metadata["gene"] = interval.metadata["ID"]
            if interval.metadata["ID"] in aa_records:
                interval.metadata["translation"] = aa_records[interval.metadata["ID"]]
        # need to capture genbank output so we can combine into a multi-genbank file
        capture_print = io.StringIO()
        seq.write(capture_print, format="genbank")
        genbank_records += capture_print.getvalue()

    # if no output then return as string, otherwise write write to output location
    if output_gbk is None:
        return genbank_records
    else:
        open(output_gbk, "w").write(genbank_records)


def generate_genbank():
    pass


@click.command(
    "generate_genbank",
    context_settings=dict(help_option_names=["-h", "--help"]),
)
@click.pass_context
def generate_genbank_cmd(
    ctx: click.Context,
):
    """
    Make A DRAM GenBank File (Not Ready)
    ___

    Using an annotations file and a set of called genes make a GenBank file with select information from dram.

    I have not looked into implementing this yet and I don't know what will be required.

    """
    print("This command requires more work to work in dram2")
