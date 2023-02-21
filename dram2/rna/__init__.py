import click

__version__ = '2.0.0'

RAW_RRNA_COLUMNS = [
    "scaffold",
    "tool_name",
    "type",
    "begin",
    "end",
    "e-value",
    "strand",
    "empty",
    "note",
]
RRNA_COLUMNS = ["fasta", "begin", "end", "strand", "type", "e-value", "note"]

def get_dups(columns):
    keep = list()
    seen = list()
    for column in columns:
        if column in seen:
            keep.append(False)
        else:
            seen.append(column)
            keep.append(True)
    return keep

def count_motifs(gene_faa, motif="(C..CH)"):
    motif_count_dict = dict()
    for seq in read_sequence(gene_faa, format="fasta"):
        motif_count_dict[seq.metadata["id"]] = len(list(seq.find_with_regex(motif)))
    return motif_count_dict


def make_trnas_interval(scaffold, row, i):
    if row["Begin"] < row["End"]:
        begin = row["Begin"]
        end = row["End"]
        strand = "+"
    else:
        begin = row["End"]
        end = row["Begin"]
        strand = "-"
    metadata = {
        "source": "tRNAscan-SE",
        "type": "tRNA",
        "score": row["Score"],
        "strand": strand,
        "phase": 0,
        "ID": "%s_tRNA_%s" % (scaffold, i),
        "codon": row["Codon"],
        "product": "tRNA-%s" % row["Type"],
    }
    if not pd.isna(row["Note"]):
        metadata["Note"] = row["Note"]
    return begin, end, metadata


def run_barrnap(fasta, fasta_name, logger, threads=10):
    raw_rrna_str = run_process(
        ["barrnap", "--threads", str(threads), fasta],
        logger,
        capture_stdout=True,
        check=False,
    )
    raw_rrna_table = pd.read_csv(
        io.StringIO(raw_rrna_str),
        skiprows=1,
        sep="\t",
        header=None,
        names=RAW_RRNA_COLUMNS,
        index_col=0,
    )
    rrna_table_rows = list()
    for gene, row in raw_rrna_table.iterrows():
        rrna_row_dict = {
            entry.split("=")[0]: entry.split("=")[1] for entry in row["note"].split(";")
        }
        rrna_table_rows.append(
            [
                fasta_name,
                row.begin,
                row.end,
                row.strand,
                rrna_row_dict["Name"].replace("_", " "),
                row["e-value"],
                rrna_row_dict.get("note", ""),
            ]
        )
    if len(raw_rrna_table) > 0:
        return pd.DataFrame(
            rrna_table_rows, index=raw_rrna_table.index, columns=RRNA_COLUMNS
        ).reset_index()
    else:
        logger.warning("No rRNAs were detected, no rrnas.tsv file will be created.")
        return None




def run_trna_scan(
    fasta,
    tmp_dir,
    fasta_name,
    logger,
    threads=10,
):
    """Run tRNAscan-SE on scaffolds and create a table of tRNAs as a separate output"""
    raw_trnas = path.join(tmp_dir, "raw_trnas.txt")
    run_process(
        ["tRNAscan-SE", "-G", "-o", raw_trnas, "--thread", str(threads), fasta],
        logger,
    )
    if path.isfile(raw_trnas) and stat(raw_trnas).st_size > 0:
        trna_frame = pd.read_csv(raw_trnas, sep="\t", skiprows=[0, 2])
        trna_frame.columns = [i.strip() for i in trna_frame.columns]
        # if begin.1 or end.1 are in trnas then drop, else drop the second begin or end
        if "Begin.1" in trna_frame.columns:
            trna_frame = trna_frame.drop(["Begin.1"], axis=1)
        if "End.1" in trna_frame.columns:
            trna_frame = trna_frame.drop(["End.1"], axis=1)
        trna_frame = trna_frame.loc[:, get_dups(trna_frame.columns)]
        trna_frame.insert(0, "fasta", fasta_name)
        return trna_frame
    else:
        logger.warning("No tRNAs were detected, no trnas.tsv file will be created.")
        return None

@click.command('pull_rrna')
@click.version_option(__version__)
def pull_rrna():
    print("This comand is comming soon")

@click.command('pull_trna')
@click.version_option(__version__)
def pull_trna():
    print("This comand is comming soon")


#trash
import click

__version__ = '2.0.0'

def get_dups(columns):
    keep = list()
    seen = list()
    for column in columns:
        if column in seen:
            keep.append(False)
        else:
            seen.append(column)
            keep.append(True)
    return keep


def run_trna_scan(
    fasta,
    tmp_dir,
    fasta_name,
    logger,
    threads=10,
):
    """Run tRNAscan-SE on scaffolds and create a table of tRNAs as a separate output"""
    raw_trnas = path.join(tmp_dir, "raw_trnas.txt")
    run_process(
        ["tRNAscan-SE", "-G", "-o", raw_trnas, "--thread", str(threads), fasta],
        logger,
    )
    if path.isfile(raw_trnas) and stat(raw_trnas).st_size > 0:
        trna_frame = pd.read_csv(raw_trnas, sep="\t", skiprows=[0, 2])
        trna_frame.columns = [i.strip() for i in trna_frame.columns]
        # if begin.1 or end.1 are in trnas then drop, else drop the second begin or end
        if "Begin.1" in trna_frame.columns:
            trna_frame = trna_frame.drop(["Begin.1"], axis=1)
        if "End.1" in trna_frame.columns:
            trna_frame = trna_frame.drop(["End.1"], axis=1)
        trna_frame = trna_frame.loc[:, get_dups(trna_frame.columns)]
        trna_frame.insert(0, "fasta", fasta_name)
        return trna_frame
    else:
        logger.warning("No tRNAs were detected, no trnas.tsv file will be created.")
        return None

@click.command('pull_rrna')
@click.version_option(__version__)
def pull_rrna():
    print("This comand is comming soon")

@click.command('pull_trna')
@click.version_option(__version__)
def pull_trna():
    print("This comand is comming soon")
