
import click
import logging
from pathlib import Path
from pkg_resources import resource_filename
from skbio.io import read as read_sequence, write as write_sequence
from dram2.db_kits.utils import Fasta, run_process

DEFAULT_MIN_CONTIG_SIZE: int = 2500
DEFAULT_PRODIGAL_MODE: str = "meta"
DEFAULT_TRANS_TABLE: str = "11"

def get_config_loc():
    return Path(resource_filename("dram2.utils", "CONFIG"))
def filter_and_call_genes(
    fastas: list[Fasta],
    logger: logging.Logger,
    min_contig_size,
    keep_tmp,
    prodigal_mode,
    trans_table: str,
):
    for fasta in fastas:
        # filter input fasta
        filtered_fasta = fasta.tmp_dir / "filtered_fasta.fa"
        filter_fasta(fasta.origin, min_contig_size, filtered_fasta)

        if filtered_fasta.stat().st_size == 0:
            logger.warning(
                f"No sequences in {fasta.name} were longer than min_contig_size"
            )
            continue
        # predict ORFs with prodigal
        # TODO: handle when prodigal returns no genes
        faa, fna, gff = run_prodigal(
            filtered_fasta,
            fasta.tmp_dir,
            logger,
            mode=prodigal_mode,
            trans_table=trans_table,
        )
        if faa.stat().st_size == 0:
            logger.warning(f"No genes were returned by Prodigal for {fasta.name}")
            continue
        if not keep_tmp:
            filtered_fasta.unlink()
        yield Fasta(fasta.name, fasta.origin, fasta.tmp_dir, faa, fna, gff, None)

def filter_fasta(fasta_loc: Path, min_len=5000, output_loc=None):
    """Removes sequences shorter than a set minimum from fasta files, outputs an object or to a file"""
    kept_seqs = (
        seq
        for seq in read_sequence(str(fasta_loc.absolute()), format="fasta")
        if len(seq) >= min_len
    )
    if output_loc is None:
        return list(kept_seqs)
    else:
        write_sequence(kept_seqs, format="fasta", into=output_loc.absolute().as_posix())


def run_prodigal(
    filtered_fasta: Path,
    tmp_dir: Path,
    logger: logging.Logger,
    mode=DEFAULT_PRODIGAL_MODE,
    trans_table: str = DEFAULT_TRANS_TABLE,
) -> tuple[Path, Path, Path]:
    """Runs the prodigal gene caller on a given fasta file, outputs resulting files to given directory"""
    faa = tmp_dir / "genes.faa"
    fna = tmp_dir / "genes.fna"
    gff = tmp_dir / "genes.gff"

    run_process(
        [
            "prodigal",
            "-i",
            filtered_fasta.resolve(),
            "-p",
            mode,
            "-g",
            trans_table,
            "-f",
            "gff",
            "-o",
            gff.resolve(),
            "-a",
            faa.resolve(),
            "-d",
            fna.resolve(),
        ],
        logger,
    )
    return faa, fna, gff

def get_fasta_name(fasta_loc:Path):
    if fasta_loc.suffix == 'gz':
        return Path(fasta_loc.stem).stem
    else:
        return fasta_loc.stem

def get_fasta_names(fasta_locs,working_dir: Path):
    fasta_names = [get_fasta_name(i) for i in fasta_locs]
    # make temporary directory
    if len(fasta_names) != len(set(fasta_names)):
        raise ValueError(
            "Genome file names must be unique. At least one name appears twice in this search."
        )
    # make tmp_dirs
    tmp_dirs = [working_dir / i for i in fasta_names]
    for i in tmp_dirs:
        i.mkdir(exist_ok=True)

    return [
        Fasta(i, j, k, None, None, None, None)
        for i, j, k in zip(fasta_names, fasta_locs, tmp_dirs)
    ]


@click.command("call")
@click.argument(
    "fasta_paths",
    type=click.Path(exists=True, path_type=Path),
    nargs=-1,
)
@click.option(
    "--prodigal_mode",
    default=DEFAULT_PRODIGAL_MODE,
    type=click.Choice(["train", "meta", "single"], case_sensitive=False),
    help="Mode of prodigal to use for gene calling. NOTE: normal or single mode require genomes which are high quality with low contamination and long contigs (average length >3 Kbp).",
)
@click.option(
    "--prodigal_trans_tables",
    type=click.Choice([str(i) for i in range(1, 26)], case_sensitive=False),
    default=DEFAULT_TRANS_TABLE,
    help="Mode of prodigal to use for gene calling. NOTE: normal or single mode require genomes which are high quality with low contamination and long contigs (average length >3 Kbp).",
)
@click.pass_context
def call_genes(
    ctx,
    fasta_paths: list[Path],
    min_contig_size=DEFAULT_MIN_CONTIG_SIZE,
    prodigal_mode=DEFAULT_PRODIGAL_MODE,
    prodigal_trans_tables=DEFAULT_TRANS_TABLE,

):
    logger = ctx.obj.get_logger()
    working_dir = ctx.obj.get_working_dir()
    # get assembly locations
    fastas = get_fasta_names(fasta_paths, working_dir)
    keep_tmp = True
    fastas = [
        i
        for i in filter_and_call_genes(
            fastas,
            logger,
            min_contig_size,
            keep_tmp,
            prodigal_mode,
            prodigal_trans_tables,
        )
    ]


    if len(fastas) < 1:
        raise ValueError("No genes found, DRAM2 can't proceed")
