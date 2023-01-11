"""
Call Genes
==========

This is a very simple tool to call the genes in a set of MAG or other set of fastas. It dose some othere things too of coarse. You may think of them as side efects to the main perpose but they are still important. It renames the scafolds, it filters the called gehes based on the minimum contig size.

"""
import click
import logging
from pathlib import Path
from pkg_resources import resource_filename
import collections
from skbio.io import read as read_sequence, write as write_sequence
from dram2.db_kits.utils import Fasta, run_process
from multiprocessing import Pool
from functools import partial
from typing import Optional
from datetime import datetime
from shutil import rmtree


DEFAULT_MIN_CONTIG_SIZE: int = 2500
DEFAULT_PRODIGAL_MODE: str = "meta"
DEFAULT_TRANS_TABLE: str = "11"

"""
import os
os.system("dram2 -o test call ./tests/data/NC_001422.fasta # pass")
os.system("dram2 -o test call ./tests/data/NC_001422.fasta # fail")
os.system("dram2 -o test call -f ./tests/data/NC_001422.fasta # pass")
"""

@click.command("call")
@click.argument(
    "fasta_paths",
    type=click.Path(exists=True, path_type=Path),
    nargs=-1,
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help= "Remove all called genes and information about them, you will only get the current set of genes from the command"
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
    force:bool = False,
):
    """
    Call Genes and organize files
    ________________________

    As previously mentoned Prodigal is one of many tools that we use in the DRAM pipline. You will notice that this function not only calls prodigal it also performs a number of checks and organizes the files.  These first steps alow us to be confindent we will not fail down the line.

    comandline expample:
    '''
    import os
    os.system("dram2 -o test call ./tests/data/NC_001422.fasta # pass")
    os.system("dram2 -o test call ./tests/data/NC_001422.fasta # fail")
    os.system("dram2 -o test call -f ./tests/data/NC_001422.fasta # pass")
    '''

    :param ctx: The context passed from click
    :param fasta_paths: A list of fasta files probably each representing a BIN from a MAG.
    :param min_contig_size: The minimum contig size for DRAM to consider calling genes on
    :param prodigal_mode: Mode of prodigal to use
    :param prodigal_trans_tables: The number of trans_tables to use for prodigal see the prodigal docs
    :raises ValueError:
    """
    logger:logging.Logger = ctx.obj.get_logger()
    output_dir:Path = ctx.obj.get_output_dir()
    keep_tmp:bool = ctx.obj.keep_tmp
    cores:int = ctx.obj.cores
    timestamp_id:str = datetime.now().strftime('%Y%m%d%H%M%S')
    working_dir = output_dir / timestamp_id
    project_config:dict = ctx.obj.get_project_config()
    # get assembly locations
    try:
        if force:
            clean_called_genes(project_config)
        if project_config.get("genes_called") is not None:
            old_names = [j[0] for i in project_config["genes_called"].values() for j in i['fastas']]
        else:
            old_names = []
        fastas_named: list[Fasta] = get_fasta_names_dirs(fasta_paths, working_dir, cores, old_names, logger)
        with Pool(cores) as p:
            fastas_called: list[Optional[Fasta]] = p.map(
                partial(
                   filter_and_call_genes,
                   logger=logger,
                   min_contig_size=min_contig_size,
                   keep_tmp=keep_tmp,
                   prodigal_mode=prodigal_mode,
                   trans_table=prodigal_trans_tables,
                ),
                fastas_named)
        fastas: list[Fasta]  = [ i for i in fastas_called if i is not None]
        new_config = {"genes_called": {
            timestamp_id: {
                "min_contig_size": min_contig_size,
                "prodigal_mode": prodigal_mode,
                "prodigal_trans_tables": prodigal_trans_tables,
                "annotated": False,
                "working_dir": working_dir.absolute().as_posix(),
                "fastas": [i.export() for i in fastas],
        }}}
        project_config.update(new_config)
        ctx.obj.set_project_config(project_config)

    except Exception as e:
        logger.error(e)
        logger.exception("Fatal error in calling genes")
        raise e

    if len(fastas) < 1:
        raise ValueError("No genes found, DRAM2 will not be able to proceed")

def clean_called_genes(project_config:dict):
    if project_config.get("genes_called") is None:
        return
    for j in project_config["genes_called"].values():
        rmtree(Path(j['working_dir']))
    del(project_config['genes_called'])


def filter_fasta(fasta_loc: Path, min_len, output_loc) -> Optional[list]:
    """
    Removes sequences shorter than a set minimum from fasta files, outputs an object or to a file

    TODO:

     - Type hint the result better

    :param fasta_loc: A fasta file path probably a BIN from a MAG.
    :param min_len: The minimum contige size for DRAM to consider calling genes on
    :param output_loc: Its the output location
    :returns:
    """
    kept_seqs = (
        seq
        for seq in read_sequence(str(fasta_loc.absolute()), format="fasta")
        if len(seq) >= min_len
    )
    if output_loc is None:
        return list(kept_seqs)
    else:
        write_sequence(kept_seqs, format="fasta", into=output_loc.absolute().as_posix())

def filter_and_call_genes(
    fasta: Fasta,
    logger: logging.Logger,
    min_contig_size,
    keep_tmp,
    prodigal_mode,
    trans_table: str,
) -> Optional[Fasta]:
    """
    :param fasta: A fasta object, probly representing a BIN from a MAG.
    :param logger: Standered python logger
    :param min_contig_size: The minimum contige size for DRAM to consider calling genes on
    :param keep_tmp: True or False keep the temp file
    :param prodigal_mode: Mode of prodigal to use
    :param trans_table: Prodigal trans table seting, look it up on Prodigal website
    :returns: A Fasta object
    """
    # filter input fasta
    filtered_fasta:Path = fasta.tmp_dir / "filtered_fasta.fa"
    filter_fasta(fasta.origin, min_contig_size, filtered_fasta)

    if filtered_fasta.stat().st_size == 0:
        logger.warning(
            f"No sequences in {fasta.name} were longer than min_contig_size"
        )
        return None
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
        return
    if not keep_tmp:
        filtered_fasta.unlink()
    return Fasta(fasta.name, fasta.origin, fasta.tmp_dir, faa, fna, gff, None)

def run_prodigal(
    filtered_fasta: Path,
    tmp_dir: Path,
    logger: logging.Logger,
    mode=DEFAULT_PRODIGAL_MODE,
    trans_table: str = DEFAULT_TRANS_TABLE,
) -> tuple[Path, Path, Path]:
    """
    Run Prodigal
    _____________


    Runs the prodigal gene caller on a given fasta file, outputs resulting files to the "tmp_dir" directory, this is usualy but not alwase temporary.

    TODO:
      - Prodigal should be multi threaded in the stable releas by the time you read this. So increse to atleast 2 threads each should help with eficient excution. That is however only an assumption and needs profiling to conferm
      - Filtering may not be totaly necessary, explore this option.


    :param filtered_fasta: The alread filterd fasta file, filtering may not be necciary
    :param tmp_dir: Output file, usualy tempory but not every time
    :param logger: Standered python logger
    :param mode: Prodigal mode, look it up
    :param trans_table: Prodigal trans table seting, look it up on prodigs webcite
    :returns: A tuple the faa path the fna path and the gff path.
    """
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
    """
    :param fasta_loc: A path to a fasta file or fasta.gz file
    :returns: The name of the fasta file as a string
    """
    if fasta_loc.suffix == 'gz':
        return Path(fasta_loc.stem).stem
    else:
        return fasta_loc.stem

def mkdir(i: str, working_dir:Path):
    (working_dir / i).mkdir(exist_ok=False, parents=True)

def get_fasta_names_dirs(fasta_paths: list[Path], working_dir: Path, cores: int, old_fasta_names: Optional[list[str]], logger: logging.Logger) -> list[Fasta]:
    """
    :param fasta_paths: A list of fasta files probly each representing a BIN from a MAG.
    :param working_dir:
    :raises ValueError: If the file names are not unique
    :returns:
    """
    with Pool(cores) as p:
        fasta_names: list[str] = p.map(get_fasta_name, fasta_paths)
        # make temporary directory
        all_fasta_names = fasta_names
        if old_fasta_names is not None:
            all_fasta_names += old_fasta_names
        duplicated = [item for item, count in collections.Counter(all_fasta_names).items() if count > 1]
        if len(duplicated) > 0:
            logger.debug(f"duplicated names: {','.join(duplicated)}")
            raise ValueError(
                f"Genome file names must be unique. There are {len(duplicated)} name/s that appear twice in this search."
            )
        # make tmp_dirs
        p.map(partial(mkdir, working_dir=working_dir), fasta_names)

    return [
        Fasta(i, j, (working_dir / i), None, None, None, None)
        for i, j in zip(fasta_names, fasta_paths)
    ]


