"""
==========
Call Genes
==========

This is a very simple tool to call the genes in a set of MAG or other set of FASTAs. It does some other things too of course. You may think of them as side effects to the main purpose, but they are still important. It renames the scaffolds, it filters the called genes based on the minimum contig size.

"""

import click
import logging
from pathlib import Path
from pkg_resources import resource_filename
from collections import Counter
from skbio.io import read as read_sequence, write as write_sequence
from multiprocessing import Pool
from functools import partial
from typing import Optional
from shutil import rmtree

from dram2.utils import Fasta, run_process, DramUseError
from dram2.utils.globals import FASTAS_CONF_TAG
from dram2.cli.context import DramContext, get_time_stamp_id, __version__

DEFAULT_MIN_CONTIG_SIZE: int = 2500
DEFAULT_PRODIGAL_MODE: str = "meta"
DEFAULT_TRANS_TABLE: str = "11"
DEFAULT_GENES_FILE: str = "genes"
GENES_RUN_TAG: str = "genes"

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
    help=(
        "Remove all called genes and information about them, you"
        " will only get the current set of genes from the command"
    ),
)
@click.option(
    "--prodigal_mode",
    default=DEFAULT_PRODIGAL_MODE,
    type=click.Choice(["train", "meta", "single"], case_sensitive=False),
    help=(
        "Mode of prodigal to use for gene calling. NOTE: normal"
        " or single mode require genomes which are high quality with"
        " low contamination and long contigs (average length >3 Kbp)."
    ),
)
@click.option(
    "--genes_dir",
    default=None,
    type=click.Path(path_type=Optional[Path]),
    help=(
        "Mode of prodigal to use for gene calling. NOTE: normal or"
        " single mode require genomes which are high quality with low"
        " contamination and long contigs (average length >3 Kbp)."
    ),
)
@click.option(
    "--prodigal_trans_tables",
    type=click.Choice([str(i) for i in range(1, 26)], case_sensitive=False),
    default=DEFAULT_TRANS_TABLE,
    help=(
        "Mode of prodigal to use for gene calling. NOTE:"
        " normal or single mode require genomes which are high quality"
        " with low contamination and long contigs"
        " (average length >3 Kbp)."
    ),
)
@click.pass_context
def call_genes_cmd(
    ctx: click.Context,
    fasta_paths: list[Path],
    genes_dir: Optional[Path],
    min_contig_size=DEFAULT_MIN_CONTIG_SIZE,
    prodigal_mode=DEFAULT_PRODIGAL_MODE,
    prodigal_trans_tables=DEFAULT_TRANS_TABLE,
    force: bool = False,
):
    """
    Call Genes and Filter Fastas
    ----------------------------

    Prodigal is one of many tools that we use in the DRAM pipeline. You will
    notice that this function not only calls Prodigal, it also performs a number
    of checks and organizes the files.  These first steps allow us to be
    confident we will not fail down the line.

    This command takes a positional argument/arguments namely FASTAs. The FASTA
    are path/paths to FASTAs representing mags or other genome collections of
    uncalled genes. This means that the use of the program will look like
    this::

       dram2 -o dram_dir call <option> /some/path/*.fasta

    or This ::

       dram2 -o dram_dir call <option> /some/path/fasta1.fasta /some/path/fasta2.fasta

    :param ctx:
    :param fasta_paths:
    :param genes_dir:
    :param min_contig_size:
    :param prodigal_mode:
    :param prodigal_trans_tables:
    :param force:
    """
    context: DramContext = ctx.obj
    logger = context.get_logger()
    output_dir: Path = context.get_output_dir()
    keep_tmp: bool = context.keep_tmp
    cores: int = context.cores
    run_id: str = get_time_stamp_id(GENES_RUN_TAG)

    if genes_dir is None:
        genes_dir = output_dir / DEFAULT_GENES_FILE

    project_config: dict = context.get_project_config()
    # get assembly locations

    try:

        new_config = call_genes(
            output_dir=output_dir,
            fasta_paths=fasta_paths,
            cores=cores,
            logger=logger,
            run_id=run_id,
            working_dir=genes_dir,
            project_config=project_config,
            keep_tmp=keep_tmp,
            min_contig_size=min_contig_size,
            prodigal_mode=prodigal_mode,
            prodigal_trans_tables=prodigal_trans_tables,
            force=force,
        )

        project_config.update(new_config)
        context.set_project_config(project_config)

    except Exception as e:

        logger.error(e)
        logger.exception("Fatal error in calling genes")

        raise e

def call_genes(
    output_dir: Path,
    fasta_paths: list[Path],
    cores: int,
    logger: logging.Logger,
    working_dir: Path,
    project_config: dict,
    run_id: str,
    keep_tmp: bool,
    min_contig_size=DEFAULT_MIN_CONTIG_SIZE,
    prodigal_mode=DEFAULT_PRODIGAL_MODE,
    prodigal_trans_tables=DEFAULT_TRANS_TABLE,
    force: bool = False,
) -> dict:
    """
    Call Genes in FASTAs With Prodigal
    --------------------

    Prodigal is one of many tools that we use in the DRAM pipeline. You will
    notice that this function not only calls prodigal it also performs a number
    of checks and organizes the files.  These first steps allow us to be
    confident we will not fail down the line.

    TODO:

      - split this into smaller functions


    :param ctx: The context passed from click
    :param fasta_paths: A list of FASTA files, probably each representing a BIN from a MAG.
    :param min_contig_size: The minimum contig size for DRAM to consider calling genes on
    :param prodigal_mode: Mode of prodigal to use, see prodigal documentation for more info
    :param prodigal_trans_tables: The number of trans tables to use for prodigal see the
        prodigal documentation for more information on this
    :raises ValueError:
    """

    # get assembly locations
    if force:
        logger.info(
            f"The force flag is being used, the old genes directories will be"
            f" fully deleted from: {working_dir}"
        )
        _ = [rmtree(x) for x in working_dir.glob("**/*") if not x.is_file()]
        clean_called_genes(output_dir, project_config, logger)
        old_names = []
    else:
        old_names = (
            []
            if "genes_called" not in project_config
            else [
                j[0]
                for i in project_config["genes_called"].values()
                for j in i["fastas"]
            ]
        )
    logger.info(f"Started calling genes for {len(fasta_paths)} fasta/s.")
    fastas_named: list[Fasta] = get_fasta_names_dirs(
        fasta_paths, working_dir, cores, old_names, logger
    )
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
            fastas_named,
        )
    logger.info("gene calling was a success, updating DRAM logs")
    fastas: list[Fasta] = [i for i in fastas_called if i is not None]
    new_config = {
        "genes_called": {
            run_id: {
                "min_contig_size": min_contig_size,
                "prodigal_mode": prodigal_mode,
                "prodigal_trans_tables": prodigal_trans_tables,
                # "annotated": False, # may use in the future
                "working_dir": working_dir.relative_to(output_dir).as_posix(),
                FASTAS_CONF_TAG: [i.name for i in fastas],
            }
        },
        FASTAS_CNF_TAG: [i.export(output_dir) for i in fastas],
    }
    if len(fastas) < 1:
        raise DramUsageError("No genes found, DRAM2 will not be able to proceed")
    return new_config




def clean_called_genes(output_dir: Path, project_config: dict, logger: logging.Logger):
    if (genes_called := project_config.get("genes_called")) is not None:
        for i in genes_called.values():
            if (
                genes_path := import_posible_path(i["working_dir"], output_dir)
            ) is not None:
                rmtree(genes_path)
    del genes_called


def filter_fasta(fasta_loc: Path, min_len, output_loc) -> Optional[list]:
    """
    Removes sequences shorter than a set minimum from FASTA files, outputs an object or to a file

    TODO:

     - Type hint the result better

    :param fasta_loc: A FASTA file path, probably a BIN from a MAG.
    :param min_len: The minimum contig size for DRAM to consider calling genes on
    :param output_loc: It's the output location
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
    :param fasta: A FASTA object, probably representing a BIN from a MAG.
    :param logger: Standard python logger
    :param min_contig_size: The minimum contig size for DRAM to consider calling genes on
    :param keep_tmp: True or False keep the temp file
    :param prodigal_mode: Mode of prodigal to use
    :param trans_table: Prodigal trans table setting, look it up on Prodigal website
    :returns: A Fasta object
    """
    # filter input fasta
    filtered_fasta: Path = fasta.tmp_dir / "filtered_fasta.fa"
    filter_fasta(fasta.origin, min_contig_size, filtered_fasta)

    if filtered_fasta.stat().st_size == 0:
        logger.warning(f"No sequences in {fasta.name} were longer than min_contig_size")
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
    -----------

    Runs the prodigal gene caller on a given FASTA file, outputs resulting files
    to the "tmp_dir" directory, this is usually but not always temporary.

    TODO:
      - Prodigal should be multithreaded in the stable release by the time you
        read this. So increasing to at least 2 threads each should help with
        efficient execution. That is however only an assumption and needs
        profiling to confirm.
      - Filtering may not be totally necessary, explore this option.


    :param filtered_fasta: The already filtered FASTA file, filtering may not be necessary
    :param tmp_dir: Output file, usually temporary but not every time
    :param logger: Standard python logger
    :param mode: Prodigal mode, look it up
    :param trans_table: Prodigal trans table setting, look it up on prodigals website
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


def get_fasta_name(fasta_loc: Path):
    """
    :param fasta_loc: A path to a fasta file or fasta.gz file
    :returns: The name of the fasta file as a string
    """
    if fasta_loc.suffix == "gz":
        return Path(fasta_loc.stem).stem
    else:
        return fasta_loc.stem


def mkdir(i: str, working_dir: Path):
    (working_dir / i).mkdir(exist_ok=False, parents=True)


def get_fasta_names_dirs(
    fasta_paths: list[Path],
    working_dir: Path,
    cores: int,
    old_fasta_names: Optional[list[str]],
    logger: logging.Logger,
) -> list[Fasta]:
    """
    :param fasta_paths: A list of FASTA files, probably each representing a BIN from a MAG.
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
        duplicated = [
            item for item, count in Counter(all_fasta_names).items() if count > 1
        ]
        if len(duplicated) > 0:
            logger.debug(f"duplicated names: {','.join(duplicated)}")
            raise DramUsageError(
                f"Genome file names must be unique. There is/are {len(duplicated)} name/s that appear twice in this search."
            )
        # make tmp_dirs
        p.map(partial(mkdir, working_dir=working_dir), fasta_names)

    return [
        Fasta(i, j, (working_dir / i), None, None, None, None)
        for i, j in zip(fasta_names, fasta_paths)
    ]
