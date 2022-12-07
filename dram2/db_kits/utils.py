"""General utils for database objects including the template class"""
import re
import subprocess
from os import path, stat
import pandas as pd
import logging
from typing import Callable
from dram2.utils.utils import run_process
HMMSCAN_ALL_COLUMNS = [
    "query_id",
    "query_ascession",
    "query_length",
    "target_id",
    "target_ascession",
    "target_length",
    "full_evalue",
    "full_score",
    "full_bias",
    "domain_number",
    "domain_count",
    "domain_cevalue",
    "domain_ievalue",
    "domain_score",
    "domain_bias",
    "target_start",
    "target_end",
    "alignment_start",
    "alignment_end",
    "query_start",
    "query_end",
    "accuracy",
    "description",
]
HMMSCAN_COLUMN_TYPES = [
    str,
    str,
    int,
    str,
    str,
    int,
    float,
    float,
    float,
    int,
    int,
    float,
    float,
    float,
    float,
    int,
    int,
    int,
    int,
    int,
    int,
    float,
    str,
]

BOUTFMT6_COLUMNS = [
    "qId",
    "tId",
    "seqIdentity",
    "alnLen",
    "mismatchCnt",
    "gapOpenCnt",
    "qStart",
    "qEnd",
    "tStart",
    "tEnd",
    "eVal",
    "bitScore",
]


def process_reciprocal_best_hits(
    forward_output_loc, reverse_output_loc, target_prefix="target"
):
    """Process the forward and reverse best hits results to find reverse best hits
    Returns the query gene, target gene, if it was a reverse best hit, % identity, bit score and e-value
    """
    forward_hits = pd.read_csv(
        forward_output_loc, sep="\t", header=None, names=BOUTFMT6_COLUMNS
    )
    forward_hits = forward_hits.set_index("qId")
    reverse_hits = pd.read_csv(
        reverse_output_loc, sep="\t", header=None, names=BOUTFMT6_COLUMNS
    )
    reverse_hits = reverse_hits.set_index("qId")

    def check_hit(row: pd.Series):
        rbh = False
        if row.tId in reverse_hits.index:
            rbh = row.name == reverse_hits.loc[row.tId].tId
        return {
            "%s_hit" % target_prefix: row.tId,
            "%s_RBH" % target_prefix: rbh,
            "%s_identity" % target_prefix: row.seqIdentity,
            "%s_bitScore" % target_prefix: row.bitScore,
            "%s_eVal" % target_prefix: row.eVal,
            "index": row.name,
        }

    hits = forward_hits.apply(check_hit, axis=1, result_type="expand")
    # NOTE these lines may not be necessary
    hits.set_index("index", drop=True, inplace=True)
    hits.index.name = None
    return hits


def get_reciprocal_best_hits(
    query_db,
    target_db,
    logger,
    output_dir=".",
    query_prefix="query",
    target_prefix="target",
    bit_score_threshold=60,
    rbh_bit_score_threshold=350,
    threads=10,
    verbose=False,
):
    """Take results from best hits and use for a reciprocal best hits search"""
    # TODO: Make it take query_target_db as a parameter
    # create subset for second search
    query_target_db_top_filt = path.join(
        output_dir,
        "%s_%s.tophit.minbitscore%s.mmsdb"
        % (query_prefix, target_prefix, bit_score_threshold),
    )  # I DON'T LIKE THIS
    query_target_db_filt_top_swapped = path.join(
        output_dir,
        "%s_%s.minbitscore%s.tophit.swapped.mmsdb"
        % (query_prefix, target_prefix, bit_score_threshold),
    )
    # swap queries and targets in results database
    run_process(
        [
            "mmseqs",
            "swapdb",
            query_target_db_top_filt,
            query_target_db_filt_top_swapped,
            "--threads",
            str(threads),
        ],
        logger,
        verbose=verbose,
    )
    target_db_filt = path.join(output_dir, "%s.filt.mmsdb" % target_prefix)
    # create a subdatabase of the target database with the best hits as well as the index of the target database
    run_process(
        [
            "mmseqs",
            "createsubdb",
            query_target_db_filt_top_swapped,
            target_db,
            target_db_filt,
        ],
        logger,
        verbose=verbose,
    )
    run_process(
        [
            "mmseqs",
            "createsubdb",
            query_target_db_filt_top_swapped,
            "%s_h" % target_db,
            "%s_h" % target_db_filt,
        ],
        logger,
        verbose=verbose,
    )

    return get_best_hits(
        target_db_filt,
        query_db,
        logger,
        output_dir,
        target_prefix,
        query_prefix,
        rbh_bit_score_threshold,
        threads,
        verbose,
    )


def do_blast_style_search(
    query_db,
    target_db,
    working_dir,
    db_handler,
    formater,
    logger,
    db_name="database",
    bit_score_threshold=60,
    rbh_bit_score_threshold=350,
    threads=10,
    verbose=False,
):
    """A convenience function to do a blast style reciprocal best hits search"""
    # Get kegg hits
    logger.info("Getting forward best hits from %s" % db_name)
    forward_hits = get_best_hits(
        query_db,
        target_db,
        logger,
        working_dir,
        "gene",
        db_name,
        bit_score_threshold,
        threads,
        verbose=verbose,
    )
    if stat(forward_hits).st_size == 0:
        return pd.DataFrame(columns=[f"{db_name}_hit"])
    logger.info("Getting reverse best hits from %s" % db_name)
    reverse_hits = get_reciprocal_best_hits(
        query_db,
        target_db,
        logger,
        working_dir,
        "gene",
        db_name,
        bit_score_threshold,
        rbh_bit_score_threshold,
        threads,
        verbose=verbose,
    )
    hits = process_reciprocal_best_hits(forward_hits, reverse_hits, db_name)
    logger.info("Getting descriptions of hits from %s" % (db_name))
    if "%s_description" % db_name in db_handler.get_database_names():
        header_dict = db_handler.get_descriptions(
            hits["%s_hit" % db_name], "%s_description" % db_name
        )
    else:
        header_dict = multigrep(
            hits["%s_hit" % db_name], "%s_h" % target_db, logger, "\x00", working_dir
        )
    hits = formater(hits, header_dict)
    return hits


def make_mmseqs_db(
    fasta_loc, output_loc, logger, create_index=True, threads=10, verbose=False
):
    """Takes a fasta file and makes a mmseqs2 database for use in blast searching and hmm searching with mmseqs2"""
    run_process(["mmseqs", "createdb", fasta_loc, output_loc], logger, verbose=verbose)
    if create_index:
        tmp_dir = path.join(path.dirname(output_loc), "tmp")
        run_process(
            ["mmseqs", "createindex", output_loc, tmp_dir, "--threads", str(threads)],
            logger,
            verbose=verbose,
        )


def run_hmmscan(
    genes_faa: str,
    db_loc: str,
    db_name: str,
    output_loc: str,
    formater: Callable,
    logger: logging.Logger,
    threads: int = 2,
    verbose: bool = False,
):
    output = path.join(output_loc, f"{db_name}_results.unprocessed.b6")
    run_process(
        ["hmmsearch", "--domtblout", output, "--cpu", str(threads), db_loc, genes_faa],
        logger,
        verbose=verbose,
    )
    # Parse hmmsearch output
    if not (path.isfile(output) and stat(output).st_size > 0):
        return pd.DataFrame()
    hits = parse_hmmsearch_domtblout(output)
    if len(hits) < 1:
        return pd.DataFrame()
    return formater(hits)


def parse_hmmsearch_domtblout(file):
    df_lines = list()
    for line in open(file):
        if not line.startswith("#"):
            line = line.split()
            line = line[:22] + [" ".join(line[22:])]
            df_lines.append(line)
    hmmsearch_frame = pd.DataFrame(df_lines, columns=HMMSCAN_ALL_COLUMNS)
    for i, column in enumerate(hmmsearch_frame.columns):
        hmmsearch_frame[column] = hmmsearch_frame[column].astype(
            HMMSCAN_COLUMN_TYPES[i]
        )
    return hmmsearch_frame


def get_best_hits(
    query_db,
    target_db,
    logger,
    output_dir=".",
    query_prefix="query",
    target_prefix="target",
    bit_score_threshold=60,
    threads=10,
    verbose=False,
):
    """Uses mmseqs2 to do a blast style search of a query db against a target db, filters to only include best hits
    Returns a file location of a blast out format 6 file with search results
    """
    # make query to target db
    tmp_dir = path.join(output_dir, "tmp")
    query_target_db = path.join(
        output_dir, "%s_%s.mmsdb" % (query_prefix, target_prefix)
    )
    run_process(
        [
            "mmseqs",
            "search",
            query_db,
            target_db,
            query_target_db,
            tmp_dir,
            "--threads",
            str(threads),
        ],
        logger,
        verbose=verbose,
    )
    # filter query to target db to only best hit
    query_target_db_top = path.join(
        output_dir, "%s_%s.tophit.mmsdb" % (query_prefix, target_prefix)
    )
    run_process(
        [
            "mmseqs",
            "filterdb",
            query_target_db,
            query_target_db_top,
            "--extract-lines",
            "1",
        ],
        logger,
        verbose=verbose,
    )
    # filter query to target db to only hits with min threshold
    query_target_db_top_filt = path.join(
        output_dir,
        "%s_%s.tophit.minbitscore%s.mmsdb"
        % (query_prefix, target_prefix, bit_score_threshold),
    )
    run_process(
        [
            "mmseqs",
            "filterdb",
            "--filter-column",
            "2",
            "--comparison-operator",
            "ge",
            "--comparison-value",
            str(bit_score_threshold),
            "--threads",
            str(threads),
            query_target_db_top,
            query_target_db_top_filt,
        ],
        logger,
        verbose=verbose,
    )
    # convert results to blast outformat 6
    forward_output_loc = path.join(
        output_dir, "%s_%s_hits.b6" % (query_prefix, target_prefix)
    )
    run_process(
        [
            "mmseqs",
            "convertalis",
            query_db,
            target_db,
            query_target_db_top_filt,
            forward_output_loc,
            "--threads",
            str(threads),
        ],
        logger,
        verbose=verbose,
    )
    return forward_output_loc

from abc import ABC

class DBKit(ABC):

    def __init__(self, name: str, formal_name: str, citation: str, db_version: str):
        self.is_dbkit: bool = True
        self.name: str = name
        self.formal_name: str = formal_name
        self.db_version: str = db_version
        self.citation: str = citation
        self.settings: dict[str:str] = {}

    @abstractmethod
    def get_descriptions(self):
        pass

    def download(self, args):
        return false

    @abstractmethod
    def process(self, args):
        pass

    @abstractmethod
    def search(self):
        pass

class FastaKit(DBKit):

    def process(
        custom_fasta_loc: str,
        custom_db_name: str,
        output_dir: str,
        logger: logging.Logger,
        threads=1,
        verbose=False,
    ):
        # if none is passed from argparse then set to tuple of len 0
        mkdir(output_dir)
    
        if custom_fasta_loc is None:
            custom_fasta_loc = ()
        if custom_db_name is None:
            custom_db_name = ()
        if len(custom_fasta_loc) != len(custom_db_name):
            raise ValueError(
                "Lengths of custom db fasta list and custom db name list must be the same."
            )
        custom_dbs = {
            custom_db_name[i]: custom_fasta_loc[i] for i in range(len(custom_db_name))
        }
        custom_db_locs = dict()
        for db_name, db_loc in custom_dbs.items():
            custom_db_loc = path.join(output_dir, "%s.custom.mmsdb" % db_name)
            make_mmseqs_db(
                db_loc, custom_db_loc, logger=logger, threads=threads, verbose=verbose
            )
            custom_db_locs[db_name] = custom_db_loc
        return custom_db_locs

    def process(
        custom_fasta_loc: str,
        custom_db_name: str,
        output_dir: str,
        logger: logging.Logger,
        threads=1,
        verbose=False,
    ):
        # if none is passed from argparse then set to tuple of len 0
        mkdir(output_dir)
    
        if custom_fasta_loc is None:
            custom_fasta_loc = ()
        if custom_db_name is None:
            custom_db_name = ()
        if len(custom_fasta_loc) != len(custom_db_name):
            raise ValueError(
                "Lengths of custom db fasta list and custom db name list must be the same."
            )
        custom_dbs = {
            custom_db_name[i]: custom_fasta_loc[i] for i in range(len(custom_db_name))
        }
        custom_db_locs = dict()
        for db_name, db_loc in custom_dbs.items():
            custom_db_loc = path.join(output_dir, "%s.custom.mmsdb" % db_name)
            make_mmseqs_db(
                db_loc, custom_db_loc, logger=logger, threads=threads, verbose=verbose
            )
            custom_db_locs[db_name] = custom_db_loc
        return custom_db_locs


class HmmKit(DBKit):
    
    def process(self, custom_hmm_loc, custom_hmm_name, logger, verbose=False):
        if custom_hmm_loc is None:
            custom_hmm_loc = ()
        if custom_hmm_name is None:
            custom_hmm_name = ()
        if len(custom_hmm_loc) != len(custom_hmm_name):
            raise ValueError(
                "Lengths of custom db hmm list and custom hmm db name list must be the same."
            )
        custom_hmm_locs = dict()
        for i in range(len(custom_hmm_name)):
            run_process(
                ["hmmpress", "-f", custom_hmm_loc[i]], logger, verbose=verbose
            )  # all are pressed just in case
            custom_hmm_locs[custom_hmm_name[i]] = custom_hmm_loc[i]
        return custom_hmm_locs

