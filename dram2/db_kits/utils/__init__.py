"""General utils for database objects including the template class"""
import re
from typing import NamedTuple
from abc import ABC, abstractmethod
from os import path, stat
from typing import Callable, Union
from typing import Optional
from functools import partial
from pathlib import Path
from dataclasses import dataclass
import logging

import pandas as pd

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


@dataclass
class Fasta:
    name: Optional[str]
    origin: Optional[Path]
    tmp_dir: Optional[Path]
    faa: Optional[Path]
    fna: Optional[Path]
    gff: Optional[Path]
    mmsdb: Optional[Path]



def run_mmseqs_profile_search(
    query_db,
    pfam_profile,
    output_loc,
    logger,
    output_prefix="mmpro_results",
    db_handler=None,
    threads=10,
):
    """Use mmseqs to run a search against pfam, currently keeping all hits and not doing any extra filtering"""
    tmp_dir = path.join(output_loc, "tmp")
    output_db = path.join(output_loc, "%s.mmsdb" % output_prefix)
    run_process(
        [
            "mmseqs",
            "search",
            query_db,
            pfam_profile,
            output_db,
            tmp_dir,
            "-k",
            "5",
            "-s",
            "7",
            "--threads",
            str(threads),
        ],
        logger,
    )
    output_loc = path.join(output_loc, "%s_output.b6" % output_prefix)
    run_process(
        ["mmseqs", "convertalis", query_db, pfam_profile, output_db, output_loc],
        logger,
    )
    pfam_results = pd.read_csv(
        output_loc, sep="\t", header=None, names=BOUTFMT6_COLUMNS
    )
    if pfam_results.shape[0] > 0:
        pfam_dict = dict()
        if db_handler is not None:
            pfam_descriptions = db_handler.get_descriptions(
                set(pfam_results.tId), "%s_description" % output_prefix
            )
        else:
            pfam_descriptions = {}
        for gene, pfam_frame in pfam_results.groupby("qId"):
            if len(pfam_descriptions) < 1:
                pfam_dict[gene] = "; ".join(pfam_frame.tId)
            else:
                pfam_dict[gene] = "; ".join(
                    [
                        "%s [%s]" % (pfam_descriptions[ascession], ascession)
                        for ascession in pfam_frame.tId
                    ]
                )
        return pd.DataFrame(pfam_dict, index=[f"{output_prefix}_hits"]).T
    else:
        return pd.DataFrame(columns=[f"{output_prefix}_hits"])


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
    )


def multigrep(search_terms, search_against, logger, split_char="\n", output="."):
    # TODO: multiprocess this over the list of search terms
    """Search a list of exact substrings against a database, takes name of mmseqs db index with _h to search against"""
    hits_file = path.join(output, "hits.txt")
    with open(hits_file, "w") as f:
        f.write("%s\n" % "\n".join(search_terms))
    results = run_process(
        ["grep", "-a", "-F", "-f", hits_file, search_against],
        logger,
        capture_stdout=True,
    )
    processed_results = [
        i.strip() for i in results.strip().split(split_char) if len(i) > 0
    ]
    # remove(hits_file)
    return {i.split()[0]: i for i in processed_results if i != ""}


def do_blast_style_search(
    query_db,
    target_db,
    working_dir,
    logger,
    db_name="database",
    bit_score_threshold=60,
    rbh_bit_score_threshold=350,
    threads=10,
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
    )
    hits = process_reciprocal_best_hits(forward_hits, reverse_hits, db_name)
    logger.info("Getting descriptions of hits from %s" % (db_name))
    # if "%s_description" % db_name in db_handler.get_database_names():
    #     header_dict = db_handler.get_descriptions(
    #         hits["%s_hit" % db_name], "%s_description" % db_name
    #     )
    # else:
    return hits


def make_mmseqs_db(
    fasta_loc,
    output_loc,
    logger,
    create_index=True,
    threads=10,
):
    """Takes a fasta file and makes a mmseqs2 database for use in blast searching and hmm searching with mmseqs2, """
    run_process(
        ["mmseqs", "createdb", fasta_loc, output_loc],
        logger,
    )
    if create_index:
        tmp_dir = path.join(path.dirname(output_loc), "tmp")
        run_process(
            ["mmseqs", "createindex", output_loc, tmp_dir, "--threads", str(threads)],
            logger,
        )


def run_hmmscan(
    genes_faa: str,
    db_loc: str,
    db_name: str,
    output_loc: str,
    formater: Callable,
    logger: logging.Logger,
    threads: int = 2,
):
    output = path.join(output_loc, f"{db_name}_results.unprocessed.b6")
    run_process(
        ["hmmsearch", "--domtblout", output, "--cpu", str(threads), db_loc, genes_faa],
        logger,
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
    query_db: Union[str, Path],
    target_db: Union[str, Path],
    logger,
    output_dir: Union[str, Path]=".",
    query_prefix="query",
    target_prefix="target",
    bit_score_threshold=60,
    threads=10,
):
    """Uses mmseqs2 to do a blast style search of a query db against a target db, filters to only include best hits
    Returns a file location of a blast out format 6 file with search results
    """
    # make query to target db
    tmp_dir = path.join(output_dir, "tmp")
    query_target_db = path.join(
        output_dir, "%s_%s.mmsdb" % (query_prefix, target_prefix)
    )
    run_process( [ "mmseqs", "search", str(query_db), str(target_db), str(query_target_db), str(tmp_dir), "--threads", str(threads), ], logger,)
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
    )
    return forward_output_loc


def process_custom_hmm_db_cutoffs(
    custom_hmm_db_cutoffs_loc, custom_hmm_db_name, logger
):
    if custom_hmm_db_cutoffs_loc is None:
        return {}
    if custom_hmm_db_name is None:
        raise ValueError(
            "You can't use the custom_hmm_db_cutoffs_loc argument without the custom_hmm_db_name and"
            " custom_hmm_db_locs aguments specified."
        )
    if len(custom_hmm_db_cutoffs_loc) != len(custom_hmm_db_name):
        logger.warning(
            f"Custom hmm cutoffs and descriptions were only provided to the first {len(custom_hmm_db_cutoffs_loc)}."
            " The rest of the custom hmms will use standard cutoffs and have no descriptions."
        )
    return {custom_hmm_db_name[i]: j for i, j in enumerate(custom_hmm_db_cutoffs_loc)}




def get_basic_description(hits, header_dict, db_name):
    """Get viral gene full descriptions based on headers (text before first space)"""
    hit_list = list()
    description = list()
    for hit in hits["%s_hit" % db_name]:
        header = header_dict[hit]
        hit_list.append(hit)
        description.append(header)
    new_df = pd.DataFrame(
        [hit_list, description],
        index=["%s_id" % db_name, "%s_hit" % db_name],
        columns=hits.index,
    )
    return pd.concat(
        [new_df.transpose(), hits.drop("%s_hit" % db_name, axis=1)], axis=1, sort=False
    )


def get_sig_row(row, evalue_lim: float = 1e-15):
    """Check if hmm match is significant, based on dbCAN described parameters"""
    tstart, tend, tlen, evalue = row[
        ["target_start", "target_end", "target_length", "full_evalue"]
    ].values
    perc_cov = (tend - tstart) / tlen
    if perc_cov >= 0.35 and evalue <= evalue_lim:
        return True
    else:
        return False


# TODO decide if we need use_hmmer_thresholds:bool=False
def generic_hmmscan_formater(
    hits: pd.DataFrame, db_name: str, hmm_info_path: Path = None, top_hit: bool = True
):
    if hmm_info_path is None:
        hmm_info = None
        hits_sig = hits[hits.apply(get_sig_row, axis=1)]
    else:
        hmm_info = pd.read_csv(hmm_info_path, sep="\t", index_col=0)
        hits_sig = sig_scores(hits, hmm_info)
    if len(hits_sig) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame()
    if top_hit:
        # Get the best hits
        hits_sig = hits_sig.sort_values("full_evalue").drop_duplicates(
            subset=["query_id"]
        )
    hits_df = hits_sig[["target_id", "query_id"]]
    hits_df.set_index("query_id", inplace=True, drop=True)
    hits_df.rename_axis(None, inplace=True)
    hits_df.columns = [f"{db_name}_id"]
    if hmm_info is not None:
        hits_df = hits_df.merge(
            hmm_info[["definition"]],
            how="left",
            left_on=f"{db_name}_id",
            right_index=True,
        )
        hits_df.rename(columns={"definition": f"{db_name}_hits"}, inplace=True)
    return hits_df


def sig_scores(hits: pd.DataFrame, score_db: pd.DataFrame) -> pd.DataFrame:
    is_sig = list()
    for i, frame in hits.groupby("target_id"):
        row = score_db.loc[i]
        if row["score_type"] == "domain":
            score = frame.domain_score
        elif row["score_type"] == "full":
            score = frame.full_score
        elif row["score_type"] == "-":
            continue
        else:
            raise ValueError(row["score_type"])
        frame = frame.loc[score.astype(float) > float(row.threshold)]
        is_sig.append(frame)
    if len(is_sig) > 0:
        return pd.concat(is_sig)
    else:
        return pd.DataFrame()

class DBKit(ABC):
    """
    DBKit Abstract Class
    ____________________

    Use this as a base class to model all other absract classess off of.

    """

    # this name will apear in lists as the
    name: str = ""
    formal_name: str = ""
    logger: logging.Logger
    working_dir: Path
    bit_score_threshold: int
    rbh_bit_score_threshold: int
    past_annotations_path: str
    kofam_use_dbcan2_thresholds: bool
    threads: int
    make_new_faa: bool
    dry: bool
    force: bool
    extra: dict
    config: dict = {}
    selectable: bool = True
    dram_db_loc: Path

    def set_universals(
        self, name: str, formal_name: str, config: dict, citation: str, db_version: str
    ):
        self.name: str = name
        self.config: dict= config
        self.is_dbkit: bool = True
        self.formal_name: str = formal_name
        self.db_version: str = db_version
        self.citation: str = citation

    def __init__(self, config: dict, args: dict):
        self.config = config
        self.set_args(**args)
        self.check_setup()

    @classmethod
    def download(cls):
        pass

    @classmethod
    def pre_process(cls):
        pass

    def set_args(
        self,
        logger: logging.Logger,
        working_dir: Path,
        output_dir: Path,
        bit_score_threshold: int,
        rbh_bit_score_threshold: int,
        past_annotations_path: str,
        kofam_use_dbcan2_thresholds: bool,
        threads: int,
        make_new_faa: bool,
        dry: bool,
        force: bool,
        extra: dict,
        db_path: Path
    ):
        self.kofam_use_dbcan2_thresholds: bool = kofam_use_dbcan2_thresholds
        self.logger: logging.Logger = logger
        self.working_dir: Path = working_dir
        self.output_dir: Path = output_dir
        self.bit_score_threshold: int = bit_score_threshold
        self.rbh_bit_score_threshold: int = rbh_bit_score_threshold
        self.past_annotations_path: str = past_annotations_path
        self.threads: int = threads
        self.make_new_faa: bool = make_new_faa
        self.dry: bool = dry
        self.force: bool = force
        self.extra: dict = extra
        self.db_path = self.setup_db_path(db_path)

    @staticmethod
    def setup_db_path(db_path:Path):
        if db_path is None:
            return db_path
        db_path = Path(db_path)
        if db_path.exists():
            return db_path
        else:
            db_path.mkdir(parents=True)


    @abstractmethod
    def check_setup(self):
        """
        This will be used to check if the database is setup. Unless you overwrite the constructor this functon will be called during annotation after the values have been stored. So you can use this to check user arguments even read in cusom arguments.
        """
        pass

    def check_on_fly_setup(self):
        '''
        TODO:
        - This should be a match but I don't feel like updating today
        -

        '''
        if self.config.get('default_db_dir') is None:
            self.config['default_db_dir'] = self.dram_db_loc
        if self.db_path is None:
            self.dram_db_loc = self.config['default_db_dir']
        if self.db_path is None and self.config.get('default_db_dir') is None:
            raise ValueError("Without a dram_db_directory defided, database can't be built on the fly")


    @abstractmethod
    def search(self):
        pass

    @abstractmethod
    def get_descriptions(self):
        pass

    @classmethod
    @abstractmethod
    def get_ids(cls, annotatons):
        pass


class FastaKit(DBKit):

    name = "custom_fasta_db"
    selectable: bool = False

    def __init__(self, name: str, loc: Path, config: dict, args: dict):
        self.set_universals(
            name, name, {}, "Custom FASTA", "No citation for custom DBs"
        )
        self.set_args(**args)
        self.config: dict = config
        self.fasta_loc: Path = loc
        self.pre_process()
        # if none is passed from argparse then set to tuple of len 0

    def pre_process(self):

        temp_dir = self.working_dir / f"{self.name}_fasta_db"
        if not temp_dir.exists():
            temp_dir.mkdir()
        self.mmsdb_target = temp_dir / f"{self.name}.mmsdb"
        make_mmseqs_db(
            self.fasta_loc,
            self.mmsdb_target.as_posix(),
            logger=self.logger,
            threads=self.threads,
        )

    def check_setup(self):
        pass

    def search(self, query_ob: Fasta) -> pd.DataFrame:
        annotatons = do_blast_style_search(
            query_ob.mmsdb,
            self.mmsdb_target,
            self.working_dir,
            self.logger,
            self.name,
            self.bit_score_threshold,
            self.rbh_bit_score_threshold,
            self.threads,
        )
        return annotatons

    def get_descriptions(self, annotatons):
        header_dict = multigrep(
            hits[f"{db_name}_hit"], f"{target_db}_h", self.logger, "\x00", self.working_dir
        )
        hits = get_basic_description(annotatons, header_dict, self.name)
        return hits

    @classmethod
    def get_ids(cls, annotatons, name):
        pass



class HmmKit(DBKit):

    name = "custom_hmm_db"
    selectable: bool = False

    def __init__(
        self, name: str, loc: Path, descriptions: Path, config: dict, args: dict
    ):
        self.set_universals(name, name, {}, "Custom hmm", "No citation for custom DBs")
        self.hmm_loc: Path = loc
        self.descriptions: Path = descriptions
        self.set_args(**args)
        self.process()
        self.config: dict = config
        self.fasta_loc: Path = loc
        self.process()

    def check_setup(self):
        pass

    def process(self):
        self.logger.info(f"Pre processing custom hmm database {self.name}")
        run_process(
            ["hmmpress", "-f", self.hmm_loc], self.logger
        )  # all are pressed just in case

    def search(self, query_ob: Fasta):
        self.logger.info(f"Annotating custom hmm database {self.name}")
        annotatons = run_hmmscan(
            genes_faa=query_ob.faa.as_posix(),
            db_loc=self.hmm_loc.as_posix(),
            db_name=self.name,
            threads=self.threads,
            output_loc=self.working_dir.as_posix(),
            formater=partial(
                generic_hmmscan_formater,
                db_name=self.name,
                hmm_info_path=self.descriptions,
                top_hit=True,
            ),
            logger=self.logger,
        )
        return annotatons

    def get_descriptions(self):
        pass

    @classmethod
    def get_ids(cls, annotatons, name):
        pass
