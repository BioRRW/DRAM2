from os import path, stat
from typing import Optional
import tarfile
from shutil import move, rmtree
from dram2.db_kits.fegenie_kit import process
from dram2.utils.utils import download_file, run_process, get_package_path, Fasta
from dram2.db_kits.utils import (
    make_mmseqs_db,
    run_hmmscan,
    get_best_hits,
    BOUTFMT6_COLUMNS,
    DBKit,
    get_sig_row,
    FILE_LOCATION_TAG,
    DRAM_DATAFOLDER_TAG,
    DBKIT_TAG,
)


from pathlib import Path
from functools import partial
import logging
import pandas as pd

VERSION = "1.0.0-beta.1"
CITATION: str = "CAMPER has no citeation and is in beta so you should not be using it."
FA_DB_KEY = "camper_fa_db"
FA_DB_CUTOFFS_KEY = "camper_fa_db_cutoffs"
HMM_KEY = "camper_hmm"
HMM_CUTOFFS_KEY = "camper_hmm_cutoffs"
DISTILLATE_KEY = "camper_distillate"


def rank_per_row(row):
    r_a = row["A_rank"]
    r_b = row["B_rank"]
    score = row["bitScore"]
    if score is None:
        return None
    if float(score) >= float(r_a):
        return "A"
    if pd.isnull(r_b):
        return None
    if float(score) >= float(r_b):
        return "B"
    return None


def blast_search_formater(hits_path, db_name, info_db, logger):
    if stat(hits_path).st_size == 0:
        return pd.DataFrame()
    hits = pd.read_csv(
        hits_path, sep="\t", header=None, names=BOUTFMT6_COLUMNS, index_col="qId"
    )
    hits = hits.merge(info_db, how="left", left_on="tId", right_index=True)
    rank_col = f"{db_name}_rank"
    hits[rank_col] = hits.apply(rank_per_row, axis=1)
    hits.dropna(subset=[rank_col], inplace=True)
    logger.debug("Getting descriptions of hits from %s" % db_name)
    hits = hits[["tId", rank_col, "bitScore", "ID_for_distillate", "definition"]]
    hits.rename(
        columns={
            "tId": f"{db_name}_hits",
            "ID_for_distillate": f"{db_name}_id",
            "bitScore": f"{db_name}_bitScore",
            "definition": f"{db_name}_definition",
        },
        inplace=True,
    )
    hits[f"{db_name}_search_type"] = "blast"
    return hits


def bitScore_per_row(row):
    if row["score_type"] == "domain":
        return row.domain_score
    elif row["score_type"] == "full":
        return row.full_score
    elif row["score_type"] == "-":
        return None
    else:
        raise ValueError("The score_type must be 'domain', 'full', or 'j")


def hmmscan_formater(
    hits: pd.DataFrame, db_name: str, hmm_info_path: str, top_hit: bool = True
):
    if hmm_info_path is None:
        hmm_info = None
        hits = hits[hits.apply(get_sig_row, axis=1)]
    else:
        hmm_info = pd.read_csv(hmm_info_path, sep="\t", index_col=0)
        hits = hits.merge(hmm_info, how="left", left_on="target_id", right_index=True)
        hits["bitScore"] = hits.apply(bitScore_per_row, axis=1)
        hits["score_rank"] = hits.apply(rank_per_row, axis=1)
        hits.dropna(subset=["score_rank"], inplace=True)
    if len(hits) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame()
    if top_hit:
        # Get the best hits
        # TODO check we want top hit
        hits = hits.sort_values("full_evalue").drop_duplicates(subset=["query_id"])
    hits.set_index("query_id", inplace=True, drop=True)
    hits.rename_axis(None, inplace=True)
    if "definition" in hits.columns:
        hits = hits[["target_id", "score_rank", "bitScore", "definition"]]
        hits.columns = [
            f"{db_name}_id",
            f"{db_name}_rank",
            f"{db_name}_bitScore",
            f"{db_name}_hits",
        ]
    else:
        hits = hits[["target_id", "score_rank", "bitScore"]]
        hits.columns = [f"{db_name}_id", f"{db_name}_rank", f"{db_name}_bitScore"]
    # Rename
    hits[f"{db_name}_search_type"] = "hmm"
    return hits


def get_minimum_bitscore(info_db):
    bit_score_threshold = min(info_db[["A_rank", "B_rank"]].min().values)
    return bit_score_threshold

def blast_search(
    query_db,
    target_db,
    working_dir,
    info_db_path,
    db_name,
    logger,
    threads,
):
    """A convenience function to do a blast style forward best hits search"""
    # Get kegg hits
    info_db = pd.read_csv(info_db_path, sep="\t", index_col=0)
    bit_score_threshold = get_minimum_bitscore(info_db)
    hits_path = get_best_hits(
        query_db,
        target_db,
        logger,
        working_dir,
        "gene",
        db_name,
        bit_score_threshold,
        threads,
    )
    return blast_search_formater(hits_path, db_name, info_db, logger)


# in the future the database will get the same input as was given in the data
class CamperKit(DBKit):

    name = "camper"
    formal_name: str = "CAMPER"
    version: str = VERSION
    citation: str = CITATION
    camper_fa_db: Path
    camper_fa_db_cutoffs: Path
    camper_hmm: Path
    camper_hmm_cutoffs: Path
    camper_distillate: Path
    search_type: str = "hmm_and_blast_style"
    has_genome_summary: bool= True

    def get_genome_summary(self) -> Path:
        genome_summary_form = self.request_config_path("genome_summary_form")
        if genome_summary_form is None:
            return get_package_path(Path("db_kits", "camper_kit", "CAMPER_distillate.tsv"))
        return genome_summary_form

    def search(self, fasta: Fasta) -> pd.DataFrame | pd.Series:
        if fasta.name is None:
            raise ValueError('A fasta file needs a name')
        self.logger.debug(f"Annotating {fasta.name} with {self.formal_name}.")
        if (
            self.camper_fa_db is None
            or self.camper_fa_db_cutoffs is None
            or self.camper_hmm is None
            or self.camper_hmm_cutoffs is None
        ):
            raise ValueError(
                "You must first load a valid config before you can search with "
                f"{self.formal_name}."
            )
        fasta_temp_dir = self.working_dir / self.name / fasta.name
        fasta_temp_dir.mkdir(parents=True, exist_ok=False)
        blast = blast_search(
            query_db=fasta.mmsdb,
            target_db=self.camper_fa_db,
            working_dir=fasta_temp_dir,
            info_db_path=Path(self.camper_fa_db_cutoffs),
            db_name=self.name,
            logger=self.logger,
            threads=self.threads,
        )
        hmm = run_hmmscan(
            genes_faa=fasta.faa.as_posix(),
            db_loc=self.camper_hmm.as_posix(),
            db_name=self.name,
            threads=self.threads,
            output_loc=fasta_temp_dir.as_posix(),
            logger=self.logger,
            formater=partial(
                hmmscan_formater,
                db_name=self.name,
                hmm_info_path=self.camper_hmm_cutoffs.as_posix(),
                top_hit=True,
            ),
        )
        if not self.keep_tmp:
            rmtree(fasta_temp_dir)
        full = pd.concat([blast, hmm])
        if len(full) < 1:
            return pd.DataFrame()

        self.check_counter_after_annotation()
        return full.groupby(full.index).apply(
            lambda x: (
                x.sort_values(
                    f"{self.name}_search_type", ascending=True
                )  # make sure hmm is first
                .sort_values(f"{self.name}_bitScore", ascending=False)
                .iloc[0]
            )
        )

    def load_dram_config(self):
        self.camper_fa_db = self.get_config_path("mmsdb")
        self.camper_hmm = self.get_config_path("hmmdb")
        self.camper_fa_db_cutoffs = self.get_config_path("mmsdb_cutoffs")
        self.camper_hmm_cutoffs = self.get_config_path("hmmdb_cutoffs")
        self.camper_distillate = self.get_config_path("genome_summary_form")

    # def setup(self):
    #     if self.db_path is None:
    #         raise ValueError(f"CAMPER needs an output location to setup the db")
    #     tarfile = self.download(self.working_dir, self.logger)
    #     # tarfile = "./CAMPER-1.0.0-beta.1.tar.gz"
    #     self.config = self.pre_process(tarfile, self.db_path, self.logger, self.threads)
    #     # raise ValueError("I have not made this yet")

    @classmethod
    def download(cls, temporary, logger, version=None):
        """
        Retrieve CAMPER release tar.gz

        This will get a tar file that is automatically generated from making a campers release on git hub.  In order to
        avoid changes in CAMPER being blindly excepted into DRAM, a new number must be put into the OPTIONS global
        variable in order to change this.

        :param temporary: Usually in the output dir
        :param verbose: TODO replace with logging setting
        :returns: Path to tar
        """
        raise NotImplementedError(
            "Your request depends on a feature that has not yet been implimented.\n",
            "More specify, CAMPER is not public at time of writing so cant be Downloaded",
        )
        if version is None:
            version = cls.version
        camper_database = path.join(temporary, f"CAMPER_{version}.tar.gz")
        # Note the 'v' in the name, GitHub wants it in the tag then it just takes it out. This could be a problem
        download_file(
            f"https://github.com/WrightonLabCSU/CAMPER/archive/refs/tags/v{version}.tar.gz",
            logger,
            camper_database,
        )
        return camper_database

    @classmethod
    def get_descriptions(self, annotation):
        return pd.DataFrame()

    @classmethod
    def setup(
        cls,
        camper_tar_gz,
        output_dir,
        logger,
        threads,
        version=None,
    ) -> dict:
        # Check if all the locations are in config
        # if not call setup
        if version is None:
            version = cls.version
        temp_dir = path.dirname(camper_tar_gz)
        tar_paths = {
            "camper_fa_db": path.join(f"CAMPER-{version}", "CAMPER_blast.faa"),
            "camper_hmm": path.join(f"CAMPER-{version}", "CAMPER.hmm"),
            "camper_fa_db_cutoffs": path.join(
                f"CAMPER-{version}", "CAMPER_blast_scores.tsv"
            ),
            "camper_distillate": path.join(
                f"CAMPER-{version}", "CAMPER_distillate.tsv"
            ),
            "camper_hmm_cutoffs": path.join(
                f"CAMPER-{version}", "CAMPER_hmm_scores.tsv"
            ),
        }

        output_dir = Path(output_dir).as_posix()
        final_paths = {
            "camper_fa_db": path.join(output_dir, "CAMPER_blast.mmsdb"),
            "camper_hmm": path.join(output_dir, "CAMPER.hmm"),
            "camper_fa_db_cutoffs": path.join(output_dir, "CAMPER_blast_scores.tsv"),
            "camper_distillate": path.join(output_dir, "CAMPER_distillate.tsv"),
            "camper_hmm_cutoffs": path.join(output_dir, "CAMPER_hmm_scores.tsv"),
        }

        with tarfile.open(camper_tar_gz) as tar:
            for v in tar_paths.values():
                tar.extract(v, temp_dir)

        # move tsv files, and hmm to location
        for i in [
            "camper_fa_db_cutoffs",
            "camper_distillate",
            "camper_hmm_cutoffs",
            "camper_hmm",
        ]:
            move(path.join(temp_dir, tar_paths[i]), final_paths[i])

        # build dbs
        make_mmseqs_db(
            path.join(temp_dir, tar_paths["camper_fa_db"]),
            final_paths["camper_fa_db"],
            logger,
            threads=threads,
        )
        run_process(
            ["hmmpress", "-f", final_paths["camper_hmm"]], logger
        )  # all are pressed just in case
        return final_paths
