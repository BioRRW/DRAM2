"""
If you actualy read this that is cool its my last day and thinging is geting hard.

any whay use this like
conda activate ./dram2_env
dram2 -d cant_hyd_db build_db --make_db cant_hyd --update
rm -r cant_hyd_search
dram2 -d cant_hyd_search annotate --use_db cant_hyd dram2/db_kits/cant_hyd_kit/data/BacMet_ExpVerified_BiocideRes_genes_SHORT.faa tests/data/NC_001422.faa
exit
continue

databases[0].search(fasta)
search_results: list[pd.DataFrame] = [j.search(fasta) for j in databases]
"""
from os import path, stat
from shutil import move, rmtree, copy
from pathlib import Path
import logging
import urllib.parse
from functools import partial

import pandas as pd

from dram2.utils import download_file, run_process, Fasta, get_package_path
from dram2.db_kits.utils import (
    make_mmseqs_db,
    run_hmmscan,
    get_best_hits,
    BOUTFMT6_COLUMNS,
    DBKit,
    get_sig_row,
    run_hmmscan,
)


VERSION: str = "1.0.0-beta.1"
CITATION: str = "CAMPER has no citeation and is in beta so you should not be using it."

HMM_KEY: str = "cant_hyd_hmmdb"
FAA_KEY: str = "cant_hyd_faa"
HMM_SCORES_KEY: str = "cant_hyd_hmmdb_cutoffs"
FAA_SCORES_KEY: str = "cant_hyd_faa_scores"
DRAM2_MODUAL_KEY: str = "cant_hyd_dram2_modual"
GENOME_SUMMARY_KEY = "genome_summary_form"
MMSDB_KEY = "mmsdb"
MMSDB_CUTOFFS_KEY = "mmsdb_cutoffs"

HMM_DEST: str = "CANT_HYD.hmm"
FAA_DEST: str = "CANT_HYD.faa"
MSDB_DEST: str = "CANT_HYD.msdb"
HMM_SCORES_IN: str = "CANT_HYD_HMM_scores.csv"
FAA_SCORES_IN: str = "CANT_HYD_BLAST_scores.csv"
HMM_SCORES_DEST: str = "CANT_HYD_HMM_scores.tsv"
FAA_SCORES_DEST: str = "CANT_HYD_BLAST_scores.tsv"
DRAM2_MODUAL_DEST: str = "engineeredsys_dram_module.tsv"


def get_minimum_bitscore(info_db):
    bit_score_threshold = min(info_db[["A_rank", "B_rank"]].min().values)
    return bit_score_threshold


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
        query_db=query_db,
        target_db=target_db,
        logger=logger,
        output_dir=working_dir,
        bit_score_threshold=bit_score_threshold,
        query_prefix="gene",
        target_prefix=db_name,
        threads=threads,
    )
    return blast_search_formater(hits_path, db_name, info_db, logger)


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


class CantHydKit(DBKit):
    """A tool implement a database in DRAM"""

    name = "cant_hyd"  # the name as used by dram all lowercase leters and _
    formal_name: str = "CANT-HYD"  # The actual name, any ascii caricter you want
    version: str = "1.0"  # the version as a string
    # formated citation in ascii
    citation: str = (
        "Khot V, Zorz J, Gittins DA, Chakraborty A, Bell E, Bautista MA, "
        "Paquette AJ, Hawley AK, Novotnik B, Hubert CRJ, Strous M and "
        "Bhatnagar S(2022) CANT-HYD: A Curated Database of "
        "Phylogeny-Derived Hidden Markov Models for Annotation of Marker "
        "Genes Involved in Hydrocarbon Degradation. Front. Microbiol. "
        "12: 764058. doi: 10.3389/fmicb.2021.764058"
    )
    search_type: str = (
        "hmm_and_blast_style"  # describe the type hmm_style or blast_style
    )
    has_genome_summary: bool = True  # if there is a geneome summary acociated say so.
    selectable: bool = True
    location_keys: list[str] = [
        HMM_KEY,
        FAA_KEY,
        HMM_SCORES_KEY,
        FAA_SCORES_KEY,
        DRAM2_MODUAL_KEY,
    ]

    def get_genome_summary(self) -> Path:
        """Get the ids from a complete annotations pandas DataFrame."""
        genome_summary = self.get_config_path(GENOME_SUMMARY_KEY)
        return genome_summary

    def search(self, fasta: Fasta) -> pd.DataFrame | pd.Series:
        """Perform a search, be that HMM, Blast or something else."""
        if fasta.faa is None:
            raise ValueError(
                f"Fasta with out called genes faa was passed to the search function for {self.formal_name}"
            )
        self.logger.info(f"Getting hits from {self.formal_name}")
        fasta_temp_dir = self.working_dir / self.name / fasta.name
        fasta_temp_dir.mkdir(parents=True, exist_ok=False)
        hmm_annotations = run_hmmscan(
            genes_faa=fasta.faa.absolute().as_posix(),
            db_loc=self.hmmdb.absolute().as_posix(),
            db_name=self.name,
            output_loc=fasta_temp_dir.absolute().as_posix(),
            threads=self.threads,
            logger=self.logger,
            formater=partial(
                hmmscan_formater,
                db_name=self.name,
                hmm_info_path=self.hmm_cutoffs,
                top_hit=True,
            ),
        )
        blast_annotatons = blast_search(
            query_db=fasta.mmsdb,
            target_db=self.mmsdb,
            working_dir=fasta_temp_dir,
            info_db_path=Path(self.mmsdb_cutoffs),
            db_name=self.name,
            logger=self.logger,
            threads=self.threads,
        )
        full = pd.concat([hmm_annotations, blast_annotatons])
        if len(full) < 1:
            return pd.DataFrame()

        return full.groupby(full.index).apply(
            lambda x: (
                x.sort_values(
                    f"{self.name}_search_type", ascending=True
                )  # make sure hmm is first
                .sort_values(f"{self.name}_bitScore", ascending=False)
                .iloc[0]
            )
        )

    # def get_descriptions(self, hits):
    #     header_dict = multigrep(
    #         hits[f"{self.name}_hit"],
    #         f"{self.mmsdb_target}_h",
    #         self.logger,
    #         "\x00",
    #         self.working_dir,
    #     )
    #     hits = blast_search_formater(hits, header_dict, self.name)
    #     return hits

    def load_dram_config(self):
        """Extract data from the larger DRAM config yaml"""
        self.mmsdb = self.get_config_path(MMSDB_KEY)
        self.hmmdb = self.get_config_path(HMM_KEY)
        self.mmsdb_cutoffs = self.get_config_path(MMSDB_CUTOFFS_KEY)
        self.hmm_cutoffs = self.get_config_path(HMM_SCORES_KEY)

    def download(self, user_locations_dict: dict[str, Path]) -> dict[str, Path]:
        """Download your raw data, return a full location dictionary."""
        if HMM_KEY not in user_locations_dict:
            hmm_path = self.working_dir / "CANT-HYD.hmm"
            download_file(
                (
                    f"https://raw.githubusercontent.com/dgittins/"
                    f"CANT-HYD-HydrocarbonBiodegradation/main/HMMs/"
                    f"{urllib.parse.quote('concatenated HMMs/')}"
                    f"CANT-HYD.hmm"
                ),
                hmm_path,
                self.logger,
            )
        else:
            hmm_path = user_locations_dict[HMM_KEY]

        # The rest we just pull locally
        if FAA_KEY in user_locations_dict:
            faa_path = user_locations_dict[FAA_KEY]
        else:
            faa_path = get_package_path(
                Path(
                    "db_kits",
                    "cant_hyd_kit",
                    "data",
                    "BacMet_ExpVerified_BiocideRes_genes_SHORT.faa",
                )
            )

        if FAA_SCORES_KEY in user_locations_dict:
            faa_scores_path = user_locations_dict[FAA_SCORES_KEY]
        else:
            faa_scores_path = get_package_path(
                Path(
                    "db_kits",
                    "cant_hyd_kit",
                    "data",
                    FAA_SCORES_IN,
                )
            )
        if HMM_SCORES_KEY in user_locations_dict:
            hmm_scores_path = user_locations_dict[HMM_SCORES_KEY]
        else:
            hmm_scores_path = get_package_path(
                Path("db_kits", "cant_hyd_kit", "data", HMM_SCORES_IN)
            )
        if DRAM2_MODUAL_KEY in user_locations_dict:
            dram2_modual_path = user_locations_dict[DRAM2_MODUAL_KEY]
        else:
            dram2_modual_path = get_package_path(
                Path("db_kits", "cant_hyd_kit", "data", "engineeredsys_dram_module.tsv")
            )

        location_dictionary: dict[str, Path] = {
            HMM_KEY: hmm_path,
            FAA_KEY: faa_path,
            HMM_SCORES_KEY: hmm_scores_path,
            FAA_SCORES_KEY: faa_scores_path,
            DRAM2_MODUAL_KEY: dram2_modual_path,
        }
        return location_dictionary

    def setup(
        self,
        location_dict: dict[str, Path],
        output_dir: Path,
    ) -> dict[str, dict[str, str]]:
        """Do whatever you need to process the data at the locations provided."""
        move(location_dict[HMM_KEY], (output_dir / HMM_DEST))
        (
            pd.read_csv(location_dict[HMM_SCORES_KEY]).to_csv(
                output_dir / HMM_SCORES_DEST, sep="\t", index=False
            )
        )
        (
            pd.read_csv(location_dict[FAA_SCORES_KEY]).to_csv(
                output_dir / FAA_SCORES_DEST, sep="\t", index=False
            )
        )
        copy(location_dict[DRAM2_MODUAL_KEY], (output_dir / DRAM2_MODUAL_DEST))
        make_mmseqs_db(
            location_dict[FAA_KEY],
            output_dir / MSDB_DEST,
            self.logger,
            threads=self.threads,
        )
        run_process(
            ["hmmpress", "-f", (output_dir / HMM_DEST).as_posix()], self.logger
        )  # all are pressed just in case
        return {
            GENOME_SUMMARY_KEY: {"location": DRAM2_MODUAL_DEST},
            MMSDB_KEY: {"location": MSDB_DEST},
            MMSDB_CUTOFFS_KEY: {"location": FAA_SCORES_DEST},
            HMM_KEY: {"location": HMM_DEST},
            HMM_SCORES_KEY: {"location": HMM_SCORES_DEST},
        }


def test_dbkit_CantHydKit_Setup():
    """Test the CANT HYD DBKIT"""
    test_dbkit = CantHydKit({}, logging.getLogger())
    test_dbkit.set_args(Path("temp2", "temp"))
    download_paths = test_dbkit.download({})
    received_output = test_dbkit.setup(download_paths, Path("temp2"))
    expected_output = {
        "genome_summary_form": {"location": "engineeredsys_dram_module.tsv"},
        "mmsdb": {"location": "CANT_HYD.msdb"},
        "mmsdb_cutoffs": {"location": "CANT_HYD_BLAST_scores.tsv"},
        "hmmdb": {"location": "CANT_HYD.hmm"},
        "hmmdb_cutoffs": {"location": "CANT_HYD_HMM_scores.tsv"},
    }
    assert received_output == expected_output


def test_dbkit_CantHydKit_search():
    """Test the CANT HYD DBKIT"""
    pass
