from os import path, stat
from shutil import move, rmtree
from pathlib import Path
import logging

import pandas as pd

from dram2.utils import download_file, run_process, Fasta, get_package_path
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


VERSION = "1.0.0-beta.1"
CITATION: str = "CAMPER has no citeation and is in beta so you should not be using it."

HMM_KEY = "cant_hyd_hmms"
FAA_KEY = "cant_hyd_faa"
HMM_SCORES_KEY = "cant_hyd_hmm_scores"
FAA_SCORES_KEY = "cant_hyd_faa_scores"
DRAM2_MODUAL_KEY = "cant_hyd_dram2_modual"


class CantHydHmmKit(DBKit):
    """A tool implement a database in DRAM"""

    name = "cant_hyd_hmm"  # the name as used by dram all lowercase leters and _
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
    location_collumns: list[str] = [
        HMM_KEY,
        FAA_KEY,
        HMM_SCORES_KEY,
        FAA_SCORES_KEY,
        DRAM2_MODUAL_KEY,
    ]

    def get_genome_summary(self) -> Path:
        """Get the ids from a complete annotations pandas DataFrame."""
        pass

    def search(self, fasta: Fasta) -> pd.DataFrame | pd.Series:
        """Perform a search, be that HMM, Blast or something else."""
        pass

    def load_dram_config(self):
        """Extract data from the larger DRAM config yaml"""
        pass

    def download(self, user_locations_dict):
        """Download your raw data, return a full location dictionary."""
        if HMMS_KEY not in user_locations_dict:
            hmm_path = self.working_dir / "CANT-HYD.hmm"
            download_file(
                (
                    "https://github.com/dgittins/CANT-HYD-HydrocarbonBiodegradation/raw/main/"
                    "HMMs/concatenated % 20HMMs/CANT-HYD.hmm"
                ),
                self.logger,
                hmm_path,
            )
        else:
            hmm_path = user_locations_dict[HMM_KEY]

        # The rest we just pull locally
        if FAA_KEY in user_locations_dict:
            faa_path = user_locations_dict[FAA_KEY]
        else:
            faa_path = (
                get_package_path(
                    Path(
                        "db_kits",
                        "cant_hyd_kit",
                        "data",
                        "BacMet_ExpVerified_BiocideRes_genes_SHORT.faa",
                    )
                ),
            )
        if FAA_SCORES_KEY in user_locations_dict:
            faa_scores_path = user_locations_dict[FAA_SCORES_KEY]
        else:
            faa_scores_path = get_package_path(
                Path(
                    "db_kits",
                    "cant_hyd_kit",
                    "data",
                    "BacMet_blast_db_biocideresist.csv",
                )
            )
        if HMM_SCORES_KEY in user_locations_dict:
            hmm_scores_path = user_locations_dict[HMM_SCORES_KEY]
        else:
            hmm_scores_path = get_package_path(
                Path("db_kits", "cant_hyd_kit", "data", "CANT-HYD_hmm_scores.csv")
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
    ) -> dict:
        """Do whatever you need to process the data at the locations provided."""
        pass
