from os import path, stat
import re
from dram2.utils import download_file, run_process, Fasta
from dram2.db_kits.utils import (
    run_hmmscan,
    DBKit,
    get_sig_row,
)


from pathlib import Path
from functools import partial
from sqlalchemy import Column, String
import logging
from typing import Optional
import pandas as pd
from shutil import move

from dram2.db_kits.utils.sql_descriptions import SQLDescriptions, BASE

KNOWN_DBCAN_NONE_DESCRIPTORS = {"GT2_Glycos_transf_2", "GT2_Glyco_tranf_2_3"}
DEFAULT_DBCAN_RELEASE = "11"
DEFAULT_DBCAN_DATE = "08062022"

VERSION = "11"
DATE = "08062022"
CITATION = (
    'Y. Yin, X. Mao, J. Yang, X. Chen, F. Mao, and Y. Xu, "dbcan'
    ": a web resource for automated carbohydrate-active enzyme an"
    'notation," Nucleic acids research, vol. 40, no. W1, pp. W44'
    "5–W451, 2012."
)
DESCRIPTION_COL = "description"
EC_COL = "subfam_ec"
# partial(
#     self.process_pfam_descriptions,
#     self.config.get("database_descriptions")["pfam_hmm"],
# ),
# partial(
#     self.process_vogdb_descriptions,
#     self.config.get("database_descriptions")["vog_annotations"],
# ),


class DbcanDescription(BASE):
    __tablename__ = "dbcan_description"

    id = Column(String(30), primary_key=True, nullable=False, index=True)

    description = Column(String(1000))
    ec = Column(String(1000))

    @property
    def serialize(self):
        return {
            "dbcan_id": self.id,
            DESCRIPTION_COL: self.description,
            EC_COL: self.ec,
        }


def find_best_dbcan_hit(genome: str, group: pd.DataFrame):
    group["perc_cov"] = group.apply(
        lambda x: (x["target_end"] - x["target_start"]) / x["target_length"], axis=1
    )
    group.sort_values("perc_cov", inplace=True)
    group.columns
    group.sort_values("full_evalue", inplace=True)
    return group.iloc[0]["target_id"]


def process_dbcan_descriptions(dbcan_fam_activities, dbcan_subfam_ec):
    def line_reader(line):
        if not line.startswith("#") and len(line.strip()) != 0:
            line = line.strip().split()
            if len(line) == 1:
                description = line[0]
            elif line[0] == line[1]:
                description = " ".join(line[1:])
            else:
                description = " ".join(line)
            return pd.DataFrame(
                {"id": line[0], DESCRIPTION_COL: description.replace("\n", " ")},
                index=[0],
            )

    with open(dbcan_fam_activities) as f:
        description_data = (pd.concat([line_reader(line) for line in f.readlines()])
                            .groupby("id")
                            .apply(lambda x: ",".join(x["description"].unique()))
                            )
        description_data = pd.DataFrame(description_data, columns=["description"]).reset_index()

    ec_data = pd.read_csv(
        dbcan_subfam_ec, sep="\t", names=["id", "id2", "ec"], comment="#"
    )[["id", "ec"]].drop_duplicates()
    ec_data = (
        pd.concat([ec_data["id"], ec_data["ec"].str.split("|", expand=True)], axis=1)
        .melt(id_vars="id", value_name="ec")
        .dropna(subset=["ec"])[["id", "ec"]]
        .groupby("id")
        .apply(lambda x: ",".join(x["ec"].unique()))
    )
    ec_data = pd.DataFrame(ec_data, columns=["ec"]).reset_index()
    data = pd.merge(description_data, ec_data, how="outer", on="id").fillna("")
    data.drop_duplicates()
    return [i.to_dict() for _, i in data.iterrows()]


def description_pull(x: str, sql_descriptions: SQLDescriptions):
    id_list = ([re.findall("^[A-Z]*[0-9]*", str(x))[0] for x in x.split("; ")],)
    id_list = [y for x in id_list for y in x if len(x) > 0]
    description_list = sql_descriptions.get_descriptions(
        id_list, DESCRIPTION_COL, KNOWN_DBCAN_NONE_DESCRIPTORS
    ).values()
    description_str = "; ".join(description_list)
    return description_str


def dbcan_hmmscan_formater(
    hits: pd.DataFrame, db_name: str, sql_descriptions: Optional[SQLDescriptions]
):
    """
    format the ouput of the dbcan database.

    Note these changes
    introduce new column for best hit from cazy database this will be the best
    hit above the already established threshold (0.35 coverage, e-18 evalue) -
    then for the distillate pull info from the best hit only

    introduce new column for corresponding EC number information from sub
    families (EC numbers are subfamily ECs)

    Make sure that ids and descriptions are separate (descriptions these are family based)

    :param hits:
    :param db_name:
    :returns:
    """
    hits_sig = hits[hits.apply(partial(get_sig_row, evalue_lim=1e-18), axis=1)]
    if len(hits_sig) == 0:
        # if nothing significant then return nothing, don't get descriptions
        return pd.DataFrame()
    hit_groups = hits_sig.groupby("query_id")
    all_hits = hit_groups.apply(
        lambda x: "; ".join(x["target_id"].apply(lambda y: y[:-4]).unique())
    )
    hits_df = pd.DataFrame(all_hits)
    hits_df.columns = [f"{db_name}_ids"]

    if sql_descriptions is not None:
        hits_df[f"{db_name}_hits"] = hits_df[f"{db_name}_ids"].apply(
            partial(description_pull, sql_descriptions=sql_descriptions)
        )
        hits_df[f"{db_name}_subfam_ec"] = hits_df[f"{db_name}_ids"].apply(
            lambda x: "; ".join(
                sql_descriptions.get_descriptions(
                    x.split("; "), "ec", KNOWN_DBCAN_NONE_DESCRIPTORS
                ).values()
            )
        )
    hits_df[f"{db_name}_best_hit"] = [find_best_dbcan_hit(*i) for i in hit_groups]
    hits_df.rename_axis(None, inplace=True)
    hits_df.columns
    return hits_df


class dbCANKit(DBKit):

    name = "dbcan"
    formal_name: str = "dbCAN"
    max_threads: int = 2
    version: str = VERSION
    citation: str = CITATION
    date: str = DATE
    hmm_db: Path
    description_db: SQLDescriptions
    location_keys: list[str] = [
        "dbcan_hmm",
        "dbcan_ec",
        "dbcan_fam",
    ]

    def load_dram_config(self):
        self.hmm_db = self.get_config_path("hmmdb")
        # self.dbcan_subfam_ec = self.get_config_path("dbcan_subfam_ec")
        # self.dbcan_fam_activities = self.get_config_path("dbcan_fam_activities")
        self.description_db = SQLDescriptions(
            self.get_config_path("description_db"),
            self.logger,
            DbcanDescription,
            self.name,
        )

    def search(self, fasta: Fasta) -> pd.DataFrame | pd.Series:

        if fasta.faa is None:
            raise ValueError(
                f"Fasta with out called genes faa was passed to the search function for {self.formal_name}"
            )
        self.logger.info("Getting hits from dbCAN")
        annotations = run_hmmscan(
            genes_faa=fasta.faa.absolute().as_posix(),
            db_loc=self.hmm_db.absolute().as_posix(),
            db_name=self.name,
            output_loc=self.working_dir.absolute().as_posix(),
            threads=self.threads,
            logger=self.logger,
            formater=partial(
                dbcan_hmmscan_formater,
                db_name=self.name,
                sql_descriptions=self.description_db,
            ),
        )
        return annotations

    def get_ids(self, annotations: pd.Series) -> list:
        main_id = f"{self.name}_best_hit"
        if main_id not in annotations:
            self.logger.debug(f"Expected {main_id} to be in annotations but not found")
        elif not pd.isna(annotations[main_id]):
            return [str(annotations[main_id]).split("_")[0]]
        return []

    def download(self, user_locations_dict: dict[str, Path]):
        """
        you know
        """
        if "dbcan_ec" in user_locations_dict:
            dbcan_subfam_ec = user_locations_dict["dbcan_ec"].as_posix()
        else:
            dbcan_subfam_ec = path.join(
                self.working_dir, f"CAZyDB.{DEFAULT_DBCAN_DATE}.fam.subfam.ec.txt"
            )
            url = (
                f"https://bcb.unl.edu/dbCAN2/download/Databases/"
                f"V{DEFAULT_DBCAN_RELEASE}/CAZyDB.{DEFAULT_DBCAN_DATE}.fam.subfam.ec.txt"
            )
            self.logger.info(f"Downloading dbCAN sub-family encumber from : {url}")
            download_file(url, dbcan_subfam_ec, self.logger)

        if "dbcan_hmm" in user_locations_dict:
            dbcan_subfam_ec = user_locations_dict["dbcan_hmm"].as_posix()
        else:
            dbcan_hmm = path.join(
                self.working_dir, f"dbCAN-HMMdb-V{DEFAULT_DBCAN_RELEASE}.txt"
            )
            link_path = f"http://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V{DEFAULT_DBCAN_RELEASE}.txt"
            self.logger.debug(f"Downloading dbCAN from: {link_path}")
            download_file(link_path, dbcan_hmm, self.logger)
        if "dbcan_fam" in user_locations_dict:
            dbcan_subfam_ec = user_locations_dict["dbcan_fam"].as_posix()
        else:
            dbcan_fam_activities = path.join(
                self.working_dir, f"CAZyDB.{DEFAULT_DBCAN_RELEASE}.fam-activities.txt"
            )
            url = (
                f"https://bcb.unl.edu/dbCAN2/download/Databases/"
                f"V{DEFAULT_DBCAN_RELEASE}/CAZyDB."
                f"{DEFAULT_DBCAN_DATE}.fam-activities.txt"
            )
            self.logger.info(f"Downloading dbCAN family activities from : {url}")
            download_file(url, dbcan_fam_activities, self.logger)
        return {
            "dbcan_hmm": Path(dbcan_hmm),
            "dbcan_ec": Path(dbcan_subfam_ec),
            "dbcan_fam": Path(dbcan_fam_activities),
        }

    def setup(
        self,
        location_dict: dict[str, Path],
        output_dir: Path,
    ) -> dict:
        hmm_path = location_dict["dbcan_hmm"]
        dbcan_ec = location_dict["dbcan_ec"]
        dbcan_fam = location_dict["dbcan_fam"]
        hmm_out = output_dir / path.basename(hmm_path)
        description_out = output_dir / "dbcan.sqlight"
        move(hmm_path, hmm_out)
        run_process(["hmmpress", "-f", hmm_out], self.logger)
        self.logger.info("dbCAN database processed")

        description_db = SQLDescriptions(
            description_out,
            self.logger,
            DbcanDescription,
            self.name,
        )
        description_db.populate_description_db(
            description_out,
            self.name,
            partial(
                process_dbcan_descriptions,
                dbcan_ec,
                dbcan_fam,
            ),
        )
        return {
            "description_db": {"location": description_out.relative_to(output_dir).as_posix()},
            "hmmdb": {"location": hmm_out.relative_to(output_dir).as_posix()},
        }
