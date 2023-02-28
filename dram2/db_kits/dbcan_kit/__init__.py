from os import path, stat
import re
import tarfile
from shutil import move, rmtree
from dram2.db_kits.fegenie_kit import process
from dram2.utils.utils import download_file, run_process
from dram2.db_kits.utils import (
    make_mmseqs_db,
    run_hmmscan,
    get_best_hits,
    BOUTFMT6_COLUMNS,
    DBKit,
    get_sig_row,
)
from dram2.utils.utils import Fasta


from pathlib import Path
from functools import partial
from sqlalchemy import Column, String
import logging
from typing import Optional
import pandas as pd

from dram2.db_kits.utils.sql_descriptions import SQLDescriptions, BASE

KNOWN_DBCAN_NONE_DESCRIPTORS = {"GT2_Glycos_transf_2", "GT2_Glyco_tranf_2_3"}

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
        description_data = pd.concat([line_reader(line) for line in f.readlines()])

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
        sql_descriptions.start_db_session()
        hits_df[f"{db_name}_hits"] = hits_df[f"{db_name}_ids"].apply(
            partial(description_pull, sql_descriptions=sql_descriptions)
        )
        hits_df[f"{db_name}_subfam_ec"] = hits_df[f"{db_name}_ids"].apply(
            lambda x: "; ".join(
                sql_descriptions.get_descriptions(
                    x.split("; "),
                    "ec", KNOWN_DBCAN_NONE_DESCRIPTORS
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
    version: str = VERSION
    citation: str = CITATION
    date: str = DATE
    hmm_db: Path
    description_db: SQLDescriptions

    def setup() -> dict:
        pass

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
        self.logger.info("{self.formal_name} looks ready to use!")

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

    def get_descriptions(self, hits):
        "fix"
        return hits

    def get_ids(self, annotations: pd.Series) -> list:
        main_id = "cazy_best_hit"
        if main_id in annotations:
            return [annotations[main_id].split("_")[0]]
        return []

    # "cazy_id": lambda x: [i.split("_")[0] for i in x.split("; ")],
    # "cazy_hits": lambda x: [
    #     f"{i[1:3]}:{i[4:-1]}" for i in re.findall(r"\(EC [\d+\.]+[\d-]\)", x)
    # ],
    # "cazy_subfam_ec": lambda x: [f"EC:{i}" for i in re.findall(r"[\d+\.]+[\d-]", x)],
