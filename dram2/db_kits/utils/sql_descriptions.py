from os import path, remove, getenv
from pkg_resources import resource_filename
import json
import gzip
import logging
from shutil import copy2
import warnings
from datetime import datetime
from functools import partial
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from abc import ABC, abstractmethod
from typing import Callable

from pathlib import Path
import pandas as pd
from typing import Optional

from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base

BASE = declarative_base()


def create_description_db(db_loc):
    engine = create_engine("sqlite:///%s" % db_loc)
    BASE.metadata.create_all(engine)


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i : i + n]


def make_header_dict_from_mmseqs_db(mmseqs_db):
    mmseqs_headers_handle = open("%s_h" % mmseqs_db, "rb")
    mmseqs_headers = mmseqs_headers_handle.read().decode(errors="ignore")
    mmseqs_headers = [
        i.strip() for i in mmseqs_headers.strip().split("\n\x00") if len(i) > 0
    ]
    mmseqs_headers_split = []
    mmseqs_ids_unique = set()
    mmseqs_ids_not_unique = set()
    # TODO this could be faster with numpy
    for i in mmseqs_headers:
        header = {"id": i.split(" ")[0], "description": i}
        if header["id"] not in mmseqs_ids_unique:
            mmseqs_headers_split += [header]
            mmseqs_ids_unique.add(header["id"])
        else:
            mmseqs_ids_not_unique.add(header["id"])
    if len(mmseqs_ids_not_unique) > 0:
        warnings.warn(
            f"There are {len(mmseqs_ids_not_unique)} non unique headers "
            f"in {mmseqs_db}! You should definitly investigate this!"
        )
    return mmseqs_headers_split


class SQLDescriptions:
    def __init__(
        self,
        description_loc: Path,
        logger: logging.Logger,
        description_class,
        db_name: str,
    ):
        self.logger = logger
        self.description_class = description_class
        self.description_loc = description_loc
        self.db_name = db_name

    def setup(self):
        pass

    def start_db_session(self):
        engine = create_engine(f"sqlite:///{self.description_loc}")
        db_session = sessionmaker(bind=engine)
        self.session = db_session()

    # functions for adding descriptions to tables
    def add_descriptions_to_database(self, description_list, clear_table=True):
        if clear_table:
            self.session.query(self.description_class).delete()
        # TODO: try batching
        self.session.bulk_save_objects(
            [self.description_class(**i) for i in description_list]
        )
        self.session.commit()
        self.session.expunge_all()

    # functions for getting descriptions from tables
    def get_description(self, annotation_id, return_ob=False):
        return (
            self.session.query(self.description_class)
            .filter_by(id=annotation_id)
            .one()
            .description
        )

    def get_descriptions(
        self,
        ids: list | pd.Series,
        description_name="description",
        none_descriptors: Optional[set] = None,
    ):
        """
        Pull The descriptions from SQLDB
        ________________________________

        the dram database
        """
        if not isinstance(ids, list):
            ids: list = ids.dropna().values
        self.start_db_session()
        descriptions = [
            des
            for chunk in divide_chunks(list(ids), 499)
            for des in self.session.query(self.description_class)
            .filter(self.description_class.id.in_(chunk))
            .all()
        ]
        self.session.close()
        if len(descriptions) == 0:
            for i in list(ids):
                if none_descriptors is None or i not in none_descriptors:
                    self.logger.warn(
                        'No descriptions were found for your id\'s. Does the id "%s" look like an id from %s'
                        % (i, self.db_name)
                    )
                break

        return {i.id: i.__dict__[description_name] for i in descriptions}

    def populate_description_db(
        self, output_loc: Path, name: str, process_function: Callable | None = None
    ):
        self.start_db_session()
        # I don't think this is needed
        if output_loc.exists():
            remove(output_loc)
        create_description_db(output_loc)

        def check_db(db_name, db_function):
            self.add_descriptions_to_database(
                db_function(), f"{db_name}_description", clear_table=True
            )
            self.logger.info(f"Description updated for the {db_name} database")

        if process_function is None:
            process_function = partial(make_header_dict_from_mmseqs_db, name)

            check_db(name, process_function)
