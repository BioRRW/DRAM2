"""General utilitys, avalible to all dram2 tools"""
import re
import subprocess
from os import path, stat
from urllib.request import urlopen, urlretrieve
from urllib.error import HTTPError
import pandas as pd
import logging
from typing import Callable
from os import path, getenv
from typing import NamedTuple
from pathlib import Path


import json

def load_config(alt_location:Path, logger: logging.Logger):
    """If all_loc is none the """
    if alt_location is not None:
        location = alt_location
    elif (envioment_location := Path(getenv('DRAM_CONFIG_LOCATION'))) is not None:
        location = envioment_location
    elif (user_location := (Path.home() / '.config' / 'dram2bio'/ 'config')).exists(): 
        location = user_location
    elif (global_location := Path("/etc", "dram2bio", "config")).exists(): 
        location = global_location
    else:
        logger.info(f"No config found any config that is created will go to {user_location}")
        config = {}
        config['config_location'] = user_location
        return config
    logger.info(f"Loading config from: {location}")
    config = json.loads(open(location).read())
    config['config_location'] = location
    return config



def download_file(url, logger, output_file=None, verbose=True):
    if verbose:
        print("downloading %s" % url)
    if output_file is None:
        return urlopen(url).read().decode("utf-8")
    else:
        try:
            urlretrieve(url, output_file)
        except HTTPError as error:
            logger.critical(f"Something went wrong with the download of the url: {url}")
            raise error


def setup_logger(logger, *log_file_paths, level=logging.INFO):
    logger.setLevel(level)
    formatter = logging.Formatter("%(asctime)s - %(message)s")
    # create console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    for log_file_path in log_file_paths:
        fh = logging.FileHandler(log_file_path)
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        # add the handlers to the logger
        logger.addHandler(fh)


def get_ids_from_annotations_by_row(data, logger):
    functions = {i: j for i, j in ID_FUNCTION_DICT.items() if i in data.columns}
    missing = [i for i in ID_FUNCTION_DICT if i not in data.columns]
    logger.info(
        "Note: the fallowing id fields "
        f"were not in the annotations file and are not being used: {missing},"
        f" but these are {list(functions.keys())}"
    )
    out = data.apply(
        lambda x: {
            i
            for k, v in functions.items()
            if not pd.isna(x[k])
            for i in v(str(x[k]))
            if not pd.isna(i)
        },
        axis=1,
    )
    return out


def get_ids_from_annotations_all(data, logger):
    data = get_ids_from_annotations_by_row(data, logger)
    data.apply(list)
    out = Counter(chain(*data.values))
    return out


def run_process(
    command,
    logger,
    shell: bool = False,
    capture_stdout: bool = True,
    save_output: str = None,
    check: bool = False,
    stop_on_error: bool = True,
    verbose: bool = False,
) -> str:
    """
    Standardization of parameters for using subprocess.run, provides verbose mode and option to run via shell
    """
    # TODO just remove check
    try:
        results = subprocess.run(
            command,
            check=check,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as error:
        logger.critical(f"The subcommand {command} experienced an error")
        if stop_on_error:
            raise error
    if results.returncode != 0:
        logger.critical(
            f"The subcommand {command} experienced an error: {results.stderr}"
        )
        logging.debug(results.stdout)
        if stop_on_error:
            raise subprocess.SubprocessError(
                f"The subcommand {' '.join(command)} experienced an error, see the log for more info."
            )

    if save_output is not None:
        with open(save_output, "w") as out:
            out.write(results.stdout)

    if capture_stdout:
        return results.stdout


def merge_files(files_to_merge, outfile, has_header=False):
    """It's in the name, if has_header assumes all files have the same header"""
    with open(outfile, "w") as outfile_handle:
        if has_header:
            outfile_handle.write(open(files_to_merge[0]).readline())
        for file in files_to_merge:
            with open(file) as f:
                if has_header:
                    _ = f.readline()
                outfile_handle.write(f.read())


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i : i + n]


def remove_prefix(text, prefix):
    if text.startswith(prefix):
        return text[len(prefix) :]
    return text  # or whatever


def remove_suffix(text, suffix):
    if text.endswith(suffix):
        return text[: -1 * len(suffix)]
    return text  # or whatever


def get_ordered_uniques(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x) or pd.isna(x))]



