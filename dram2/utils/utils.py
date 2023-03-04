"""General utilitys, avalible to all dram2 tools"""
import subprocess
from dataclasses import dataclass
from urllib.request import urlopen, urlretrieve
from urllib.error import HTTPError
import pandas as pd
import logging

# from os import getenv
from typing import Optional
from pathlib import Path

# import json


# def load_config(alt_location:Path, logger: logging.Logger):
#     """If all_loc is none the """
#     if alt_location is not None:
#         location = alt_location
#     elif (envioment_location := Path(getenv('DRAM_CONFIG_LOCATION'))) is not None:
#         location = envioment_location
#     elif (user_location := (Path.home() / '.config' / 'dram2bio'/ 'config')).exists():
#         location = user_location
#     elif (global_location := Path("/etc", "dram2bio", "config")).exists():
#         location = global_location
#     else:
#         logger.info(f"No config found any config that is created will go to {user_location}")
#         config = {}
#         config['config_location'] = user_location
#         return config
#     logger.info(f"Loading config from: {location}")
#     config = json.loads(open(location).read())
#     config['config_location'] = location
#     return config


class DramUsageError(Exception):
    "Raised when dram is not used corectly, usally it means you are missing a step"
    pass


def get_package_path(local_path: Path):
    """
    Locate the package data or non python files

    :param local_path:
    :returns:
    """
    abs_snake_path = Path(__file__).parent.parent.absolute() / local_path
    return abs_snake_path


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


def get_annotation_ids_by_row(data, logger):
    # functions = {i: j for i, j in ID_FUNCTION_DICT.items() if i in data.columns}
    # missing = [i for i in ID_FUNCTION_DICT if i not in data.columns]
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


def get_all_annotation_ids(data, logger):
    data = get_ids_from_annotations_by_row(data, logger)
    data.apply(list)
    out = Counter(chain(*data.values))
    return out


def export_posible_path(
    path: Optional[Path], relative_path: Optional[Path] = None
) -> Optional[str]:
    if path is None:
        return None
    out_path = path.absolute()
    if relative_path is not None and relative_path in out_path.parents:
        out_path = out_path.relative_to(relative_path)
    return out_path.as_posix()


@dataclass
class Fasta:
    name: Optional[str]
    origin: Optional[Path]
    tmp_dir: Optional[Path]
    faa: Optional[Path]
    fna: Optional[Path]
    gff: Optional[Path]
    mmsdb: Optional[Path]

    def export(self, output_dir):
        return (
            self.name,
            export_posible_path(self.origin),
            export_posible_path(self.tmp_dir, output_dir),
            export_posible_path(self.faa, output_dir),
            export_posible_path(self.fna, output_dir),
            export_posible_path(self.gff, output_dir),
            export_posible_path(self.mmsdb, output_dir),
        )

    @classmethod
    def import_strings(
        cls,
        relative_path: Path,
        name: str,
        origin: str,
        tmp_dir: str,
        faa: str,
        fna: str,
        gff: str,
        mmsdb: str,
    ):
        ob = cls(
            name,
            import_posible_path(origin),
            import_posible_path(tmp_dir, relative_path),
            import_posible_path(faa, relative_path),
            import_posible_path(fna, relative_path),
            import_posible_path(gff, relative_path),
            import_posible_path(mmsdb, relative_path),
        )
        return ob


def import_posible_path(
    path: Optional[str], relative_path: Optional[Path] = None
) -> Optional[Path]:
    if path is None:
        return None
    out_path = Path(path)
    if relative_path is None:
        return out_path.absolute()
    return (relative_path / out_path).absolute()


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
