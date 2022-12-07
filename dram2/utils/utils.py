"""General utilitys, avalible to all dram2 tools"""
import re
import subprocess
from os import path, stat
from urllib.request import urlopen, urlretrieve
from urllib.error import HTTPError
import pandas as pd
import logging
from typing import Callable



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


def get_basic_description(hits, header_dict, db_name="viral"):
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
        verbose=False,
    )
    processed_results = [
        i.strip() for i in results.strip().split(split_char) if len(i) > 0
    ]
    # remove(hits_file)
    return {i.split()[0]: i for i in processed_results if i != ""}


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
    hits: pd.DataFrame, db_name: str, hmm_info_path: str = None, top_hit: bool = True
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
