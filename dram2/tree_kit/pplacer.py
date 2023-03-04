"""
=================
Tools For PPlacer
=================
"""
import os
import json
import logging
import pandas as pd
from tempfile import TemporaryDirectory
from dram2.utils.utils import run_process

QUERY_TMP_NAME = "query"


class DramTree:
    def __init__(
        self,
        name,
        pplacer_profile: str,
        reference_seq: str,
        gene_mapping_path: str,
        color_mapping_path: str,
        target_dbs: list[set["str"]],
        target_ids: list["str"],
    ):
        self.name = name
        self.pplacer_profile = pplacer_profile
        self.reference_seq = reference_seq
        self.gene_mapping = pd.read_csv(gene_mapping_path, sep="\t", index_col=0)
        self.color_mapping = pd.read_csv(color_mapping_path, sep="\t", index_col=0)
        self.mapping = pd.merge(
            self.gene_mapping,
            self.color_mapping,
            left_on="call",
            right_index=True,
            how="outer",
        )
        # tree.mapping['call'].unique()
        self.target_ids = set(target_ids)
        self.target_dbs = target_dbs

    def set_logger(self, logger: logging.Logger):
        self.logger = logger

    def _align(self, fasta, working_dir):
        """ """
        self.logger.info("Aligning sequences with MAFFT")
        # The file at the end of the output path has the extension 'fasta'!
        output_path = os.path.join(working_dir, f"{QUERY_TMP_NAME}.fasta")
        if os.path.exists(output_path):
            self.logger.warn(
                f"The alignment file {output_path} already "
                "exists, so the program is being called twice in "
                "the same working directory."
            )
        run_process(
            ["mafft", "--add", fasta, "--reorder", self.reference_seq],
            self.logger,
            capture_stdout=False,
            save_output=output_path,
        )

        self.logger.info("Sequences aligned")
        return output_path

    def pplacer_place_sequences(self, query_faa: str, working_dir: str, threads: int):
        self.logger.info("Sequence placement with pplacer started")
        output_path = os.path.join(working_dir, f"{QUERY_TMP_NAME}.jplace")
        run_process(
            [
                "pplacer",
                "-p",
                "-c",
                self.pplacer_profile,
                self._align(query_faa, working_dir),
                "--out-dir",
                working_dir,
                "-j",
                str(threads),
            ],
            self.logger,
        )
        self.logger.info("Sequence placement with pplacer finished")

        return output_path


"""
One day we may want to be able to make these profiles

taxit_comand= ['taxit', 'create', '-l', '16s_rRNA', '-P', 'nxr_nar.refpkg',
 '--aln-fasta', '$NXR_NAR_FAA', '--tree-stats',
 './RAxML_info.nxr_nar_raxml', '--tree-file', './nxr_nar.tre']
"""
