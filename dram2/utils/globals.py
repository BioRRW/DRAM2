"""Store global variables, in this modual by default"""

"""defaults for prodigal"""
from pathlib import Path

PRODIGAL_MODE_DFT = "meta"
PRODIGAL_MODE_CHOICES = ["train", "meta", "single"]
PRODIGAL_TRANS_TABLE_CHOICES = [str(i) for i in range(1, 26)]
PRODIGAL_TRANS_TABLE_DFT = 11

MIN_CONTIG_SIZE_DFT = 2500
BIT_SCORE_THRESHOLD_DFT = 60
RBH_BIT_SCORE_THRESHOLD_DFT = 350

GENOMES_PER_PRODUCT = 1000

# all piplines
DEFAULT_FORCE: bool = False
DEFAULT_OUTPUT_DIR: Path = Path(".")
