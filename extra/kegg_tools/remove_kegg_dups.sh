#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --mem=40gb
#SBATCH --time=12-00:00:00
#SBATCH --job-name=mk_new_kegg
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=rmflynn@colostate.edu
#SBATCH --partition=debug
#SBATCH --output=remove_kegg_dups_%j.out

eval "$(conda shell.bash hook)"
source /opt/Miniconda2/miniconda2/bin/activate scripts

python ./remove_kegg_dups.py
