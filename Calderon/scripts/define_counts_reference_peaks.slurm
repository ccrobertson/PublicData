#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --account rich_immunochip_impute
#SBATCH --partition=standard
#SBATCH --mem=10000


module load gcc/7.1.0
module load R/3.5.1

Rscript ${scripts}/define_counts_reference_peaks.R
