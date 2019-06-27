#!/bin/bash
#SBATCH --mem='20000'
#SBATCH --cpus-per-task='1'
#SBATCH --time='2-00:00:00'
#SBATCH --partition='standard'
#SBATCH --account rich_immunochip_impute

file=$1
module load htslib

#compress file
date 
echo ${file}
bgzip -f ${file}
wait
date
