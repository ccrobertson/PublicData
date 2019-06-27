#!/bin/bash
#SBATCH --output='fetch_data.log'
#SBATCH --mem='8000'
#SBATCH --cpus-per-task='1'
#SBATCH --time='2-00:00:00'
#SBATCH --partition='standard'
#SBATCH --account rich_immunochip_impute

#load sratoolkit
module load sratoolkit/2.9.1

#change default path
mkdir ${HOME}/.ncbi
echo "/repository/user/main/public/root = \"${PWD}\"" > ${HOME}/.ncbi/user-settings.mkfg

#SRA Project page SRP156496
cd ${WORK}
prefetch $(<SraAccList.txt)
