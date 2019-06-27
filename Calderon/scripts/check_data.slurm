#!/bin/bash
#SBATCH --output='check_data.log'
#SBATCH --mem='200000'
#SBATCH --cpus-per-task='1'
#SBATCH --time='2-00:00:00'
#SBATCH --partition='standard'
#SBATCH --account rich_immunochip_impute

cd ${WORK}
for file in `ls sra/*.sra`; do 
	vdb-validate ${file}
done
