#!/bin/bash
#SBATCH --output='convert_to_fastq.log'
#SBATCH --mem='200000'
#SBATCH --cpus-per-task='1'
#SBATCH --time='2-00:00:00'
#SBATCH --partition='standard'
#SBATCH --account rich_immunochip_impute


#load sratoolkit
module load sratoolkit/2.9.1
cd ${WORK}

#convert to fastq
for file in `ls sra/*.sra`; do
	date 
	echo ${file}
	fasterq-dump ${file} --outdir fastq/ --verbose --force
	wait
	date
	echo -e " "
done 