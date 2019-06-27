#!/bin/bash


cd ${WORK}/fastq

#compress file
for file in `ls *.fastq`; do 
	date
	echo ${file}
	sbatch ${SCRIPTS}/gzip_fastq.slurm ${file} 
	echo -e " "
done 