#!/bin/bash

nickname=$1


cd ${WORK}


ls -1 fastq/ | grep "fastq.gz" | sed 's/_[0-9].fastq.gz//g' | sort | uniq > sample_list

cat sample_list | awk 'BEGIN {OFS=","; print "sample_name","protocol","organism","read1","read2","read_type"} {print $1, "ATAC", "human", "R1","R2", "paired"}' > ${nickname}.csv

sed "s/NICKNAME/${nickname}/g" ${SCRIPTS}/template.yaml > ${nickname}.yaml