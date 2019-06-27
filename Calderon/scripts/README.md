

### Setup
```bash
module load gcc/7.1.0
module load R/3.5.1
export WORK=${PUBLIC_DATA}/Calderon/atac_data
export scripts=${PUBLIC_DATA}/Calderon/scripts
export analysis=${WORK}/analysis
cd ${WORK}
```



### Download and prepare data
Fetch data. Requires list of sra accession numbers (e.g. 'SRR7650803') in file called SraAccList.txt
```bash
sbatch ${scripts}/fetch_data.slurm
```

Make sure data download was complete
```bash
sbatch ${scripts}/check_data.slurm
grep -v 'ok\|consistent' check_data.log
```

Convert to fastq
```bash
sbatch ${scripts}/convert_to_fastq.slurm
```

Compress fastq files
```bash
bash ${scripts}/gzip_fastq.sh > gzip_fastq2.log 2>&1
```

Prep config files for pepatac
```bash
bash ${scripts}/prep_config_files.sh calderon_atac
```

### Run pepatac
```bash
looper run calderon_atac.yaml --dry-run --compute singularity_slurm
looper run calderon_atac.yaml --compute singularity_slurm
```

### Map sra ids to meaningful names
New naming scheme is `{cell_type}_{donor}_{replicate}`
  ```bash
  mkdir -p ${WORK}/results_combined
  mkdir -p ${WORK}/results_combined/bigwigs
  mkdir -p ${WORK}/results_combined/peaks
  awk 'BEGIN {a[$12"_"$13]=0} NR>1 {a[$12"_"$13]++; print $9".sra", $12"_"$13"_"a[$12"_"$13]}' SraRunTable.txt > sra_to_sampleid
  while read old new; do
    oldbw=${WORK}/output/results_pipeline/${old}/aligned_hg38_exact/${old}_exact.bw
    newbw=${WORK}/results_combined/bigwigs/${new}_hg38_exact.bw
    oldpk=${WORK}/output/results_pipeline/${old}/peak_calling_hg38/${old}_peaks.narrowPeak
    newpk=${WORK}/results_combined/peaks/${new}_hg38_peaks.narrowPeak
    echo $old "-->" $new
    ln -s ${oldbw} ${newbw}
    ln -s ${oldpk} ${newpk}
  done < sra_to_sampleid
  ```

### Analysis 
  ```bash
  bedops --merge `ls ${WORK}/results_combined/peaks/*narrowPeak | tr '@\n' ' '` > ${analysis}/merged_peaks.bed
  ls ${WORK}/results_combined/bigwigs/*_exact.bw | tr '@\n' ' ' > ${analysis}/paths_to_bigWigs.txt
  sbatch ${scripts}/define_counts_reference_peaks.slurm
  Rscript ${scripts}/calderon_analysis.R
  ```
