library(stringr)
#library(tidyverse)
library(DESeq2)
library(lattice)
library(bigWig)

setwd(paste0(Sys.getenv("WORK"),"/analysis"))
#source Mike's functions
#source("../scripts/mikes_functions.R")

#read reference peaks
merged_peaks = read.table("merged_peaks.bed")
row.names(merged_peaks) <- apply(merged_peaks, MARGIN=1, FUN= function(x) {paste(str_trim(as.character(x)), collapse="_")} )


#read bigWig paths into a vector
paths_to_bigWigs = scan("paths_to_bigWigs.txt", what="character")
#reassign vector names to sample ids
names(paths_to_bigWigs) <- sapply(paths_to_bigWigs, function(a) {gsub(".*/([^/]+)_exact.bw", "\\1", a)}, USE.NAMES=FALSE)


#get counts for reference peaks
get.raw.counts.interval = function(bigWigFile, ranges) {
	loaded.bw = load.bigWig(bigWigFile)
	counts = bed.region.bpQuery.bigWig(bw=loaded.bw, bed=ranges)
}
counts.df = sapply(paths_to_bigWigs, get.raw.counts.interval, merged_peaks)
row.names(counts.df) <- row.names(merged_peaks)


#save count matrix
saveRDS(counts.df, file="counts_merged_peaks.rds")
