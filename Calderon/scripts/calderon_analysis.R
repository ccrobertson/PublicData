library(pheatmap)
library(ggplot2)
library(edgeR)
library(GenomicRanges)

setwd(Sys.getenv('analysis'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ANALYZE CELL-TYPE SPECIFICITY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
YfilteredNorm = readRDS(file="counts_filtered_ranknormalized.rds")
sampleDF = readRDS(file="sampleDF.rds")

### rank peaks
runAOV = function(y, x) {
  summary(aov(y~x))[[1]][1,]
}
getAOVResults = function(Y,x) {
  aovRes = apply(Y, MARGIN=1, FUN=runAOV, x=x)
  aovResMat = do.call("rbind", aovRes)
  names(aovResultsMat) <- c("Df","SumSq","MeanSq","F","pvalue")
}


#by broadest cell category
aovRes.category = getAOVResults(YfilteredNorm, sampleDF[,"category"])
saveRDS(aovRes.category, file="aov_results_peak_by_category.rds")

#by T cell subtype
t_cells = row.names(sampleDF[sampleDF$category=="T",])
aovRes.T_cell_types = getAOVResults(YfilteredNorm[,t_cells], sampleDF[t_cells,"type"])
saveRDS(aovRes.T_cell_types, file="aov_results_peak_by_T_subtype.rds")

#aovRes.T_cell_types = apply(YfilteredNorm[,t_cells], MARGIN=1, FUN=getAOV, x=sampleDF[t_cells,"type"])
#aovResMat.T_cell_types = do.call("rbind", aovRes.T_cell_types
#names(aovResMat.T_cell_types) <- c("Df","SumSq","MeanSq","F","pvalue")
#saveRDS(aovResMat.T_cell_types, file="aov_results_peak_by_T_subtype.rds")

#by B cell subtype
b_cells = row.names(sampleDF[sampleDF$category=="B",])
aovRes.B_cell_types = getAOVResults(YfilteredNorm[,b_cells], sampleDF[b_cells,"type"])
saveRDS(aovRes.B_cell_types, file="aov_results_peak_by_B_subtype.rds")


pdf("calderon_analysis.pdf")
par(mfrow=c(2,2))
hist(aovRes.category[,"Pr(>F)"])
hist(aovRes.T_cell_types[,"Pr(>F)"])
hist(aovRes.B_cell_types[,"Pr(>F)"])
dev.off()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# VISUALIZE PEAKS AND OVERLAP WITH T1D CREDIBLE SNPS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
YfilteredNorm = readRDS(file="counts_filtered_ranknormalized.rds")
sampleDF = readRDS(file="sampleDF.rds")
aovResultsMat = readRDS(file="aov_results_peak_by_category.rds")
head(aovResultsMat[order(aovResultsMat$F, decreasing=TRUE),], n=100)

### get region info
getRegion = function(x){
  s = strsplit(x,split="_")[[1]]
  n = length(s)
  c(s[1], s[(n-1)],s[n])
}
regDF = data.frame(do.call("rbind",lapply(row.names(aovResultsMat),getRegion)))
row.names(regDF) = row.names(aovResultsMat)
names(regDF) <- c("chromosome", "start","end")
regDF$region = row.names(regDF)
regDF = merge(regDF, aovResultsMat, by="row.names")
row.names(regDF) = regDF$region

#get credible variants
#credible_variants = scan("mega_credible_variants_ppsum_0.5.txt", what="character")
credDF_all = read.table("/nv/vol185/MEGA/release4/IMPUTED_TOPMED/fine_mapping/REDO/all_creds.txt", header=TRUE)
credDF = credDF_all[credDF_all$ppsum>0.8,]

#create genomic ranges objects
atacGR = GRanges( seqnames = regDF$chromosome, ranges = IRanges(start=as.numeric(regDF$start), end=as.numeric(regDF$end)), mcols = regDF)
credGR = GRanges( seqnames = paste0("chr",credDF$chromosome), ranges = IRanges(start=credDF$position, end=credDF$position+1), mcols = credDF)

sum(countOverlaps(atacGR, credGR))
credGR_overlap = subsetByOverlaps(credGR, atacGR)
atacGR_overlap = subsetByOverlaps(atacGR, credGR)

credDF$in_atac_peak = credDF$MarkerName %in% mcols(credGR_overlap)$mcols.MarkerName
regDF$in_cred_set = regDF$region %in% mcols(atacGR_overlap)$mcols.region

credDF_overlap = credDF[credDF$in_atac_peak==TRUE,]
length(unique(credDF_overlap$MarkerName))
length(unique(credDF_overlap$tag))
length(unique(credDF$tag))
regDF_overlap = regDF[regDF$in_cred_set==TRUE,]
length(unique(regDF_overlap$region))
wilcox.test(regDF$F ~ regDF$in_cred_set)$p.value

pdf("atacGR_overlap.pdf")
#ggplot(regDF) + geom_density(aes(x=F))
regDF$T1D_credible = factor(regDF$in_cred_set, levels=c(TRUE, FALSE), labels=c("Overlapping", "Non-overlapping"))
ggplot(regDF) + geom_density(aes(x=F, colour=T1D_credible, fill=T1D_credible), alpha=0.5) + xlab("ANOVA F statistic for cell-specificity")
#regDF$F ~ regDF$in_cred_set)
dev.off()

pdf("atacGR_heatmap.pdf")
top100peaks = row.names(head(aovResultsMat[order(aovResultsMat$F, decreasing=TRUE),], n=100))
top100peaksMat = YfilteredNorm[top100peaks,]
pheatmap(top100peaksMat, cluster_rows=TRUE, show_rownames=FALSE,show_colnames=FALSE,
         cluster_cols=TRUE, annotation_col=sampleDF[,c("donor","category","condition")],treeheight_row = 0, treeheight_col = 0)
pheatmap(YfilteredNorm[regDF_overlap$region,], cluster_rows=TRUE, show_rownames=TRUE,show_colnames=FALSE,
        cluster_cols=TRUE, annotation_col=sampleDF[,c("donor","category","condition")], treeheight_row = 0, treeheight_col = 0, fontsize_row=2)
dev.off()

pdf("atacGR_T_cells_heatmap.pdf")
top100peaksTcells = row.names(head(aovResMat.T_cell_types[order(aovResMat.T_cell_types$F, decreasing=TRUE),], n=100))
top100peaksTcellsMat = YfilteredNorm[top100peaksTcells,]
pheatmap(top100peaksTcellsMat, cluster_rows=TRUE, show_rownames=FALSE,show_colnames=FALSE,
        cluster_cols=TRUE, annotation_col=sampleDF[,c("donor","type","condition")],treeheight_row = 0, treeheight_col = 0)
dev.off()
#summary(mcols(atacGR)$mcols.F)
#summary(mcols(atacGR_overlap)$mcols.F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RUN SIMULATION OF IMMUNOCHIP OVERLAP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load("/nv/vol185/MEGA/release4/IMPUTED_TOPMED/meta_analysis/resultsmeta_5pcs_vcf_5pc_diff.RData")
sites = res[,c("MarkerName","chromosome","position")]
atacGR = GRanges( seqnames = regDF$chromosome, ranges = IRanges(start=as.numeric(regDF$start), end=as.numeric(regDF$end)), mcols = regDF)

getSimulatedRes = function(seed, atacGR, atacDF) {

  #define genomic ranges object based on random set of snps
  set.seed(seed)
  #randomDF = sites[base::sample(x=dim(sites)[1], size=dim(credDF)[1]),]
  randomDF = sites[base::sample(x=dim(sites)[1], size=length(unique(credDF$tag))),]
  randomGR = GRanges( seqnames = paste0("chr",randomDF$chromosome), ranges = IRanges(start=randomDF$position, end=randomDF$position+1), mcols = randomDF)

  #get snps and peaks that overlap
  randomGR_overlap = subsetByOverlaps(randomGR, atacGR)
  atacGR_overlap = subsetByOverlaps(atacGR, randomGR)

  #create indicators for snps overlapping peaks (or peaks overlapping snps)
  randomDF$in_atac_peak = randomDF$MarkerName %in% mcols(randomGR_overlap)$mcols.MarkerName
  atacDF$in_random_set = atacDF$region %in% mcols(atacGR_overlap)$mcols.region

  #test if F statistic is different between overlapping and non overlapping regions
  wilcox_out = wilcox.test(atacDF$F ~ atacDF$in_random_set)
  overlap = sum(countOverlaps(atacGR, randomGR))
  W = wilcox_out$statistic
  pval = wilcox_out$p.value

  return(c(overlap, W, pval))
  #creat data.frames with only overlap
  #randomDF_overlap = randomDF[randomDF$in_atac_peak==TRUE,]
  #regionsDF_overlap = regionsDF[regionsDF$in_random_set==TRUE,]
}

simulation = list()
for (i in 1:1000) {
  cat(i,"\n")
  simulation[[i]] = getSimulatedRes(i, atacGR=atacGR, atacDF=regDF)
}
simulationDF = as.data.frame(do.call("rbind", simulation))
names(simulationDF) <- c("overlap","W","pval")


tagDF = credDF[credDF$snp%in%credDF$tag,]

pdf('simulations.pdf')
#ggplot() + geom_density(aes(x=-log10(simulationDF$pval))) + geom_vline(xintercept=-log10(wilcox.test(regDF$F ~ regDF$in_cred_set)$p.value))
ggplot() + geom_density(aes(x=simulationDF$overlap)) + geom_vline(xintercept=23)
dev.off()

# Summary of overlap
# Twenty-three of 74 index snps (tag snps) overlap an accessible chromatin site.
# To determine how this number compares to what we could expect by chance due to the immune-enrichment of the Immunochip,
# we sampled 74 independent variants from the imputed immunochip genotypes 1000 times.
# The median number of variants overlapping accessible chromatin peaks in immune cell types
# was [insert]. Only # of 1000 simulations had more than 15 variants overlapping open chromatin.
# In contrast, we saw no difference in frequency of overlap with muscle, adipose, or __ cell accessible sites between T1D index SNPs
# randomly samples SNP sets.

# Summary of cell-specificity of overlap
