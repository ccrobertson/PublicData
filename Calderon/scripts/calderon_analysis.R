library(pheatmap)
library(ggplot2)
library(edgeR)
library(GenomicRanges)

setwd(Sys.getenv('analysis'))

### get peaks
countMat = readRDS("counts_merged_peaks.rds")

### create smaple info data.frame
parseSampleId = function(x) {
  v = unlist(strsplit(x, split="_"))
  n = length(v)
  c = paste(v[seq(1,(n-3))], collapse="_")
  return(c(c, v[seq((n-2),n)]))
}
sampleDF = data.frame(t(sapply(colnames(countMat),parseSampleId)))
names(sampleDF) <- c("type","donor","condition","build")
sampleDF$sample_id = row.names(sampleDF)
sampleDF$category = NA
sampleDF$category[sampleDF$type %in% c("Myeloid_DCs","pDCs")] <- "DC"
sampleDF$category[sampleDF$type %in% c("Monocytes")] <- "Monocyte"
sampleDF$category[sampleDF$type %in% c("Immature_NK","Mature_NK","Memory_NK")] <- "NK"
sampleDF$category[sampleDF$type %in% c("Bulk_B","Mem_B","Plasmablasts","Naive_B")] <- "B"
sampleDF$category[sampleDF$type %in% c("CD8pos_T","Central_memory_CD8pos_T","Effector_CD4pos_T", "Effector_memory_CD8pos_T", "Follicular_T_Helper","Gamma_delta_T","Memory_Teffs","Memory_Tregs","Naive_CD8_T","Naive_Teffs","Naive_Tregs","Regulatory_T","Th1_precursors","Th17_precursors","Th2_precursors")] <- "T"
head(sampleDF)


### generate normalized countmatrix
#create edgeR object
y_raw <- DGEList(counts=countMat, samples=sampleDF)

#filter for low cpm (note, eventhough y contains raw counts, this function filters based on cpm)
min.count=10
keep.exprs <- filterByExpr(y=y_raw, group=y_raw$samples$type, min.count=min.count)
y <- y_raw[keep.exprs,, keep.lib.sizes=FALSE]

#calculate normalization factors
y <- calcNormFactors(y, method="TMM")

#rank normalize TMM counts
rankNormalize = function(x) {
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}
Y = cpm(y)
Ynorm = t(apply(Y, MARGIN=1, FUN=rankNormalize))

#kmeans clustering QC (should update this to be based on TMM not raw counts)
y_kmeans = t(apply(Ynorm, MARGIN=1, FUN=function(x) {kmeans(x, centers=2)$cluster}))
y_kmeans_summary = data.frame(cluster1=rowSums(y_kmeans==1), cluster2=rowSums(y_kmeans==2))
keep2 = row.names(y_kmeans_summary)[y_kmeans_summary$cluster1>1 & y_kmeans_summary$cluster2>1]
y_filtered <- y[keep2, , keep.lib.sizes=FALSE]
Yfiltered = cpm(y_filtered)
YfilteredNorm = t(apply(Yfiltered, MARGIN=1, FUN=rankNormalize))


### rank peaks by cell category specificity
getAOV = function(x) {
  summary(aov(x~sampleDF$category))[[1]][1,]
}
aovResults = apply(YfilteredNorm, MARGIN=1, FUN=getAOV)
aovResultsMat = do.call("rbind", aovResults)
saveRDS(aovResultsMat, file="aov_results_peak_by_category.rds")
pdf("calderon_analysis.pdf")
hist(aovResultsMat[,"Pr(>F)"])
dev.off()

names(aovResultsMat) <- c("Df","SumSq","MeanSq","F","pvalue")
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

#summary(mcols(atacGR)$mcols.F)
#summary(mcols(atacGR_overlap)$mcols.F)
