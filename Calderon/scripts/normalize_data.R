library(edgeR)
library(limma)
library(ggplot2)
library(Rtsne)

setwd(Sys.getenv('analysis'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

sampleDF$subcategory_cdtype = NA
sampleDF$subcategory_cdtype[sampleDF$type %in% c("Effector_CD4pos_T","Follicular_T_Helper","Memory_Teffs","Memory_Tregs","Naive_Teffs","Naive_Tregs","Regulatory_T","Th1_precursors","Th17_precursors","Th2_precursors")] <- "CD4"
sampleDF$subcategory_cdtype[sampleDF$type %in% c("CD8pos_T","Central_memory_CD8pos_T","Effector_memory_CD8pos_T","Naive_CD8_T")] <- "CD8"

sampleDF$subcategory_memory = NA
sampleDF$subcategory_memory[sampleDF$type %in% c("Memory_NK","Mem_B","Central_memory_CD8pos_T","Effector_memory_CD8pos_T","Memory_Teffs","Memory_Tregs")] <- "Memory"
sampleDF$subcategory_memory[sampleDF$type %in% c("Naive_B","Naive_CD8_T","Naive_Teffs","Naive_Tregs")] <- "Naive"

saveRDS(sampleDF, file="sampleDF.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filter low count peaks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#create edgeR object
Y_raw <- DGEList(counts=countMat, samples=sampleDF)

#filter for low cpm (note, eventhough y contains raw counts, this function filters based on cpm)
min.count=10
keep.exprs <- filterByExpr(y=Y_raw, group=Y_raw$samples$type, min.count=min.count)
Y_filtered <- Y_raw[keep.exprs,, keep.lib.sizes=FALSE]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# normalize
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#calculate normalization factors
Y_filtered <- calcNormFactors(Y_filtered, method="TMM")
Y_cpm = cpm(Y_filtered)
Y_lcpm = cpm(Y_filtered, log=TRUE)

#rank normalize TMM counts
rankNormalize = function(x) {
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}
Yranknorm = t(apply(Y_cpm, MARGIN=1, FUN=rankNormalize))

#voom normalize
mm = model.matrix(~0 + sampleDF$type + sampleDF$condition)
Yvoomnorm = voom(Y_cpm, mm)

pdf('rank_vs_voom_normalization.pdf')
i = 1
plot(Y_cpm[i,],Y_filtered$counts[i,], main=rownames(Y_cpm)[i])
plot(Y_cpm[i,],Y_lcpm[i,], main=rownames(Y_cpm)[i])
plot(Y_cpm[i,],Yranknorm[i,], main=rownames(Y_cpm)[i])
plot(Y_cpm[i,], Yvoomnorm$E[i,], main=rownames(Y_cpm)[i])
i = 100000
plot(Y_cpm[i,],Y_filtered$counts[i,], main=rownames(Y_cpm)[i])
plot(Y_cpm[i,],Y_lcpm[i,], main=rownames(Y_cpm)[i])
plot(Y_cpm[i,],Yranknorm[i,], main=rownames(Y_cpm)[i])
plot(Y_cpm[i,], Yvoomnorm$E[i,], main=rownames(Y_cpm)[i])
dev.off()


saveRDS(Yvoomnorm, file="counts_voomnormalized.rds")
saveRDS(Yranknorm, file="counts_ranknormalized.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# remove batch effect
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Yvoomnorm = readRDS(file="counts_voomnormalized.rds")
Yranknorm = readRDS(file="counts_ranknormalized.rds")
sampleDF = readRDS(file="sampleDF.rds")

Yvoomnorm_nobatch = removeBatchEffect(Yvoomnorm$E, batch=sampleDF$donor)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# exploratory plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
runMDS = function(Y, main) {
  mds = plotMDS(Y, plot=FALSE)
  mds_df = data.frame(x=mds$x, y=mds$y, donor=sampleDF$donor, type=sampleDF$type, category=sampleDF$category, condition=sampleDF$condition)
  p1 = ggplot(mds_df) + geom_point(aes(x=x,y=y,colour=as.factor(donor))) + ggtitle(main)
  p2 = ggplot(mds_df) + geom_point(aes(x=x,y=y,colour=as.factor(category))) + ggtitle(main) + facet_grid(.~condition)
  print(p1)
  print(p2)
}

runTSNE = function(Y, main) {
  tsne = Rtsne(t(Y))
  tsne_df = data.frame(tsne$Y,donor=sampleDF$donor, type=sampleDF$type, category=sampleDF$category, condition=sampleDF$condition)
  p1 = ggplot(tsne_df) + geom_point(aes(x=X1,y=X2,colour=as.factor(donor))) + ggtitle(main)
  p2 = ggplot(tsne_df) + geom_point(aes(x=X1,y=X2,colour=as.factor(category))) + ggtitle(main) + facet_grid(.~condition)
  print(p1)
  print(p2)
}

#get top 20% of variable peaks
variances=NULL
medians=NULL
for (i in 1:dim(Yvoomnorm$E)[1]) {
  variances[i] = var(Yvoomnorm$E[i,])
  medians[i] = median(Yvoomnorm$E[i,])
}
top = variances > quantile(variances,probs=0.90)

pdf("mds_plots_voom.pdf")
runMDS(Yvoomnorm$E[top,], main="Before Batch Removal (MDS)")
runMDS(Yvoomnorm_nobatch[top,], main="After Batch Removal (MDS)")
runTSNE(Yvoomnorm$E[top,], main="Before Batch Removal (tSNE)")
runTSNE(Yvoomnorm_nobatch[top,], main="After Batch Removal (tSNE)")
dev.off()

pdf("mds_plots_ranknorm.pdf")
runTSNE(Yranknorm[top,], main="Before Batch Removal (tSNE)")
dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# kmeans QC (remove peaks with outliers)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kmeansQC = function(Y) {
  #input: Y is a matrix (each row is a peak, each column is a sample)
  y_kmeans = t(apply(Y, MARGIN=1, FUN=function(x) {kmeans(x, centers=2)$cluster}))
  y_kmeans_summary = data.frame(cluster1=rowSums(y_kmeans==1), cluster2=rowSums(y_kmeans==2))
  keep = row.names(y_kmeans_summary)[y_kmeans_summary$cluster1>1 & y_kmeans_summary$cluster2>1]
  return(Y[keep,])
}

Yranknorm_clean = kmeansQC(Yranknorm)
Yvoomnorm_clean = kmeansQC(Yvoomnorm_nobatch)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### save filtered and normalized data
saveRDS(Yranknorm_clean, file="counts_filtered_ranknormalized.rds")
saveRDS(Yvoomnorm_clean, file="counts_filtered_voomnormalized.rds")
