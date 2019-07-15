
# Install WGCNA
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("GO.db", "preprocessCore", "impute"))
#install.packages(c("matrixStats", "Hmisc", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
#biocLite(c("GO.db", "KEGG.db", "topGO", "org.Hs.sgd.db","org.Hs.eg.db", "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
#install.packages("WGCNA")

# Setup
library(WGCNA)
library(pheatmap)
setwd(Sys.getenv('analysis'))
enableWGCNAThreads()

# Read in data
Y = readRDS(file="counts_filtered_voomnormalized.rds")
sampleDF = readRDS(file="sampleDF.rds")

# Try kmeans clustering
kclusters = list()
kvals = seq(5,25)
for (i in 1:length(kvals)) {
  k = kvals[i]
  cat("k = ", k,"\n")
  kclusters[[i]] = kmeans(Y, centers=k)
}
kclusters = kclusters[1:7]
fit = sapply(kclusters, function(m) m[["betweenss"]]/m[["totss"]])
clust = data.frame(sapply(kclusters, function(m) paste0("cluster",m[["cluster"]])))
names(clust) = paste0("k_",kvals)
row.names(clust) <- names(kclusters[[1]]$cluster)

#get top 20% of variable peaks
variances=NULL
medians=NULL
for (i in 1:dim(Y)[1]) {
  variances[i] = var(Y[i,])
  medians[i] = median(Y[i,])
}
top = variances > quantile(variances,probs=0.999)

pdf("wgcna_atac_kmeans.pdf")
plot(kvals, fit)
pheatmap(Y[top,],
  cluster_rows=TRUE, cluster_cols=TRUE,
  show_rownames=FALSE,show_colnames=FALSE,
  annotation_row=clust[,c("k_6","k_10")],
  annotation_col=sampleDF[,c("donor","category","condition")],
  treeheight_row = 0, treeheight_col = 0)
dev.off()





# Transpose
Y_subset = Y[sample(1:dim(Y)[1],1000),]
datExpr = as.data.frame(t(Y_subset))
datTraits = data.frame(type=factor(sampleDF$type), donor = factor(sampleDF$donor), condition=factor(sampleDF$condition))
datTraitsNumeric = sapply(datTraits, MARGIN=2, FUN=as.numeric)

# Cluster samples
sampleTree = hclust(dist(datExpr), method = "average")
pdf("wgcna_atac_sample_cluster.pdf")
#plotDendroAndColors(sampleTree, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
traitColors = labels2colors(datTraitsNumeric)
plotDendroAndColors(sampleTree, colors = traitColors, main = "Sample dendrogram and trait heatmap", cex.dendroLabels=0.2)
dev.off()

# Create network
type = "Bulk_B"
samples = sampleDF$sample_id[sampleDF$type==type]

donor = "1001"
samples = sampleDF$sample_id[sampleDF$donor==donor && sampleDF$condition==1]

datExpr_sub = datExpr[samples,]
datTrait_sub = datTraits[samples,]

#choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#plot scale-free topology fit index as a function of the soft-thresholding power
sft = pickSoftThreshold(datExpr_sub, powerVector = powers, verbose = 5)
pdf("wgcna_atac_softthresh.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") # this line corresponds to using an R^2 cut-off of h
#plot mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
