
# Setup
library(pheatmap)
setwd(Sys.getenv('analysis'))

# Read in data
Y = readRDS(file="counts_filtered_voomnormalized.rds")
sampleDF = readRDS(file="sampleDF.rds")


#get top 20% of variable peaks
variances=NULL
medians=NULL
for (i in 1:dim(Y)[1]) {
  variances[i] = var(Y[i,])
  #medians[i] = median(Y[i,])
}
top = row.names(Y)[variances > quantile(variances,probs=0.999)]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# kmeans clustering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
kclusters = list()
kvals = seq(5,10, by=2)
for (i in 1:length(kvals)) {
  k = kvals[i]
  cat("k = ", k,"\n")
  kclusters[[i]] = kmeans(Y[top,], centers=k)
}
fit = sapply(kclusters, function(m) m[["betweenss"]]/m[["totss"]])
clust = data.frame(sapply(kclusters, function(m) paste0("cluster",m[["cluster"]])))
names(clust) = paste0("k_",kvals)
row.names(clust) <- names(kclusters[[1]]$cluster)

pdf("wgcna_atac_kmeans.pdf")
plot(kvals, fit)
pheatmap(Y[top,],
  cluster_rows=TRUE, cluster_cols=TRUE,
  show_rownames=FALSE,show_colnames=FALSE,
  annotation_row=clust[,c("k_5","k_7","k_9")],
  annotation_col=sampleDF[,c("donor","category","condition")],
  treeheight_row = 0, treeheight_col = 0)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Hierarchical clustering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tree = hclust(dist(Y[top,]))
plot(tree)
hclusters10 = cutree(tree, k = 10)[top]
hclusters25 = cutree(tree, k = 25)[top]

pheatmap(Y[top,],
  cluster_rows=TRUE, cluster_cols=TRUE,
  show_rownames=FALSE,show_colnames=FALSE,
  #annotation_row=data.frame(h10=as.factor(hclusters10),h25=as.factor(hclusters25)),
  annotation_row=data.frame(h10=as.factor(hclusters10)),
  annotation_col=sampleDF[,c("category","condition")],
  treeheight_row = 0, treeheight_col = 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sub-clustering
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rest = sampleDF$sample_id[sampleDF$condition==1]
stim = sampleDF$sample_id[sampleDF$condition==2]

Y_rest = Y[,rest]
Y_stim = Y[,stim]



### Define subgroups
categories = unique(sampleDF$category)
subcategories1 = unique(sampleDF$subcategory_cdtype[!is.na(sampleDF$subcategory_cdtype)])
subcategories2 = unique(sampleDF$subcategory_memory[!is.na(sampleDF$subcategory_memory)])
conditions = unique(sampleDF$condition)
count=0
hsubclusters = list()
for (i in 1:length(categories)) {
  samples_to_include = with(sampleDF, sampleDF[which(category==categories[i]),"sample_id"])
  if (length(samples_to_include)>0) {
    count = count + 1
    label = paste0(categories[i])
    cat(count,":",label,":",length(samples_to_include),"\n")
    hsubclusters[[label]] = list()
    hsubclusters[[label]]$samples = samples_to_include
  }
  for (j in 1:length(conditions)) {
    samples_to_include = with(sampleDF, sampleDF[which(category==categories[i] & condition==conditions[j]),"sample_id"])
    if (length(samples_to_include)>0) {
      count = count + 1
      label = paste0(categories[i],"-",conditions[j])
      cat(count,":",label,":",length(samples_to_include),"\n")
      hsubclusters[[label]] = list()
      hsubclusters[[label]]$samples = samples_to_include
    }
  }
}

for (i in 1:length(subcategories1)) {
  samples_to_include = with(sampleDF, sampleDF[which(subcategory_cdtype==subcategories1[i]),"sample_id"])
  if (length(samples_to_include)>0) {
    count = count + 1
    label = paste0(subcategories1[i])
    cat(count,":",label,":",length(samples_to_include),"\n")
    hsubclusters[[label]] = list()
    hsubclusters[[label]]$samples = samples_to_include
  }
  for (j in 1:length(conditions)) {
    samples_to_include = with(sampleDF, sampleDF[which(subcategory_cdtype==subcategories1[i] & condition==conditions[j]),"sample_id"])
    if (length(samples_to_include)>0) {
      count = count + 1
      label = paste0(subcategories1[i],"-",conditions[j])
      cat(count,":",label,":",length(samples_to_include),"\n")
      hsubclusters[[label]] = list()
      hsubclusters[[label]]$samples = samples_to_include
    }
  }
}

# for (i in 1:length(subcategories2)) {
#   samples_to_include = with(sampleDF, sampleDF[which(subcategory_memory==subcategories2[i]),"sample_id"])
#   if (length(samples_to_include)>0) {
#     count = count + 1
#     label = paste0(subcategories2[i])
#     cat(count,":",label,":",length(samples_to_include),"\n")
#     hsubclusters[[label]] = list()
#     hsubclusters[[label]]$samples = samples_to_include
#   }
#   for (j in 1:length(conditions)) {
#     samples_to_include = with(sampleDF, sampleDF[which(subcategory_memory==subcategories2[i] & condition==conditions[j]),"sample_id"])
#     if (length(samples_to_include)>0) {
#       count = count + 1
#       label = paste0(subcategories2[i],"-",conditions[j])
#       cat(count,":",label,":",length(samples_to_include),"\n")
#       hsubclusters[[label]] = list()
#       hsubclusters[[label]]$samples = samples_to_include
#     }
#   }
# }


### For each subgroup, calculate peak variances across samples
#hsubclusters = lapply(hsubclusters, "[[<-", "test", NULL)
for (i in 1:length(hsubclusters)) {
  cat(names(hsubclusters)[i],"\n")
  hsubclusters[[i]][["variances"]] = apply(X=Y[,hsubclusters[[i]][["samples"]]], MARGIN=1, FUN=var)
}
saveRDS(hsubclusters, file="hierarchical_sublustering_results.rds")


### For each subgroup, cluster and plot most variable peaks
getClustersAndPlots = function(label, ngenes=250) {
  cat(label,"\n")
  top_peaks = names(sort(-hsubclusters[[label]][["variances"]]))[1:ngenes]
  samples_to_include = hsubclusters[[label]][["samples"]]
  Y_sub = Y[top_peaks,samples_to_include]
  sampleDF_sub = sampleDF[samples_to_include,]

  tree = hclust(dist(Y_sub))
  hclusters5 = cutree(tree, k = 5)[top_peaks]
  hclusters10 = cutree(tree, k = 10)[top_peaks]
  plot(tree, main=label)

  print(pheatmap(Y_sub,
    cluster_rows=TRUE, cluster_cols=TRUE,
    show_rownames=FALSE,show_colnames=TRUE,
    annotation_row=data.frame(h5=as.factor(hclusters5), h10=as.factor(hclusters10), row.names=top_peaks),
    annotation_col=data.frame(type=sampleDF_sub$type, condition=sampleDF_sub$condition, row.names=sampleDF_sub$sample_id),
    treeheight_row = 0, treeheight_col = 0, main=label))
}

pdf("subclustering.pdf")
for (i in 1:length(hsubclusters)) {
  label=names(hsubclusters)[i]
  getClustersAndPlots(label)
}
dev.off()

pdf("subclustering_cd4.pdf")
#remove outlier
hsubclusters[["CD4-1"]]$samples = hsubclusters[["CD4-1"]]$samples[!hsubclusters[["CD4-1"]]$samples=="Th17_precursors_1003_1_hg38"]
#recalculate variances without outlier
hsubclusters[["CD4-1"]][["variances"]] = apply(X=Y[,hsubclusters[["CD4-1"]][["samples"]]], MARGIN=1, FUN=var)
#replot
getClustersAndPlots("CD4-1", ngenes=50)
dev.off()






saveRDS(hsubclusters, file="hierarchical_sublustering_results.rds")
