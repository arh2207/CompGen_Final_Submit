# install github, etc.
###INITIAL PIPELINE
setwd('/Users/andrew/Library/CloudStorage/OneDrive-Personal/Columbia/2022-1 Spring/CBMFW4761/APCluster_CompGen')
#setwd('/Volumes/External NVMe (MacOS)/OneDrive/Columbia/2022-1 Spring/CBMFW4761/APCluster_CompGen')
#setwd('/Users/Andrew/Desktop/APCluster_CompGen')
#setwd('/Users/kellyjones/Desktop/CBMF4671_APCluster')
set.seed(1000)
library(stringr)
library(apcluster)
library(igraph)
library(PISCES)
library(Seurat)
library(Matrix.utils)
library(data.table)
library(tidyverse)
library(APTest)
source('Rscripts/Gen_Funcs.R')
tissue <- 'pancreas'
counts_matrix <- readRDS(paste('tissue_counts/', tissue, '_counts.rds', sep = ''))

#######
# Quality Control
#######

# 1.1) Baron et al. removed mitochondrial and ribosomal genes, and had some filtering of data
# Seurat object from pancreas counts
pancreas_seurat <- CreateSeuratObject(counts = counts_matrix, project = paste('panc_', tissue, sep = ''), min.cells = 3, min.features = 200)

# Quality control
# Redundant due to 1.1, but needed for QCPlots
temp <- GeneNameConvert(counts_matrix, species='human', from='gn', to='ensg')
mitochondrial_genes <- intersect(mitochondrial_genes$hum.ensg, rownames(temp))
pancreas_seurat[["percent.mt"]] <- PercentageFeatureSet(pancreas_seurat, features = mitochondrial_genes)
quality_control <- QCPlots(pancreas_seurat)
ggsave(file=paste('tissue_analysis/', tissue, '/', 'figures', '/', tissue, '_qc.jpg', sep = ''),
     width = 8, height = 4, units = "in", dpi = 2000)
dev.off()
print(quality_control)

pancreas_seurat <- subset(pancreas_seurat, subset = nCount_RNA < 10000)
# Normalizes and scales data, leaving multiple transformations, namely cpm and sct
pancreas_seurat <- SCTransform(pancreas_seurat, verbose = FALSE)
# Standard feature transformations
pancreas_seurat <- RunPCA(pancreas_seurat, verbose=FALSE, nPCs = 50)
pancreas_seurat <- RunUMAP(pancreas_seurat, dims = 1:50, verbose = FALSE)
saveRDS(pancreas_seurat, file=paste('local/', tissue, '/', tissue, '_post-QC.rds', sep=''))

#######
# Louvain Clustering
#######
# Generate knn graph
pancreas_seurat <- readRDS(file=paste('local/', tissue, '/', tissue, '_post-QC.rds', sep=''))
pancreas_seurat <- FindNeighbors(pancreas_seurat, dims = 1:30, verbose = FALSE)
sim_mat <- as.matrix(pancreas_seurat$SCT_snn)
saveRDS(sim_mat, file=paste('tissue_analysis/',tissue,'/sim_mat_pca.rds', sep=""))
# Clustering with Louvain 
pancreas_seurat <- FindClusters(pancreas_seurat, verbose = FALSE)
# Displays number of clusters and count of cells in each
table(Idents(pancreas_seurat))
louvain_cluster <- pancreas_seurat[['seurat_clusters']][,1]
names(louvain_cluster) <- rownames(pancreas_seurat[['seurat_clusters']])
louvain_cluster <- as.matrix(louvain_cluster)
# louvain.rds is named vector with assignments as entries
saveRDS(louvain_cluster, file = paste('tissue_analysis/', tissue,'/',tissue,'_louvain.rds', sep=''))

#######
# K-Means Clustering
#######
cell.embeddings <- pancreas_seurat[['pca']]@cell.embeddings
kmeans.obj <- kmeans(cell.embeddings, 14)
# kmeans.rds is a named vector with clustering assignments
saveRDS(kmeans.obj$cluster, file = paste('tissue_analysis/', tissue,'/',tissue,'_kmeans.rds', sep=''))

#######
# AP Clustering first-pass.  Calls APTest, a function built for another class,
# but the present analysis is unique.
#######
# Commented out due to runtime
#apresult.obj <- APTest(pancreas_seurat, PCA = TRUE, npcs = 30, spearman = FALSE, pearson = FALSE, negmanhattan = TRUE, negeuclidean = FALSE)
#This took a long time to generate, so we save the whole object
saveRDS(apresult.obj, file=paste('tissue_analysis/', tissue,'/',tissue,'_firstpass-apres_obj.rds', sep=''))
#the best cluster is clearly the non-pca cpm neg manhattan distance with minimum on the diagonal
apresult.best <- apresult.obj$`NON-PCA`$cpm.manhattan_Min
#function I wrote to vectorize into the named vector data structure
ap_cluster <- APResultToVec(apresult.best)
# ap.rds is a named vector with assignments
saveRDS(ap_cluster, file = paste('tissue_analysis/', tissue,'/',tissue,'_ap.rds', sep=''))


#######
# First-Pass Analysis Objects Loading
#######
assignments <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_assignments.rds', sep=''))
kmeans <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_kmeans.rds', sep=''))
louvain <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_louvain.rds', sep=''))
ap <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_ap.rds', sep=''))
#singleR analysis was removed from experiment.R.  Not useful in this analysis.
#singleR <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_singleR.rds', sep=''))
pancreas_seurat <- readRDS(paste('local/', tissue, '/', tissue, '_post-QC.rds', sep=''))

#######
# First-Pass Analysis
#######
#Compare NMI values
FindNMI(assignments, louvain) #0.75429
FindNMI(assignments, kmeans) #0.6898
FindNMI(assignments, ap) #0.3503
#FindNMI(assignments, singleR) #0.3139

#Saving plots of the different clusters using tSNE reduction
pancreas_seurat <- RunTSNE(pancreas_seurat)
plot <- plotClusters(pancreas_seurat, assignments, 'Assignments')
print(plot)
ggsave(file=paste('tissue_analysis', '/', tissue, '/', 'figures', '/', tissue, '_assignments_plot.png',sep=''))
plot <- plotClusters(pancreas_seurat, louvain, 'Louvain')
print(plot)
ggsave(file=paste('tissue_analysis', '/', tissue, '/', 'figures', '/', tissue,'_louvain_plot.png',sep=''))
plot <- plotClusters(pancreas_seurat, kmeans, 'K-Means')
print(plot)
ggsave(file=paste('tissue_analysis', '/', tissue, '/', 'figures', '/', tissue, '_kmeans_plot.png',sep=''))
plot <- plotClusters(pancreas_seurat, ap, 'Affinity Propagation')
print(plot)
ggsave(file=paste('tissue_analysis', '/', tissue, '/', 'figures', '/', tissue, '_ap_plot.png',sep=''))
#plot <- plotClusters(pancreas_seurat, singleR, 'Single R')
#print(plot)
#ggsave(file=paste('tissue_analysis', '/', tissue, '/', 'figures', '/', tissue, '_singleR_plot.png',sep=''))

#######
# Variance explained by PCs
#######
graphPCvariance(pancreas_seurat)
ggsave(filename = paste('tissue_analysis/', tissue, '/', 'figures', '/', tissue, '_pca-var.jpg', sep=''), dpi = 500, height = 7, width = 10)

#######
# First-pass parameter space exploration of AP clustering parameter space
#######
pancreas_seurat <- readRDS(paste('local/', tissue, '/', tissue, '_post-QC.rds', sep=''))
#####  LOAD BELOW  ##### APparamspace.bigresult <- APTest(pancreas_seurat, skipnonpca = TRUE, npcs = 11, iterate = TRUE)
#####  LOAD BELOW  ##### saveRDS(APparamspace.bigresult, paste('tissue_analysis/', tissue, '/', tissue, '_paramspace-11PCs_PCit-apres_obj.rds', sep=''))
#####  LOAD BELOW  ##### APparamspace2.bigresult <- APTest(pancreas_seurat, PCA = FALSE, iterate = FALSE) #basic, on no PCs
#####  LOAD BELOW  ##### saveRDS(APparamspace2.bigresult, paste('tissue_analysis/', tissue, '/', tissue, '_paramspace-nopca_apres-obj.rds', sep=''))

#######
# Extracting the best solution using correlation and negative distance metrics
#######
pancreas_seurat <- RunTSNE(pancreas_seurat)
APparamspace.bigresult <- readRDS(paste('tissue_analysis/', tissue, '/', tissue, '_paramspace-11PCs_PCit-apres_obj.rds', sep=''))
APparamspace2.bigresult <- readRDS(paste('tissue_analysis/', tissue, '/', tissue, '_paramspace-nopca_apres-obj.rds', sep=''))

# Due to a very unique quirk in the APTest analysis, I must remove a result from APparamspace.bigresult.  One-off.
# I found these by using SHUTUP=FALSE and seeing which nmis would not compute, and removing them.
APparamspace.bigresult[["11PCit_manhattan_1=step"]][[6]][["manhattan_7PCs_Max"]] <- NULL
APparamspace.bigresult[["11PCit_manhattan_1=step"]][[8]][["manhattan_9PCs_Max"]] <- NULL
APparamspace2.bigresult[["NON-PCA"]][["cpm.spearman_Max"]] <- NULL

apnmis <- scoreAPTest(APparamspace.bigresult, assignments)
apnmis <- append(apnmis, scoreAPTest(APparamspace2.bigresult, assignments))

# We have 107 nmi results... wow!  Time to sort and then plot the best one
apnmis <- sort(apnmis, decreasing = TRUE)

# The bestresult from the first-pass analysis; NMI = 0.642
print(apnmis[1])
# This must be extrapolated manually
best.ap <- APResultToVec(APparamspace.bigresult[["11-PCS"]][["pca.pearson_Min"]])
plot <- plotClusters(pancreas_seurat, best.ap, paste("AP Pearson (First-pass)", sep="")) 
#plot <- plotClusters(pancreas_seurat, best.ap, paste("Best AP Result (nmi = ", round(apnmis[1], digits=2), ")", "\n", names(apnmis)[1], sep=""))
print(plot)
ggsave(file=paste('tissue_analysis', '/', tissue, '/', 'figures', '/', tissue, '_best-APTest.png',sep=''))

# Maybe we look for a good result that also has a reasonable number of clusters.
# A manual search of the clusters in the *.bigresults shows this example with
best.ap2 <- APResultToVec(APparamspace2.bigresult[["NON-PCA"]][["cpm.manhattan_Min"]])
plot <- plotClusters(pancreas_seurat, best.ap2, paste("AP Result (nmi = ", round(apnmis["cpm.manhattan_Min"], digits=2), ")", "\n", "cpm.manhattan_Min", sep=""))
print(plot)
ggsave(file=paste('tissue_analysis', '/', tissue, '/', 'figures', '/', tissue, '_best-APTest2.png',sep=''))

#######
# Make KNN-similarity matrix
#######
pancreas_seurat <- readRDS(paste('local/', tissue, '/', tissue, '_post-QC.rds', sep=''))
assignments <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_assignments.rds', sep=''))
# runtime START 5/6 @ 2:04am
# very large with similarity matrices saved
jacresult1 <- APJaccardSimKRange(pancreas_seurat, transformation = 'scale', target.clust.vect = assignments)
saveRDS(jacresult, file=paste('local/', tissue, '/', tissue, '_jacresult1.rds', sep=''))
# remove similarity matrices
jacresult1$jacc.name <- NULL
saveRDS(jacresult1, file=paste('tissue_analysis/', tissue,'/',tissue,'_jacresult1-pruned.rds', sep=''))

jacresult2 <- APJaccardSimKRange(pancreas_seurat, transformation = 'log', target.clust.vect = assignments)
saveRDS(jacresult2, file=paste('tissue_analysis/', tissue,'/',tissue,'_jacresult2-pruned.rds', sep=''))

jacresult3 <- APJaccardSimKRange(pancreas_seurat, transformation = 'stretch', target.clust.vect = assignments)
saveRDS(jacresult3, file=paste('tissue_analysis/', tissue,'/',tissue,'_jacresult3-pruned.rds', sep=''))

jacresult4 <- APJaccardSimKRange(pancreas_seurat, transformation = 'tan', target.clust.vect = assignments)
saveRDS(jacresult4, file=paste('tissue_analysis/', tissue,'/',tissue,'_jacresult4-pruned.rds', sep=''))

jacresult5 <- APJaccardSimKRange(pancreas_seurat, transformation = 'tanh', target.clust.vect = assignments)
saveRDS(jacresult5, file=paste('tissue_analysis/', tissue,'/',tissue,'_jacresult5-pruned.rds', sep=''))

#######
# Analyze KNN-similarity AP Cluster results
#######
#jacresult1 <- readRDS(file=paste('tissue_analysis/', tissue,'/',tissue,'_jacresult1-pruned.rds', sep=''))
#jacresult2 <- readRDS(file=paste('tissue_analysis/', tissue,'/',tissue,'_jacresult2-pruned.rds', sep=''))
#jacresult3 <- readRDS(file=paste('tissue_analysis/', tissue,'/',tissue,'_jacresult3-pruned.rds', sep=''))
#jacresult4 <- readRDS(file=paste('tissue_analysis/', tissue,'/',tissue,'_jacresult4-pruned.rds', sep=''))
#jacresult5 <- readRDS(file=paste('tissue_analysis/', tissue,'/',tissue,'_jacresult5-pruned.rds', sep=''))

# We notice that tangent performed the best.  Now, I modify Gen_Funcs to see if we can
# Transform the data more substantially using tangent; and, to iterate past k=50
# As we have not hit a local maximum.  We make sure to use 'tan-sharp' on k=50 first
# To make sure that we get at least the same NMI score as the 'tan'
jacbest1 <- APJaccardSimKRange(pancreas_seurat, transformation = 'tan-sharp', target.clust.vect = assignments, kmin = 50, kstep = 5, kmax = 100)
saveRDS(jacbest1, file=paste('tissue_analysis/', tissue,'/',tissue,'_jacbest1-pruned.rds', sep=''))

jacbest2 <- APJaccardSimKRange(pancreas_seurat, transformation = 'tan-sharp', target.clust.vect = assignments, kmin = 100, kstep = 10, kmax = 200)
saveRDS(jacbest2, file=paste('tissue_analysis/', tissue,'/',tissue,'_jacbest2-pruned.rds', sep=''))

#NMI keeps increasing...
jacbest3 <- APJaccardSimKRange(pancreas_seurat, transformation = 'tan-sharp', target.clust.vect = assignments, kmin = 200, kstep = 10, kmax = 300)
saveRDS(jacbest3, file=paste('tissue_analysis/', tissue,'/',tissue,'_jacbest3-pruned.rds', sep=''))

# Optimal parametrization: tan-sharp with k=280 neighbors, 0.695 NMI

#jacbest1 <- readRDS(file=paste('tissue_analysis/', tissue,'/',tissue,'_jacbest1-pruned.rds', sep=''))
#jacbest2 <- readRDS(file=paste('tissue_analysis/', tissue,'/',tissue,'_jacbest2-pruned.rds', sep=''))
jacbest3 <- readRDS(file=paste('tissue_analysis/', tissue,'/',tissue,'_jacbest3-pruned.rds', sep=''))
jaccard <- jacbest3$solutions$`k=280`
saveRDS(jaccard, file = paste('tissue_analysis/', tissue,'/',tissue,'_jaccard.rds', sep=''))

#######
# Graph best jaccard similarity result
#######
pancreas_seurat <- readRDS(paste('local/', tissue, '/', tissue, '_post-QC.rds', sep=''))
pancreas_seurat <- RunTSNE(pancreas_seurat)
jaccard <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_jaccard.rds', sep=''))
plotClusters(pancreas_seurat, jaccard, 'AP Jaccard Tangent, K=280')
ggsave(file=paste('tissue_analysis', '/', tissue, '/', 'figures', '/', tissue, '_jaccard_plot.png',sep=''))
