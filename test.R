### Sample code to run parts of our pipeline, test various functions
# You must set your working directory as the base directory!  Example below:
setwd('/Users/andrew/Library/CloudStorage/OneDrive-Personal/Columbia/2022-1 Spring/CBMFW4761/CompGen_Final_Submit')
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

# Loads the counts matrix, which was curated into a *.rds file using Rscripts/curatepancreas.R
counts_matrix <- readRDS(paste('tissue_counts/', tissue, '_counts.rds', sep = ''))

#######
# Quality Control example, uses PISCES package
#######

# 1.1) Baron et al. removed mitochondrial and ribosomal genes, and had some filtering of data
# Seurat object from pancreas counts
pancreas_seurat <- CreateSeuratObject(counts = counts_matrix, project = paste('panc_', tissue, sep = ''), min.cells = 3, min.features = 200)

# Quality control
# Redundant due to 1.1, but needed for QCPlots
temp <- GeneNameConvert(counts_matrix, species='human', from='gn', to='ensg')
mitochondrial_genes <- intersect(mt.genes$hum.ensg, rownames(temp))
pancreas_seurat[["percent.mt"]] <- PercentageFeatureSet(pancreas_seurat, features = mitochondrial_genes)
quality_control <- QCPlots(pancreas_seurat)
ggsave(file=paste('tissue_analysis/', tissue, '/', 'figures', '/', tissue, '_qc.jpg', sep = ''),
     width = 8, height = 4, units = "in", dpi = 2000)
dev.off()

# As we can see, mt features were filtered in the original dataset!  (MT% 0)
print(quality_control)

# Filter abnormally high read depth cells
pancreas_seurat <- subset(pancreas_seurat, subset = nCount_RNA < 10000)

# Normalizes and scales data, leaving multiple transformations, namely cpm and sct
pancreas_seurat <- SCTransform(pancreas_seurat, verbose = FALSE)

# Standard feature transformations
pancreas_seurat <- RunPCA(pancreas_seurat, verbose=FALSE, nPCs = 50)
pancreas_seurat <- RunUMAP(pancreas_seurat, dims = 1:50, verbose = FALSE)

#######
# Louvain Clustering example
#######
# Generate knn graph
pancreas_seurat <- FindNeighbors(pancreas_seurat, dims = 1:30, verbose = FALSE)
sim_mat <- as.matrix(pancreas_seurat$SCT_snn)
# Clustering with Louvain 
pancreas_seurat <- FindClusters(pancreas_seurat, verbose = FALSE)
# Displays number of clusters and count of cells in each
table(Idents(pancreas_seurat))
louvain_cluster <- pancreas_seurat[['seurat_clusters']][,1]
names(louvain_cluster) <- rownames(pancreas_seurat[['seurat_clusters']])
louvain_cluster <- as.matrix(louvain_cluster)
# louvain.rds is named vector with assignments as entries

#######
# Load some cluster solutions and score NMIs
#######
assignments <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_assignments.rds', sep=''))
kmeans <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_kmeans.rds', sep=''))
louvain <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_louvain.rds', sep=''))
ap <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_ap.rds', sep=''))

#Compare NMI values
FindNMI(assignments, louvain) #0.75429
FindNMI(assignments, kmeans) #0.6898
FindNMI(assignments, ap) #0.3503

# tSNE space plotting.  Must call Seurat RunTSNE first
pancreas_seurat <- RunTSNE(pancreas_seurat)

# The assignments
plot <- plotClusters(pancreas_seurat, assignments, 'Assignments')
print(plot)

# The Louvain solution
plot <- plotClusters(pancreas_seurat, louvain, 'Louvain')
print(plot)

# A first-pass AP solution
plot <- plotClusters(pancreas_seurat, ap, 'Affinity Propagation')
print(plot)

# One example of running APJaccardSimKRange, runtime is not feasible for test file.
# These results below are garbage results, and the ones we use in the real analysis are 
# read from .rds files below.  This is simply to showcase a function.  Inspect
# the object by accessing its elements or using the environment object inspector
# in an environment like RStudio
#jacresult <- APJaccardSimKRange(pancreas_seurat, transformation = 'default', target.clust.vect = assignments,
#                                kmin = 5, kmax = 15, kstep = 5)

# Examples of running APJaccardSimKRange deleted here, outputs are in tissue_analysis/pancreas/.
# See Rscripts/experiment.R to see the full exploration of the Jaccard index

# The best result from our paper is here, commented out, but you can read it the line below
#jacbest <- APJaccardSimKRange(pancreas_seurat, transformation = 'tan-sharp', target.clust.vect = assignments, kmin = 200, kstep = 10, kmax = 300)
jacbest <- readRDS(file=paste('tissue_analysis/', tissue,'/',tissue,'_jacbest3-pruned.rds', sep=''))

# Optimal parametrization: tan-sharp with k=280 neighbors, 0.695 NMI
jaccard <- jacbest$solutions$`k=280`

# And now we have the best AP solution!
# Here is the NMI: 0.694.  You can also see the NMI in an object viewer by opening
# 'jacbest' and seeing how the NMI changes based on k
FindNMI(assignments, jaccard)
plotClusters(pancreas_seurat, jaccard, 'AP Jaccard Tangent, K=280')
