#setwd('.')
setwd('/Users/andrew/Library/CloudStorage/OneDrive-Personal/Columbia/2022-1 Spring/CBMFW4761/APCluster_CompGen')
set.seed(1000)
library(data.table)
library(stringr)
library(Seurat)

tissue <- 'pancreas'

assignments <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_assignments.rds', sep=''))
louvain <-readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_louvain.rds', sep=''))
singler <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_singleR.rds', sep=''))
kmeans <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_kmeans.rds', sep=''))
ap <- readRDS(file = paste('tissue_analysis/', tissue,'/',tissue,'_ap.rds', sep=''))

kmeans <- as.vector(kmeans)
assignments <- as.vector(assignments)
louvain <- as.vector(t(louvain))

louvain.vs.assignments.table <- table(assignments, louvain)
colnames(louvain.vs.assignments.table) <- paste('Louvain', colnames(louvain.vs.assignments.table), sep = '.')
rownames(louvain.vs.assignments.table) <- paste(rownames(louvain.vs.assignments.table), sep = '.')
write.table(louvain.vs.assignments.table, sep = ',', quote = FALSE, col.names=NA,
            file = paste('tissue_analysis/', tissue, '/', tissue, '_louvain-assignments-table.csv', sep = ''))

kmeans.vs.assignments.table <- table(assignments, kmeans)
colnames(kmeans.vs.assignments.table) <- paste('kmeans', colnames(kmeans.vs.assignments.table), sep = '.')
rownames(kmeans.vs.assignments.table) <- paste(rownames(kmeans.vs.assignments.table), sep = '.')
write.table(kmeans.vs.assignments.table, sep = ',', quote = FALSE, col.names=NA,
            file = paste('tissue_analysis/', tissue, '/', tissue, '_kmeans-assignments-table.csv', sep = ''))

ap.vs.assignments.table <- table(assignments, ap)
colnames(ap.vs.assignments.table) <- paste('ap', colnames(ap.vs.assignments.table), sep = '.')
rownames(ap.vs.assignments.table) <- paste(rownames(ap.vs.assignments.table), sep = '.')
write.table(ap.vs.assignments.table, sep = ',', quote = FALSE, col.names=NA,
            file = paste('tissue_analysis/', tissue, '/', tissue, '_ap-assignments-table.csv', sep = ''))

