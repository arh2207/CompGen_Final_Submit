## Script Andrew used to curate pancreas data 
setwd('/Users/andrew/Library/CloudStorage/OneDrive-Personal/Columbia/2022-1 Spring/CBMFW4761/APCluster_CompGen')
example <- readRDS('tissue_counts/lung_counts.rds')

p1 <- read.csv('Source/pancreas/GSM2230757_human1_umifm_counts.csv')
p2 <- read.csv('Source/pancreas/GSM2230758_human2_umifm_counts.csv')
p3 <- read.csv('Source/pancreas/GSM2230759_human3_umifm_counts.csv')
p4 <- read.csv('Source/pancreas/GSM2230760_human4_umifm_counts.csv')

# stack data
stacked <- rbind(p1, p2, p3, p4)

# rename cells
for (i in 1:nrow(stacked)) {
  stacked[i,"X"] = paste("pancreas.", i, sep='')
}

# rename columns
names(stacked)[names(stacked) == 'X'] <- 'cell'
names(stacked)[names(stacked) == 'assigned_cluster'] <- 'type'

# create assignment vector
assignments <- stacked[,c(1,3)]
assignments <- as.matrix(assignments)
assignment.vec <- assignments[,2]
#NAMES and not ROWNAMES (referring to matrix issue earlier)
names(assignment.vec) <- assignments[,1]
saveRDS(assignment.vec, file = "tissue_analysis/pancreas/pancreas_assignments.rds")

# create data vector, like 'example' above
output <- as.matrix(stacked[,4:ncol(stacked)])
rownames(output) <- stacked[1:nrow(stacked),1]
output <- as.matrix(t(output))
saveRDS(output, file="tissue_counts/pancreas_counts.rds")