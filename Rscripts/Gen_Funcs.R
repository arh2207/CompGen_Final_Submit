library(APTest)
library(Seurat)
library(apcluster)
library(Matrix.utils)
library(stringr)
library(PISCES)
library(data.table)
library(tidyverse)
library(igraph)

#' A function that computes the Normalized Mutual Information  (NMI) between two 
#' cluster vectors.  
#' 
#' @param clusters1 a (named) vector with elements being assignment labels
#' @param clusters2 a (named) vector with elements being assignment labels
#' @return A numerical value on the range [0,1] corresponding to the NMI
#' @export
FindNMI <- function(clusters1, clusters2){
  mat1 <- as.matrix(clusters1)
  mat2 <- as.matrix(clusters2)
  for(x in list(mat1, mat2)){
    x[,1] <- as.numeric(as.factor(x[,1]))
  }
  combined <- join.Matrix(mat1, mat2, by.x=rownames(mat1), by.y=rownames(mat2), type='on')
  nmi <- compare(combined[,1], combined[,2], method='nmi')
  return(nmi)
}

#' (Saving plots of the different clusters using tSNE reduction)
#' 
#' @param object ...
#' @param clusters ...
#' @param plotName ...
#' @return ....
#' @export
plotClusters <- function(object, clusters, plotName){
  clusters <- as.matrix(clusters)
  rownames(clusters) <- substr(rownames(clusters), start=10, stop=14)
  clusters <- clusters[order(as.integer(rownames(clusters))),]
  clusters <- as.numeric(as.factor(clusters))
  temp <- Idents(object)
  Idents(object) <- clusters
  plot <- DimPlot(object, reduction='tsne', group.by='ident', label=TRUE, label.size=3, repel=TRUE)
  plot <- plot+labs(title=plotName)
  plot <- plot+theme(aspect.ratio=20/30, legend.position='none')
  Idents(object) <-
  return(plot)
}

# Calculate NMI for all vectors in a bigresult (returned from APTest), returns
# a named vector with the nmi values and the result from the bigresult
scoreAPTest <- function(bigresult, assignments, SHUTUP = TRUE) {
  names <- c()
  nmis <- c()
  for (majorclass in 1:length(bigresult)) {
    majorclass.name <- names(bigresult)[majorclass]
    print(majorclass.name)
    for (minorclass in 1:length(bigresult[[majorclass]])) {
      minorclass.name <- names(bigresult[[majorclass]])[minorclass]
      if (is.null(minorclass.name)) {
        # minorclass.name is null for the iteration,we have the PC iteration result
        for (itclass in 1:length(bigresult[[majorclass]][[minorclass]])) {
          itclass.name <- names(bigresult[[majorclass]][[minorclass]])[itclass]
          itclass.vect <- APResultToVec(bigresult[[majorclass]][[minorclass]][[itclass]])
          itclass.nmi <- FindNMI(assignments, itclass.vect)
          if (!SHUTUP) {
            print(paste(itclass.name, "COMPUTED"))
          }
          names <- append(names, itclass.name)
          nmis <- append(nmis, itclass.nmi)
        }
      } else {
        # for all other results, we have the actual cluster vector
        minorclass.vect <- APResultToVec(bigresult[[majorclass]][[minorclass]])
        minorclass.nmi <- FindNMI(assignments, minorclass.vect)
        if (!SHUTUP) {
          print(paste(minorclass.name, "COMPUTED"))
        }
        names <- append(names, minorclass.name)
        nmis <- append(nmis, minorclass.nmi)
      }
    }
  }
  names(nmis) <- names
  return(nmis)
}

#Create Jaccard Similarity matrix with specified number of neighbors and transformation
FindJaccardMatrix <- function(object, k = 5, transformation = 'default'){
  object <- FindNeighbors(object, k = k, verbose=FALSE, reduction='pca', dims=1:30)
  mat <- as.matrix(object$SCT_snn)
  if (is.null(transformation)){
    return(mat)
  } else if(transformation == 'log'){
    mat[mat==0] <- -25
    mat[mat!=-25] <- log(mat[mat!=-25])
    return(mat)
  } else if (transformation == 'inverse') {
    # this function attempts to squish all similarities close to 0 to the same 
    # negative score (-9), while accentuating scores close to 1 (100)
    mat <- 1/(1.01-mat) - 10
  } else if (transformation == 'stretch') {
    # this function simply applies a linear stretch of all data presumed in [0,1]
    # e.g. values of 0 will be -500, values of 1 will be 500
    mat <- (mat-0.5) * 1000
  } else if (transformation == 'tan') {
    # This one may seem confusing.  I am trying to use the tangent function on
    # [0,1] to accentuate both very high and very low values (stretch toward
    # positive and negative infinity).  However, I need to apply a squeeze 
    # (multiplying mat) to make this occur in a range of [0,1] instead of pi,
    # then apply a shift (make this occur over [0,1] instead of [-.5,.5]), and
    # I multiply by 0.99999 because we can't actually have values tend toward positive
    # infinity or negative infinity, because there does exist 0 and 1, but we 
    # want them very high.
    mat <- tan(0.9999*pi*mat - 0.9999*pi*0.5)
  } else if (transformation == 'tanh') {
    # This example will apply the hyperbolic tangent function (sigmoid may
    # also have been used) to essentially flatten values both close to 0
    # and close to 1 (opposite ov tangent above), and scale values between.
    # Then, there is a linear stretch so that all values are on [-1000,1000]
    mat <- 1000 * tanh(5*mat-2.5)
  } else if (transformation == 'tan-sharp') {
    mat <- tan(0.999999*pi*mat - 0.999999*pi*0.5)
  } else {
    return(mat)
  }
}

#Create Jaccard Similarity matrices using the specified 'transformation', and then
#iterate over the number of neighbors k used.  APCluster based on that result.
# Wraps FindJaccardMatrix.  If a target
#cluster vector is provided, it will score each vector in the list and return that as well.
#Inspired by PISCES::LouvainKRange
# Also keeps a list of similarity matrices, so you can compare how they affected the result
APJaccardSimKRange <- function(data.obj, transformation = 'default', APq = 0, 
                             kmin = 5, kmax = 50, kstep = 5, target.clust.vect= NULL,
                             save.simmats = FALSE) {
  if (is.null(target.clust.vect)) {
    has.target <- FALSE
  } else {
    has.target <- TRUE
  }
  cluster.list <- list()
  if (save.simmats) {
    simmat.list <- list()
  }
  if (has.target) {
    nmi.list <- list()
  }
  k <- kmin
  # Iterate until kmax, stepwise (kstep = 5 by default)
  while (k <= kmax) {
    print(paste("Generating Jaccard similarity matrix with transformation: ", transformation, "; k = ", k, "...", sep = ''))
    # Generate cluster and score based on NMI
    sim.mat <- FindJaccardMatrix(data.obj, k, transformation)
    clust.vec <- APResultToVec(apcluster(sim.mat, q=APq))
    if (has.target) {
      nmi.score <- FindNMI(target.clust.vect, clust.vec)
    }  
    
    # Append results to output list
    listindex <- paste('k=', k, sep = '')
    cluster.list[[listindex]] <- clust.vec
    if (save.simmats) {
      simmat.list[[listindex]] <- sim.mat
    }
    if (has.target) {
      nmi.list[[listindex]] <- nmi.score
    }  
    
    k <- k + kstep
  }
  # Return list of clustering solutions
  jacc.name <- paste('jaccard sim; ', transformation, sep = '')
  if (has.target) {
    if (save.simmats) {
      return(list('solutions' = cluster.list, 'nmis' = nmi.list, jacc.name = simmat.list))
    } else {
      return(list('solutions' = cluster.list, 'nmis' = nmi.list))
    }
  } else {
    if (save.simmats) {
      return(list('solutions' = cluster.list, jacc.name = simmat.list))
    } else {
      return(list('solutions' = cluster.list))
    }
  }
}
