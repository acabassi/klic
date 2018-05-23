#' Cluster-Of-Clusters Analysis
#'
#' This function allows to do Cluster-Of-Clusters-Analysis on a binary matrix where each column is a clustering of the data,
#' each row corresponds to a data point and the element in position (i,j) is equal to 1 if data point i belongs to cluster
#' j, 0 otherwise.
#' @param data N X C data matrix, where C is the total number of clusters considered.
#' @param K number of clusters.
#' @param B number of iterations of the Consensus Cluster step.
#' @param pItem proportion of items sampled at each iteration of the Consensus Cluster step.
#' @param hclustMethod method to be used by the hclust function for the hierarchical clustering step.
#' @return The output is a consensus matrix, that is a symmetric matrix where the element in position (i,j) corresponds to
#' the proportion of times that items i and j have been clustered together.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references Monti, S., Tamayo, P., Mesirov, J. and Golub, T., 2003. Consensus clustering: a resampling-based method for
#' class discovery and visualization of gene expression microarray data. Machine learning, 52(1-2), pp.91-118.
#' @references The Cancer Genome Atlas, 2012. Comprehensive molecular portraits of human breast tumours. Nature,
#' 487(7407), pp.61–70.
#' @export

coca = function(dataset, K, B = 1000, pItem = 0.8, hclustMethod = 'average'){

  # # Step 1. Build matrix of clusters
  # do.call(rbind,matrix.list)

  # Step 1. Compute the consensus matrix
  consensusMatrix <- consensusCluster(dataset, K, B, pItem)

  # Step 2. Use hierarchical clustering on the consensus matrix to find the clustering
  distances <- as.dist(1-consensusMatrix)
  hClustering <- hclust(distances, method = hclustMethod)
  clusterLabels <- cutree(hClustering, K)

}
