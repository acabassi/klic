#' Cophenetic correlation coefficient
#'
#' Compute the cophenetic correlation coefficient of a kernel matrix, which is
#' a measure of howfaithfully hierarchical clustering would preserve the
#' pairwise distances between the original data points.
#' @param kernelMatrix kernel matrix.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references Sokal, R.R. and Rohlf, F.J., 1962. The comparison of dendrograms by objective methods. Taxon, 11(2), pp.33-40.
#' @examples
#' # Load kernel matrix
#' km <- as.matrix(read.csv(system.file("extdata",
#' "kernel-matrix1.csv", package = "klic"), row.names = 1))
#' # Compute cophenetic correlation
#' ccc <- copheneticCorrelation(km)
#' print(ccc)
#' @export

copheneticCorrelation = function(kernelMatrix){

  # Convert the similarities into distances
  KMasDist   <- stats::as.dist(1 - kernelMatrix)

  # Perform hierarchical clustering on the distances derived from the kernel matrix
  hclustered  <- stats::hclust(KMasDist)

  # Compute the cophenetic distances of the hierarchical clustering
  cophDist <- stats::cophenetic(hclustered)

  # Compute the correlation between the cophenetic distances and the actual distances
  cophCor  <- stats::cor(KMasDist, cophDist)

}

