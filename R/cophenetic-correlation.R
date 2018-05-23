#' Cophenetic correlation coefficient
#'
#' Compute the cophenetic correlation coefficient of a kernel matrix.
#' @param kernelMatrix kernel matrix.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references Sokal, R.R. and Rohlf, F.J., 1962. The comparison of dendrograms by objective methods. Taxon, 11(2), pp.33-40.
#' @export

copheneticCorrelation = function(kernelMatrix){
  
  KMasDist   <- as.dist(1 - kernelMatrix)
  hclustered  <- hclust(KMasDist)
  cophDist <- cophenetic(hclustered)
  cophCor  <- cor(KMasDist, cophDist)

}

