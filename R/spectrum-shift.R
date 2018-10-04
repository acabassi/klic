#' Spectrum shift
#'
#' Make a symmetric matrix positive semi-definite
#' @param kernelMatrix symmetric matrix
#' @param coeff Coefficient by which the minimum eigenvalue is multiplied when shifting the eigenvalues, in order to avoid numeric problems. Default is 1.2.
#' @param verbose Boolean flag: if TRUE, information about the shift is printed to screen. Default is FALSE.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @examples
#' # Load one dataset with 300 observations, 2 variables, 6 clusters
#' data <- as.matrix(read.csv(system.file("extdata", "dataset1.csv", package = "klic"), row.names = 1))
#' # Compute consensus clustering with K=6 clusters
#' cm <- consensusCluster(data, 6)
#' # Shift eigenvalues of the matrix by a constant: (min eigenvalue) * (coeff)
#' km <- spectrumShift(cm, coeff = 1.05)
#' @export

spectrumShift = function(kernelMatrix, coeff = 1.2, verbose = FALSE){

  if(!isSymmetric(kernelMatrix))
    stop("The kernel matrix must be symmetric!")

  N <- dim(kernelMatrix)[1]

  # Get smallest eigenvalue
  min_eig <- eigen(kernelMatrix)$values[N]

  if(min_eig < 0){

    if(verbose) cat("The smallest eigenvalue is negative:", min_eig,"\n")

    kernelMatrix <- kernelMatrix + diag(dim(kernelMatrix)[1])*abs(min_eig*coeff)
    kernelMatrix <- kernelMatrix/kernelMatrix[1,1] # rescales

    if(verbose) cat("Shifting by a coefficient:  ", abs(min_eig*coeff), "\n")

    new_min_eig <- eigen(kernelMatrix)$values[N]

    if(verbose) cat("The smallest eigenvalue is now: ", new_min_eig, "\n")
  }

  kernelMatrix
}

