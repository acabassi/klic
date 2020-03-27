#' Kernel k-means
#'
#' Perform the training step of kernel k-means.
#' @param K Kernel matrix.
#' @param parameters A list containing the number of clusters
#' \code{number_count}.
#' @return This function returns a list containing:
#' \item{clustering}{the cluster labels for each element (i.e. row/column) of
#' the kernel matrix.}
#' \item{objective}{the value of the objective function for the given
#' clustering.}
#' \item{parameters}{same parameters as in the input.}
#' @author Mehmet Gonen
#' @references Gonen, M. and Margolin, A.A., 2014. Localized data fusion for
#' kernel k-means clustering with application to cancer biology. In Advances in
#' Neural Information Processing Systems (pp. 1305-1313).
#' @examples
#' \dontrun{clustering}
#' # Load one dataset with 100 observations, 2 variables, 4 clusters
#' data <- as.matrix(read.csv(system.file("extdata", "dataset1.csv",
#' package = "klic"), row.names = 1))
#' # Compute consensus clustering with K=4 clusters
#' cm <- coca::consensusCluster(data, 4)
#' # Shift eigenvalues of the matrix by a constant: (min eigenvalue) * (coeff)
#' km <- spectrumShift(cm, coeff = 1.05)
#' # Initalize the parameters of the algorithm
#' parameters <- list()
#' # Set the number of clusters
#' parameters$cluster_count <- 4
#' # Perform training
#' state <- kkmeans(km, parameters)
#' # Display the clustering
#' print(state$clustering)
#' }
#' @export
#'
kkmeans <- function(K, parameters) {
    state <- list()
    state$time <- system.time({
        H <- eigen(K, symmetric = TRUE)$vectors[, 1:parameters$cluster_count]
        objective <- sum(diag(t(H) %*% K %*% H)) - sum(diag(K))
        H_normalized <- H/matrix(sqrt(rowSums(H^2, 2)),
                                 nrow(H),
                                 parameters$cluster_count,
                                 byrow = FALSE)
        H_normalized[sqrt(rowSums(H^2, 2)) == 0, ] <- 0

        set.seed(NULL)
        state$clustering <- stats::kmeans(H_normalized,
                                          centers = parameters$cluster_count,
                                          iter.max = 1000,
                                          nstart = 10)$cluster
        state$objective <- objective
        state$parameters <- parameters
    })
    state
}
