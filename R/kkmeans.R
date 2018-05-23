#' Kernel k-means
#'
#' Perform the training step of kernel k-means
#' @param K kernel matrix
#' @param parameters list of parameters
#' @author Mehmet Gonen
#' @references Gönen, M. and Margolin, A.A., 2014. Localized data fusion for kernel k-means clustering with application to cancer biology. In Advances in Neural Information Processing Systems (pp. 1305-1313).
#' @export
#'
kkmeansTrain <- function(K, parameters) {
  state <- list()
  state$time <- system.time({
    H <- eigen(K, symmetric = TRUE)$vectors[, 1:parameters$cluster_count]
    objective <- sum(diag(t(H) %*% K %*% H)) - sum(diag(K))
    H_normalized <- H / matrix(sqrt(rowSums(H^2, 2)), nrow(H), parameters$cluster_count, byrow = FALSE)

    set.seed(NULL)
    state$clustering <- kmeans(H_normalized, centers = parameters$cluster_count, iter.max = 1000, nstart = 10)$cluster
    state$objective <- objective
    state$parameters <- parameters
  })
  state
}
