#' Localised multiple kernel k-means
#'
#' Perform the training step of the localised multiple kernel k-means.
#' @param Km An array of size N x N x M containing M different N x N kernel
#' matrices.
#' @param parameters A list of parameters containing the desired number of
#' clusters, \code{cluster_count}, and the number of iterations of the
#' algorithm to be run, \code{iteration_count}.
#' @param verbose Boolean flag. If TRUE, at each iteration the iteration number
#' is printed. Default is FALSE.
#' @return This function returns a list containing:
#' \item{clustering}{the cluster labels for each element (i.e. row/column) of
#' the kernel matrix.}
#' \item{objective}{the value of the objective function for the given
#' clustering.}
#' \item{parameters}{same parameters as in the input.}
#' \item{Theta}{N x M matrix of weights, each row corresponds to an observation
#' and each column to one of the kernels.}
#' @author Mehmet Gonen
#' @references Gonen, M. and Margolin, A.A., 2014. Localized data fusion for
#' kernel k-means clustering with application to cancer biology. In Advances in
#' Neural Information Processing Systems (pp. 1305-1313).
#' @examples
#' if(requireNamespace("Rmosek", quietly = TRUE) &&
#' (!is.null(utils::packageDescription("Rmosek")$Configured.MSK_VERSION))){
#'
#' # Initialise 100 x 100 x 3 array containing M kernel matrices
#' # representing three different types of similarities between 100 data points
#' km <- array(NA, c(100, 100, 3))
#' # Load kernel matrices
#' km[,,1] <- as.matrix(read.csv(system.file('extdata',
#' 'kernel_matrix1.csv', package = 'klic'), row.names = 1))
#' km[,,2] <- as.matrix(read.csv(system.file('extdata',
#' 'kernel_matrix2.csv', package = 'klic'), row.names = 1))
#' km[,,3] <- as.matrix(read.csv(system.file('extdata',
#' 'kernel_matrix3.csv', package = 'klic'), row.names = 1))
#'
#' # Initalize the parameters of the algorithm
#' parameters <- list()
#' # Set the number of clusters
#' parameters$cluster_count <- 4
#' # Set the number of iterations
#' parameters$iteration_count <- 10
#'
#' # Perform training
#' state <- lmkkmeans(km, parameters)
#'
#' # Display the clustering
#' print(state$clustering)
#' # Display the kernel weights
#' print(state$Theta)
#' }
#' @export

lmkkmeans <- function(Km, parameters, verbose = FALSE) {
    state <- list()
    state$time <- system.time({
        N <- dim(Km)[2]
        P <- dim(Km)[3]

        # Initialise weight matrix assigning equal weights to each object in
        # each kernel
        Theta <- matrix(1/P, N, P)

        # Initialise weighted kernel matrix
        K_Theta <- matrix(0, nrow(Km), ncol(Km))
        for (m in 1:P) {
            K_Theta <- K_Theta + (Theta[, m, drop = FALSE] %*%
                                      t(Theta[, m, drop = FALSE])) * Km[, , m]
        }

        # Initialise vector of objective functions
        objective <- rep(0, parameters$iteration_count)

        for (iter in 1:parameters$iteration_count) {

            if (verbose)
                print(sprintf("running iteration %d...", iter))
            H <- eigen(K_Theta,
                       symmetric = TRUE)$vectors[, 1:parameters$cluster_count]
            HHT <- H %*% t(H)

            Q <- matrix(0, N * P, N * P)
            for (m in 1:P) {
                start_index <- (m - 1) * N + 1
                end_index <- m * N
                Q[start_index:end_index, start_index:end_index] <-
                    diag(1, N, N) * Km[, , m] - HHT * Km[, , m]
            }

            ### Solve QP problem ###
            problem <- list()
            # problem$sense: Objective sense: e.g. 'min' or 'max'
            problem$sense <- "min"
            # problem$c: Objective coefficient array
            problem$c <- rep(0, N * P)
            # problem$A: Constraint sparse matrix
            problem$A <- Matrix::Matrix(rep(diag(1, N, N), P),
                                        nrow = N,
                                        ncol = N * P, sparse = TRUE)
            # problem$bc: Lower and upper constraint bounds
            problem$bc <- rbind(blc = rep(1, N), buc = rep(1, N))
            # problem$bx: Lower and upper variable bounds
            problem$bx <- rbind(blx = rep(0, N * P), bux = rep(1, N * P))
            # problem$qobj: Quadratic convex optimization
            I <- matrix(1:(N * P), N * P, N * P, byrow = FALSE)
            J <- matrix(1:(N * P), N * P, N * P, byrow = TRUE)
            problem$qobj <- list(i = I[lower.tri(I, diag = TRUE)],
                                 j = J[lower.tri(J, diag = TRUE)],
                                 v = Q[lower.tri(Q, diag = TRUE)])
            opts <- list()
            # opts$verbose: Output logging verbosity
            opts$verbose <- 0

            # Solve QP problem
            result <- Rmosek::mosek(problem, opts)

            # Extract Theta and put it in matrix form
            Theta <- matrix(result$sol$itr$xx, N, P, byrow = FALSE)

            # Update weighted kernel
            K_Theta <- matrix(0, nrow(Km), ncol(Km))
            for (m in 1:P) {
                K_Theta <-
                    K_Theta + (Theta[, m, drop = FALSE] %*%
                                   t(Theta[, m, drop = FALSE])) * Km[, , m]
            }

            # Update objective function
            objective[iter] <- sum(
                diag(t(H) %*% K_Theta %*% H)) - sum(diag(K_Theta))

            diff <- objective[iter]-objective[iter-1]
            if(iter > 2 && abs(diff) < abs(1e-04*objective[iter-1])) break
        }
        H_normalized <- H/matrix(sqrt(rowSums(H^2, 2)),
                                 nrow(H),
                                 parameters$cluster_count,
                                 byrow = FALSE)

        state$clustering <- stats::kmeans(H_normalized,
                                          centers = parameters$cluster_count,
                                          iter.max = 1000, nstart = 10)$cluster
        state$objective <- objective[1:iter]
        state$parameters <- parameters
        state$Theta <- Theta
    })
    state
}
