#' Localised multiple kernel k-means
#'
#' Perform the training step of the localised multiple kernel k-means.
#' @param Km Array of size N X N X M containing M  different N x N kernel
#' matrices.
#' @param parameters A list of parameters containing the desired number of
#' clusters, \code{cluster_count}, and the number of iterations of the
#' algorithm to be run, \code{iteration_count}.
#' @param missing Matrix of size N X M containing missingness indicators, i.e.
#' missing[i,j] = 1 (or = TRUE) if observation \code{i} is missing in dataset
#' \code{j}, missing[i,j] = 0 (or = FALSE).
#' @param verbose Boolean flag. If TRUE, at each iteration the iteration number
#' is printed. Defaults to FALSE.
#' @author Mehmet Gönen, Alessandra Cabassi
#' @references Gönen, M. and Margolin, A.A., 2014. Localized data fusion for
#' kernel k-means clustering with application to cancer biology. In Advances in
#' Neural Information Processing Systems (pp. 1305-1313).
#' @examples
#' \dontrun{
#' # Intialise 100 x 100 x 3 array containing M kernel matrices
#' # representing three different types of similarities between 100 data points
#' km <- array(NA, c(100, 100, 3))
#' # Load kernel matrices
#' km[,,1] <- as.matrix(read.csv(system.file('extdata',
#' 'kernel_matrix1.csv', package = 'klic'), row.names = 1))
#' km[,,2] <- as.matrix(read.csv(system.file('extdata',
#' 'kernel_matrix2.csv', package = 'klic'), row.names = 1))
#' km[,,3] <- as.matrix(read.csv(system.file('extdata',
#' 'kernel_matrix3.csv', package = 'klic'), row.names = 1))
#' # Introduce some missing data
#' km[76:80, , 1] <- NA
#' km[, 76:80, 1] <- NA
#'
#' # Define missingness indicators
#' missing <- matrix(FALSE, 100, 3)
#' missing[76:80,1] <- TRUE
#'
#' # Initalize the parameters of the algorithm
#' parameters <- list()
#' # Set the number of clusters
#' parameters$cluster_count <- 4
#' # Set the number of iterations
#' parameters$iteration_count <- 10
#'
#' # Perform training
#' state <- lmkkmeans_missingData(km, parameters, missing)
#'
#' # Display the clustering
#' print(state$clustering)
#' # Display the kernel weights
#' print(state$Theta)
#' }
#' @export

lmkkmeans_missingData <-
    function(Km,
             parameters,
             missing = NULL,
             verbose = FALSE) {

    state <- list()

    N <- dim(Km)[2]
    P <- dim(Km)[3]

    if (!is.null(missing)) {
        avail <- abs(1 - missing)
    } else {
        avail <- matrix(1, N, P)
    }

    # Initialise weight matrix assigning equal weights to each object in each
    # kernel
    Theta <- matrix(NA, N, P)
    for (i in 1:N) {
        Theta[i, ] <- 1/sum(avail[i, ])
    }
    Theta <- Theta * avail  # Set to zero the weights of missing observations

    # Initialise weighted kernel matrix
    K_Theta <- matrix(0, nrow(Km), ncol(Km))
    for (m in 1:P) {
        avail_m <- avail[, m] > 0
        K_Theta[avail_m, avail_m] <- K_Theta[avail_m, avail_m] +
            (Theta[avail_m, m, drop = FALSE] %*%
                 t(Theta[avail_m, m, drop = FALSE])) * Km[avail_m, avail_m, m]
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

            avail_m <- avail[, m] > 0
            start_index <- (m - 1) * N + 1
            end_index <- m * N

            Q[(start_index:end_index)[avail_m],
              (start_index:end_index)[avail_m]] <-
                diag(1, sum(avail_m), sum(avail_m)) *
                Km[avail_m, avail_m, m] - HHT[avail_m, avail_m] *
                Km[avail_m, avail_m, m]
            # See Gönen & Margolin 2014 NIPS, page 5
        }

        avail_vec <- as.logical(as.vector(avail))
        sum_avail_vec <- sum(avail_vec)
        Q <- Q[avail_vec, avail_vec]

        ### Solve QP problem ###
        problem <- list()
        # problem$sense: Objective sense: e.g. 'min' or 'max'
        problem$sense <- "min"
        # problem$c: Objective coefficient array
        problem$c <- rep(0, sum_avail_vec)
        # problem$A: Constraint sparse matrix
        A <- Matrix::Matrix(rep(diag(1, N, N), P),
                            nrow = N,
                            ncol = N * P,
                            sparse = TRUE)
        problem$A <- A[, avail_vec, drop = F]
        # problem$bc: Lower and upper constraint bounds
        problem$bc <- rbind(blc = rep(1, N), buc = rep(1, N))
        # problem$bx: Lower and upper variable bounds
        problem$bx <- rbind(blx = rep(0, sum_avail_vec),
                            bux = rep(1, sum_avail_vec))
        # problem$qobj: Quadratic convex optimization
        I <-
            matrix(1:sum_avail_vec, sum_avail_vec, sum_avail_vec, byrow = FALSE)
        J <-
            matrix(1:sum_avail_vec, sum_avail_vec, sum_avail_vec, byrow = TRUE)
        problem$qobj <-
            list(i = I[lower.tri(I, diag = TRUE)],
                 j = J[lower.tri(J, diag = TRUE)],
                 v = Q[lower.tri(Q, diag = TRUE)])
        opts <- list()
        # opts$verbose: Output logging verbosity
        opts$verbose <- 0

        # Solve QP problem
        result <- Rmosek::mosek(problem, opts)

        # Extract Theta and put it in matrix form
        Theta <- matrix(0, N, P)
        count <- 0
        for (i in 1:P) {
            avail_i <- which(avail[, i] == 1)
            startt <- (count + 1)
            endt <- count + sum(avail[, i])
            Theta[avail_i, i] <- result$sol$itr$xx[startt:endt]
            count <- count + sum(avail[, i])
        }

        # Update weighted kernel
        K_Theta <- matrix(0, nrow(Km), ncol(Km))
        for (m in 1:P) {
            avail_m <- avail[, m] > 0
            K_Theta[avail_m, avail_m] <-
                K_Theta[avail_m, avail_m] +
                (Theta[avail_m, m, drop = FALSE] %*%
                     t(Theta[avail_m, m, drop = FALSE])) *
                Km[avail_m, avail_m, m]
        }


        # Update objective function
        objective[iter] <-
            sum(diag(t(H) %*% K_Theta %*% H)) - sum(diag(K_Theta))
    }

    normalize <- which(rowSums(H^2, 2) > .Machine$double.eps)
    H_normalized <- matrix(0, N, parameters$cluster_count)
    H_normalized[normalize, ] <- H[normalize, ]/matrix(sqrt(rowSums(H[normalize,
        ]^2, 2)), nrow(H[normalize, ]), parameters$cluster_count, byrow = FALSE)

    state$clustering <-
        stats::kmeans(
            H_normalized,
            centers = parameters$cluster_count,
            iter.max = 1000,
            nstart = 10
        )$cluster
    state$objective <- objective
    state$parameters <- parameters
    state$Theta <- Theta

    state
}
