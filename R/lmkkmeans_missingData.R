#' Localised multiple kernel k-means
#'
#' Perform the training step of the localised multiple kernel k-means
#' @param Km Array of size N X N X M containing M kernel matrices.
#' @param parameters List of parameters, containing `iteration_count`, `cluster_count`.
#' @param missing Matrix of size N X M containing missingness indicators. missing[i,j]=1 if
#' observation `i` is missing in dataset `j`
#' @param verbose Boolean flag. If TRUE, at each iteration the iteration number is printed.
#' Defaults to FALSE.
#' @author Mehmet Gonen, Alessandra Cabassi
#' @references Gonen, M. and Margolin, A.A., 2014. Localized data fusion for kernel k-means
#' clustering with application to cancer biology. In Advances in Neural Information Processing
#' Systems (pp. 1305-1313).
#' @examples
#' # Intialise 300 x 300 x 3 array containing M kernel matrices
#' # representing three different types of similarities
#' # between 300 data points
#' km <- array(NA, c(300, 300, 3))
#' # Load kernel matrices
#' km[,,1] <- as.matrix(read.csv(system.file("extdata",
#'                                          "kernel-matrix1.csv", package = "klic"), row.names = 1))
#' km[,,2] <- as.matrix(read.csv(system.file("extdata",
#'                                          "kernel-matrix2.csv", package = "klic"), row.names = 1))
#' km[,,3] <- as.matrix(read.csv(system.file("extdata",
#'                                          "kernel-matrix3.csv", package = "klic"), row.names = 1))
#' # Introduce some missing data
#' km[76:100, , 1] <- NA
#' km[, 76:100, 1] <- NA
#'
#' # Define missingness indicators
#' missing <- matrix(FALSE, 300, 3)
#' missing[76:100,1] <- TRUE
#'
#' #initalize the parameters of the algorithm
#' parameters <- list()
#' #set the number of clusters
#' parameters$cluster_count <- 2
#' # set the number of iterations
#' parameters$iteration_count <- 10
#' #perform training
#' state <- lmkkmeans_missingData(km, parameters, missing)
#' #display the clustering
#' print(state$clustering)
#' #display the kernel weights
#' print(state$Theta)
#' @export

lmkkmeans_missingData <- function(Km, parameters, missing = NULL, verbose = FALSE){

    # Assumiamo che `missing` sia una matrice N X P in cui missing[i,j]=1 se il dato `i`
    # manca nel dataset `j`

    state <- list()

    N <- dim(Km)[2]
    P <- dim(Km)[3]

    if(!is.null(missing)){
        avail <- abs(1 - missing)
    }else{
        avail <- matrix(1, N, P)
    }

    # Initialise weight matrix assigning equal weights to each object in each kernel
    Theta <- matrix(NA, N, P)
    for(i in 1:N){
        Theta[i,] <- 1/sum(avail[i,])
    }
    Theta <- Theta*avail # Set to zero the weights of missing observations

    # Initialise weighted kernel matrix
    K_Theta <- matrix(0, nrow(Km), ncol(Km))
    for (m in 1:P) {
        avail_m <- avail[,m]>0
        K_Theta[avail_m, avail_m] <- K_Theta[avail_m, avail_m] +
            (Theta[avail_m, m, drop = FALSE] %*% t(Theta[avail_m, m, drop = FALSE] )) *
            Km[avail_m, avail_m, m]
    }

    # Initialise vector of objective functions
    objective <- rep(0, parameters$iteration_count)

    for (iter in 1:parameters$iteration_count) {

        if(verbose) print(sprintf("running iteration %d...", iter))
        H <- eigen(K_Theta, symmetric = TRUE)$vectors[, 1:parameters$cluster_count]
        HHT <- H %*% t(H)

        Q <- matrix(0, N * P, N * P)
        for (m in 1:P) {
            start_index <- (m - 1) * N + 1
            end_index <- m * N
            Q[start_index:end_index, start_index:end_index] <- diag(1, N, N) * Km[,,m] -
                HHT * Km[,,m]
            # See Gonen & Margolin 2014 NIPS, page 5
        }

        avail_vec <- as.logical(as.vector(avail))
        sum_avail_vec <- sum(avail_vec)
        Q <- Q[avail_vec, avail_vec]

        ### Solve QP problem ###
        problem <- list()
        # problem$sense: Objective sense: e.g. "min" or "max"
        problem$sense <- "min"
        # problem$c: Objective coefficient array
        problem$c <- rep(0, sum_avail_vec)
        # problem$A: Constraint sparse matrix
        A <- Matrix::Matrix(rep(diag(1, N, N), P), nrow = N, ncol = N * P, sparse = TRUE)
        problem$A <- A[,avail_vec,drop=F]
        # problem$bc: Lower and upper constraint bounds
        problem$bc <- rbind(blc = rep(1, N), buc = rep(1, N))
        # problem$bx: Lower and upper variable bounds
        problem$bx <- rbind(blx = rep(0, sum_avail_vec), bux = rep(1, sum_avail_vec))
        # problem$qobj: Quadratic convex optimization
        I <- matrix(1:sum_avail_vec, sum_avail_vec, sum_avail_vec, byrow = FALSE)
        J <- matrix(1:sum_avail_vec, sum_avail_vec, sum_avail_vec, byrow = TRUE)
        problem$qobj <- list(i = I[lower.tri(I, diag = TRUE)],
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
        for(i in 1:P){
            avail_i <- which(avail[,i]==1)
            startt <- (count+1)
            endt <- count+sum(avail[,i])
            Theta[avail_i,i] <- result$sol$itr$xx[startt:endt]
            count <- count + sum(avail[,i])
        }

        # Update weighted kernel
        K_Theta <- matrix(0, nrow(Km), ncol(Km))
        for (m in 1:P) {
            K_Theta <- K_Theta + (Theta[,m,drop = FALSE] %*% t(Theta[,m,drop = FALSE])) * Km[,,m]
        }

        # Update objective function
        objective[iter] <- sum(diag(t(H) %*% K_Theta %*% H)) - sum(diag(K_Theta))
    }

    normalize <- which(rowSums(H^2, 2)>.Machine$double.eps)
    H_normalized <- matrix(0, N, parameters$cluster_count)
    H_normalized[normalize,] <- H[normalize,] / matrix(sqrt(rowSums(H[normalize,]^2, 2)), nrow(H[normalize,]), parameters$cluster_count, byrow = FALSE)

    set.seed(NULL)
    state$clustering <- stats::kmeans(H_normalized, centers = parameters$cluster_count, iter.max = 1000, nstart = 10)$cluster
    state$objective <- objective
    state$parameters <- parameters
    state$Theta <- Theta

    state
}

