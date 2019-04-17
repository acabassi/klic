#' Kernel learning integrative clustering
#'
#' This function allows to perform Kernel Learning Integrative Clustering on `M` data sets relative
#' to the same observations. The similarities between the observations in each data set
#' are summarised into `M` different kernels, that are then fed into a kernel k-means clustering
#' algorithm. The output is a clustering of the observations that takes into account all the
#' available data types and a set of weights that sum up to one, indicating how much each data set
#' contributed to the kernel k-means clustering.
#'
#' @param data List of M datasets, each of size N X P_m, m = 1, ..., M.
#' @param M number of datasets.
#' @param individualK Vector containing the number of clusters in each dataset. Default is NULL.
#' If the number of clusters is not provided, then all the possible values between 2 and individualMaxK
#' are considered and the best value is chosen for each dataset by maximising the silhouette.
#' @param individualMaxK Maximum number of clusters considered for the individual data. Default is 6.
#' @param individualClAlgorithm Clustering algorithm used for clustering of each dataset individually
#' if is required to find the best number of clusters.
#' @param globalK Number of global clusters. Default is NULL. If the number of clusters is not provided,
#' then all the possible values between 2 and globalMaxK are considered and the best value is chosen by
#' maximising the silhouette.
#' @param globalMaxK Maximum number of clusters considered for the final clustering. Default is 6.
#' @param B Number of iterations for consensus clustering. Default is 1000.
#' @param C Maximum number of iterations for localised kernel k-means. Default is 100.
#' @param scale Boolean. If TRUE, each dataset is scaled such that each column has zero mean and
#' unitary variance.
#' @param savePlots Boolean. If TRUE, a plot of the silhouette is saved in the working folder.
#' Default is FALSE.
#' @param fileName If savePlots is TRUE, this is the name of the png file.
#' @param verbose Boolean. Default is TRUE.
#' @param annotations Data frame containing annotations for final plot.
#' @param ccClMethods The i-th element of this vector goes into the `clMethod`` argument of consensusCluster()
#' for the i-th dataset.
#' @param ccDistHCs The i-th element of this vector goes into the `dist` argument of consensusCluster()
#' for the i-th dataset.
#' @param widestGap Boolean. If TRUE, compute also widest gap index to choose best number
#' of clusters. Default is FALSE.
#' @param dunns Boolean. If TRUE, compute also Dunn's index to choose best number
#' of clusters. Default is FALSE.
#' @param dunn2s Boolean. If TRUE, compute also alternative Dunn's index to choose best number
#' of clusters. Default is FALSE.
#' @return The function returns `consensusMatrices`, which is an array containing one consensus matrix
#' per data set, `weights`, that is a vector containing the weights assigned by the kernel k-means
#' algorithm to each consensu matrix, `globalClusterLabels`, a vector containing the cluster labels
#' of the observations, according to kernel k-means clustering done on the kernel matrices. If the number
#' of clusters is not provided, the function also returns `bestK`, the best number of clusters between
#' 2 and `maxIndividualK` for each data set. Similarly, if the final number of clusters is not
#' provided, the best number of clusters for the final (global) clustering `globalK` is also returned.
#' This too is chosen so as to maximise the silhouette index.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @references Cabassi, A. and Kirk, P. D. W. (2019). Multiple kernel learning for
#' integrative consensus clustering of genomic datasets. arXiv preprint. arXiv:1904.07701.
#' @examples
#' # Load synthetic data
#' data1 <- as.matrix(read.csv(system.file("extdata",
#' "dataset1.csv", package = "klic"), row.names = 1))
#' data2 <- as.matrix(read.csv(system.file("extdata",
#' "dataset2.csv", package = "klic"), row.names = 1))
#' data3 <- as.matrix(read.csv(system.file("extdata",
#' "dataset3.csv", package = "klic"), row.names = 1))
#' data <- list(data1, data2, data3)
#'
#' # Perform clustering with KLIC assuming to know the
#' # number of clusters in each individual dataset and in
#' # the final clustering
#' klicOutput <- klic(data, 3, individualK = c(6,6,6),
#' globalK = 6, B = 50, C = 5)
#'
#' # Extract cluster labels
#' klic_labels <- klicOutput$globalClusterLabels
#'
#' cluster_labels <- as.matrix(read.csv(system.file("extdata",
#' "cluster_labels.csv", package = "klic"), row.names = 1))
#' # Compute ARI
#' ari <- mclust::adjustedRandIndex(klic_labels, cluster_labels)
#' @export
#'
klic = function(data, M, individualK = NULL, individualMaxK = 6,
                individualClAlgorithm = "kkmeans",
                globalK = NULL, globalMaxK = 6,
                B = 1000, C = 100, scale = FALSE, savePlots = FALSE, fileName = "klic",
                verbose = TRUE, annotations = NULL,
                ccClMethods = 'km', ccDistHCs = 'euclidean',
                widestGap = FALSE, dunns = FALSE, dunn2s = FALSE){

    ### Data check ###
    N = dim(data[[1]])[1]
    for(i in 1:M){
        if(dim(data[[i]])[1]!=N) stop("All datasets must have the same number of rows.")
    }
    if(verbose)
        print(paste("All datasets contain the same number of observations, ", N,".", sep = ""))
        print(paste("We assume that the observations are the same in each dataset and that they are in the same order."))

    # Check values of K
    if(individualMaxK == 2){
        individualK = 2
        warning('Since individualMaxK = 2, individualK is automatically set to 2.')
    }

    if(globalMaxK == 2){
        globalK = 2
        warning('Since globalMaxK = 2, globalK is automatically set to 2.')
    }

    # Initialise empty list for output
    output = list()

    ### Consensus clustering ###
    CM = array(NA, c(N, N, M))

    # If individual numbers of clusters are not known
    if(is.null(individualK)){

        if(verbose)
            print("*** Choosing the number of clusters for each dataset ***")

        # Initialise empty vector for best number of clusters in each dataset
        output$bestK <- rep(NA, M)
        # Initialise empty consensus matrices for all possible numbers of clusters
        tempCM = array(NA, c(N, N, individualMaxK-1))
        # Initialise empty cluster labels for all possible numbers of clusters
        clLabels = matrix(NA, individualMaxK-1, N)

        # For each dataset
        for(i in 1:M){

            if(verbose){
                print(paste("Dataset", i, sep = " "))
                pb = txtProgressBar(min = 0, max = individualMaxK-1, style = 3) # create progress bar
            }

            if(scale){
                # Scale the data such that each column has zero mean and unitary variance
                dataset_i = scale(data[[i]])
            }else{
                dataset_i = data[[i]]
            }

            ccClMethod_i = ccClMethods[i]
            ccDistHC_i = ccDistHCs[i]

            # For each possible number of clusters
            for(j in 2:individualMaxK){

                # Compute consensus matrix
                tempCM[,,j-1] <- coca::consensusCluster(dataset_i, j, B,
                                                        clMethod = ccClMethod_i, dist = ccDistHC_i)
                # Make consensus matrix positive definite
                tempCM[,,j-1] <- spectrumShift(tempCM[,,j-1])

                # If the chosen clustering algorithm is kernel k-means
                if(individualClAlgorithm == "kkmeans"){

                    # Initialise parameters for kernel k-means
                    parameters_kkmeans <- list()
                    # Set number of clusters for kernel k-means
                    parameters_kkmeans$cluster_count <- j
                    # Train kernel k-means
                    kkm <- kkmeans(tempCM[,,j-1], parameters_kkmeans)
                    # Extract cluster labels
                    clLabels[j-1,] <- kkm$clustering

                # If the chosen clustering algorithm is hierarchical clustering
                }else if(individualClAlgorithm == 'hclust'){

                    # Find clusters through hiearchical clustering
                    hCl <- stats::hclust(stats::as.dist(1-tempCM[,,j-1]), method = "average")
                    # Extract cluster labels
                    clLabels[j-1,] <-  stats::cutree(hCl, j)

                }else if(individualClAlgorithm == 'pam'){

                    # Find clusters through partitioning around medoids
                    clLabels[j-1,] <- cluster::pam(stats::as.dist(1-tempCM[,,j-1]))$clustering

                }else{ # If the value of individualClAlgorithm is not any of the above
                    stop("If there are factors among the covariates, individualClAlgorithm must be 'kkmeans', 'hclust', or 'pam'")
                }
                if(verbose) setTxtProgressBar(pb, j-1)
            }
            if(verbose) close(pb)

            # Find the number of clusters that maximises the silhouette
            maxSil <- coca::maximiseSilhouette(tempCM, clLabels, maxK = individualMaxK, savePNG = savePlots,
                                               fileName = paste("KLIC_dataset",i,sep=""),
                                               widestGap = widestGap,
                                               dunns = dunns, dunn2s = dunn2s)
            # If there is more than one, choose smallest number of clusters among the ones that maximise the silhouette
            bestK <- output$bestK[i] <- maxSil$K[1]

            if(verbose)
                print(paste("K =", bestK, sep = " "))

            # For dataset i, retain the consensus matrix corresponding to the smallest number of clusters for which the silhouette is maximised
            CM[,,i] <- tempCM[,,bestK-1]

            # Save plot of similarity matrix
            if(savePlots){

                if(!dir.exists("consensus-matrix")) dir.create("consensus-matrix", showWarnings = FALSE)
                fileName_i = paste('consensus-matrix/',fileName, "_CM", i, ".png", sep = "")

                clLab_i_bestK <- as.factor(clLabels[bestK-1,])
                names(clLab_i_bestK) <- as.character(1:N)
                CM_i_bestK <- as.matrix(CM[,,i])
                rownames(CM_i_bestK) <- colnames(CM_i_bestK) <- names(clLab_i_bestK)

                plotSimilarityMatrix2(CM_i_bestK, y = as.data.frame(clLab_i_bestK),
                                      fileName = fileName_i, save = TRUE)
            }
        }

    }else{ # If individual numbers of clusters are known

        if(verbose){
            print("*** Generating similarity matrices ***")
            pb = txtProgressBar(min = 0, max = M, style = 3) # create progress bar
        }

        # For each dataset
        for(i in 1:M){

            if(scale){

            # Scale the data such that each column has zero mean and unitary variance
            dataset_i = scale(data[[i]])

        }
            else{

            dataset_i = data[[i]]
        }

            # Compute consensus matrix
            CM[,,i] <- coca::consensusCluster(dataset_i, individualK[i], B)

            # Make consensus matrix positive definite
            CM[,,i] <- spectrumShift(CM[,,i])

            # Save plot of similarity matrix
            if(savePlots){
                if(!dir.exists("consensus-matrix")) dir.create("consensus-matrix", showWarnings = FALSE)
                fileName_i = paste('consensus-matrix/',fileName, "_CM", i, ".png", sep = "")

                CM_i <- CM[,,i]
                rownames(CM_i) <- colnames(CM_i) <- as.character(1:N)
                plotSimilarityMatrix2(CM_i, clr = TRUE, clc = TRUE,
                                      fileName = fileName_i, save = TRUE)
            }

            if(verbose)
                setTxtProgressBar(pb, i)
        }

        if(verbose) close(pb)
    }

    # Save individual consensus matrices for each dataset
    output$consensusMatrices <- CM

    ### Localised kernel k-means ###

    # Initialise parameters for localised multiple kernel k-means
    parameters <- list()
    # Set maximum number of iterations
    parameters$iteration_count <- C

    # If the global number of clusters is not known
    if(is.null(globalK)){

        if(verbose){
            print("*** Choosing the number of clusters for the global clustering ***")
            pb = txtProgressBar(min = 0, max = globalMaxK-1, style = 3) # create progress bar
        }

        # Initialise empty kernel matrix for each possible number of clusters
        KM <- array(0, c(N, N, globalMaxK-1))
        # Initialise empty cluster labels for each possible number of clusters
        clLabels <- array(NA, c(globalMaxK-1, N))
        # Initialise empty sets of weights for each possible number of clusters
        weights <- array(NA, c(N, M, globalMaxK - 1))

        # For every possible number of clusters
        for(i in 2:globalMaxK){
            # Set the number of clusters K
            parameters$cluster_count <- i
            # Train localised multiple kernel k-means
            lmkkm <- lmkkmeans(CM, parameters)
            # Save the weights
            weights[,,i-1] <- lmkkm$Theta
            # Save the combined kernel matrix
            for(j in 1:M){
                KM[,,i-1] <- KM[,,i-1] + (lmkkm$Theta[,j]%*%t(lmkkm$Theta[,j]))*CM[,,j]
            }
            # Save the cluster labels
            clLabels[i-1,] <- lmkkm$clustering

            if(verbose) setTxtProgressBar(pb, i-1)
        }
        if(verbose) close(pb)

        # Find the number of clusters that maximises the silhouette
        maxSil <- coca::maximiseSilhouette(KM, clLabels, maxK = globalMaxK,
                                           savePNG = TRUE, fileName = "KLIC_global",
                                           widestGap = widestGap,
                                           dunns = dunns, dunn2s = dunn2s)
        globalK <- output$globalK <- maxSil$K

        if(verbose)
            print(paste("Global K =", globalK, sep = " "))

        # Save chosen cluster labels
        output$globalClusterLabels <- clLabels[globalK-1,]
        # Save chosen weights
        output$weights <- weights[,,globalK-1]
        # Save chosen consensus matrix
        output$weightedKM <- KM[,,globalK-1]

    }else{ # If the global number of clusters is known

        if(verbose)
            print("*** Finding the global clustering ***")

        # Set number of clusters for localised multiple kernel k-means
        parameters$cluster_count <- globalK
        # Train localised multiple kernel k-means
        lmkkm <- lmkkmeans(CM, parameters)
        # Initialise empty weighted matrix
        weightedKM <- matrix(0, N, N)
        # Compute weighted matrix
        for(j in 1:M){
            weightedKM <- weightedKM + (lmkkm$Theta[,j]%*%t(lmkkm$Theta[,j]))*CM[,,j]
        }
        # Save cluster labels
        output$globalClusterLabels <-lmkkm$clustering
        # Save weights
        output$weights <- lmkkm$Theta
        # Save weighted kernel matrix
        output$weightedKM <- weightedKM
    }

    # Save plot of similarity matrix
    if(savePlots){

        if(!dir.exists("consensus-matrix")) dir.create("consensus-matrix", showWarnings = FALSE)
        fileName_i = paste('consensus-matrix/',fileName, "_weightedCM", i, ".png", sep = "")

        glClLab <- as.factor(output$globalClusterLabels)
        names(glClLab) <- as.character(1:N)
        weightedKM <- output$weightedKM
        rownames(weightedKM) <- colnames(weightedKM) <- names(glClLab)

        allAnnotations <- as.data.frame(glClLab)

        if(!is.null(annotations))
            allAnnotations <- cbind(allAnnotations, annotations)

        plotSimilarityMatrix2(weightedKM, y = allAnnotations, fileName = fileName_i,
                              save = TRUE)
    }

    return(output)
}
