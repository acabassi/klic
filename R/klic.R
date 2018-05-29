#' Kernel learning integrative clustering
#'
#'Perform kernel learning integrative clustering
#' @param data List of M datasets, each of size N X P_m, m = 1, ..., M.
#' @param M number of data.
#' @param individualK Vector containing the number of clusters in each dataset. Default is NULL. If the number of clusters is not provided, then all the possible values between 2 and individualMaxK are considered and the best value is chosen for each dataset by maximising the silhouette.
#' @param individualMaxK Maximum number of clusters considered for the individual data. Default is 6.
#' @param globalK Number of global clusters. Default is NULL. If the number of clusters is not provided, then all the possible values between 2 and globalMaxK are considered and the best value is chosen by maximising the silhouette.
#' @param globalMaxK Maximum number of clusters considered for the final clustering. Default is 6.
#' @param B Number of iterations for consensus clustering. Default is 1000.
#' @param C Maximum number of iterations for localised kernel k-means. Default is 100.
#' @param savePlots If TRUE, a plot of the silhouette is saved in the working folder. Default is FALSE.
#' @param fileName If savePlots is TRUE, this is the name of the png file.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @export
klic = function(data, M, individualK = NULL, individualMaxK = 6, globalK = NULL, globalMaxK = 6,
                B = 1000, C = 100, savePlots = FALSE, fileName = "klic"){

  ## Data check
  N = dim(data[[1]])[1]
  for(i in 1:M){
    if(dim(data[[i]])[1]!=N) stop("All datasets must have the same number of rows.")
  }

  ## Initialise empty list for output
  output = list()

  ## Consensus clustering
  CM = array(NA, c(N, N, M))
  if(is.null(individualK)){
      output$bestK <- rep(NA, M)
      tempCM = array(NA, c(N, N, individualMaxK-1))
      clLabels = matrix(NA, individualMaxK-1, N)
      for(i in 1:M){
        scaledDataset = scale(data[[i]])
        for(j in 2:individualMaxK){
          tempCM[,,j-1] <- consensusCluster(scaledDataset, j, B)
          tempCM[,,j-1] <- spectrumShift(tempCM[,,j-1])
          hCl <- hclust(as.dist(1-tempCM[,,j-1]), method = "average")
          clLabels[j-1,] <-  cutree(hCl, j)
        }
        maxSil <- maximiseSilhouette(tempCM, clLabels, individualMaxK)
        bestK <- output$bestK[i] <- maxSil$k[1] # choose smallest number of clusters
        CM[,,i] <- tempCM[,,bestK]
        if(savePlots){
          plotSimilarityMatrix(CM[,,i],
                               file_name = cat(fileName, "_CM", i, "_K", bestK, ".png", sep = ""),
                               savePNG = TRUE)
        }
    }
  }else{
    for(i in 1:M){
        scaledDataset = scale(data[[i]])
        CM[,,i] <- consensusCluster(scaledDataset, individualK[i], B)
        CM[,,i] <- spectrumShift(CM[,,i])
        if(savePlots){
          plotSimilarityMatrix(CM[,,i],
                               file_name = cat(fileName, "_CM", i, "_K", individualK[i], ".png", sep = ""),
                               savePNG = TRUE)
        }
      }
    }

    ## Localised kernel k-means
    parameters <- list()
    parameters$iteration_count <- C # set the maximum number of iterations
    if(is.null(globalK)){
      KM <- array(0, c(N, N, globalMaxK-1))
      clLabels <- array(NA, c(globalMaxK-1, N))
      weights <- array(NA, c(N, M, globalMaxK - 1))
      for(i in 2:globalMaxK){
        parameters$cluster_count <- i # set the number of clusters K
        lmkkm <- lmkkmeansTrain(CM, parameters)
        weights[,,i-1] <- lmkkm$Theta
        for(j in 1:M){
          KM[,,i-1] <- KM[,,i-1] + (lmkkm$Theta[,j]%*%t(lmkkm$Theta[,j]))*CM[,,j]
        }
      clLabels[i-1,] <- lmkkm$clustering
      }
      maxSil <- maximiseSilhouette(KM, clLabels, maxK = globalMaxK)
      globalK <- output$globalK <- maxSil$k
      output$globalClusterLabels <- clLabels[globalK-1,]
      output$weights <- weights[,,globalK-1]
      output$weightedKM <- KM[,,globalK-1]

    }else{
      parameters$cluster_count <- globalK
      lmkkm <- lmkkmeansTrain(CM, parameters)
      weightedKM <- matrix(0, N, N)
      for(j in 1:M){
        weightedKM <- weightedKM + (lmkkm$Theta[,j]%*%t(lmkkm$Theta[,j]))*CM[,,j]
      }
      output$globalClusterLabels <-lmkkm$clustering
      output$weights <- lmkkm$Theta
      output$weightedKM <- weightedKM
    }
    if(savePlots){
      plotSimilarityMatrix(weightedKM, clusLabels = globalClusterLabels, savePNG = TRUE,
                           file_name = cat(fileName, "_weightedCM", i , "_K", globalK, ".png", sep = ""),
                           savePNG = TRUE)
    }

    output$consensusMatrices <- CM

  return(output)
}
