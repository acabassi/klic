#' Kernel learning integrative clustering
#'
#'Perform kernel learning integrative clustering
#' @param data List of M datasets, each of size N X P_m, m = 1, ..., M.
#' @param M number of datasets.
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

  ### Data check ###
  N = dim(data[[1]])[1]
  for(i in 1:M){
    if(dim(data[[i]])[1]!=N) stop("All datasets must have the same number of rows.")
  }

  # Initialise empty list for output
  output = list()

  ### Consensus clustering ###
  CM = array(NA, c(N, N, M))
  
  # If individual numbers of clusters are not known
  if(is.null(individualK)){ 
    
      # Initialise empty vector for best number of clusters in each dataset
      output$bestK <- rep(NA, M) 
      # Initialise empty consensus matrices for all possible numbers of clusters
      tempCM = array(NA, c(N, N, individualMaxK-1)) 
      # Initialise empty cluster labels for all possible numbers of clusters
      clLabels = matrix(NA, individualMaxK-1, N) 
      
      # For each dataset 
      for(i in 1:M){ 
        
        # Scale the data such that each column has zero mean and unitary variance
        scaledDataset = scale(data[[i]]) 
        
        # For each possible number of clusters
        for(j in 2:individualMaxK){ 
          
          # Compute consensus matrix
          tempCM[,,j-1] <- consensusCluster(scaledDataset, j, B) 
          # Make consensus matrix positive definite
          tempCM[,,j-1] <- spectrumShift(tempCM[,,j-1]) 
          # Find clusters through hiearchical clustering
          hCl <- hclust(as.dist(1-tempCM[,,j-1]), method = "average") 
          # Extract cluster labels
          clLabels[j-1,] <-  cutree(hCl, j) 
        }
        # Find the number of clusters that maximises the silhouette
        maxSil <- maximiseSilhouette(tempCM, clLabels, individualMaxK) 
        # If there is more than one, choose smallest number of clusters among the ones that maximise the silhouette
        bestK <- output$bestK[i] <- maxSil$k[1] 
        # For dataset i, retain the consensus matrix corresponding to the smallest number of clusters for which the silhouette is maximised
        CM[,,i] <- tempCM[,,bestK] 
        
        # Save plot of similarity matrix
        if(savePlots){
          plotSimilarityMatrix(CM[,,i],
                               file_name = cat(fileName, "_CM", i, "_K", bestK, ".png", sep = ""),
                               savePNG = TRUE)
        }
      }
      
      }else{ # If individual numbers of clusters are known
        
      # For each dataset
      for(i in 1:M){
        
        # Scale the data such that each column has zero mean and unitary variance
        scaledDataset = scale(data[[i]])  
        # Compute consensus matrix
        CM[,,i] <- consensusCluster(scaledDataset, individualK[i], B) 
        # Make consensus matrix positive definite
        CM[,,i] <- spectrumShift(CM[,,i]) 
        
        # Save plot of similarity matrix
        if(savePlots){
          plotSimilarityMatrix(CM[,,i],
                               file_name = cat(fileName, "_CM", i, "_K", individualK[i], ".png", sep = ""),
                               savePNG = TRUE)
        }
      }
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
        lmkkm <- lmkkmeansTrain(CM, parameters) 
        # Save the weights
        weights[,,i-1] <- lmkkm$Theta 
        # Save the combined kernel matrix
        for(j in 1:M){ 
          KM[,,i-1] <- KM[,,i-1] + (lmkkm$Theta[,j]%*%t(lmkkm$Theta[,j]))*CM[,,j]
        }
        # Save the cluster labels
        clLabels[i-1,] <- lmkkm$clustering 
      }
      
      # Find the number of clusters that maximises the silhouette
      maxSil <- maximiseSilhouette(KM, clLabels, maxK = globalMaxK) 
      globalK <- output$globalK <- maxSil$k
      # Save chosen cluster labels
      output$globalClusterLabels <- clLabels[globalK-1,] 
      # Save chosen weights
      output$weights <- weights[,,globalK-1] 
      # Save chosen consensus matrix
      output$weightedKM <- KM[,,globalK-1] 

    }else{ # If the global number of clusters is known
      
      # Set number of clusters for localised multiple kernel k-means
      parameters$cluster_count <- globalK
      # Train localised multiple kernel k-means
      lmkkm <- lmkkmeansTrain(CM, parameters) 
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
    
    # Save plot of kernel matrix 
    if(savePlots){
      plotSimilarityMatrix(weightedKM, clusLabels = globalClusterLabels, savePNG = TRUE,
                           file_name = cat(fileName, "_weightedCM", i , "_K", globalK, ".png", sep = ""),
                           savePNG = TRUE)
    }

  return(output)
}
