## ----consensus_cluster, fig.show='hold', warning=FALSE, cache=TRUE------------
## Compute co-clustering matrices for each dataset
library(klic)
CM <- array(NA, c(N, N, n_datasets))
for(i in 1: n_datasets){
  # Scale the columns to have zero mean and unitary variance
  scaledData <- scale(data[[i]]) 
  # Use consensus clustering to find the consensus matrix of each dataset
  CM[,,i] <- coca::consensusCluster(scaledData, K = 4, B = 50)
}
## Plot consensus matrix of one of the datasets
trueLab <- as.factor(trueLab)
names(trueLab) <- as.character(1:N)
CM3 <- as.matrix(CM[,,3])
rownames(CM3) <- colnames(CM3) <- names(trueLab)
plotSimilarityMatrix(CM3, y = as.data.frame(trueLab))

## ----lmk_kmeans, fig.show='hold', warning=FALSE, cache=TRUE-------------------
## Perform localised kernel k-means on the consensus matrices
library(Matrix)
library(Rmosek)
parameters <- list()
parameters$cluster_count <- 4 # set the number of clusters K
parameters$iteration_count <- 100 # set the maximum number of iterations
lmkkm <- lmkkmeans(CM, parameters)

## ----maximise_silhouette, fig.show='hold', warning=FALSE, cache=TRUE----------
## Find the value of k that maximises the silhouette

# Initialise array of kernel matrices 
maxK = 6
KM <- array(0, c(N, N, maxK-1))
clLabels <- array(NA, c(maxK-1, N))

parameters <- list()
parameters$iteration_count <- 100 # set the maximum number of iterations

for(i in 2:maxK){
  
  # Use kernel k-means with K=i to find weights and cluster labels
  parameters$cluster_count <- i # set the number of clusters K
  lmkkm <- lmkkmeans(CM, parameters)
  
  # Compute weighted matrix
  for(j in 1:dim(CM)[3]){
    KM[,,i-1] <- KM[,,i-1] + (lmkkm$Theta[,j]%*%t(lmkkm$Theta[,j]))*CM[,,j]
  }
  
  # Save cluster labels
  clLabels[i-1,] <- lmkkm$clustering 
}

# Find value of K that maximises silhouette
maxSil <- coca::maximiseSilhouette(KM, clLabels, maxK = 4)
maxSil$k

## ----klic, fig.show='hold', warning=FALSE, cache=TRUE-------------------------
klic <- klic(data, M = n_datasets, individualK = c(4, 4, 4))
klic$globalK

## ----k_kmeans, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE------
## Set parameters of the kernel k-means algorithm
parameters <- list()
parameters$cluster_count <- 4
parameters$iteration_count <- 100
## Run kernel k-means
kkm <- kkmeans(CM[,,3], parameters)
## Compare clustering to the true labels
clusterLabels <- kkm$clustering
adjustedRandIndex(trueLab, lmkkm$clustering) 

## ----cophenetic_correlation, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE----
## Compute cophenetic correlation coefficient for each consensus matrix
ccc <- rep(NA, n_datasets)
for(i in 1:n_datasets){
  ccc[i] <- copheneticCorrelation(CM[,,i])
}
ccc

