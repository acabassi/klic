## ----generate data, fig.show='hold', warning=FALSE, cache=TRUE-----------
## Load synthetic data
data1 <- as.matrix(read.csv(system.file("extdata", "dataset1.csv", package = "klic"), row.names = 1))
data2 <- as.matrix(read.csv(system.file("extdata","dataset2.csv", package = "klic"), row.names = 1))
data3 <- as.matrix(read.csv(system.file("extdata", "dataset3.csv", package = "klic"), row.names = 1))
data <- list(data1, data2, data3)
n_datasets <- 3
N <- dim(data[[1]])[1]

## ----consensus_cluster, fig.show='hold', warning=FALSE, cache=TRUE-------
## Compute co-clustering matrices for each dataset
library(klic)
CM <- array(NA, c(N, N, n_datasets))
for(i in 1: n_datasets){
  # Scale the columns to have zero mean and unitary variance
  scaledData <- scale(data[[i]]) 
  # Use consensus clustering to find the consensus matrix of each dataset
  CM[,,i] <- coca::consensusCluster(scaledData, K = 6, B = 50)
}
## Plot consensus matrix of one of the datasets
plotSimilarityMatrix(CM[,,3])

## ----spectrum_shift, fig.show='hold', warning=FALSE, cache=TRUE----------
## Check if consensus matrices are PSD and shift eigenvalues if needed.
for(i in 1: n_datasets){
  CM[,,i] <- spectrumShift(CM[,,i], verbose = FALSE)
}
## Plot updated consensus matrix of one of the datasets
plotSimilarityMatrix(CM[,,3])

## ----lmk_kmeans, fig.show='hold', warning=FALSE, cache=TRUE--------------
## Perform localised kernel k-means on the consensus matrices
library(Matrix)
library(Rmosek)
parameters <- list()
parameters$cluster_count <- 6 # set the number of clusters K
parameters$iteration_count <- 100 # set the maximum number of iterations
lmkkm <- lmkkmeans(CM, parameters)

## ----ari, fig.show='hold', warning=FALSE, cache=TRUE---------------------
## Compare clustering found with KLIC to the true one
ones <- rep(1, N/parameters$cluster_count)
true_labels <- c(ones, ones*2, ones*3, ones*4, ones*5, ones*6)
library(mclust, verbose = FALSE)
adjustedRandIndex(true_labels, lmkkm$clustering) 

## ----maximise_silhouette, fig.show='hold', warning=FALSE, cache=TRUE-----
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
maxSil <- coca::maximiseSilhouette(KM, clLabels, maxK = 6)
maxSil$k

## ----klic, fig.show='hold', warning=FALSE, cache=TRUE--------------------
klic <- klic(data, M = n_datasets, individualK = c(6,6,6))
klic$globalK

## ----k_kmeans, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE----
## Set parameters of the kernel k-means algorithm
parameters <- list()
parameters$cluster_count <- 6
parameters$iteration_count <- 100
## Run kernel k-means
kkm <- kkmeans(CM[,,3], parameters)
## Compare clustering to the true labels
clusterLabels <- kkm$clustering
adjustedRandIndex(true_labels, lmkkm$clustering) 

## ----cophenetic_correlation, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE----
## Compute cophenetic correlation coefficient for each consensus matrix
ccc <- rep(NA, n_datasets)
for(i in 1:n_datasets){
  ccc[i] <- copheneticCorrelation(CM[,,i])
}
ccc

