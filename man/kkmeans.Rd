% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kkmeans.R
\name{kkmeans}
\alias{kkmeans}
\title{Kernel k-means}
\usage{
kkmeans(K, parameters, seed = NULL)
}
\arguments{
\item{K}{Kernel matrix.}

\item{parameters}{A list containing the number of clusters
\code{number_count}.}

\item{seed}{The seed used inside the \code{kmeans} function. By default it is
NULL, which means that it is randomly set by the \code{set.seed} function.}
}
\value{
This function returns a list containing:
\item{clustering}{the cluster labels for each element (i.e. row/column) of
the kernel matrix.}
\item{objective}{the value of the objective function for the given
clustering.}
\item{parameters}{same parameters as in the input.}
}
\description{
Perform the training step of kernel k-means.
}
\examples{
# Load one dataset with 100 observations, 2 variables, 4 clusters
data <- as.matrix(read.csv(system.file("extdata", "dataset1.csv",
package = "klic"), row.names = 1))
# Compute consensus clustering with K=4 clusters
cm <- coca::consensusCluster(data, 4)
# Shift eigenvalues of the matrix by a constant: (min eigenvalue) * (coeff)
km <- spectrumShift(cm, coeff = 1.05)
# Initalize the parameters of the algorithm
parameters <- list()
# Set the number of clusters
parameters$cluster_count <- 4
# Perform training
state <- kkmeans(km, parameters)
# Display the clustering
print(state$clustering)
}
\references{
Gonen, M. and Margolin, A.A., 2014. Localized data fusion for
kernel k-means clustering with application to cancer biology. In Advances in
Neural Information Processing Systems (pp. 1305-1313).
}
\author{
Mehmet Gonen
}
