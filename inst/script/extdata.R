n_elements_per_cluster <- 25
n_clusters <- 4

# Vector of true cluster labels
true_labels <- c(rep(1, n_elements_per_cluster),
                 rep(2, n_elements_per_cluster),
                 rep(3, n_elements_per_cluster),
                 rep(4, n_elements_per_cluster))
write.csv(true_labels, "inst/extdata/cluster_labels.csv")

# Datasets
set.seed(1)
dataset1 <- t(mvtnorm::rmvnorm(2, true_labels*3))
write.csv(dataset1, "inst/extdata/dataset1.csv")
set.seed(2)
dataset2 <- t(mvtnorm::rmvnorm(2, true_labels*3))
write.csv(dataset2, "inst/extdata/dataset2.csv")
set.seed(3)
dataset3 <- t(mvtnorm::rmvnorm(2, true_labels*3))
write.csv(dataset3, "inst/extdata/dataset3.csv")

# Consensus matrices
consensus_matrix1 <- coca::consensusCluster(dataset1, K = n_clusters)
write.csv(consensus_matrix1, "inst/extdata/consensus_matrix1.csv")
consensus_matrix2 <- coca::consensusCluster(dataset2, K = n_clusters)
write.csv(consensus_matrix2, "inst/extdata/consensus_matrix2.csv")
consensus_matrix3 <- coca::consensusCluster(dataset3, K = n_clusters)
write.csv(consensus_matrix3, "inst/extdata/consensus_matrix3.csv")

# Kernel matrices
shift = 0.1
temp_kernel <- consensus_matrix1 + diag(dim(consensus_matrix1)[1]) * shift
kernel_matrix1 <- temp_kernel/temp_kernel[1, 1]  # Rescale
write.csv(kernel_matrix1, "inst/extdata/kernel_matrix1.csv")
temp_kernel <- consensus_matrix2 + diag(dim(consensus_matrix2)[1]) * shift
kernel_matrix2 <- temp_kernel/temp_kernel[1, 1]  # Rescale
write.csv(kernel_matrix2, "inst/extdata/kernel_matrix2.csv")
temp_kernel <- consensus_matrix3 + diag(dim(consensus_matrix3)[1]) * shift
kernel_matrix3 <- temp_kernel/temp_kernel[1, 1]  # Rescale
write.csv(kernel_matrix3, "inst/extdata/kernel_matrix3.csv")
