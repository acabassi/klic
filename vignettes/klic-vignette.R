## ----cophenetic_correlation, fig.show='hold', message=FALSE, warning=FALSE, cache=TRUE----
## Compute cophenetic correlation coefficient for each consensus matrix
ccc <- rep(NA, n_datasets)
for(i in 1:n_datasets){
  ccc[i] <- copheneticCorrelation(CM[,,i])
}
ccc

