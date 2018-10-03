#' Plot silhouette
#'
#' Plot average silhouette.
#' @param sil vector of the average silhouette for K from 2 to some value maxK
#' @param fileName name of the png file
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
plotSilhouette = function(sil, fileName){

  fileName = paste(fileName, "png", sep = ".")
  grDevices::png(fileName, width = 400, height = 400)
  graphics::plot(sil, xlab = "Number of clusters", ylab = "Average silhouette")
  grDevices::dev.off()

}

#' Choose K that maximises the silhouette from a set kernel matrices and clusterings
#'
#' @param kernelMatrix N X N X (maxK-1) array of kernel matrices.
#' @param clLabels (maxK-1) X N matrix containing the clusterings obtained for different values of K
#' @param maxK Maximum number of clusters considered.
#' @param savePNG If TRUE, a plot of the silhouette is saved in the working folder. Defaults to FALSE.
#' @param fileName If savePNG is TRUE, this is the name of the png file.
#' @author Alessandra Cabassi \email{ac2051@cam.ac.uk}
#' @export
maximiseSilhouette = function(kernelMatrix, clLabels, maxK,
                              savePNG = FALSE, fileName = "silhouette"){

  # initialise vector of average silhouette
  sil <- rep(5, maxK-1)

  for(i in 2:maxK){
    # use the kernel matrix as distance matrix
    DM <- stats::as.dist(1 - kernelMatrix[,,i-1], diag = FALSE, upper = FALSE)

    # calculate the average silhouette over all the points
    sil[i-1]<- mean(cluster::silhouette(clLabels[i-1,], DM)[,'sil_width'])
  }
  if(savePNG) plotSilhouette(sil, paste('silhouette/', fileName, '.png', sep = ''))

  return(list(silhouette = sil, k = which.max(sil)+1 ))
  # the "+1" is there because the silhouette is stored in a vector where element i corresponds to number of clusters equal to i+1
}
