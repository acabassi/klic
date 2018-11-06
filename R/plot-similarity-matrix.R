#' Plot similarity matrix
#'
#' Plot similarity matrix. It is possible to plot a side vector that corresponds to a
#' response variable Y and to order the rows (and columns) according to some clustering
#' structure specified by the variable clusLabels.
#' @param X Matrix
#' @param y Response
#' @param clusLabels Cluster labels
#' @param colX Colours for the matrix
#' @param colY Colours for the response
#' @param myLegend Vector of strings with the names of the variables
#' @param file_name Name of the file
#' @param savePNG Boolean flag: if TRUE, the plot is saved as a png file.
#' @param saveTIKZ Boolean flag: if TRUE, the plot is saved as a tikz file.
#' @param saveEPS Boolean flag: if TRUE, the plot is saved as an eps file.
#' @param semiSupervised Boolean flag: if TRUE, the response is plotted next to the matrix.
#' @param scale Used as input for the parameter "scale" of the gplot::heatmap.2() function.
#' Can be either "none" or "columns".
#' @param labRow Vector of row labels, default is NA.
#' @param labCol Vector of column labels, default is NA.
#' @param dendro If 'both', plot dendrogram on rows and columns, if 'none' no dendrograms are shown. Default is 'none'.
#' @examples
#' # Load one dataset with 300 observations, 2 variables, 6 clusters
#' data <- as.matrix(read.csv(system.file("extdata", "dataset1.csv",
#' package = "klic"), row.names = 1))
#' # Compute consensus clustering with K=6 clusters
#' cm <- coca::consensusCluster(data, 6)
#' # Plot consensus (similarity) matrix
#' plotSimilarityMatrix(cm)
#' @export
#'
plotSimilarityMatrix = function(X, y = NULL, clusLabels = NULL, colX = NULL, colY = NULL,
                                myLegend = NULL, file_name = "myheatmap.png", savePNG = FALSE,
                                saveTIKZ = FALSE, saveEPS = FALSE, semiSupervised = FALSE,
                                scale = "none", labRow = NA, labCol = NA, dendro = "none"){

  ### Set the colours for the similarity matrix
  if(is.null(colX)){
    colX= grDevices::colorRampPalette(c(grDevices::rgb(232/255, 119/255, 34/255),
                                        grDevices::rgb(1,1,1),
                                        grDevices::rgb(0/255, 114/255, 206/255))) # white - blue
  }else if(colX == "dark"){
    colX= grDevices::colorRampPalette(c(grDevices::rgb(190/255, 77/255, 0/255),
                                        grDevices::rgb(1,1,1),
                                        grDevices::rgb(0/255, 60/255, 113/255))) # white - blue
  }else if(colX == "default"){
    colX = "heat.colors"
  }else if(colX == "weights"){
    colX = grDevices::colorRampPalette(c(grDevices::rgb(213/255, 0/255, 50/255),
                                         grDevices::rgb(1,1,1),
                                         grDevices::rgb(0/255, 176/255, 185/255)))
  }

  ### Set the options to save the plot in a file
  if(savePNG){
    grDevices::png(file_name, width = 600, height = 600)
  }else if(saveTIKZ){
    tikzDevice::tikz(file_name, width = 3, height = 3)
  }else if(saveEPS){
    grDevices::setEPS()
    grDevices::postscript(file_name, width = 3, height = 3)
  }

  ### If the input data include a response variable...
  if(semiSupervised){

    if(!is.integer(clusLabels))
      stop("Cluster labels must be integers.")

    n_clusters = length(table(clusLabels))
    riordina = NULL
    for (i in 1:n_clusters){
      riordina = c(riordina, which(clusLabels==i))
    }
    X_new = X[riordina,riordina]
    y_new = y[riordina]

    y_new_int <- rep(0, length(y_new))
    # Replace values in y_new with integers
    count <- 0
    for(i in unique(y_new)){
        count <- count + 1
        y_new_int[which(y_new==i)] = count
    }

    ### LEGEND MUST BE REORDERED TOO ###

    if(is.null(colY)){
      colY = wesanderson::wes_palette("Royal2")
      colY = colY[c(5,4,3,1,2)]
      colY = c(colY, grDevices::rgb(108/255, 172/255, 228/255))
    }

    if(dim(table(clusLabels))>length(colY)) {
      cl <- grDevices::colors(distinct = TRUE)
      set.seed(151) # to set random generator seed
      colY <- sample(cl, dim(table(clusLabels)))
      #stop("There are more classes than colours available.")
    }

    rowseparators = cumsum(table(clusLabels))

    if(dendro == 'both'){
      Rowvv = TRUE
      Colvv = TRUE
    }else{
      Rowvv = FALSE
      Colvv = FALSE
    }

    gplots::heatmap.2(X_new,
              rowsep = rowseparators,
              sepwidth = c(0.1,0.1),
              sepcolor = "black",
              scale = scale,
              dendrogram = dendro,
              Rowv = Rowvv,
              Colv = Colvv,
              RowSideColors = colY[y_new_int],
              trace = 'none',
              col = colX,
              labRow = labRow,
              labCol = labCol,
              key = FALSE
    )
    if(!is.null(myLegend)) graphics::legend("bottomleft", legend = myLegend, col = colY,
                                            lty = 1, lwd = 4, bty ="n")
  }else{

    if(dendro == 'both'){
      Rowvv = TRUE
      Colvv = TRUE
    }else{
    ### in the unsupervised setting...
      Rowvv = FALSE
      Colvv = FALSE
    }

    gplots::heatmap.2(X,
              # rowsep = rowseparators,
              # sepwidth = c(0.1,0.1),
              # sepcolor = "black",
              scale = scale,
              dendrogram = dendro,
              Rowv = Rowvv,
              Colv = Colvv,
              # RowSideColors = colY[y_new],
              trace = 'none',
              col = colX,
              labRow = labRow,
              labCol = labCol,
              # key = FALSE,
              symkey = FALSE
    )
    if(!is.null(myLegend)) graphics::legend("bottomleft", legend = myLegend, col = colY,
                                            lty = 1, lwd = 4, bty ="n")
  }

  if(savePNG || saveTIKZ) grDevices::dev.off()

}
