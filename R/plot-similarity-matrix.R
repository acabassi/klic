#' Plot similarity matrix
#'
#' Plot similarity matrix. It is possible to plot a side vector that corresponds to a
#' response variable Y and to order the rows (and columns) according to some clustering
#' structure specified by the variable clusLabels.
#' @param X matrix
#' @param y response
#' @param clusLabels cluster labels
#' @param colX colours for the matrix
#' @param colY colours for the response
#' @param myLegend vector of strings with the names of the variables
#' @param file_name name of the file
#' @param savePNG boolean flag. if true, saves the plot as a png file
#' @param saveTIKZ boolean flag. if true saves the plot as a tikz file
#' @param semiSupervised boolean flag. if true, the response is plotted next to the matrix.
#' @param scale parameter scale of heatmap.2. can be "none" or "columns"
#' @param labRow vector of row labels. default NA.
#' @param labCol vector of column labels. default NA.
#' @export
#'
plotSimilarityMatrix = function(X, y = NULL, clusLabels = NULL, colX = NULL, colY = NULL, myLegend = NULL,
                                file_name = "myheatmap.png", savePNG = FALSE, saveTIKZ = FALSE, saveEPS = FALSE, semiSupervised = FALSE,
                                scale = "none", labRow = NA, labCol = NA, dendro = 'none'){

  if(!require(gplots)) require(gplots)

  ### Set the colours for the similarity matrix
  if(is.null(colX)){
    colX= colorRampPalette(c(rgb(232/255, 119/255, 34/255), rgb(1,1,1), rgb(0/255, 114/255, 206/255))) # white - blue
  }else if(colX == "dark"){
    colX= colorRampPalette(c(rgb(190/255, 77/255, 0/255), rgb(1,1,1), rgb(0/255, 60/255, 113/255))) # white - blue
  }else if(colX == "default"){
    colX = "heat.colors"
  }else if(colX == "weights"){
    colX = colorRampPalette(c(rgb(213/255, 0/255, 50/255), rgb(1,1,1), rgb(0/255, 176/255, 185/255)))
  }

  ### Set the options to save the plot in a file
  if(savePNG){
    png(file_name, width = 600, height = 600)
  }else if(saveTIKZ){
    library(tikzDevice)
    tikz(file_name, width = 3, height = 3)
  }else if(saveEPS){
    setEPS()
    postscript(file_name, width = 3, height = 3)
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

    if(is.null(colY)){
      library(wesanderson)
      colY = wes_palette("Royal2")
      colY = colY[c(5,4,3,1,2)]
      colY = c(colY, rgb(108/255, 172/255, 228/255))
    }

    if(dim(table(clusLabels))>length(colY)) {
      cl <- colors(distinct = TRUE)
      set.seed(15887) # to set random generator seed
      colY <- sample(cl, dim(table(clusLabels)))#stop("There are more classes than colours available.")
    }

    rowseparators = cumsum(table(clusLabels))

    if(dendro == 'both'){
      Rowvv = TRUE
      Colvv = TRUE
    }else{
      Rowvv = FALSE
      Colvv = FALSE
    }

    heatmap.2(X_new,
              rowsep = rowseparators,
              sepwidth = c(0.1,0.1),
              sepcolor = "black",
              scale = scale,
              dendrogram = dendro,
              Rowv = Rowvv,
              Colv = Colvv,
              RowSideColors = colY[y_new],
              trace = 'none',
              col = colX,
              labRow = labRow,
              labCol = labCol,
              key = FALSE#,
              # lmat = rbind(c(0,3),c(2,1),c(0,4)),
              # lwid = c(1,8),
              # lhei = c(1,8,1)
              # lmat=rbind(c(2),c(3),c(1),c(4)),
              # lhei=c(1,1,10,1),
              # lwid=c(10)
    )
    if(!is.null(myLegend)) legend("bottomleft", legend = myLegend, col = colY, lty = 1, lwd = 4, bty ="n")
  }else{

    if(dendro == 'both'){
      Rowvv = TRUE
      Colvv = TRUE
    }else{
    ### in the unsupervised setting...
      Rowvv = FALSE
      Colvv = FALSE
    }

    heatmap.2(X,
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
              symkey = FALSE#,
              # lmat = rbind(c(0,3),c(2,1),c(0,4)),
              # lwid = c(1,8),
              # lhei = c(1,8,1)
              # lmat=rbind(c(2),c(3),c(1),c(4)),
              # lhei=c(1,1,10,1),
              # lwid=c(10)
    )
    if(!is.null(myLegend)) legend("bottomleft", legend = myLegend, col = colY, lty = 1, lwd = 4, bty ="n")
  }

  if(savePNG || saveTIKZ) dev.off()

}
