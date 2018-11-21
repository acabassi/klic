#' Plot similarity matrix with pheatmap
#'
#' @param X Similarity matrix.
#' @param y Vector
#' @param clusLabels Cluster labels
#' @param colX Colours for the matrix
#' @param colY Colours for the response
#' @param myLegend Vector of strings with the names of the variables
#' @param fileName Name of pdf file
#' @param save Boolean flag: if TRUE, the plot is saved as a png file.
#' @param semiSupervised Boolean flag: if TRUE, the response is plotted next to the matrix.
#' @param scale Used as input for the parameter "scale" of the gplot::heatmap.2() function.
#' Can be either "none" or "columns".
#' @param showObsNames Boolean. If TRUE, observation names are shown in the plot. Default is FALSE.
#' @param clr Boolean. If TRUE, rows are ordered by hierarchical clustering. Default is FALSE.
#' @param clc Boolean. If TRUE, columns are ordered by hierarchical clustering. Default is FALSE.
#' @param plotWidth Plot width. Default is 500.
#' @param plotHeight Plot height. Default is 450.
#' @examples
#' # Load one dataset with 300 observations, 2 variables, 6 clusters
#' data <- as.matrix(read.csv(system.file("extdata", "dataset1.csv",
#' package = "klic"), row.names = 1))
#' cluster_labels <- as.matrix(read.csv(system.file("extdata",
#' "cluster_labels.csv", package = "klic"), row.names = 1))
#' # Compute consensus clustering with K=6 clusters
#' cm <- coca::consensusCluster(data, 6)
#' # Plot consensus (similarity) matrix
#' plotSimilarityMatrix2(cm)
#' # Plot consensus (similarity) matrix with response
#' names(cluster_labels) <- as.character(1:300)
#' rownames(cm) <- names(cluster_labels)
#' plotSimilarityMatrix2(cm, y = cluster_labels)
#' @export

plotSimilarityMatrix2 = function(X, y = NULL, clusLabels = NULL, colX = NULL, colY = NULL,
                                 myLegend = NULL, fileName = "posteriorSimilarityMatrix.pdf", save = FALSE,
                                 semiSupervised = FALSE,
                                 scale = "none",  showObsNames = FALSE,
                                 clr = FALSE, clc = FALSE,
                                 plotWidth = 500, plotHeight = 450){

    if(!is.null(y)){
        # Check if the rownames correspond to the ones in the similarity matrix
        check <- sum(1-rownames(X)%in%row.names(y))
        if(check==1) stop('X and y must have the same row names.')
    }

    if(!is.null(clusLabels)){

        if(!is.integer(clusLabels))
            stop("Cluster labels must be integers.")

        n_clusters <- length(table(clusLabels))
        riordina <- NULL
        for (i in 1:n_clusters){
            print(i)
            print(riordina)
            print(which(clusLabels==i))
            riordina <- c(riordina, which(clusLabels==i))
        }
        print(riordina)
        X <- X[riordina,riordina]
        y <- y[riordina,]
        y <- as.data.frame(y)
        print(dim(X))
    }

    if(save) grDevices::png(fileName, width = plotWidth, height = plotHeight)

    if(!is.null(y)){
        pheatmap::pheatmap(X, legend = TRUE,
                           color =  c("white", (RColorBrewer::brewer.pal(n = 6, name = "PuBu"))),
                           cluster_rows = clr,
                           cluster_cols = clc,
                           annotation_col = y,
                           show_rownames = showObsNames,
                           show_colnames = showObsNames,
                           drop_levels = FALSE, na_col = "seashell2")
    }else{
        pheatmap::pheatmap(X, legend = TRUE,
                           color =  c("white", (RColorBrewer::brewer.pal(n = 6, name = "PuBu"))),
                           cluster_rows = clr,
                           cluster_cols = clc,
                           show_rownames = showObsNames,
                           show_colnames = showObsNames,
                           drop_levels = FALSE, na_col = "seashell2")
    }


    if(save){
        grDevices::dev.off()
        warning('After saving a pheatmap plot to file, you sometimes have to repeat the
                `dev.off()` command in order to shut down the plotting device completely.')
    }

}
