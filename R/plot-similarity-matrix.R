#' Plot similarity matrix with pheatmap
#'
#' @param X Similarity matrix.
#' @param y Vector
#' @param clusLabels Cluster labels
#' @param colX Colours for the matrix
#' @param colY Colours for the response
#' @param myLegend Vector of strings with the names of the variables
#' @param savePNG Boolean: if TRUE, the plot is saved as a png file. Default is
#' FALSE.
#' @param fileName If "savePNG" is TRUE, this is the string containing the name
#' of the output file. Can be used to specify the folder path too. Default is
#' "posteriorSimilarityMatrix". The extension ".png" is automatically added to
#' this string.
#' @param semiSupervised Boolean flag: if TRUE, the response is plotted next to
#' the matrix.
#' @param showObsNames Boolean. If TRUE, observation names are shown in the
#' plot. Default is FALSE.
#' @param clr Boolean. If TRUE, rows are ordered by hierarchical clustering.
#' Default is FALSE.
#' @param clc Boolean. If TRUE, columns are ordered by hierarchical clustering.
#' Default is FALSE.
#' @param plotWidth Plot width. Default is 500.
#' @param plotHeight Plot height. Default is 450.
#' @author Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk}
#' @examples
#' # Load one dataset with 100 observations, 2 variables, 4 clusters
#' data <- as.matrix(read.csv(system.file('extdata', 'dataset1.csv',
#' package = 'klic'), row.names = 1))
#' # Load cluster labels
#' cluster_labels <- as.matrix(read.csv(system.file('extdata',
#' 'cluster_labels.csv', package = 'klic'), row.names = 1))
#'
#' # Compute consensus clustering with K=4 clusters
#' cm <- coca::consensusCluster(data, 4)
#'
#' # Plot consensus (similarity) matrix
#' plotSimilarityMatrix(cm)
#'
#' # Plot consensus (similarity) matrix with response
#' names(cluster_labels) <- as.character(1:100)
#' rownames(cm) <- names(cluster_labels)
#' plotSimilarityMatrix(cm, y = cluster_labels)
#' @export


plotSimilarityMatrix = function(X,
                                y = NULL,
                                clusLabels = NULL,
                                colX = NULL,
                                colY = NULL,
                                myLegend = NULL,
                                fileName = "posteriorSimilarityMatrix",
                                savePNG = FALSE,
                                semiSupervised = FALSE,
                                showObsNames = FALSE,
                                clr = FALSE,
                                clc = FALSE,
                                plotWidth = 500,
                                plotHeight = 450) {


    if (!is.null(y)) {
        # Check if the rownames correspond to the ones in the similarity matrix
        check <- sum(1 - rownames(X) %in% row.names(y))
        if (check == 1)
            stop("X and y must have the same row names.")
    }

    if (!is.null(clusLabels)) {

        if (!is.integer(clusLabels))
            stop("Cluster labels must be integers.")

        n_clusters <- length(table(clusLabels))
        riordina <- NULL
        for (i in 1:n_clusters) {
            riordina <- c(riordina, which(clusLabels == i))
        }

        X <- X[riordina, riordina]
        y <- y[riordina, ]
        y <- as.data.frame(y)

    }

    if (savePNG)
        grDevices::png(paste(fileName, ".png", sep = ""),
                       width = plotWidth,
                       height = plotHeight)


    if (!is.null(y)) {
        pheatmap::pheatmap(
            X,
            legend = TRUE,
            color = c("white", (
                RColorBrewer::brewer.pal(n = 6,
                                         name = "PuBu")
            )),
            cluster_rows = clr,
            cluster_cols = clc,
            annotation_col = y,
            show_rownames = showObsNames,
            show_colnames = showObsNames,
            drop_levels = FALSE,
            na_col = "seashell2"
        )
    } else {
        pheatmap::pheatmap(
            X,
            legend = TRUE,
            color = c("white", (
                RColorBrewer::brewer.pal(n = 6,
                                         name = "PuBu")
            )),
            cluster_rows = clr,
            cluster_cols = clc,
            show_rownames = showObsNames,
            show_colnames = showObsNames,
            drop_levels = FALSE,
            na_col = "seashell2"
        )
    }


    if (savePNG)
        grDevices::dev.off()

}
