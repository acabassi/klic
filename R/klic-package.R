#' Kernel Learning Integrative Clustering
#'
#' Kernel Learning Integrative Clustering (KLIC) is an algorithm that allows to
#' combine multiple kernels, each representing a different measure of the
#' similarity between a set of observations. The contribution of each kernel on
#' the final clustering is weighted according to the amount of information
#' carried by it. As well as providing the functions required to perform the
#' kernel-based clustering, this package also allows the user to simply give the
#' data as input: the kernels are then built using consensus clustering.
#' Different strategies to choose the best number of clusters are also
#' available. For further details please see Cabassi and Kirk (2019)
#' <arXiv:1904.07701>.
#'
#' @author
#' Alessandra Cabassi \email{alessandra.cabassi@mrc-bsu.cam.ac.uk},
#' Paul D. W. Kirk,
#' Mehmet Gonen
#' @references Cabassi, A. and Kirk, P. D. W. (2019). Multiple kernel learning
#' for integrative consensus clustering of genomic datasets. arXiv preprint.
#' arXiv:1904.07701.
#'
#' Gonen, M. and Margolin, A.A. (2014). Localized data fusion for
#' kernel k-means clustering with application to cancer biology. In Advances in
#' Neural Information Processing Systems (pp. 1305-1313).
"_PACKAGE"
#> [1] "_PACKAGE"
