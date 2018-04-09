#' Convex Analysis of Mixtures
#'
#' This function performs a fully unsupervised computational deconvolution method
#'     to identify marker genes that define each of the multiple subpopulations
#'     and estimate the proportions of these subpopulations in the mixed samples
#'     as well as their respective expression profiles.
#' @param data Matrix of mixture expression profiles.
#'     Data frame will be internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     data should be in non-log linear space with non-negative numerical values (i.e. >= 0).
#'     Missing values are not supported.
#' @param K The candidate subpopulation number(s), e.g. K = 2:8.
#' @param corner.strategy The method to find corners of convex hull.
#'     1: minimum sum of margin-of-errors; 2: minimum sum of reconstruction errors.
#'     The default is 2.
#' @param dim.rdc Reduced data dimension; should be not less than maximum candidate K.
#' @param thres.low The percentage of genes the user wants to remove with lowest norm.
#'     The range should be between 0 and 1. The default is 0.05.
#' @param thres.high The percentage of genes the user wants to remove with highest norm.
#'     The range should be between 0 and 1. The default is 0.95.
#' @param cluster.method The method to do clustering.
#'     The default "K-Means" will use \code{\link{kmeans}}.
#'     The alternative "apcluster" will use \code{\link[apcluster]{apclusterK}}.
#' @param cluster.num The number of clusters; should be much larger than K.
#'     The default is 50.
#' @param MG.num.thres The clusters with the gene number smaller than
#'     MG.num.thres will be treated as outliers.
#'     The default is 20.
#' @param lof.thres Remove local outlier using \code{\link[DMwR]{lofactor}}.
#'     MG.num.thres is used as the number of neighbours in the calculation of
#'     the local outlier factors.
#'     The default value 0.02 will remove top 2\% local outliers.
#'     Zero value will disable lof.
#' @param cores The number of system cores for parallel computing.
#'     If not provided, the default back-end is used.
#'     Zero value will disable parallel computing.
#' @param seed For reproducibility, the seed of the random number generator for k-Means.
#' @details This function includes three necessary steps to decompose a matrix
#' of mixture expression prefiles: data preprocessing, marker gene cluster
#' search, and matrix decomposition. They are implemented in
#' \code{\link{CAMPrep}}, \code{\link{CAMMGCluster}}, and
#' \code{\link{CAMASest}}, seperately.
#' More details can be find in the help document of each function.
#'
#' For \code{\link{CAM}} function, you needs to specify the range of possible
#' subpopulation numbers and the percentage of low/high-expressed molecueles to
#' be removed. Typically, 30\% ~ 50\% low-expressed genes can be removed from
#' gene expression data. Much less low-expressed proteins are removed,
#' e.g. 0\% ~ 10\%, due to a limited number of proteins in proteomics data.
#' The removel of high-expressed molecules has much less impact on results,
#' and usually set to be 0\% ~ 10\%.
#' @return A list containing the following components:
#' \item{PrepResult}{An object of class "CAMPrepObj" containing
#' data preprocessing results from \code{\link{CAMPrep}} function.}
#' \item{MGResult}{An object of class "CAMMGObj" containg
#' marker gene detection results from \code{\link{CAMMGCluster}} function.}
#' \item{ASestResult}{An object of class "CAMASObj" containing
#' estimated proportions, subpopution-specific expressions and mdl values
#' from \code{\link{CAMASest}} function.}
#' @export
#' @examples
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' #CAM with known subpopulation number
#' rCAM <- CAM(data, K = 3, dim.rdc= 3, thres.low = 0.30, thres.high = 0.95)
#' #A larger dim.rdc can improve performance but increase time complexity
#'
#' \donttest{
#' #CAM with a range of subpopulation number
#' rCAM <- CAM(data, K = 2:5, dim.rdc= 10, thres.low = 0.30, thres.high = 0.95)
#'
#' #Use "apcluster" to aggregate gene vectors in CAM
#' rCAM <- CAM(data, K = 2:5, dim.rdc= 10, thres.low = 0.30, thres.high = 0.95,
#'             cluster.method = 'apcluster')
#' }

CAM <- function(data, K = NULL, corner.strategy = 2, dim.rdc = 10,
                thres.low = 0.05, thres.high = 0.95,
                cluster.method = c('K-Means', 'apcluster'),
                cluster.num = 50, MG.num.thres = 20, lof.thres = 0.02,
                cores = NULL, seed = NULL){
    if (is.null(K)) {
        stop("K is missing")
    }
    if (!is.numeric(K)) {
        stop("K is not numeric")
    }
    if (class(data) == "data.frame") {
        data <- as.matrix(data)
    } else if (class(data) != "matrix") {
        stop("Only matrix and data frame are supported for expression data!")
    }
    if (sum(is.na(data)) > 0) {
        stop("Data with missing values are not supported!")
    }
    if (sum(data<0) > 0) {
        stop("Only non-negative data are supported!")
    }
    if (is.null(rownames(data))) {
        rownames(data) <- seq_len(nrow(data))
    }
    if (dim.rdc < max(K)) {
        warning("dim.rdc is less than max(K)!")
    }

    coreParam <- NULL
    if (length(K) > 1 && (is.null(cores) || cores > 0)) {
        registered()
        if (is.null(cores)) {
            coreParam <- SnowParam(workers = length(K), type = "SOCK")
        } else if (cores > 0) {
            coreParam <- SnowParam(workers = cores, type = "SOCK")
        }
    }


    ################ Preprocessing ################
    message('Preprocessing\n')
    PrepResult <- CAMPrep(data, dim.rdc, thres.low, thres.high, cluster.method,
                          cluster.num, MG.num.thres, lof.thres, seed)

    ################ Marker Gene Selection ################
    message('Marker Gene Selection\n')
    MGResult<-list()

    if (is.null(coreParam)) {
        MGResult <- lapply(K, CAMMGCluster, PrepResult)
    } else {
        MGResult <- bplapply(K, CAMMGCluster, PrepResult,
                            BPPARAM = coreParam)
    }
    names(MGResult) <- as.character(K)

    ################ A and S Matrix Estimation ################
    message('A and S Matrix Estimation\n')
    if (is.null(coreParam)) {
        ASestResult <- lapply(MGResult, CAMASest, PrepResult, data,
                            corner.strategy)
    } else {
        ASestResult <- bplapply(MGResult, CAMASest, PrepResult, data,
                                corner.strategy, BPPARAM = coreParam)
    }
    names(ASestResult) <- as.character(K)

    return(list(PrepResult=PrepResult,
                MGResult=MGResult,
                ASestResult=ASestResult))
}
