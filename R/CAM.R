#' Convex Analysis of Mixtures
#'
#' This function performs a fully unsupervised computational deconvolution
#'     to identify marker genes that define each of the multiple
#'     subpopulations, and estimate the proportions of these subpopulations in
#'     the mixture tissues as well as their respective expression profiles.
#' @param data Matrix of mixture expression profiles.
#'     Data frame, SummarizedExperiment or ExpressionSet object will be
#'     internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     Data should be in non-log linear space with non-negative numerical values
#'     (i.e. >= 0). Missing values are not supported.
#'     All-zero rows will be removed internally.
#' @param K The candidate subpopulation number(s), e.g. K = 2:8.
#' @param corner.strategy The method to find corners of convex hull.
#'     1: minimum sum of margin-of-errors; 2: minimum sum of reconstruction
#'     errors. The default is 2.
#' @param dim.rdc Reduced data dimension;
#'     should be not less than maximum candidate K.
#' @param thres.low The lower bound of percentage of genes to keep for CAM
#'     with ranked norm. The value should be between 0 and 1.
#'     The default is 0.05.
#' @param thres.high The higher bound of percentage of genes to keep for CAM
#'     with ranked norm. The value should be between 0 and 1.
#'     The default is 0.95.
#' @param cluster.method The method to do clustering.
#'     The default "K-Means" will use \code{\link{kmeans}}.
#'     The alternative "apcluster" will use
#'     \code{\link[apcluster]{apclusterK-methods}}.
#' @param cluster.num The number of clusters; should be much larger than K.
#'     The default is 50.
#' @param MG.num.thres The clusters with the gene number smaller than
#'     MG.num.thres will be treated as outliers.
#'     The default is 20.
#' @param lof.thres Remove local outlier using \code{\link[DMwR]{lofactor}}.
#'     MG.num.thres is used as the number of neighbors in the calculation of
#'     the local outlier factors.
#'     The default value 0.02 will remove top 2\% local outliers.
#'     Zero value will disable lof.
#' @param quick.select  The number of candidate corners kept after quickhull
#'     and SFFS greedy search. If Null, only quickhull is applied.
#'     The default is 20. If this value is larger than the number of candidate
#'     corners after quickhull, greedy search will also not be applied.
#' @param sample.weight Vector of sample weights. If NULL, all samples have
#'     the same weights. The length should be the same as sample numbers.
#'     All values should be positive.
#' @param generalNMF If TRUE, the decomposed proportion matrix has no sum-to-one
#'     constraint for each row. The default is FALSE.
#'     TRUE value brings two changes: (1) Without assuming samples are
#'     normalized, the first principal component will not forced to be along
#'     c(1,1,..,1) but a standard PCA will be applied during preprocessing.
#'     (2) Without sum-to-one constraint for each row, the scale ambiguity of
#'     each column vector in proportion matrix will not be removed.
#' @param cores The number of system cores for parallel computing.
#'     If not provided, one core for each element in K will be invoked.
#'     Zero value will disable parallel computing.
#' @details This function includes three necessary steps to decompose a matrix
#' of mixture expression profiles: data preprocessing, marker gene cluster
#' search, and matrix decomposition. They are implemented in
#' \code{\link{CAMPrep}}, \code{\link{CAMMGCluster}} and
#' \code{\link{CAMASest}}, separately.
#' More details can be found in the help document of each function.
#'
#' For this function, you needs to specify the range of possible
#' subpopulation numbers and the percentage of low/high-expressed genes to
#' be removed. Typically, 30\% ~ 50\% low-expressed genes can be removed from
#' gene expression data. The removal of high-expressed genes has much less
#' impact on results, and usually set to be 0\% ~ 10\%.
#'
#' This function can also analyze other molecular expression data, such as
#' proteomics data. Much less low-expressed proteins need to be removed,
#' e.g. 0\% ~ 10\%, due to a limited number of proteins without missing values.

#' @return An object of class "\code{\link{CAMObj}}" containing the following
#' components:
#' \item{PrepResult}{An object of class "\code{\link{CAMPrepObj}}" containing
#' data preprocessing results from \code{\link{CAMPrep}} function.}
#' \item{MGResult}{A list of "\code{\link{CAMMGObj}}" objects containing
#' marker gene detection results from \code{\link{CAMMGCluster}} function for
#' each K value.}
#' \item{ASestResult}{A list of "\code{\link{CAMASObj}}" objects containing
#' estimated proportions, subpopulation-specific expressions and mdl values
#' from \code{\link{CAMASest}} function for each K value.}
#' @export
#' @examples
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' #set seed to generate reproducible results
#' set.seed(111)
#'
#' #CAM with known subpopulation number
#' rCAM <- CAM(data, K = 3, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)
#' #Larger dim.rdc can improve performance but increase time complexity
#'
#' \dontrun{
#' #CAM with a range of subpopulation number
#' rCAM <- CAM(data, K = 2:5, dim.rdc = 10, thres.low = 0.30, thres.high = 0.95)
#'
#' #Use "apcluster" to aggregate gene vectors in CAM
#' rCAM <- CAM(data, K = 2:5, dim.rdc = 10, thres.low = 0.30, thres.high = 0.95,
#'             cluster.method = 'apcluster')
#'
#' #CAM with quick selection to reduce time complexity
#' rCAM <- CAM(data, K = 3, dim.rdc = 10, thres.low = 0.30, thres.high = 0.95,
#'         quick.select = 20)
#'
#' #CAM with different sample weights (e.g. adjusted based on sample quality)
#' rCAM <- CAM(data, K = 3, dim.rdc = 5, thres.low = 0.30, thres.high = 0.95,
#'         sample.weight = c(rep(10,11),rep(1,10)))
#'
#' #CAM for general NMF (no sum-to-one contraint for proportion matrix)
#' rCAM <- CAM(data, K = 3, dim.rdc = 5, thres.low = 0.30, thres.high = 0.95,
#'        generalNMF = TRUE)
#' }
CAM <- function(data, K = NULL, corner.strategy = 2, dim.rdc = 10,
                thres.low = 0.05, thres.high = 0.95,
                cluster.method = c('K-Means', 'apcluster'),
                cluster.num = 50, MG.num.thres = 20, lof.thres = 0.02,
                quick.select = NULL, sample.weight = NULL,
                generalNMF = FALSE, cores = NULL){
    if (is.null(K)) {
        stop("K is missing")
    }
    if (!is.numeric(K)) {
        stop("K is not numeric")
    }
    if (is(data, "data.frame")) {
        data <- as.matrix(data)
    } else if (is(data, "SummarizedExperiment")) {
        data <- SummarizedExperiment::assay(data)
    } else if (is(data, "ExpressionSet")) {
        data <- Biobase::exprs(data)
    } else if (is(data, "matrix") == FALSE) {
        stop("Only matrix, data frame, SummarizedExperiment and ExpressionSet
            object are supported for expression data!")
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

    data <- data[rowSums(data) > 0,]

    ################ Preprocessing ################
    message('Preprocessing\n')
    PrepResult <- CAMPrep(data, dim.rdc, thres.low, thres.high, cluster.method,
                        cluster.num, MG.num.thres, lof.thres, quick.select,
                        sample.weight, generalNMF)

    ################ Marker Gene Selection ################
    coreParam <- NULL
    if (length(K) > 1 && (is.null(cores) || cores > 0)) {
        registered()
        if (is.null(cores)) {
            coreParam <- SnowParam(workers = length(K), type = "SOCK")
        } else if (cores > 0) {
            coreParam <- SnowParam(workers = cores, type = "SOCK")
        }
    }

    message('Marker Gene Selection\n')
    MGResult<-list()

    if (is.null(coreParam)) {
        MGResult <- lapply(K, CAMMGCluster, PrepResult, generalNMF)
    } else {
        MGResult <- bplapply(K, CAMMGCluster, PrepResult, generalNMF,
                            BPPARAM = coreParam)
    }
    names(MGResult) <- as.character(K)

    ################ A and S Matrix Estimation ################
    message('A and S Matrix Estimation\n')
    if (is.null(coreParam)) {
        ASestResult <- lapply(MGResult, CAMASest, PrepResult, data,
                            corner.strategy, generalNMF)
    } else {
        ASestResult <- bplapply(MGResult, CAMASest, PrepResult, data,
                                corner.strategy, generalNMF,
                                BPPARAM = coreParam)
    }
    names(ASestResult) <- as.character(K)

    return(new("CAMObj",PrepResult=PrepResult,
                MGResult=MGResult,
                ASestResult=ASestResult))
}
