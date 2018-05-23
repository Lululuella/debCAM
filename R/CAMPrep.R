#' Data preprocessing for CAM
#'
#' This function perform preprocessing for CAM, including norm-based filtering,
#' dimension deduction, perspective projection, local outlier removal and
#' aggregation of gene expression vectors by clustering.
#' @param data Matrix of mixture expression profiles.
#'     Data frame, SummarizedExperiment or ExpressionSet object will be
#'     internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     Data should be in non-log linear space with non-negative numerical values
#'     (i.e. >= 0). Missing values are not supported.
#'     All-zero rows will be removed internally.
#' @param dim.rdc Reduced data dimension; should be not less than maximum
#'     candidate K.
#' @param thres.low The lower bound of percentage of genes to keep for CAM
#'     with ranked norm. The value should be between 0 and 1.
#'     The default is 0.05.
#' @param thres.high The higer bound of percentage of genes to keep for CAM
#'     with ranked norm. The value should be between 0 and 1.
#'     The default is 0.95.
#' @param cluster.method The method to do clustering. The default "K-Means" will
#'     use \code{\link{kmeans}} fucntion. The alternative "apcluster" will use
#'     \code{\link[apcluster]{apclusterK-methods}}.
#' @param cluster.num The number of clusters; should be much larger than K.
#'     The default is 50.
#' @param MG.num.thres The clusters with the gene number smaller than
#'     MG.num.thres will be treated as outliers. The default is 20.
#' @param lof.thres Remove local outlier using \code{\link[DMwR]{lofactor}}
#'     function. MG.num.thres is used as the number of neighbours in the
#'     calculation of the local outlier factors.
#'     The default value 0.02 will remove top 2\% local outliers.
#'     Zero value will disable lof.
#' @details This function is used internally by \code{\link{CAM}} function to
#' preprocess data, or used when you want to perform CAM step by step.
#'
#' Low/high-expressed genes are filtered by their L2-norm ranks.
#' Dimension reduction is slightly different from PCA.
#' The first loading vector is forced to be c(1,1,...,1) with unit norm
#' normaliztion. The remaining are eignvectors from PCA in the space orthogonal
#' to the first vector.
#' Perspective projection is to project dimension-reduced gene expression
#' vectors to the hyperplane orthogonal to c(1,0,...,0), i.e., the first axis
#' in the new coordinate system.
#' local outlier removal is optional to exclude outliers in simplex formed
#' after perspective projection.
#' Finallly, gene expression vectors are aggregated by clustering
#' to further reduce the impact of noise/outler and help improve the efficiency
#' of simplex corner detction.
#' @return An object of class "\code{\link{CAMPrepObj}}" containing the
#' following components:
#' \item{Valid}{logical vector to indicate the genes left after filtering.}
#' \item{Xprep}{Preprocessed data matrix.}
#' \item{Xproj}{Preprocessed data matrix after perspective projection.}
#' \item{W}{The matrix whose rows are loading vectors.}
#' \item{cluster}{cluster results including two vectors.
#'     The first indicates the cluster to which each gene is allocated.
#'     The second is the number of genes in each cluster.}
#' \item{c.outlier}{The clusters with the gene number smaller than
#'     MG.num.thres.}
#' \item{centers}{The centers of candidate corner clusters (candidate clusters
#'     containing marker genes).}
#' @export
#' @examples
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' #set seed to generate reproduciable results
#' set.seed(111)
#'
#' #preprocess data
#' rPrep <- CAMPrep(data, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)
CAMPrep <- function(data, dim.rdc = 10, thres.low = 0.05, thres.high = 0.95,
                    cluster.method = c('K-Means', 'apcluster'),
                    cluster.num = 50, MG.num.thres = 20, lof.thres = 0){
    if (is(data, "data.frame")) {
        data <- as.matrix(data)
    } else if (is(data, "SummarizedExperiment")) {
        data <- assay(data)
    } else if (is(data, "ExpressionSet")) {
        data <- exprs(data)
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
    if (thres.low >= thres.high) {
        stop("thres.low must be smaller than thres.high!")
    }
    if (is.null(rownames(data))) {
        rownames(data) <- seq_len(nrow(data))
    }

    data <- data[rowSums(data) > 0,]
    X <- t(data)

    ######### data preprocessing ############
    sigNorm <- apply(X, 2, function(x) norm(matrix(x),"F") )
    Valid <- sigNorm > quantile(sigNorm, thres.low) &
        sigNorm <= quantile(sigNorm, thres.high)
    X <- X[,Valid]

    Lorg <- dim(X)[1]
    dataSize<- dim(X)[2]
    Xorg <- X

    ################ PCA dimension deduction #################
    L <- min(dim.rdc, Lorg)
    Xmean <- rowMeans(Xorg)
    Xmeanrm <- Xorg - Xmean
    P1 <- rep(1,Lorg) / sqrt(Lorg)
    C1 <- P1 %*% Xmeanrm
    Xmeanrm <- Xmeanrm - matrix(rep(C1, 1, each=Lorg), Lorg) *
        matrix(P1, Lorg, dataSize)
    r <- eigen(Xmeanrm %*% t(Xmeanrm))
    Ppca <- t(r$vectors[,seq_len(L-1)])
    X <- rbind(P1,Ppca) %*% Xorg
    weightMatrix <- rbind(P1, Ppca)

    ################ perspective projection #################
    Xcenter <- c(1, rep(0, L-1))
    denom <- Xcenter %*% X
    denom <- denom[rep(1,L),]
    Xproj <- X / denom

    ################ local outlier removal #################
    if(lof.thres > 0){
        lof.factor <- DMwR::lofactor(t(Xproj), MG.num.thres)
        lof.outlier <- lof.factor > quantile(lof.factor, 1-lof.thres)
        Valid[Valid==TRUE][lof.outlier]<- FALSE
        X <- X[,!lof.outlier]
        Xproj <- Xproj[,!lof.outlier]
    }

    ################ clustering #################
    if (length(cluster.method) > 1) cluster.method <- cluster.method[1]
    if (cluster.method == 'K-Means') {
        clusterRes <- kmeans(t(Xproj), cluster.num, iter.max=100)
        #repeat K-means and use the best one
        for (i in seq_len(50)){
            tmp <- kmeans(t(Xproj), cluster.num, iter.max=100)
            if (clusterRes$tot.withinss > tmp$tot.withinss){
                clusterRes <- tmp
            }
        }
        cluster <- list(cluster = clusterRes$cluster, size = clusterRes$size)
    } else {
        if (cluster.method != 'apcluster') {
            stop("Only K-Means and apcluster are supported.")
        }
        clusterRes <- apcluster::apclusterK(apcluster::negDistMat(r=2),
                                            t(Xproj),  K=cluster.num)
        clusterSize <- unlist(lapply(clusterRes@clusters, length))
        clusterIdx <- rep(0, length(clusterSize))
        for(i in seq_along(clusterSize)) {
            clusterIdx[clusterRes@clusters[[i]]] <- i
        }
        cluster <- list(cluster = clusterIdx, size = clusterSize)
    }


    cluster.valid <- which(cluster$size >= MG.num.thres)
    c.outlier <- which(cluster$size < MG.num.thres)
    message('outlier cluster number: ',sum(cluster$size < MG.num.thres),"\n")

    medcenters <- vapply(cluster.valid, function(x)
        pcaPP::l1median(t(Xproj[,cluster$cluster==x])), numeric(L))
    colnames(medcenters) <- cluster.valid

    ################ quickhull #################
    J <- length(cluster.valid)
    corner <- seq_len(J)
    convex <- geometry::convhulln(rbind(t(medcenters),0), options = "QbB")
    corner <- unique(c(convex))
    corner <- corner[-which(corner == (J+1))]  # throw away the origin point
    message('convex hull cluster number: ',length(corner),"\n")


    return(new("CAMPrepObj", Valid=Valid, Xprep=X, Xproj=Xproj, W=weightMatrix,
                    cluster=cluster, c.outlier=c.outlier,
                    centers=medcenters[,corner]))
}
