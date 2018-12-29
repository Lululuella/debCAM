#' Re-estimate A, S matrix
#'
#' This function re-estimates proportion and expression matrix iteratively
#' by Alternating Least Square (ALS) method. The initial values are from
#' markers or known proportion matrix or known expression matrix,
#' @param data A data set that will be internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     data should be in non-log linear space with non-negative numerical values
#'     (i.e. >= 0).
#'     Missing values are not supported.
#'     All-zero rows will be removed internally.
#' @param MGlist A list of vectors, each of which contains CAM-detected markers
#'     and/or prior markers for one subpopulation.
#' @param A Initial proportion matrix. If NULL, it will be estimated from
#'     initial expression matrix. If initial expression matrix is also NULL, it
#'     will be estimated from \code{MGlist} using \code{\link{AfromMarkers}}.
#' @param S Initial expression matrix. If NULL, it will be estimated from
#'     initial proportion matrix.
#' @param generalNMF If TRUE, the decomposed proportion matrix has no sum-to-one
#'     constraint for each row. Without this constraint, the scale ambiguity of
#'     each column vector in proportion matrix will not be removed.
#'     The default is FALSE.
#' @param maxIter maximum number of iterations for Alternating Least
#'     Square (ALS) method. The default is 2. If zero, ALS is not applied.
#' @param methy Should be TRUE when dealing with methylation data, whose
#'     expression levels are confined between 0 and 1.
#' @details If only markers are provided, they are used to estimate initial
#' proportion matrix and then expression matrix. If proportion matrix or
#' expression matrix is provided, it will be treated as initial matrix to
#' estimate the other one. Then Alternating Least Square (ALS) method is
#' applied to estimate two matrix alternatively. Note only markers'
#' squared errors will be counted in ALS, which facilitates (1) faster
#' computational running time and (2) a greater signal-to-noise ratio owing to
#' markers' discriminatory power.
#'
#' This function can be used to refine CAM estimation or perform supervised
#' deconvolution. Note that allowing too many iterations may bring the risk of
#' a significant deviation from initial values.
#' @return A list containing the following components:
#' \item{Aest}{Proportion matrix after re-estimation and possible refinement.}
#' \item{Sest}{expression matrix after re-estimation and possible refinement.}
#' \item{mse}{Mean squared error (i.e. mean of reconstruction errors)
#' for input markers}
#' @export
#' @examples
#' #obtain data and run CAM
#' data(ratMix3)
#' data <- ratMix3$X
#' rCAM <- CAM(data, K = 3, dim.rdc= 3, thres.low = 0.30, thres.high = 0.95)
#' #obtain marker genes detected by CAM for estimating A
#' MGlist <- MGsforA(rCAM, K = 3)
#'
#' #Re-estimation based on marker list
#' rre <- redoASest(data, MGlist, maxIter = 10)
#' Aest <- rre$Aest #re-estimated A matrix
#' Sest <- rre$Sest #re-estimated S matrix
#'
#' #Re-estimation with initial A matrix
#' rre <- redoASest(data, MGlist, A=ratMix3$A, maxIter = 10)
#'
#' #Re-estimation with initial S matrix
#' rre <- redoASest(data, MGlist, S=ratMix3$S, maxIter = 10)
redoASest <- function(data, MGlist, A = NULL, S = NULL,
                    generalNMF = FALSE, maxIter = 2, methy = FALSE) {
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

    data <- data[rowSums(data) > 0,]
    MGlist <- lapply(MGlist, function(x) intersect(x,rownames(data)))
    X <- data[unlist(MGlist),]

    if (is.null(A) && is.null(S)) {
        A0 <- AfromMarkers(data, MGlist, !generalNMF)
        estimateA <- FALSE
    } else if (is.null(S)) {
        A0 <- A
        estimateA <- FALSE
    } else {
        S0 <- S[unlist(MGlist),]
        estimateA <- TRUE
    }

    err <- c()
    iter0 <- -1
    while(iter0 < maxIter*2) {
        if (estimateA) {
            if (!generalNMF) {
                A0 <- t(apply(rbind(1e+5,X), 2, function(x)
                    coef(nnls::nnls(rbind(1e+5,S0), x))))
            } else {
                A0 <- t(apply(X, 2, function(x) coef(nnls::nnls(S0, x))))
            }
            estimateA <- FALSE
        } else {
            S0 <- t(.fcnnls(A0, t(X))$coef)
            if (methy) {
                S0[S0>1] <- 1
            }
            estimateA <- TRUE
        }

        err0 <- mean((X - S0%*%t(A0))^2)
        iter0 <- iter0 + 1
        err <- c(err, err0)
    }

    S0 <- t(.fcnnls(A0, t(data))$coef)
    if (methy) {
        S0[S0>1] <- 1
    }

    return(list(Aest=A0, Sest=S0, mse=err))
}



