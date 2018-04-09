#' MG cluster detection for CAM
#'
#' This function finds corner clusters as MG clusters (clusters containing marker genes).
#' @param K The candidate subpopulation number.
#' @param PrepResult An object of class "CAMPrepObj" from \code{\link{CAMPrep}} function.
#' @param nComb The number of possible combinations of clusters as corner clusters.
#'     Within these possible combinations ranked by margin errors,
#'     We can further select the best one based on reconstruction errors.
#'     The default is 200.
#' @details This function is used internally by \code{\link{CAM}} function to
#' detect clusters containing marker genes,
#' or used when you want to perfrom CAM step by step.
#'
#' This function provides two solutions. The first is the combination of clusters
#'     yielding the minimum sum margin-of-errors for cluster ceneters. In the second,
#'     nComb possible combinations are selected by ranking sum margin-of-errors
#'     for cluster centers. Then the best one is selected based on
#'     reconstruction errors of all data points in original space.
#' @return An object of class "CAMMGObj" containing the following components:
#' \item{idx}{Two numbers which are two solutions' ranks by margin error.}
#' \item{corner}{The indexes of clusters as detected corners. Each row is a solution.}
#' \item{error}{Two rows. The first row is margin errors for nComb possible combinations.
#'     The second row is reconstruction errors for nComb possible combinations.}
#' @export
#' @examples
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' #preprocess data
#' rPrep <- CAMPrep(data, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)
#'
#' #Marker gene cluster detection with a fixed K = 3
#' rMGC <- CAMMGCluster(3, rPrep)
CAMMGCluster <- function(K, PrepResult, nComb = 200) {
    X <- PrepResult$centers
    if (K > nrow(X) || K > ncol(X)){
        warning("Return NULL for K = ", K, " larger than ", max(dim(X)))
        return(NULL)
    }
    if (K == ncol(X)){
        return(list(idx=c(1,1),
                    corner=matrix(as.integer(colnames(X)),1)[c(1,1),]))
    }

    M <- nrow(X)
    topconv <- cornerSort(X, K, nComb)
    idx <- topconv$idx
    error1 <- topconv$error
    nComb <- length(error1)


    Xall <- PrepResult$Xprep[,!(PrepResult$cluster$cluster %in%
                                    PrepResult$c.outlier)]

    errCalcu <- function (p, idx, X, Xall, W) {
        A <- X[,idx[,p]]
        #scale <- as.vector(coef(nnls(A,matrix(rowSums(W),M,1))))
        scale <- c(.fcnnls(A, matrix(rowSums(W),M,1))$coef)
        scale[scale<1e-10] <- 0.01/(sqrt(colSums(A^2)))[scale<1e-10]
        A<-A%*%diag(scale)
        #sum(apply(Xall,2, function(x) (nnls(A,x))$deviance))
        sum((Xall - A%*%(.fcnnls(A, Xall)$coef))^2)
    }
    error2 <- unlist(lapply(seq_len(nComb), errCalcu,
                            idx, X, Xall, PrepResult$W))

    idx1 <- which.min(error1)
    idx2 <- which.min(error2)
    ind1 <- as.integer(colnames(X[,idx[,idx1]]))
    ind2 <- as.integer(colnames(X[,idx[,idx2]]))


    structure(list(idx=c(idx1,idx2), corner=rbind(ind1,ind2),
                   error=cbind(error1,error2)), class = "CAMMGObj")
}


