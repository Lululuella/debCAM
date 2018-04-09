#' A estimation from marker gene list
#'
#' This function estimate proportion matrix based on marker gene list and observed
#' mixture expression data.
#' @param data A data set that will be internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     data should be in non-log linear space with non-negative numerical values (i.e. >= 0).
#'     Missing values are not supported.
#' @param MGlist A list of vectors, each of which contains known markers
#'     and/or CAM-detected markers for one subpopulation.
#' @param scaleRecover If TRUE, scale ambiguity of each column vector
#'     in A matrix is removed based on sum-to-one constraint on each row.
#' @details On the basis of the expression levels of subpopulation-specific marker genes
#'     detected by CAM or from literatures, the relative proportions
#'     of constituent subpopulations are estimated by spatial median using \code{\link[pcaPP]{l1median}}.
#'     Scale ambiguity is removed optionally.
#' @return returns the estiamted proportion matrix.
#' @export
#' @examples
#' #obtain data and marker genes
#' data(ratMix3)
#' S <- ratMix3$S
#' pMGstat <- MGstatistic(S, c("Liver","Brain","Lung"))
#' pMGlist.FC <- lapply(c("Liver","Brain","Lung"), function(x)
#'     rownames(pMGstat)[pMGstat$idx == x & pMGstat$OVE.FC > 10])
#'
#' #estimate A matrix from markers
#' Aest <- AfromMarkers(ratMix3$X, pMGlist.FC)
AfromMarkers <- function(data, MGlist, scaleRecover = TRUE){
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

    Xproj <- t(data/rowSums(data))

    XprojMean <- matrix(0, ncol(data), length(MGlist))
    for(i in seq_along(MGlist)){
        mg <- MGlist[[i]]
        mg <- mg[rowSums(data[mg,]) > 0]
        if(length(mg)==0) XprojMean[,i] <- 0
        else if(length(mg)==1) XprojMean[,i] <- Xproj[,mg]
        else XprojMean[,i] <- pcaPP::l1median(t(Xproj[,mg]))
    }

    A <- XprojMean
    if (scaleRecover == TRUE) {
        A1 <- A[,colSums(A)>0]
        scale.pca <- c(.fcnnls(A1, matrix(1,nrow(A1),1))$coef)
        #scale.pca <- as.vector(coef(nnls(A1,matrix(1,nrow(A1),1))))
        scale.pca[scale.pca<1e-10] <- 0.01/(sqrt(colSums(A1^2)))[scale.pca<1e-10]
        scale.pca[scale.pca==0] <- 0.0001
        A1 <- A1%*%diag(scale.pca)
        A[,colSums(A)>0] <- A1
        A <- A/rowSums(A)
    }

    A
}
