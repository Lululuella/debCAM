#' Proportion matrix estimation from marker genes
#'
#' This function estimates proportion matrix (A matrix) from observed mixture
#' expression data based on marker genes.
#' @param data A data set that will be internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     data should be in non-log linear space with non-negative numerical values
#'     (i.e. >= 0).
#'     Missing values are not supported.
#' @param MGlist A list of vectors, each of which contains known markers
#'     and/or CAM-detected markers for one subpopulation.
#' @param scaleRecover If TRUE, scale ambiguity of each column vector
#'     in A matrix is removed based on sum-to-one constraint on each row.
#' @details With the expression levels of subpopulation-specific
#' marker genes, the relative proportions of constituent subpopulations are
#' estimated by spatial median using \code{\link[pcaPP]{l1median}}.
#' Marker genes could be from unsupervised/supervised detection or
#' from literatures.
#' Scale ambiguity is optionally removed based on sum-to-one constraint of rows.
#' @return Return the estiamted proportion matrix (A matrix).
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
    } else if (class(data) == "SummarizedExperiment") {
        data <- assay(data)
    } else if (class(data) == "ExpressionSet") {
        data <- exprs(data)
    } else if (class(data) != "matrix") {
        stop("Only matrix, data frame and SummarizedExperiment object are
             supported for expression data!")
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
        scale.pca[scale.pca<1e-10] <-
            0.01/(sqrt(colSums(A1^2)))[scale.pca<1e-10]
        scale.pca[scale.pca==0] <- 0.0001
        A1 <- A1%*%diag(scale.pca)
        A[,colSums(A)>0] <- A1
        A <- A/rowSums(A)
    }

    A
}
