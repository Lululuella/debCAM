#' Statistics for ranking marker genes
#'
#' This function computes One-Versus-Everyone Fold Change (OVE-FC)
#' from subpopulation-specific expression profiles.
#' Bootstrappping is optional.
#' @param data A data set that will be internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     Data should be in non-log linear space with non-negative numerical values
#'     (i.e. >= 0). Missing values are not supported.
#' @param A When data are mixture expression profiles,
#'     A is estimated proportion matix or prior proportion matrix.
#'     When data are pure expression profiles, A is a phenotype vector to
#'     indicate which subpopulation each sample belongs to.
#' @param boot.alpha Alpha for bootstrapped OVE-FC confidence interval.
#'     The default is 0.05.
#' @param nboot The number of boots.
#' @param cores The number of system cores for parallel computing.
#'     If not provided, the default back-end is used.
#' @param seed For reproducibility, the seed of the random number generator for
#'     boot sampling.
#' @details This function calculates OVE-FC and bootstrapped OVE-FC which can be
#' used to identify markers from all genes.
#' @return A data frame containing the following components:
#' \item{idx}{Numbers or phenotypes indicating which subpopulation each gene
#' could be a marker for. If A is a proportion matrix without column name,
#' numbers are returned. Otherwize, phenotypes.}
#' \item{OVE.FC}{One-versus-Everyone fold change (OVE-FC)}
#' \item{OVE.FC.alpha}{lower confidence bound of bootstrapped OVE-FC at alpha
#'     level.}
#' @export
#' @examples
#' #data are mixture expression profiles, A is proportion matrix
#' data(ratMix3)
#' MGstat <- MGstatistic(ratMix3$X, ratMix3$A)
#' \donttest{
#' MGstat <- MGstatistic(ratMix3$X, ratMix3$A, boot.alpha = 0.05) #enable boot
#' }
#'
#' #data are pure expression profiles without replicates
#' MGstat <- MGstatistic(ratMix3$S) #boot is not applicable
#' \donttest{
#' #data are pure expression profiles with phenotypes
#' S <- matrix(rgamma(3000,0.1,0.1), 1000, 3)
#' S <- S[, c(1,1,1,2,2,2,3,3,3,3)] + rnorm(1000*10, 0, 0.5)
#' MGstat <- MGstatistic(S, c(1,1,1,2,2,2,3,3,3,3), boot.alpha = 0.05)
#' }
MGstatistic <- function(data, A = NULL, boot.alpha = NULL, nboot = 1000,
                        cores = NULL, seed = NULL) {
    if (class(data) == "data.frame") {
        data <- as.matrix(data)
    } else if (class(data) == "SummarizedExperiment") {
        data <- assay(data)
    } else if (class(data) == "ExpressionSet") {
        data <- exprs(data)
    } else if (class(data) != "matrix") {
        stop("Only matrix, data frame, SummarizedExperiment and ExpressionSet
            object are supported for expression data!")
    }
    if (sum(is.na(data)) > 0) {
        stop("Data with missing values are not supported!")
    }

    M <- ncol(data)

    if (is.null(A)) {
        A <- diag(M)
        if (!is.null(boot.alpha)) {
            warning("A is missing so that each sample is treated as a single
                    subpopulation and Bootstrapping is not applicable!")
            boot.alpha <- NULL
        }
    } else if (class(A) != "matrix") {
        if (M != length(A)) {
            stop("Sample size in data matrix must be the same with the length
                of A vector!")
        }

        K <- length(unique(A))
        label <- factor(A)
        A <- diag(K)[as.numeric(label),]
        colnames(A) <- levels(label)

        if (!is.null(boot.alpha)) {
            if (!is.null(seed)) {
                set.seed(seed)
            }
            withinGroupSampleBoot <- function(group) {
                sidx <- seq_along(group)
                for (g in unique(group)) {
                    sidx[group == g] <-
                        sample(which(group == g), replace = TRUE)
                }
                sidx
            }
            sampleId <- vapply(seq_len(nboot), function(x)
                withinGroupSampleBoot(label), integer(M))
        }
    } else if (M != nrow(A)) {
        stop("Sample size in data matrix and A matrix should be the same!")
    } else {
        if (!is.null(boot.alpha)) {
            if (!is.null(seed)) {
                set.seed(seed)
            }
            sampleId <- vapply(seq_len(nboot), function(x)
                sample(M, replace = TRUE), integer(M))
        }
    }

    data <- t(data)
    #S <- apply(data, 1, function(x) coef(nnls(A, x)))
    S <- .fcnnls(A, data)$coef
    idx <-  apply(S, 2, which.max)
    OVE.FC <- apply(S, 2, function(x) max(x)/max(x[-which.max(x)]))
    OVE.FC[is.na(OVE.FC)] <- 1

    if(!is.null(boot.alpha)){
        if (boot.alpha > 1 || boot.alpha < 0) {
            stop("boot.alpha should be in range (0,1)")
        }

        if (is.null(cores) || cores > 0) {
            registered()
        }
        if (is.null(cores)) {
            param <- bpparam()
        } else if (cores > 0) {
            param <- SnowParam(workers = cores, type = "SOCK")
        }

        bootFC <- function(p, sampleId, data, A, idx) {
            #S.boot <- apply(data[,sampleId[,p]], 1, function(x)
            #coef(nnls(A[sampleId[,p],], x)))
            S.boot <- tryCatch(.fcnnls(A[sampleId[,p],],
                data[sampleId[,p],])$coef, error = function(e) e)
            if (any(class(S.boot) == "error")) {
                return(rep(NA, ncol(data)))
            }
            OVE.FC.boot <- unlist(lapply(seq_len(ncol(S.boot)), function(x)
                S.boot[idx[x],x]/max(S.boot[-idx[x],x])))
            OVE.FC.boot[is.na(OVE.FC.boot)] <- 1
            return(OVE.FC.boot)
        }
        if (is.null(cores) || cores > 0) {
            S.boots <- do.call("cbind", bplapply(seq_len(nboot),
                bootFC, sampleId, data, A, idx, BPPARAM = param))
        } else {
            S.boots <- do.call("cbind", lapply(seq_len(nboot),
                bootFC, sampleId, data, A, idx))
        }
        OVE.FC.alpha <- apply(S.boots, 1, quantile, boot.alpha, TRUE)
    }


    if (!is.null(colnames(A))) {
        idx <- colnames(A)[idx]
    }

    if (is.null(boot.alpha)) {
        return(data.frame(idx = idx, OVE.FC = OVE.FC))
    } else {
        return(data.frame(idx = idx, OVE.FC = OVE.FC,
                            OVE.FC.alpha = OVE.FC.alpha))
    }
}
