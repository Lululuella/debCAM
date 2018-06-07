#' Deconvoluted matrix accessors
#'
#' Accessors to proportion matrix and subpopulation-specific expression matrix
#' estimated by CAM.
#' @param x a \code{\link{CAMObj}} object or a \code{\link{CAMASObj}} object.
#' @param ... additional argument list.
#' @return Estimated A matrix or S matrix.
#' @export
#' @rdname AS-accessor
#' @name AS-accessor
#' @aliases Amat
#' @examples
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' rCAM <- CAM(data, K = 3, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)
#' Aest <- Amat(rCAM, 3)
#' Sest <- Smat(rCAM, 3)
#'
#' Aest <- Amat(slot(rCAM, "ASestResult")[[1]])
#' Sest <- Smat(slot(rCAM, "ASestResult")[[1]])
setGeneric(name="Amat", def=function(x, ...) standardGeneric("Amat"))

#' @export
#' @rdname AS-accessor
#' @name AS-accessor
#' @aliases Smat
setGeneric(name="Smat", def=function(x, ...) standardGeneric("Smat"))


#' @param k subpopulation number
#' @param usingPCA If TRUE, A matrix is estimated by transforming
#' dimension-reduced A matrix back to original space. Otherwise, A matrix is
#' directly estimated in original data space. The default is TRUE.
#' @rdname AS-accessor
#' @aliases Amat
setMethod("Amat", signature(x="CAMObj"),
    function(x, k, usingPCA = TRUE) {
        m <- which(names(x@ASestResult) == as.character(k))
        if (usingPCA) return(x@ASestResult[[m]]@Aest)
        else return(x@ASestResult[[m]]@AestO)
    }
)

#' @rdname AS-accessor
#' @aliases Amat
setMethod("Amat", signature(x="CAMASObj"),
    function(x, usingPCA = TRUE) {
        if (usingPCA) return(x@Aest)
        else return(x@AestO)
    }
)

#' @rdname AS-accessor
#' @aliases Smat
setMethod("Smat", signature(x="CAMObj"),
    function(x, k, usingPCA = TRUE) {
        m <- which(names(x@ASestResult) == as.character(k))
        if (usingPCA) return(x@ASestResult[[m]]@Sest)
        else return(x@ASestResult[[m]]@SestO)
    }
)

#' @rdname AS-accessor
#' @aliases Smat
setMethod("Smat", signature(x="CAMASObj"),
    function(x, usingPCA = TRUE) {
        if (usingPCA) return(x@Sest)
        else return(x@SestO)
    }
)

#' Dimension-reduction loading matrix accessor
#'
#' Accessor to Dimension-reduction loading matrix.
#' @param x a \code{\link{CAMObj}} object or a \code{\link{CAMPrepObj}} object
#' @param ... additional argument list.
#' @return The matrix whose rows are loading vectors for dimension reduction.
#' @export
#' @examples
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' rCAM <- CAM(data, K = 3, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)
#' W <- PCAmat(rCAM)
#' W <- PCAmat(slot(rCAM, "PrepResult"))
setGeneric(name="PCAmat", def=function(x, ...) standardGeneric("PCAmat"))

#' @export
#' @rdname PCAmat
setMethod("PCAmat", signature(x="CAMObj"),
    function(x) {
        x@PrepResult@W
    }
)

#' @export
#' @rdname PCAmat
setMethod("PCAmat", signature(x="CAMPrepObj"),
    function(x) {
        x@W
    }
)
