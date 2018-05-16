#' Marker genes detected by CAM for estimating A
#'
#' This function returns marker genes detected by CAM for estimating A.
#' @param CAMResult Result from \code{\link{CAM}}.
#' @param K The candidate subpopulation number.
#' @param PrepResult An object of class "CAMPrepObj" from \code{\link{CAMPrep}}.
#' @param MGResult An object of class "CAMMGObj" from
#'     \code{\link{CAMMGCluster}}.
#' @param corner.strategy The method to detect corner clusters.
#'     1: minimum sum of margin-of-errors; 2: minimum sum of reconstruction
#'     errors. The default is 2.
#' @details This function needs to specify CAMResult and K, or PrepResult and
#' MGResult. The returned marker genes are those used by CAM for estimating A.
#' To obtain a more complete marker gene list, please refer to
#' \code{\link{MGstatistic}}.
#' @return A list of vectors, each of which contains marker genes for one
#' subpopulation.
#' @export
#' @examples
#' #obtain data and run CAM
#' data(ratMix3)
#' data <- ratMix3$X
#' rCAM <- CAM(data, K = 3, dim.rdc= 3, thres.low = 0.30, thres.high = 0.95)
#' #obtain marker genes detected by CAM for estimating A
#' MGlist <- MGsforA(rCAM, K = 3)
#'
#' #obtain data and run CAM step by step
#' rPrep <- CAMPrep(data, dim.rdc= 3, thres.low = 0.30, thres.high = 0.95)
#' rMGC <- CAMMGCluster(3, rPrep)
#' #obtain marker genes detected by CAM for estimating A
#' MGlist <- MGsforA(PrepResult = rPrep, MGResult = rMGC)
MGsforA <- function(CAMResult = NULL, K = NULL,
                    PrepResult = NULL, MGResult = NULL, corner.strategy = 2) {
    if (is.null(PrepResult) || is.null(MGResult)) {
        if (is.null(CAMResult)) {
            stop("Please provide CAMResult or PrepResult and MGResult!")
        }
        if (is.null(K)) {
            stop("Please provide K!")
        }

        PrepResult <- CAMResult@PrepResult
        MGResult <- CAMResult@MGResult[[which(names(CAMResult@MGResult) ==
                                        as.character(K))]]
    }

    MGlist <- lapply(as.character(MGResult@corner[corner.strategy,]),
        function(x) colnames(PrepResult@Xproj)[PrepResult@cluster$cluster == x])
    MGlist
}



