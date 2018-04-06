#' Candidate combinations as corners
#'
#' Given a set of data points, return possible combinations of data points as corners.
#' These combinations are selected by ranking margin errors.
#' @param X    A matrix of data. Each column is a data point.
#' @param K    The number of corner points.
#' @param nComb    The number of possible combinations of data points as corners.
#' @return A list containing the following components:
#' \item{idx}{A matrix to show the indexes of data points in combinations to construct a convex hull.
#' Each column is one combination.}
#' \item{error}{A vector of margin errors for each combination}
cornerSort <- function(X, K, nComb){
    corner.detect <- .jnew("CornerDetectTopN", .jarray(X, dispatch=TRUE), as.integer(K), as.integer(nComb))
    .jcall(corner.detect, "Z", "search")
    idx <- .jcall(corner.detect, "[[I", "getTopNConv")
    idx <- vapply(idx, .jevalArray, integer(K))
    error <- .jcall(corner.detect, "[D", "getTopNConvErr")
    return(list(idx=idx,error=error))
}
