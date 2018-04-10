#' The plot of scatter simplex
#'
#' This function shows simplex
#' @param data A data set that will be internally coerced into a matrix.
#'     Each row is a gene and each column is a sample.
#'     data should be in non-log linear space with non-negative numerical values
#'     (i.e. >= 0). Missing values are not supported.
#' @param A Prior/Estimated proporiton matrix.
#' @param MGlist A list of vectors, each of which contains known markers
#'     and/or CAM-detected markers for one subpopulation.
#' @param corner.order The order to show simplex corners counterclockwise.
#' @param data.col The color for data points. The default is "gray".
#' @param corner.col The color for corner points. The default is "red".
#' @param ... All other arguments are passed to the plotting command.
#' @details This function can show the scatter simplex and deteced marker genes
#' in a 2D plot. The corners in the high-dimenisonal simplex will still locate
#' at extreme points of low-dimensional simplex. These corners will follow the
#' order set by \code{corner.order} to display in the plot counterclockwise.
#' @return NULL
#' @export
#' @examples
#' #obtain data, A matrix, marker genes
#' data(ratMix3)
#' data <- ratMix3$X
#' A <- ratMix3$A
#' pMGstat <- MGstatistic(ratMix3$S, c("Liver","Brain","Lung"))
#' pMGlist.FC <- lapply(c("Liver","Brain","Lung"), function(x)
#'     rownames(pMGstat)[pMGstat$idx == x & pMGstat$OVE.FC > 10])
#'
#' #plot simplex for data
#' simplexplot(data, A)
#' simplexplot(data, A, MGlist = pMGlist.FC) #Color marker genes in simplex plot
#'
#' #set differnt corner order and colors
#' simplexplot(data, A, MGlist = pMGlist.FC, corner.order = c(2,1,3),
#'             data.col = "blue", corner.col = c("red","orange","green"))
simplexplot <- function(data, A, MGlist = NULL, corner.order = NULL,
                        data.col = 'gray', corner.col = 'red', ...){
    if (class(data) == "data.frame") {
        data <- as.matrix(data)
    } else if (class(data) != "matrix") {
        stop("Only matrix and data frame are supported for expression data!")
    }
    if (sum(data<0) > 0) {
        stop("Only non-negative data are supported!")
    }
    if (nrow(A) != ncol(data)) {
        stop("The number of samples in data and in A matrix
            should be the same!")
    }
    K <- ncol(A)
    if (is.null(corner.order)) {
        corner.order <- seq_len(K)
    } else if (sum(!(corner.order %in% seq_len(K))) > 0) {
        stop("corner.order must be a permutation of 1 ~", K)
    }

    if(length(corner.col) == 1) {
        corner.col <- rep(corner.col, K)
    } else if (length(corner.col) != K) {
        stop("corner.col needs", K, "colors")
    }

    Xproj <- t(data / rowSums(data))
    A <- A / rep(colSums(A), 1, each = nrow(A))
    PS <- t(matrix(c(cos((seq_len(K) - 1) * 2 * pi / K),
                    sin((seq_len(K) - 1) * 2 * pi / K)), K))
    tmp <- PS %*% pseudoinverse(A[,corner.order])
    tmp[1,] <- tmp[1,] / c(sqrt(tmp[1,] %*% tmp[1,]))
    tmp[2,] <- tmp[2,] - c(tmp[2,] %*% tmp[1,]) * tmp[1,]
    tmp[2,] <- tmp[2,] / c(sqrt(tmp[2,] %*% tmp[2,]))
    Xp <- tmp %*% Xproj

    plot(Xp[1,], Xp[2,], col = data.col, cex = 0.8, pch = 1,
        xlab = NA, ylab = NA, asp = 1, ...)

    for(i in seq_along(MGlist)){
        points(Xp[1,MGlist[[i]]], Xp[2,MGlist[[i]]],
            cex=1.2, col=corner.col[i], pch = i-1 )
    }
}
