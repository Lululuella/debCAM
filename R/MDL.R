#' Minimum Description Length
#'
#' This function obtains minimum description length (mdl) values for each
#' candidate subpopulation number.
#' @param CAMResult Result from  \code{\link{CAM}} function.
#' @param mdl.method Approach to calculate mdl values; should be 1, 2, or 3.
#' The default is 3.
#' @details This function extracts minimum description length (mdl) values from
#' the result of \code{\link{CAM}} function, which contains mdl values form
#' three approaches for each candidate subpopulation number.
#' For more details about three approaches, refer to \code{\link{CAMASest}}.
#'
#' mdl is code length of data under the model plus
#' code length of model. Both mdl and the first term about data are returned.
#' @return An object of class "MDLObj" containing the following components:
#' \item{K}{The candidate subpopulation numbers.}
#' \item{datalengths}{For each model with a certain subpopulation number,
#' code length of data under the model.}
#' \item{mdls}{mdl values for each model with a certain subpopulation number.}
#' @export
#' @examples
#' \donttest{
#' #obtain data
#' data(ratMix3)
#' data <- ratMix3$X
#'
#' #Analysis by CAM
#' rCAM <- CAM(data, K = 2:5, thres.low = 0.30, thres.high = 0.95)
#'
#' #extract mdl values
#' MDL(rCAM)
#' MDL(rCAM, 1)
#' MDL(rCAM, 2)
#'
#' #plot MDL curves
#' plot(MDL(rCAM))
#' plot(MDL(rCAM), data.term = TRUE) #with data length curve
#' }
MDL <- function(CAMResult, mdl.method = 3) {
    valid <- unlist(lapply(CAMResult$ASestResult, function(x) !is.null(x)))
    valid <- which(valid)
    K <- as.integer(names(CAMResult$ASestResult))[valid]
    datalengths <- unlist(lapply(valid, function(x)
        CAMResult$ASestResult[[x]]$datalength[mdl.method]))
    mdls <- unlist(lapply(valid, function(x)
        CAMResult$ASestResult[[x]]$mdl[mdl.method]))
    structure(list(K = K, datalengths = datalengths, mdls=mdls),
                    class = "MDLObj")
}

#' @param x An object of class "MDLObj" from \code{\link{MDL}}..
#' @param data.term If ture, plot data term (code lenght of data under model).
#' @param ... All other arguments are passed to the plotting command.
#' @export
#' @rdname MDL
plot.MDLObj <- function(x, data.term = FALSE, ...){
    if (data.term) {
        plot(x$K, x$mdls, xlab = 'number of sources', ylab = '', type='l',
            ylim = range(c(x$datalengths, x$mdls)), col= 'blue', xaxt = "n",
            ...)
        points(x$K, x$datalengths, type = 'l', lty = 2, col = 'red')
        axis(1, at = min(x$K) : max(x$K))
        legend("topright", cex=1.5, inset=.01, c("MDL","data term"),
            lty=c(1,2), col = c('blue', 'red'))
    } else {
        plot(x$K, x$mdl, xlab = 'number of sources', ylab = '',
            main = 'MDL Curve', col='blue', type='l', xaxt = "n", ...)
        axis(1, at = min(x$K) : max(x$K))
    }
}
