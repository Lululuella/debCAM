#' Sequential Forward Floating Selection (SFFS) to approximate convex hull
#'
#' This function detects the corners of convex hull by greedy search. It will
#' be used to reduce the number of candidate corners and thus reduce the time
#' complexity in the further exhaustive search by \code{\link{cornerSort}}.
#' @param Xall Data in a matrix. Each column is a data point.
#'     The cost is computed based on all data points.
#' @param Aall Candiatie corners in a matrix. Each column is a candidate corner.
#' @param Kmax The target number of corners to be selected.
#' @param deltaK The extra number of corners that need to be searched.
#'     The default is 8 and will be truncated based on all the available
#'     corners.
#'     SFFS runs until a corner set of cardinality (\code{Kmax} + \code{deltak})
#'     is obtained. The set of cardinality Kmax might be improved during
#'     backtracking from extra corners.
#' @details The Sequential Floating Forward selection (SFFS) is one of greedy
#' search methods for feature selection. With sum of margin-of-errors
#' as cost function and candidate corners as features, SFFS is used to select
#' best \code{Kmax} corners to form an approximated hull. The best subset of
#' candidate corners is initialized as the empty set and at each step a new
#' corner is added. After that, the algorithm searches for corner that can be
#' removed from the best subset until the cost function does not decrease.
#' @return A list with length (\code{Kmax} + \code{deltak}). Each component is
#' a vector of the corner indices of the SFFS-selected subset with
#' certain cardinality.
#' @export
#' @examples
#' data <- matrix(c(0.1,0.2,1.0,0.0,0.0,0.5,0.3,
#'                  0.1,0.7,0.0,1.0,0.0,0.5,0.3,
#'                  0.8,0.1,0.0,0.0,1.0,0.0,0.4), nrow =3, byrow = TRUE)
#' rsffs <- sffsHull(data, data, 3)
#' rsffs <- sffsHull(data, data[,1:5], 3)

sffsHull <- function(Xall, Aall, Kmax, deltaK = 8){
    deltaK <- min(deltaK, ncol(Aall) - Kmax)

    costfunc <- function(combset) {
        A <- Aall[,combset,drop=FALSE]
        fit.err <- apply(Xall, 2, function(x) (nnls::nnls(A,x))$deviance)
        sum(fit.err)
    }

    rmworst <- function(combset) {
        fit.err <- unlist(
            lapply(seq_along(combset), function(x) costfunc(combset[-x]))
        )
        list(worstelem = combset[which.min(fit.err)], cost = min(fit.err))
    }

    allset <- seq_len(ncol(Aall))
    sffsset <- rep(list(NA), Kmax + deltaK)
    sffscost <- rep(Inf, Kmax + deltaK)
    convset <-  c()

    while(length(convset) < Kmax + deltaK){
        outset <- setdiff(allset, convset)
        if (length(outset) == 0) break

        #Step 1 (Inclusion)
        fit.err<- unlist(lapply(outset, function(x) costfunc(c(convset,x))))
        newelem <- outset[which.min(fit.err)]
        convset <- c(convset, newelem)
        newcost <- min(fit.err)

        if(length(convset) > 2){
            #Step 2 (Conditional exclusion)
            rmresult <- rmworst(convset)
            if(newelem != rmresult$worstelem){
                while(rmresult$cost < sffscost[length(convset)-1]){
                    convset <- setdiff(convset, rmresult$worstelem)
                    newcost <- rmresult$cost
                    if(length(convset) <= 2) break
                    rmresult <- rmworst(convset)
                }
            }
        }

        sffsset[[length(convset)]] <- convset
        sffscost[length(convset)] <- newcost
    }

    return(sffsset)
}
