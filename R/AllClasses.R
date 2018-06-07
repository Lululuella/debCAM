#' Class "CAMPrepObj"
#'
#' An S4 class for storing data preprocessing results.
#'
#' @slot Valid logical vector to indicate the genes left after filtering.
#' @slot Xprep Preprocessed data matrix.
#' @slot Xproj Preprocessed data matrix after perspective projection.
#' @slot W The matrix whose rows are loading vectors.
#' @slot cluster cluster results including two vectors.
#'     The first indicates the cluster to which each gene is allocated.
#'     The second is the number of genes in each cluster.
#' @slot c.outlier The clusters with the gene number smaller than
#'     MG.num.thres.
#' @slot centers The centers of candidate corner clusters (candidate clusters
#'     containing marker genes).
CAMPrepObj <- setClass("CAMPrepObj",
    slots = list(
        Valid = "vector",
        Xprep = "matrix",
        Xproj = "matrix",
        W = "matrix",
        cluster = "list",
        c.outlier = "vector",
        centers = "matrix")
)


#' Class "CAMMGObj"
#'
#' An S4 class for storing marker gene detection results.
#'
#' @slot idx Two numbers which are two solutions' ranks by sum of
#'     margin-of-error.
#' @slot corner The indexes of clusters as detected corners. Each row is a
#'     solution.
#' @slot error Two rows. The first row is sum of margin-of-errors for nComb
#'     possible combinations. The second row is reconstruction errors for
#'     nComb possible combinations.
CAMMGObj <- setClass("CAMMGObj",
    slots = list(
        idx = "vector",
        corner = "matrix",
        error = "matrix")
)


#' Class "CAMASObj"
#'
#' An S4 class for storing estimated proportions, subpopulation-specific
#' expressions and mdl values.
#' The mdl values are calculated in three approaches:
#' (1) based on data and A matrix in dimension-reduced space; (2) based on
#' original data with A matrix estimated by transforming dimension-reduced
#' A matrix back to original space;
#' (3) based on original data with A directly estimated in original space.
#' A and S matrix in original space estimated from the latter two approaches are
#' returned. mdl is the sum of two terms: code length of data under the model
#' and code length of model. Both mdl value and the first term (code length
#' of data) will be returned.
#'
#' @slot Aest Estimated proportion matrix from Approach 2.
#' @slot Sest Estimated subpopulation-specific expression matrix from
#'     Approach 2.
#' @slot Aest.proj Estimated proportion matrix from Approach 2, before
#'     removing scale ambiguity.
#' @slot Ascale The estimated scales to remove scale ambiguity
#'     of each column vector in Aest. Sum-to-one constraint on each row of
#'     Aest is used for scale estimation.
#' @slot AestO Estimated proportion matrix from Approach 3.
#' @slot SestO Estimated subpopulation-specific expression matrix from
#'     Approach 3.
#' @slot AestO.proj Estimated proportion matrix from Approach 3, before
#'     removing scale ambiguity.
#' @slot AscaleO The estimated scales to remove scale ambiguity
#'     of each column vector in AestO. Sum-to-one constraint on each row of
#'     AestO is used for scale estimation.
#' @slot datalength Three values for code length of data. The first is
#'     calculated based on dimension-reduced data. The second and third are
#'     based on the original data.
#' @slot mdl Three mdl values. The first is calculated based on
#'     dimension-reduced data. The second and third are based on the original
#'     data.
CAMASObj <- setClass("CAMASObj",
    slots = list(
        Aest = "matrix",
        Sest = "matrix",
        Aest.proj = "matrix",
        Ascale = "vector",
        AestO = "matrix",
        SestO = "matrix",
        AestO.proj = "matrix",
        AscaleO = "vector",
        datalength = "vector",
        mdl = "vector")
)

#' Class "CAMObj"
#'
#' An S4 class for storing results of CAM.
#'
#' @slot PrepResult An object of class "\code{\link{CAMPrepObj}}" storing data
#'     preprocessing results from \code{\link{CAMPrep}} function.
#' @slot MGResult A list of "\code{\link{CAMMGObj}}" objects
#'     storing marker gene detection
#'     results from \code{\link{CAMMGCluster}} function for each candidate
#'     subpopulation number.
#' @slot ASestResult A list of "\code{\link{CAMASObj}}" objects storing
#'     estimated proportions, subpopulation-specific expressions and mdl values
#'     from \code{\link{CAMASest}} function for each candidate
#'     subpopulation number.
CAMObj <- setClass("CAMObj",
    slots = list(
        PrepResult = "CAMPrepObj",
        MGResult = "list",
        ASestResult = "list")
)

#' Class "MDLObj"
#'
#' An S4 class for storing mdl values.
#'
#' @slot K The candidate subpopulation numbers.
#' @slot datalengths For each model with a certain subpopulation number,
#' code length of data under the model.
#' @slot mdls mdl value for each model with a certain subpopulation number.
MDLObj <- setClass("MDLObj",
    slots = list(
        K = "vector",
        datalengths = "vector",
        mdls = "vector")
)
