#' debCAM: A package for fully unsupervised deconvolution of complex tissues.
#'
#' The core function in this package is \code{\link{CAM}} which achieves fully
#' unsupervised deconvolution on mixture expression profiles.
#' Each step in \code{\link{CAM}} can also be performed separately by
#' \code{\link{CAMPrep}}, \code{\link{CAMMGCluster}}
#' and \code{\link{CAMASest}} in a more flexible workflow.
#' \code{\link{MGstatistic}} can help extract a complete marker list from CAM
#' results. \code{\link{MDL}} can help decide the underlying subpopulation
#' number. With other functions, e.g. \code{\link{AfromMarkers}} and
#' \code{\link{MGstatistic}}, this package can also perform supervised
#' deconvolution based on prior knowledge of molecular markers,
#' subpopulation-specific expression matrix (S) or proportion matrix (A).
#' Semi-supervised deconvolution can be achieved by combining molecular markers
#' from CAM and from prior knowledge to analyze mixture expressions.
#'
#' @references Wang, N., Hoffman, E. P., Chen, L., Chen, L., Zhang, Z.,
#' Liu, C., … Wang, Y. (2016). Mathematical modelling of transcriptional
#' heterogeneity identifies novel markers and subpopulations in complex
#' tissues. Scientific Reports, 6, 18909. http://doi.org/10.1038/srep18909
#'
#' @docType package
#' @name debCAM-package
#' @aliases debCAM
#'
#' @import methods
#' @import rJava
#' @import BiocParallel
#' @import stats
#' @import graphics
#' @importFrom NMF .fcnnls
#' @importFrom corpcor pseudoinverse
NULL
