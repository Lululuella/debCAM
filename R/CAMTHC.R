#' CAMTHC: A package for Tissue Heterogeneity Characterization.
#'
#' The CAMTHC package provides basic functions to perform unsupervised
#' deconvolution on mixture expression profiles by CAM and some auxiliary
#' functions to help understand the subpopulation-specific results.
#' CAMTHC also implements functions to perform supervised deconlution based on
#' prior knowledge of molecular markers, S matrix or A matrix. Semi-suprvised
#' deconvolution can also be achieved by combining molecular markers from CAM
#' and from prior knowledge to analyze mixture expressions.
#'
#' @docType package
#' @name CAMTHC-package
#' @aliases CAMTHC
#'
#' @import rJava
#' @import BiocParallel
#' @import stats
#' @import graphics
#' @importFrom NMF .fcnnls
#' @importFrom corpcor pseudoinverse
NULL
