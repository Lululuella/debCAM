# CAMTHC

The goal of CAMTHC is to perform comoputational deconvolution for tissue heterogeneity characterization. It provides basic functions for unsupervised deconvolution on mixture expression profiles by CAM and some auxiliary functions to help understand the subpopulation-specific results. CAMTHC also implements functions to perform supervised deconlution based on prior knowledge of molecular markers, S matrix or A matrix. Semi-suprvised deconvolution can also be achieved by combining molecular markers from CAM and from prior knowledge to analyze mixture expressions.

## Installation

You can install the released version of CAMTHC following the instructions CAMTHC's bioconductor page

## Example

This is a basic example which shows you how to decompose a data matrix by CAM:

``` r
## specify the range of underlying subpopulation number
rCAM <- CAM(data, K = 2:5)
```

