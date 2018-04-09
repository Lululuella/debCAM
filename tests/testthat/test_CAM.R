context("Valid CAM")
test_that("Return NULL for inapplicable K", {
    K <- 3
    S <- matrix(2 ^ rnorm(1000*K, 5, 2.5), ncol = K)
    # A <- matrix(c(0.2,0.4,0.6,0.8,
    #               0.8,0.6,0.4,0.2), ncol = K)
    A <- matrix(c(0.1,0.3,0.5,0.7,0.9,
                  0.2,0.4,0.5,0.1,0.0,
                  0.7,0.3,0.0,0.2,0.1), ncol = K)
    X <- log2(A%*%t(S))  + rnorm(1000*nrow(A), 0 , 0.2)
    X <- t (2 ^ X)

    expect_warning(
        CAMResult <- CAM(X, K = 2:4, dim.rdc = 3,thres.low = 0.3,
                         thres.high = 0.9, cluster.num = 20, MG.num.thres = 5,
                         lof.thres = 0.02, cores = 0)
    )

    expect_true(is.null(CAMResult$MGResult[[3]]))
    expect_true(is.null(CAMResult$ASestResult[[3]]))

    MDLres <- MDL(CAMResult)
    expect_length(MDLres$mdls, 2)
})
