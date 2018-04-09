context("corner_sort")

test_that("corners are detected correctly", {
    X <- matrix(c(0.1,0.2,1.0,0.0,0.0,0.5,0.3,
                  0.1,0.7,0.0,1.0,0.0,0.5,0.3,
                  0.8,0.1,0.0,0.0,1.0,0.0,0.4), nrow =3, byrow = TRUE)
    topconv <- cornerSort(X, 3, 1)
    expect_setequal(c(topconv$idx), c(3,4,5))
})
