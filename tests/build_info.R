library(BayesLinReg)

information <- blm_build_info()
stopifnot(
  identical(names(information), c("eigen_blas", "dense_kernel", "blas")),
  is.logical(information$eigen_blas),
  length(information$eigen_blas) == 1L,
  information$dense_kernel %in% c("Eigen", "external BLAS"),
  is.character(information$blas),
  length(information$blas) == 1L,
  identical(
    information$dense_kernel,
    if (information$eigen_blas) "external BLAS" else "Eigen"
  )
)

expected_eigen_blas <- Sys.getenv(
  "BAYESLINREG_EXPECT_EIGEN_BLAS", unset = ""
)
if (nzchar(expected_eigen_blas)) {
  stopifnot(
    expected_eigen_blas %in% c("yes", "no"),
    identical(information$eigen_blas, expected_eigen_blas == "yes")
  )
}
