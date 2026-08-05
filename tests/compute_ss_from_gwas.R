library(BayesLinReg)

set.seed(901)
n <- 80L
X <- scale(matrix(rnorm(n * 5L), nrow = n), center = TRUE, scale = FALSE)
y <- drop(scale(
  X %*% c(0.5, -0.25, 0, 0.75, 0) + rnorm(n),
  center = TRUE, scale = FALSE
))
predictor_names <- paste0("rs", seq_len(ncol(X)))
colnames(X) <- predictor_names

XtX <- crossprod(X)
Xty <- drop(crossprod(X, y))
yty <- sum(y^2)
diagonal <- diag(XtX)
beta <- Xty / diagonal
se <- sqrt((yty - beta^2 * diagonal) / ((n - 2) * diagonal))
names(beta) <- names(se) <- predictor_names
LD <- stats::cov2cor(XtX)

# With in-sample LD and response variance, the SuSiE-RSS conversion recovers
# the original centered cross-products and coefficient scale.
original <- compute_ss_from_gwas(
  beta, se, LD, n, response_var = stats::var(y)
)
stopifnot(
  identical(names(original), c(
    "XtX", "n", "Xty", "yty", "X_means", "y_mean",
    "reference_response_var"
  )),
  isTRUE(all.equal(original$XtX, XtX, tolerance = 1e-10)),
  isTRUE(all.equal(unname(original$Xty), unname(Xty), tolerance = 1e-10)),
  isTRUE(all.equal(original$yty, yty, tolerance = 1e-10)),
  identical(original$X_means, stats::setNames(numeric(ncol(X)),
                                               predictor_names)),
  identical(original$y_mean, 0),
  identical(original$reference_response_var, stats::var(y))
)

# Without response variance, the same inputs recover cross-products for
# sample-standardized predictors and response.
standardized <- compute_ss_from_gwas(beta, se, LD, n)
standardized_X <- scale(X)
standardized_y <- drop(scale(y))
stopifnot(
  isTRUE(all.equal(
    standardized$XtX, crossprod(standardized_X), tolerance = 1e-10
  )),
  isTRUE(all.equal(
    unname(standardized$Xty),
    unname(drop(crossprod(standardized_X, standardized_y))),
    tolerance = 1e-10
  )),
  isTRUE(all.equal(standardized$yty, n - 1, tolerance = 1e-12)),
  identical(standardized$reference_response_var, 1)
)

# Sparse LD remains sparse and gives the same sufficient statistics.
sparse <- compute_ss_from_gwas(
  beta, se, Matrix::Matrix(LD, sparse = TRUE), n,
  response_var = stats::var(y)
)
stopifnot(
  inherits(sparse$XtX, "sparseMatrix"),
  isTRUE(all.equal(as.matrix(sparse$XtX), XtX, tolerance = 1e-10)),
  isTRUE(all.equal(sparse$Xty, original$Xty, tolerance = 1e-10))
)

# Eigen output contains the complete, unfiltered decomposition.
eigen_output <- compute_ss_from_gwas(
  beta, se, LD, n, response_var = stats::var(y), output = "eigen"
)
reconstructed <- eigen_output$XtX_eigenvectors_raw %*%
  (eigen_output$XtX_eigenvalues_raw *
     t(eigen_output$XtX_eigenvectors_raw))
stopifnot(
  ncol(eigen_output$XtX_eigenvectors_raw) == ncol(X),
  length(eigen_output$XtX_eigenvalues_raw) == ncol(X),
  isTRUE(all.equal(unname(reconstructed), unname(XtX), tolerance = 1e-10)),
  eigen_output$XtX_eigenvalue_tolerance > 0
)

# Indefinite correlation input is not repaired or filtered; it produces a
# warning and every raw eigenpair is returned.
indefinite_LD <- matrix(c(
  1, 0.9, 0.9,
  0.9, 1, -0.9,
  0.9, -0.9, 1
), nrow = 3L)
dimnames(indefinite_LD) <- list(predictor_names[1:3], predictor_names[1:3])
warning_text <- NULL
indefinite <- withCallingHandlers(
  compute_ss_from_gwas(
    beta[1:3], se[1:3], indefinite_LD, n, output = "eigen"
  ),
  warning = function(condition) {
    warning_text <<- conditionMessage(condition)
    invokeRestart("muffleWarning")
  }
)
stopifnot(
  length(indefinite$XtX_eigenvalues_raw) == 3L,
  any(indefinite$XtX_eigenvalues_raw < 0),
  grepl("returned unchanged", warning_text, fixed = TRUE)
)

# Basic input validation rejects malformed conversions without attempting
# broader GWAS/LD mismatch diagnostics.
invalid_calls <- list(
  function() compute_ss_from_gwas(beta[-1], se, LD, n),
  function() compute_ss_from_gwas(beta, replace(se, 1, 0), LD, n),
  function() compute_ss_from_gwas(beta, se, LD[-1, -1], n),
  function() compute_ss_from_gwas(beta, se, LD + diag(0.01, ncol(LD)), n),
  function() compute_ss_from_gwas(beta, se, LD, 2),
  function() compute_ss_from_gwas(
    beta, se, LD, n, scale = "original"
  ),
  function() compute_ss_from_gwas(
    beta, se, LD, n, response_var = 1, scale = "standardized"
  ),
  function() compute_ss_from_gwas(
    beta, se, LD, n, residual_df = n
  )
)
stopifnot(all(vapply(
  invalid_calls,
  function(call) inherits(try(call(), silent = TRUE), "try-error"),
  logical(1)
)))
