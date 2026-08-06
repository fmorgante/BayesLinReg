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

# The default sufficient-statistics result feeds directly into a model that
# does not use expected-PVE calibration without passing the unused reference
# response variance.
gwas_fit <- blm_ss(
  n = original$n,
  XtX = original$XtX,
  Xty = original$Xty,
  yty = original$yty,
  X_means = original$X_means,
  y_mean = original$y_mean,
  ETA = list(model = "SpikeMultiSlab"),
  residual_shape = 2,
  residual_scale = original$reference_response_var,
  iterations = 30,
  burnin = 10,
  seed = 902
)
stopifnot(inherits(gwas_fit, "blm_fit"))

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

# A list of LD matrices represents exact zero LD between blocks and produces
# the same reconstruction as an explicitly block-diagonal matrix. The blocks
# need not align with the prior blocks used by blm_ss().
ld_indices <- list(region_a = 1:2, region_b = 3:5)
LD_blocks <- lapply(ld_indices, function(indices) LD[indices, indices])
block_diagonal_LD <- matrix(0, ncol(X), ncol(X))
for (indices in ld_indices) {
  block_diagonal_LD[indices, indices] <- LD[indices, indices]
}
dimnames(block_diagonal_LD) <- list(predictor_names, predictor_names)
block_original <- compute_ss_from_gwas(
  beta, se, LD_blocks, n, response_var = stats::var(y)
)
matrix_original <- compute_ss_from_gwas(
  beta, se, block_diagonal_LD, n, response_var = stats::var(y)
)
block_standardized <- compute_ss_from_gwas(beta, se, LD_blocks, n)
matrix_standardized <- compute_ss_from_gwas(beta, se, block_diagonal_LD, n)
stopifnot(
  identical(names(block_original$XtX), names(LD_blocks)),
  identical(unlist(lapply(block_original$XtX, colnames), use.names = FALSE),
            predictor_names),
  all(vapply(seq_along(ld_indices), function(i) {
    indices <- ld_indices[[i]]
    isTRUE(all.equal(
      block_original$XtX[[i]], matrix_original$XtX[indices, indices],
      tolerance = 1e-10
    ))
  }, logical(1))),
  all(vapply(seq_along(ld_indices), function(i) {
    indices <- ld_indices[[i]]
    isTRUE(all.equal(
      block_standardized$XtX[[i]],
      matrix_standardized$XtX[indices, indices], tolerance = 1e-10
    ))
  }, logical(1))),
  isTRUE(all.equal(block_original$Xty, matrix_original$Xty,
                   tolerance = 1e-10)),
  isTRUE(all.equal(block_standardized$Xty, matrix_standardized$Xty,
                   tolerance = 1e-10))
)

# Upper- and lower-triangular symmetric sparse LD blocks remain dsCMatrix
# objects after reconstruction instead of being expanded to both triangles.
symmetric_LD_blocks <- Map(function(block, uplo) {
  methods::as(
    Matrix::forceSymmetric(Matrix::Matrix(block, sparse = TRUE), uplo = uplo),
    "dsCMatrix"
  )
}, LD_blocks, c("U", "L"))
names(symmetric_LD_blocks) <- names(LD_blocks)
symmetric_block_original <- compute_ss_from_gwas(
  beta, se, symmetric_LD_blocks, n, response_var = stats::var(y)
)
stopifnot(
  all(vapply(symmetric_block_original$XtX, inherits, logical(1),
             "dsCMatrix")),
  identical(vapply(symmetric_block_original$XtX, slot, character(1), "uplo"),
            c(region_a = "U", region_b = "L")),
  all(vapply(seq_along(ld_indices), function(i) {
    indices <- ld_indices[[i]]
    isTRUE(all.equal(
      as.matrix(symmetric_block_original$XtX[[i]]),
      matrix_original$XtX[indices, indices], tolerance = 1e-10
    ))
  }, logical(1)))
)

# Direct fits from matrix and list LD follow the same chain, including when
# prior blocks cross LD-block boundaries.
crossing_eta <- list(
  normal = list(indices = c("rs1", "rs3", "rs5"), model = "Normal"),
  selection = list(indices = c("rs2", "rs4"), model = "SpikeMultiSlab")
)
fit_arguments <- list(
  n = n, Xty = matrix_original$Xty, yty = matrix_original$yty,
  X_means = matrix_original$X_means, y_mean = matrix_original$y_mean,
  ETA = crossing_eta, residual_var = 1,
  iterations = 60, burnin = 20, seed = 903
)
matrix_LD_fit <- do.call(
  blm_ss, c(list(XtX = matrix_original$XtX), fit_arguments)
)
block_LD_fit <- do.call(
  blm_ss, c(list(XtX = block_original$XtX), fit_arguments)
)
stopifnot(
  isTRUE(all.equal(matrix_LD_fit$ETA, block_LD_fit$ETA,
                   tolerance = 1e-10)),
  isTRUE(all.equal(matrix_LD_fit$intercept_mean,
                   block_LD_fit$intercept_mean, tolerance = 1e-10)),
  identical(block_LD_fit$XtX_representation, "block_diagonal")
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
  ),
  function() compute_ss_from_gwas(beta, se, list(), n),
  function() compute_ss_from_gwas(beta, se, LD_blocks[1], n),
  function() compute_ss_from_gwas(
    beta, se, list(LD_blocks[[1]], unname(LD_blocks[[2]])), n
  ),
  function() compute_ss_from_gwas(
    beta, se, LD_blocks, n, output = "eigen"
  ),
  function() compute_ss_from_gwas(
    beta, se, list(LD_blocks[[1]], LD_blocks[[2]][, 1, drop = FALSE]), n
  )
)
stopifnot(all(vapply(
  invalid_calls,
  function(call) inherits(try(call(), silent = TRUE), "try-error"),
  logical(1)
)))
