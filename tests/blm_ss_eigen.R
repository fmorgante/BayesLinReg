library(BayesLinReg)

blm_public <- BayesLinReg::blm
blm <- function(..., store_samples = TRUE, store_coefficient_cov = TRUE) {
  blm_public(..., store_samples = store_samples,
             store_coefficient_cov = store_coefficient_cov)
}
blm_ss_public <- BayesLinReg::blm_ss
blm_ss <- function(..., store_samples = TRUE, store_coefficient_cov = TRUE) {
  blm_ss_public(..., store_samples = store_samples,
                store_coefficient_cov = store_coefficient_cov)
}
blm_ss_eigen_public <- BayesLinReg::blm_ss_eigen
blm_ss_eigen <- function(..., store_samples = TRUE,
                         store_coefficient_cov = TRUE) {
  blm_ss_eigen_public(..., store_samples = store_samples,
                      store_coefficient_cov = store_coefficient_cov)
}

set.seed(601)
n <- 45
p <- 70
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("x", seq_len(p))
y <- drop(
  0.75 + X[, 1:5, drop = FALSE] %*% c(1.5, -1, 0.5, 0, -0.25) +
    rnorm(n, sd = 0.8)
)
XtX <- crossprod(X)
Xty <- drop(crossprod(X, y))
yty <- sum(y^2)
X_means <- colMeans(X)
y_mean <- mean(y)
centered_X <- sweep(X, 2L, X_means, FUN = "-")
decomposition <- eigen(crossprod(centered_X), symmetric = TRUE)
tolerance <- sqrt(.Machine$double.eps) * max(decomposition$values)
keep <- decomposition$values > tolerance
eigenvectors <- decomposition$vectors[, keep, drop = FALSE]
eigenvalues <- decomposition$values[keep]
rownames(eigenvectors) <- colnames(X)

# Internal preprocessing matches the original matrix expressions without
# creating squared, transposed, or swept p-by-q intermediates.
centered_Xty <- Xty - n * X_means * y_mean
prepared <- BayesLinReg:::prepare_eigen_statistics_cpp(
  eigenvectors, eigenvalues, centered_Xty
)
expected_projection <- drop(crossprod(eigenvectors, centered_Xty))
test_order <- as.integer(c(3L, 1L, 2L, seq.int(4L, p)))
test_scale <- seq(0.8, 1.2, length.out = p)
prepared_factor <- BayesLinReg:::build_scaled_eigen_factor_cpp(
  eigenvectors, eigenvalues, test_scale, test_order
)
expected_factor <- sweep(
  sqrt(eigenvalues) * t(eigenvectors[test_order, , drop = FALSE]),
  2L, test_scale, FUN = "/"
)
stopifnot(
  isTRUE(all.equal(
    prepared$transformed_response,
    expected_projection / sqrt(eigenvalues),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    prepared$projected_crossproduct,
    unname(drop(eigenvectors %*% expected_projection)),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    prepared$diagonal,
    unname(drop(eigenvectors^2 %*% eigenvalues)),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    prepared_factor, unname(expected_factor), tolerance = 1e-12
  ))
)

# Exact eigen sufficient statistics reproduce dense sufficient-statistic draws
# in a rank-deficient p > n problem for every coefficient prior.
for (model in c(
  "Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab"
)) {
  eta <- list(model = model)
  set.seed(602)
  dense_fit <- blm_ss(
    n, XtX, Xty, ETA = eta, yty = yty,
    X_means = X_means, y_mean = y_mean,
    residual_shape = 3, residual_scale = 1.5,
    iterations = 90, burnin = 30,
    store_coefficient_cov = FALSE
  )
  set.seed(602)
  eigen_fit <- blm_ss_eigen(
    n,
    XtX_eigenvectors = eigenvectors,
    XtX_eigenvalues = eigenvalues,
    XtX_prop_var = 1,
    Xty = Xty,
    ETA = eta,
    yty = yty,
    X_means = X_means,
    y_mean = y_mean,
    residual_shape = 3,
    residual_scale = 1.5,
    iterations = 90,
    burnin = 30,

    store_coefficient_cov = FALSE,
    check_eigenvectors = TRUE
  )
  eigen_diagnostics <- eigen_fit[c(
    "XtX_eigen_rank", "XtX_prop_var", "XtX_approximate",
    "residual_sse_offset"
  )]
  stopifnot(
    identical(dense_fit$likelihood_df, as.integer(n - 1)),
    identical(eigen_fit$likelihood_df, as.integer(n - 1))
  )
  eigen_fit[c(
    "XtX_eigen_rank", "XtX_prop_var", "XtX_approximate",
    "residual_sse_offset"
  )] <- NULL
  stopifnot(
    isTRUE(all.equal(dense_fit, eigen_fit, tolerance = 1e-8)),
    identical(eigen_diagnostics$XtX_eigen_rank, sum(keep)),
    identical(eigen_diagnostics$XtX_prop_var, 1),
    identical(eigen_diagnostics$XtX_approximate, FALSE),
    eigen_diagnostics$residual_sse_offset >= 0
  )
}

# Reordered multiple blocks and expected-PVE residual calibration also retain
# exact dense-path behavior.
eta_dense <- list(
  first = list(
    indices = c("x3", "x1", "x2"), model = "SpikeSlab",
    expected_pve = 0.2
  ),
  second = list(
    indices = paste0("x", 4:p), model = "SpikeMultiSlab",
    expected_pve = 0.3
  )
)
set.seed(603)
pve_dense <- blm_ss(
  n, XtX, Xty, ETA = eta_dense, yty = yty,
  X_means = X_means, y_mean = y_mean,
  residual_shape = 4,
  iterations = 80, burnin = 30,
  store_samples = FALSE, store_coefficient_cov = FALSE
)
set.seed(603)
pve_eigen <- blm_ss_eigen(
  n, eigenvectors, eigenvalues, 1, Xty, ETA = eta_dense, yty = yty,
  X_means = X_means, y_mean = y_mean,
  residual_shape = 4,
  iterations = 80, burnin = 30,
  store_samples = FALSE, store_coefficient_cov = FALSE
)
pve_eigen[c(
  "XtX_eigen_rank", "XtX_prop_var", "XtX_approximate",
  "residual_sse_offset"
)] <- NULL
stopifnot(isTRUE(all.equal(pve_dense, pve_eigen, tolerance = 1e-8)))

# Truncation is explicit in the result and produces a valid approximate fit.
truncated_rank <- 12L
truncated_prop <- sum(eigenvalues[seq_len(truncated_rank)]) /
  sum(eigenvalues)
truncated_fit <- blm_ss_eigen(
  n,
  eigenvectors[, seq_len(truncated_rank), drop = FALSE],
  eigenvalues[seq_len(truncated_rank)],
  truncated_prop,
  Xty,
  ETA = list(model = "Normal"),
  yty = yty,
  X_means = X_means,
  y_mean = y_mean,
  residual_var = 1,
  iterations = 70,
  burnin = 20,

  check_eigenvectors = TRUE
)
stopifnot(
  identical(truncated_fit$XtX_eigen_rank, truncated_rank),
  identical(truncated_fit$XtX_prop_var, truncated_prop),
  identical(truncated_fit$XtX_approximate, TRUE),
  all(is.finite(truncated_fit$ETA$ETA1$coefficient_mean))
)

# A fixed residual variance allows omission of yty.
fixed_fit <- blm_ss_eigen(
  n, eigenvectors, eigenvalues, 1, Xty,
  ETA = list(model = "Normal"),
  X_means = X_means, y_mean = y_mean,
  residual_var = 1,
  iterations = 60, burnin = 20
)
stopifnot(all(fixed_fit$residual_var_samples == 1))

# Invalid eigen representations and incompatible sufficient statistics fail.
invalid_calls <- list(
  function() blm_ss_eigen(
    n, eigenvectors, eigenvalues, 0, Xty,
    ETA = list(model = "Normal"), residual_var = 1
  ),
  function() blm_ss_eigen(
    n, eigenvectors, c(eigenvalues[-1], 0), 1, Xty,
    ETA = list(model = "Normal"), residual_var = 1
  ),
  function() blm_ss_eigen(
    n, eigenvectors[, -1, drop = FALSE], eigenvalues, 1, Xty,
    ETA = list(model = "Normal"), residual_var = 1
  ),
  function() blm_ss_eigen(
    n, eigenvectors * 2, eigenvalues, 1, Xty,
    ETA = list(model = "Normal"), residual_var = 1,
    check_eigenvectors = TRUE
  ),
  function() blm_ss_eigen(
    n, eigenvectors, eigenvalues, 1,
    stats::setNames(as.numeric(Xty), rev(rownames(eigenvectors))),
    ETA = list(model = "Normal"), residual_var = 1
  ),
  function() blm_ss_eigen(
    n, eigenvectors, eigenvalues, 1, Xty,
    ETA = list(model = "Normal"),
    X_means = stats::setNames(X_means, rev(rownames(eigenvectors))),
    y_mean = y_mean, residual_var = 1
  ),
  function() blm_ss_eigen(
    n, eigenvectors, eigenvalues, 1, Xty + seq_along(Xty),
    ETA = list(model = "Normal"), residual_var = 1
  ),
  function() blm_ss_eigen(
    n, eigenvectors, eigenvalues, 1, Xty,
    ETA = list(model = "Normal"), residual_shape = 2
  )
)
stopifnot(all(vapply(
  invalid_calls,
  function(call) inherits(try(call(), silent = TRUE), "try-error"),
  logical(1)
)))

# Lists of centered eigen blocks use the dedicated low-rank sampler without
# constructing either a global pseudo-design or block Gram matrices.
set.seed(606)
eigen_block_sizes <- c(first = 4L, second = 5L)
centered_gram_blocks <- lapply(eigen_block_sizes, function(size) {
  candidate <- matrix(rnorm(size^2), nrow = size)
  crossprod(candidate) + diag(size)
})
eigen_block_decompositions <- lapply(
  centered_gram_blocks, eigen, symmetric = TRUE
)
eigenvector_blocks <- lapply(
  eigen_block_decompositions, `[[`, "vectors"
)
eigenvalue_blocks <- lapply(
  eigen_block_decompositions, `[[`, "values"
)
block_predictor_names <- paste0("b", seq_len(sum(eigen_block_sizes)))
block_ends <- cumsum(eigen_block_sizes)
block_starts <- block_ends - eigen_block_sizes + 1L
for (block in seq_along(eigenvector_blocks)) {
  rownames(eigenvector_blocks[[block]]) <- block_predictor_names[
    seq.int(block_starts[block], block_ends[block])
  ]
}
block_centered_gram <- as.matrix(Matrix::bdiag(centered_gram_blocks))
dimnames(block_centered_gram) <- list(block_predictor_names, block_predictor_names)
block_means <- seq_along(block_predictor_names) / 50
block_y_mean <- 0.3
block_centered_Xty <- rnorm(length(block_predictor_names))
block_Xty <- block_centered_Xty + n * block_means * block_y_mean
block_XtX <- block_centered_gram + n * tcrossprod(block_means)
block_yty <- 200 + n * block_y_mean^2
block_eta <- list(
  crossing = list(
    indices = c("b1", "b5", "b2"), model = "Normal"
  ),
  remaining = list(
    indices = setdiff(block_predictor_names, c("b1", "b5", "b2")),
    model = "SpikeSlab"
  )
)
block_fit_arguments <- list(
  n = n, Xty = block_Xty, ETA = block_eta, yty = block_yty,
  X_means = block_means, y_mean = block_y_mean,
  residual_shape = 3, residual_scale = 1.5,
  iterations = 80, burnin = 30,
  store_coefficient_cov = FALSE, compute_pve = TRUE
)
set.seed(607)
block_dense_fit <- do.call(
  blm_ss, c(list(XtX = block_XtX), block_fit_arguments)
)
set.seed(607)
block_eigen_fit <- do.call(
  blm_ss_eigen,
  c(
    list(
      XtX_eigenvectors = eigenvector_blocks,
      XtX_eigenvalues = eigenvalue_blocks,
      XtX_prop_var = c(first = 1, second = 1),
      check_eigenvectors = TRUE
    ),
    block_fit_arguments
  )
)
block_diagnostics <- block_eigen_fit[c(
  "XtX_eigen_rank", "XtX_prop_var", "XtX_approximate",
  "residual_sse_offset", "XtX_representation", "XtX_number_of_blocks",
  "XtX_block_sizes", "XtX_eigen_rank_by_block",
  "XtX_prop_var_by_block", "XtX_cross_block_assumption", "nthreads"
)]
block_eigen_fit[names(block_diagnostics)] <- NULL
stopifnot(
  isTRUE(all.equal(block_dense_fit, block_eigen_fit, tolerance = 1e-8)),
  identical(block_diagnostics$XtX_eigen_rank, sum(eigen_block_sizes)),
  identical(block_diagnostics$XtX_prop_var, 1),
  identical(block_diagnostics$XtX_approximate, FALSE),
  identical(block_diagnostics$XtX_representation, "block_diagonal_eigen"),
  identical(block_diagnostics$XtX_number_of_blocks, 2L),
  identical(
    unname(block_diagnostics$XtX_block_sizes), unname(eigen_block_sizes)
  ),
  identical(
    unname(block_diagnostics$XtX_eigen_rank_by_block),
    unname(eigen_block_sizes)
  ),
  identical(unname(block_diagnostics$XtX_prop_var_by_block), c(1, 1)),
  identical(block_diagnostics$XtX_cross_block_assumption, "zero"),
  identical(block_diagnostics$nthreads, 1L)
)

# Threaded block-eigen sweeps cover every prior, learned residual variance,
# PVE, and online summaries. They are reproducible across runs and thread
# counts and accept nonzero means because the eigenpairs represent centered
# XtX.
threaded_eta <- list(
  normal = list(indices = c("b1", "b5"), model = "Normal"),
  selection = list(indices = c("b2", "b6"), model = "SpikeSlab"),
  shrinkage = list(indices = c("b3", "b7"), model = "GlobalLocal"),
  multi = list(indices = c("b4", "b8", "b9"), model = "SpikeMultiSlab")
)
threaded_arguments <- c(
  list(
    XtX_eigenvectors = eigenvector_blocks,
    XtX_eigenvalues = eigenvalue_blocks,
    XtX_prop_var = 1,
    nthreads = 2L
  ),
  block_fit_arguments
)
threaded_arguments$ETA <- threaded_eta
threaded_arguments$store_samples <- FALSE
set.seed(607)
threaded_eigen_fit <- do.call(blm_ss_eigen, threaded_arguments)
set.seed(607)
threaded_eigen_repeat <- do.call(blm_ss_eigen, threaded_arguments)
threaded_arguments$nthreads <- 3L
set.seed(607)
threaded_eigen_three <- do.call(blm_ss_eigen, threaded_arguments)
stopifnot(
  identical(threaded_eigen_fit, threaded_eigen_repeat),
  identical(threaded_eigen_fit$nthreads, 2L),
  identical(threaded_eigen_three$nthreads, 3L),
  identical(threaded_eigen_fit$ETA, threaded_eigen_three$ETA),
  identical(
    threaded_eigen_fit$residual_var_mean,
    threaded_eigen_three$residual_var_mean
  ),
  identical(threaded_eigen_fit$total_pve_mean,
            threaded_eigen_three$total_pve_mean),
  identical(threaded_eigen_fit$cross_block_pve_mean,
            threaded_eigen_three$cross_block_pve_mean),
  is.null(threaded_eigen_fit$residual_var_samples),
  threaded_eigen_fit$residual_var_mean > 0,
  is.finite(threaded_eigen_fit$total_pve_mean)
)

# Per-block truncation metadata includes a trace-weighted aggregate fraction.
truncated_block_values <- eigenvalue_blocks
truncated_block_vectors <- eigenvector_blocks
truncated_block_values[[1L]] <- eigenvalue_blocks[[1L]][-length(eigenvalue_blocks[[1L]])]
truncated_block_vectors[[1L]] <- eigenvector_blocks[[1L]][
  , -ncol(eigenvector_blocks[[1L]]), drop = FALSE
]
block_prop <- c(
  first = sum(truncated_block_values[[1L]]) / sum(eigenvalue_blocks[[1L]]),
  second = 1
)
truncated_block_fit <- blm_ss_eigen(
  n, truncated_block_vectors, truncated_block_values, block_prop,
  block_Xty, block_eta, block_yty, block_means, block_y_mean,
  residual_var = 1, iterations = 50, burnin = 20
)
expected_aggregate_prop <-
  sum(vapply(truncated_block_values, sum, numeric(1))) /
  sum(vapply(truncated_block_values, sum, numeric(1)) / block_prop)
stopifnot(
  identical(truncated_block_fit$XtX_approximate, TRUE),
  isTRUE(all.equal(truncated_block_fit$XtX_prop_var, expected_aggregate_prop)),
  identical(truncated_block_fit$XtX_prop_var_by_block, block_prop)
)

unnamed_eigenvector_blocks <- lapply(eigenvector_blocks, function(value) {
  rownames(value) <- NULL
  value
})
unnamed_block_fit <- blm_ss_eigen(
  n, unnamed_eigenvector_blocks, eigenvalue_blocks, 1,
  stats::setNames(block_Xty, block_predictor_names),
  ETA = list(model = "Normal"), residual_var = 1,
  iterations = 20, burnin = 10
)
stopifnot(identical(
  names(unnamed_block_fit$ETA[[1L]]$coefficient_mean),
  block_predictor_names
))

partially_named_eigenvectors <- eigenvector_blocks
rownames(partially_named_eigenvectors[[2L]]) <- NULL
duplicated_eigenvector_names <- eigenvector_blocks
rownames(duplicated_eigenvector_names[[2L]])[1L] <-
  rownames(duplicated_eigenvector_names[[1L]])[1L]
block_invalid_calls <- list(
  function() blm_ss_eigen(
    n, eigenvector_blocks, eigenvalues, 1, block_Xty, block_eta,
    residual_var = 1
  ),
  function() blm_ss_eigen(
    n, eigenvector_blocks, eigenvalue_blocks, c(1, 1, 1), block_Xty,
    block_eta, residual_var = 1
  ),
  function() blm_ss_eigen(
    n, eigenvectors, eigenvalues, 1, Xty, list(model = "Normal"),
    residual_var = 1, nthreads = 2
  ),
  function() blm_ss_eigen(
    n, eigenvector_blocks, eigenvalue_blocks, 1, block_Xty, block_eta,
    residual_var = 1, nthreads = 2, nchains = 2
  ),
  function() blm_ss_eigen(
    n, partially_named_eigenvectors, eigenvalue_blocks, 1, block_Xty,
    block_eta, residual_var = 1
  ),
  function() blm_ss_eigen(
    n, duplicated_eigenvector_names, eigenvalue_blocks, 1, block_Xty,
    block_eta, residual_var = 1
  ),
  function() blm_ss_eigen(
    n, eigenvector_blocks, eigenvalue_blocks, 1,
    stats::setNames(block_Xty, rev(block_predictor_names)), block_eta,
    residual_var = 1
  )
)
stopifnot(all(vapply(
  block_invalid_calls,
  function(call) inherits(try(call(), silent = TRUE), "try-error"),
  logical(1)
)))
