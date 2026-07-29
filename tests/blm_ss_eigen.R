library(BayesLinReg)

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

# Exact eigen sufficient statistics reproduce dense sufficient-statistic draws
# in a rank-deficient p > n problem for every coefficient prior.
for (model in c(
  "Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab"
)) {
  eta <- list(model = model)
  dense_fit <- blm_ss(
    n, XtX, Xty, ETA = eta, yty = yty,
    X_means = X_means, y_mean = y_mean,
    residual_shape = 3, residual_scale = 1.5,
    iterations = 90, burnin = 30, seed = 602,
    store_coefficient_cov = FALSE
  )
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
    seed = 602,
    store_coefficient_cov = FALSE,
    check_eigenvectors = TRUE
  )
  eigen_diagnostics <- eigen_fit[c(
    "XtX_eigen_rank", "XtX_prop_var", "XtX_approximate",
    "residual_sse_offset"
  )]
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
pve_dense <- blm_ss(
  n, XtX, Xty, ETA = eta_dense, yty = yty,
  X_means = X_means, y_mean = y_mean,
  residual_shape = 4,
  iterations = 80, burnin = 30, seed = 603,
  store_samples = FALSE, store_coefficient_cov = FALSE
)
pve_eigen <- blm_ss_eigen(
  n, eigenvectors, eigenvalues, 1, Xty, ETA = eta_dense, yty = yty,
  X_means = X_means, y_mean = y_mean,
  residual_shape = 4,
  iterations = 80, burnin = 30, seed = 603,
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
  seed = 604,
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
  iterations = 60, burnin = 20, seed = 605
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
