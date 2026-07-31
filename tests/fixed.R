library(BayesLinReg)

set.seed(701)
n <- 80L
X <- cbind(
  exposure = rnorm(n),
  treatment = rbinom(n, 1, 0.45),
  marker1 = rnorm(n),
  marker2 = rnorm(n)
)
y <- drop(1.5 + X %*% c(1.2, -0.8, 0.4, 0) + rnorm(n, sd = 0.7))
fixed_eta <- list(
  covariates = list(X = X[, 1:2], model = "Fixed"),
  markers = list(X = X[, 3:4], model = "Normal", var_scale = 2)
)

# R and Rcpp implement the same flat-prior coefficient update.
fixed_r <- blm(
  y, fixed_eta, residual_var = 0.49,
  iterations = 300, burnin = 100, seed = 702, version = "R"
)
fixed_rcpp <- blm(
  y, fixed_eta, residual_var = 0.49,
  iterations = 300, burnin = 100, seed = 702, version = "Rcpp",
  compute_pve = TRUE
)
stopifnot(
  identical(fixed_rcpp$ETA$covariates$model, "Fixed"),
  is.null(fixed_rcpp$ETA$covariates$var_shape),
  is.finite(fixed_rcpp$ETA$covariates$pve_mean),
  isTRUE(all.equal(
    fixed_r$ETA$covariates$coefficient_samples,
    fixed_rcpp$ETA$covariates$coefficient_samples,
    tolerance = 1e-12
  )),
  max(abs(
    fixed_rcpp$ETA$covariates$coefficient_mean -
      stats::coef(stats::lm(y ~ X))[-1][1:2]
  )) < 0.25
)

# Learning residual variance uses the usual explicitly sampled-coefficient
# conditional and requires an independently specified residual prior.
fixed_learned <- blm(
  y, list(X = X[, 1:2], model = "Fixed"),
  residual_shape = 3, residual_scale = 1,
  iterations = 250, burnin = 100, seed = 703,
  store_samples = FALSE
)
stopifnot(
  fixed_learned$residual_var_mean > 0,
  all(is.finite(fixed_learned$ETA$ETA1$coefficient_mean))
)

# Expected-PVE calibration remains available to penalized blocks when the
# residual prior is supplied independently.
fixed_pve <- blm(
  y,
  list(
    covariates = list(X = X[, 1:2], model = "Fixed"),
    markers = list(
      X = X[, 3:4], model = "GlobalLocal", expected_nonzero = 0.5,
      expected_pve = 0.2
    )
  ),
  residual_shape = 3, residual_scale = 1,
  iterations = 80, burnin = 30, seed = 704
)
stopifnot(
  identical(
    fixed_pve$ETA$markers$global_scale_calibration, "expected_pve"
  ),
  identical(fixed_pve$ETA$markers$expected_pve, 0.2)
)

# Raw and sufficient-statistic entry points agree under the same seed.
XtX <- crossprod(X)
Xty <- drop(crossprod(X, y))
yty <- sum(y^2)
ss_eta <- list(
  covariates = list(indices = 1:2, model = "Fixed"),
  markers = list(indices = 3:4, model = "Normal", var_scale = 2)
)
fixed_ss <- blm_ss(
  n, XtX, Xty, ss_eta, yty = yty,
  X_means = colMeans(X), y_mean = mean(y), residual_var = 0.49,
  iterations = 300, burnin = 100, seed = 702
)
fixed_sparse <- blm_ss(
  n, Matrix::Matrix(XtX, sparse = TRUE), Xty, ss_eta, yty = yty,
  X_means = colMeans(X), y_mean = mean(y), residual_var = 0.49,
  iterations = 300, burnin = 100, seed = 702
)
centered_X <- sweep(X, 2, colMeans(X), FUN = "-")
decomposition <- eigen(crossprod(centered_X), symmetric = TRUE)
keep <- decomposition$values > 1e-10
fixed_eigen <- blm_ss_eigen(
  n,
  decomposition$vectors[, keep, drop = FALSE],
  decomposition$values[keep], 1, Xty, ss_eta, yty = yty,
  X_means = colMeans(X), y_mean = mean(y), residual_var = 0.49,
  iterations = 300, burnin = 100, seed = 702
)
expected <- fixed_rcpp$ETA$covariates$coefficient_mean
stopifnot(
  isTRUE(all.equal(fixed_ss$ETA$covariates$coefficient_mean, expected,
                   tolerance = 1e-10)),
  isTRUE(all.equal(fixed_sparse$ETA$covariates$coefficient_mean, expected,
                   tolerance = 1e-10)),
  isTRUE(all.equal(fixed_eigen$ETA$covariates$coefficient_mean, expected,
                   tolerance = 1e-10))
)

# Flat-prior fields and jointly rank-deficient fixed predictors are rejected.
duplicate_X <- cbind(first = X[, 1], duplicate = X[, 1])
invalid_calls <- list(
  function() blm(
    y, list(X = X[, 1:2], model = "Fixed", expected_pve = 0.2),
    residual_var = 1
  ),
  function() blm(
    y, list(X = X[, 1:2], model = "Fixed"), residual_shape = 2
  ),
  function() blm(
    y,
    list(
      first = list(X = duplicate_X[, 1, drop = FALSE], model = "Fixed"),
      second = list(X = duplicate_X[, 2, drop = FALSE], model = "Fixed")
    ),
    residual_var = 1
  ),
  function() blm_ss(
    n, crossprod(duplicate_X), drop(crossprod(duplicate_X, y)),
    ETA = list(model = "Fixed"), yty = yty,
    X_means = colMeans(duplicate_X), y_mean = mean(y), residual_var = 1
  ),
  function() blm_ss_eigen(
    n,
    decomposition$vectors[, 1, drop = FALSE],
    decomposition$values[1],
    decomposition$values[1] / sum(decomposition$values),
    Xty,
    ETA = list(
      fixed = list(indices = 1:2, model = "Fixed"),
      markers = list(indices = 3:4, model = "Normal")
    ),
    X_means = colMeans(X), y_mean = mean(y), residual_var = 1
  )
)
stopifnot(all(vapply(
  invalid_calls,
  function(call) inherits(try(call(), silent = TRUE), "try-error"),
  logical(1)
)))
