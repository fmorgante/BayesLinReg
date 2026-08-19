library(BayesLinReg)

# Prior calibration uses likelihood_df for signal variance, while predictor
# standardization continues to use sample standard deviations.
n <- 9L
likelihood_df <- 3L
predictor_names <- paste0("x", 1:3)
gram_diagonal <- c(8, 18, 32)
XtX <- diag(gram_diagonal)
dimnames(XtX) <- list(predictor_names, predictor_names)
Xty <- stats::setNames(numeric(length(predictor_names)), predictor_names)
decomposition <- eigen(XtX, symmetric = TRUE)
rownames(decomposition$vectors) <- predictor_names

reference_response_var <- 2.5
expected_pve <- 0.24
var_shape <- 3
specifications <- list(
  Normal = list(model = "Normal", var_shape = var_shape,
                expected_pve = expected_pve),
  SpikeSlab = list(model = "SpikeSlab", var_shape = var_shape,
                    pi = c(a = 1, b = 3), expected_pve = expected_pve),
  SpikeMultiSlab = list(model = "SpikeMultiSlab", var_shape = var_shape,
    gamma = c(0, 0.25, 1.5), alpha = c(5, 3, 2),
    expected_pve = expected_pve)
)
prior_variance_factors <- c(
  Normal = 1,
  SpikeSlab = 1 / 4,
  SpikeMultiSlab = sum(c(5, 3, 2) / 10 * c(0, 0.25, 1.5))
)

fit_direct <- function(specification) {
  suppressWarnings(blm_ss(
    n = n, XtX = XtX, Xty = Xty, ETA = specification,
    residual_var = 1, reference_response_var = reference_response_var,
    likelihood_df = likelihood_df, iterations = 4, burnin = 2,
    version = "R"
  ))
}

fit_eigen <- function(specification) {
  suppressWarnings(blm_ss_eigen(
    n = n, XtX_eigenvectors = decomposition$vectors,
    XtX_eigenvalues = decomposition$values, XtX_prop_var = 1,
    Xty = Xty, ETA = specification, residual_var = 1,
    reference_response_var = reference_response_var,
    likelihood_df = likelihood_df, iterations = 4, burnin = 2,
    check_eigenvectors = TRUE
  ))
}

for (standardize in c(TRUE, FALSE)) {
  design_energy <- if (standardize) {
    length(gram_diagonal) * (n - 1) / likelihood_df
  } else {
    sum(gram_diagonal) / likelihood_df
  }
  for (model in names(specifications)) {
    specification <- specifications[[model]]
    specification$standardize <- standardize
    expected_scale <- (var_shape - 1) * expected_pve *
      reference_response_var /
      (design_energy * prior_variance_factors[[model]])

    set.seed(1901)
    direct <- fit_direct(specification)
    set.seed(1901)
    eigen_fit <- fit_eigen(specification)
    direct_block <- direct$ETA$ETA1
    eigen_block <- eigen_fit$ETA$ETA1
    implied_signal_mean <- design_energy * prior_variance_factors[[model]] *
      direct_block$var_scale / (var_shape - 1)

    stopifnot(
      isTRUE(all.equal(direct_block$var_scale, expected_scale,
                       tolerance = 1e-14)),
      isTRUE(all.equal(eigen_block$var_scale, expected_scale,
                       tolerance = 1e-14)),
      isTRUE(all.equal(direct_block$var_scale, eigen_block$var_scale,
                       tolerance = 1e-14)),
      isTRUE(all.equal(implied_signal_mean,
                       expected_pve * reference_response_var,
                       tolerance = 1e-14)),
      identical(direct_block$reference_response_var, reference_response_var),
      identical(eigen_block$reference_response_var, reference_response_var),
      identical(direct$likelihood_df, likelihood_df),
      identical(eigen_fit$likelihood_df, likelihood_df)
    )
  }
}

# Without means or an explicit likelihood dimension, blm_ss() fits no
# intercept and defaults to d_L = n.
default_specification <- specifications$Normal
default_specification$standardize <- TRUE
set.seed(1902)
default_fit <- suppressWarnings(blm_ss(
  n = n, XtX = XtX, Xty = Xty, ETA = default_specification,
  residual_var = 1, reference_response_var = reference_response_var,
  iterations = 4, burnin = 2, version = "R"
))
default_design_energy <- length(gram_diagonal) * (n - 1) / n
default_expected_scale <- (var_shape - 1) * expected_pve *
  reference_response_var / default_design_energy
stopifnot(
  identical(default_fit$likelihood_df, n),
  isTRUE(all.equal(default_fit$ETA$ETA1$var_scale, default_expected_scale,
                   tolerance = 1e-14))
)

# The GlobalLocal expected-sparsity heuristic deliberately retains n, rather
# than likelihood_df, in its prior scale. This is distinct from the
# moment-based expected-PVE calibrations above.
global_local_specification <- list(
  model = "GlobalLocal", expected_nonzero = 1, expected_pve = expected_pve
)
expected_global_scale <- 1 / (length(gram_diagonal) - 1) * sqrt(
  (1 - expected_pve) * reference_response_var / n
)
set.seed(1903)
global_local_direct <- fit_direct(global_local_specification)
set.seed(1903)
global_local_eigen <- fit_eigen(global_local_specification)
stopifnot(
  isTRUE(all.equal(
    global_local_direct$ETA$ETA1$global_scale, expected_global_scale,
    tolerance = 1e-14
  )),
  isTRUE(all.equal(
    global_local_eigen$ETA$ETA1$global_scale, expected_global_scale,
    tolerance = 1e-14
  )),
  identical(
    global_local_direct$ETA$ETA1$global_scale_calibration, "expected_pve"
  ),
  identical(
    global_local_eigen$ETA$ETA1$global_scale_calibration, "expected_pve"
  )
)
