library(blm)

x <- cbind(
  first = 1:5,
  second = c(0, 1, 0, 1, 0)
)
y <- 1 + 2 * x[, "first"] - 3 * x[, "second"]

known_fit <- multiple_blm(
  y = y,
  X = x,
  prior_var = 10,
  residual_var = 1,
  standardize = FALSE
)
expected_cov <- diag(c(1 / 10.1, 1 / 1.3))
dimnames(expected_cov) <- list(colnames(x), colnames(x))
expected_mean <- c(first = 20 / 10.1, second = -3.6 / 1.3)
expected_intercept <- mean(y) - sum(colMeans(x) * expected_mean)
expected_intercept_var <- 1 / 5 +
  drop(crossprod(colMeans(x), expected_cov %*% colMeans(x)))

stopifnot(
  identical(
    names(known_fit),
    c(
      "coefficient_mean", "coefficient_cov",
      "intercept_mean", "intercept_var"
    )
  ),
  isTRUE(all.equal(known_fit$coefficient_mean, expected_mean)),
  isTRUE(all.equal(known_fit$coefficient_cov, expected_cov)),
  isTRUE(all.equal(known_fit$intercept_mean, expected_intercept)),
  isTRUE(all.equal(known_fit$intercept_var, expected_intercept_var))
)

# A scalar prior variance and a repeated vector are equivalent.
vector_prior_fit <- multiple_blm(
  y, x, c(10, 10), residual_var = 1, standardize = FALSE
)
stopifnot(isTRUE(all.equal(known_fit, vector_prior_fit)))

# Data-frame inputs retain their predictor names.
data_frame_fit <- multiple_blm(
  y, as.data.frame(x), 10, residual_var = 1, standardize = FALSE
)
stopifnot(isTRUE(all.equal(known_fit, data_frame_fit)))

# With one predictor, multiple_blm agrees with simple_blm.
simple_y <- 1 + 2 * x[, "first"]
simple_fit <- simple_blm(simple_y, x[, "first"], 10, residual_var = 1)
one_predictor_fit <- multiple_blm(
  simple_y, x[, "first", drop = FALSE], 10, residual_var = 1,
  standardize = FALSE
)

# Standardized fits return coefficients and their covariance on the original
# predictor scale.
predictor_sd <- apply(x, 2L, stats::sd)
working_x <- sweep(x, 2L, predictor_sd, FUN = "/")
manual_standardized_fit <- multiple_blm(
  y, working_x, 10, residual_var = 1, standardize = FALSE
)
automatic_standardized_fit <- multiple_blm(y, x, 10, residual_var = 1)
stopifnot(
  isTRUE(all.equal(
    automatic_standardized_fit$coefficient_mean,
    manual_standardized_fit$coefficient_mean / predictor_sd
  )),
  isTRUE(all.equal(
    automatic_standardized_fit$coefficient_cov,
    manual_standardized_fit$coefficient_cov /
      outer(predictor_sd, predictor_sd)
  )),
  isTRUE(all.equal(
    automatic_standardized_fit$intercept_mean,
    manual_standardized_fit$intercept_mean
  )),
  isTRUE(all.equal(
    automatic_standardized_fit$intercept_var,
    manual_standardized_fit$intercept_var
  ))
)
stopifnot(
  isTRUE(all.equal(
    unname(one_predictor_fit$coefficient_mean),
    simple_fit$slope_mean
  )),
  isTRUE(all.equal(
    drop(one_predictor_fit$coefficient_cov),
    simple_fit$slope_var
  )),
  isTRUE(all.equal(one_predictor_fit$intercept_mean, simple_fit$intercept_mean)),
  isTRUE(all.equal(one_predictor_fit$intercept_var, simple_fit$intercept_var))
)

learned_fit <- multiple_blm(
  y = y,
  X = x,
  prior_var = 10,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1500,
  burnin = 500,
  thin = 2,
  seed = 42
)
stopifnot(
  identical(
    names(learned_fit),
    c(
      "coefficient_mean", "coefficient_cov", "intercept_mean",
      "intercept_var", "residual_var_mean", "residual_var_var",
      "coefficient_samples", "intercept_samples", "residual_var_samples"
    )
  ),
  identical(dim(learned_fit$coefficient_samples), c(500L, 2L)),
  identical(colnames(learned_fit$coefficient_samples), colnames(x)),
  length(learned_fit$intercept_samples) == 500,
  length(learned_fit$residual_var_samples) == 500,
  all(learned_fit$residual_var_samples > 0),
  isTRUE(all.equal(
    learned_fit$coefficient_mean,
    colMeans(learned_fit$coefficient_samples)
  )),
  isTRUE(all.equal(
    learned_fit$coefficient_cov,
    cov(learned_fit$coefficient_samples)
  )),
  identical(
    learned_fit$residual_var_mean,
    mean(learned_fit$residual_var_samples)
  )
)

# The seed makes MCMC output reproducible.
repeated_fit <- multiple_blm(
  y, x, 10,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1500,
  burnin = 500,
  thin = 2,
  seed = 42
)
stopifnot(isTRUE(all.equal(learned_fit, repeated_fit)))

# Chain combination retains the origin of every draw.
mock_chain <- list(
  coefficient_samples = matrix(1:4, ncol = 2),
  intercept_samples = c(1, 2),
  residual_var_samples = c(3, 4)
)
mock_combined <- blm:::.combine_blm_chains(
  list(mock_chain, mock_chain),
  use_spike_slab = FALSE
)
stopifnot(
  identical(mock_combined$chain_id, c(1L, 1L, 2L, 2L)),
  identical(dim(mock_combined$coefficient_samples), c(4L, 2L))
)

mock_global_local_chain <- c(
  mock_chain,
  list(
    local_var_samples = matrix(1:4, ncol = 2),
    tau_sq_samples = c(1, 2)
  )
)
mock_global_local_combined <- blm:::.combine_blm_chains(
  list(mock_global_local_chain, mock_global_local_chain),
  coefficient_prior = "global_local"
)
stopifnot(
  identical(dim(mock_global_local_combined$local_var_samples), c(4L, 2L)),
  identical(mock_global_local_combined$tau_sq_samples, c(1, 2, 1, 2))
)

# Socket-based multisession tests are enabled explicitly outside restricted
# package-check environments.
if (identical(Sys.getenv("BLM_TEST_FUTURE"), "true")) {
  parallel_fit <- multiple_blm(
    y, x, 10,
    residual_shape = 2,
    residual_scale = 1,
    iterations = 100,
    burnin = 20,
    thin = 2,
    seed = 123,
    nchains = 2,
    verbose = TRUE
  )
  repeated_parallel_fit <- multiple_blm(
    y, x, 10,
    residual_shape = 2,
    residual_scale = 1,
    iterations = 100,
    burnin = 20,
    thin = 2,
    seed = 123,
    nchains = 2,
    verbose = TRUE
  )
  stopifnot(
    isTRUE(all.equal(parallel_fit, repeated_parallel_fit)),
    parallel_fit$nchains == 2L,
    identical(parallel_fit$chain_id, rep.int(1:2, c(40L, 40L))),
    identical(dim(parallel_fit$coefficient_samples), c(80L, 2L)),
    !identical(
      parallel_fit$coefficient_samples[parallel_fit$chain_id == 1L, ],
      parallel_fit$coefficient_samples[parallel_fit$chain_id == 2L, ]
    ),
    isTRUE(all.equal(
      parallel_fit$coefficient_mean,
      colMeans(parallel_fit$coefficient_samples)
    ))
  )
}

# The Rcpp and R implementations target the same posterior.
r_fit <- multiple_blm(
  y, x, 10,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1500,
  burnin = 500,
  thin = 2,
  seed = 42,
  version = "R"
)
stopifnot(
  identical(dim(r_fit$coefficient_samples), c(500L, 2L)),
  max(abs(r_fit$coefficient_mean - learned_fit$coefficient_mean)) < 0.2,
  abs(r_fit$residual_var_mean - learned_fit$residual_var_mean) < 0.2
)

# Component-wise conditioning is exercised with correlated predictors.
correlated_X <- cbind(
  first = 1:8,
  second = c(2, 1, 4, 3, 6, 5, 8, 7)
)
correlated_y <- 1 + 2 * correlated_X[, 1] - correlated_X[, 2]
correlated_rcpp <- multiple_blm(
  correlated_y, correlated_X, 10,
  residual_shape = 2, residual_scale = 1,
  iterations = 600, burnin = 200, seed = 91,
  version = "Rcpp"
)
correlated_r <- multiple_blm(
  correlated_y, correlated_X, 10,
  residual_shape = 2, residual_scale = 1,
  iterations = 600, burnin = 200, seed = 91,
  version = "R"
)
stopifnot(
  max(abs(
    correlated_rcpp$coefficient_mean - correlated_r$coefficient_mean
  )) < 0.05,
  abs(correlated_rcpp$residual_var_mean -
    correlated_r$residual_var_mean) < 0.05
)

# A spike-and-slab prior learns predictor inclusion and the shared pi.
selection_n <- 60
selection_X <- cbind(
  signal = seq(-2, 2, length.out = selection_n),
  noise1 = sin(seq_len(selection_n)),
  noise2 = cos(seq_len(selection_n)),
  noise3 = rep(c(-1, 0, 1), 20)
)
selection_y <- 1 + 2.5 * selection_X[, "signal"] +
  0.25 * rep(c(-1, 1), 30)
selection_rcpp <- multiple_blm(
  selection_y, selection_X, 10,
  coefficient_prior = "spike_slab",
  pi_alpha = 1,
  pi_beta = 1,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1200,
  burnin = 400,
  thin = 2,
  seed = 77,
  version = "Rcpp"
)
stopifnot(
  identical(
    tail(names(selection_rcpp), 5),
    c(
      "inclusion_probability", "pi_mean", "pi_var",
      "inclusion_samples", "pi_samples"
    )
  ),
  identical(dim(selection_rcpp$inclusion_samples), c(400L, 4L)),
  all(selection_rcpp$inclusion_samples %in% 0:1),
  all(
    selection_rcpp$coefficient_samples[
      selection_rcpp$inclusion_samples == 0
    ] == 0
  ),
  selection_rcpp$inclusion_probability["signal"] > 0.9,
  max(selection_rcpp$inclusion_probability[-1]) < 0.3,
  all(selection_rcpp$pi_samples > 0 & selection_rcpp$pi_samples < 1),
  identical(selection_rcpp$pi_mean, mean(selection_rcpp$pi_samples)),
  identical(selection_rcpp$pi_var, var(selection_rcpp$pi_samples))
)

selection_r <- multiple_blm(
  selection_y, selection_X, 10,
  coefficient_prior = "spike_slab",
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1200,
  burnin = 400,
  thin = 2,
  seed = 77,
  version = "R"
)
stopifnot(
  max(abs(
    selection_rcpp$inclusion_probability -
      selection_r$inclusion_probability
  )) < 0.1,
  abs(selection_rcpp$pi_mean - selection_r$pi_mean) < 0.1
)

# The default global-local prior is Strawderman-Berger and learns both local
# variances and the global tau squared. Coefficients remain on their original
# scale after predictor standardization.
global_local_n <- 80
global_local_X <- cbind(
  signal = seq(-20, 20, length.out = global_local_n),
  noise1 = sin(seq_len(global_local_n)) / 10,
  noise2 = cos(seq_len(global_local_n)) * 5,
  noise3 = rep(c(-2, -1, 1, 2), 20)
)
global_local_y <- 0.75 + 0.4 * global_local_X[, "signal"] +
  0.15 * rep(c(-1, 1), 40)
global_local_rcpp <- multiple_blm(
  global_local_y, global_local_X,
  coefficient_prior = "global_local",
  global_scale = 1,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1600,
  burnin = 600,
  thin = 2,
  seed = 109,
  version = "Rcpp"
)
stopifnot(
  identical(dim(global_local_rcpp$coefficient_samples), c(500L, 4L)),
  identical(dim(global_local_rcpp$local_var_samples), c(500L, 4L)),
  identical(
    names(global_local_rcpp$local_var_mean),
    colnames(global_local_X)
  ),
  identical(
    names(global_local_rcpp$local_var_var),
    colnames(global_local_X)
  ),
  all(is.finite(global_local_rcpp$local_var_samples)),
  all(global_local_rcpp$local_var_samples > 0),
  all(is.finite(global_local_rcpp$tau_sq_samples)),
  all(global_local_rcpp$tau_sq_samples > 0),
  abs(global_local_rcpp$coefficient_mean["signal"] - 0.4) < 0.05,
  max(abs(global_local_rcpp$coefficient_mean[-1])) < 0.15,
  identical(
    global_local_rcpp$tau_sq_mean,
    mean(global_local_rcpp$tau_sq_samples)
  ),
  identical(global_local_rcpp$local_shape, c(a = 1, b = 0.5))
)

global_local_r <- multiple_blm(
  global_local_y, global_local_X,
  coefficient_prior = "global_local",
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1600,
  burnin = 600,
  thin = 2,
  seed = 109,
  version = "R"
)
stopifnot(
  max(abs(
    global_local_rcpp$coefficient_mean - global_local_r$coefficient_mean
  )) < 0.1,
  abs(global_local_rcpp$residual_var_mean -
    global_local_r$residual_var_mean) < 0.1
)

# The horseshoe is obtained using beta-prime shapes a = b = 1/2.
horseshoe_fit <- multiple_blm(
  global_local_y, global_local_X,
  coefficient_prior = "global_local",
  local_shape = c(a = 0.5, b = 0.5),
  residual_shape = 2,
  residual_scale = 1,
  iterations = 600,
  burnin = 200,
  seed = 111
)
stopifnot(
  identical(horseshoe_fit$local_shape, c(a = 0.5, b = 0.5)),
  all(horseshoe_fit$local_var_samples > 0),
  all(horseshoe_fit$tau_sq_samples > 0)
)

# A known residual variance still uses Gibbs sampling for global-local priors.
fixed_global_local <- multiple_blm(
  global_local_y, global_local_X,
  coefficient_prior = "global_local",
  residual_var = 0.25,
  iterations = 400,
  burnin = 200,
  seed = 110
)
stopifnot(
  identical(dim(fixed_global_local$coefficient_samples), c(200L, 4L)),
  all(fixed_global_local$residual_var_samples == 0.25),
  identical(fixed_global_local$residual_var_mean, 0.25),
  identical(fixed_global_local$residual_var_var, 0)
)

# The R and Rcpp GIG entry points match known moments and remain finite for
# challenging parameter combinations.
gig_lambda <- 1
gig_chi <- 2
gig_psi <- 3
gig_argument <- sqrt(gig_chi * gig_psi)
gig_expected_mean <- sqrt(gig_chi / gig_psi) *
  besselK(gig_argument, gig_lambda + 1) /
  besselK(gig_argument, gig_lambda)
set.seed(301)
gig_r <- blm:::.draw_gig(20000L, gig_lambda, gig_chi, gig_psi)
set.seed(302)
gig_rcpp <- blm:::draw_gig_rcpp_cpp(
  20000L, gig_lambda, gig_chi, gig_psi
)
stopifnot(
  abs(mean(gig_r) - gig_expected_mean) / gig_expected_mean < 0.03,
  abs(mean(gig_rcpp) - gig_expected_mean) / gig_expected_mean < 0.03,
  all(is.finite(blm:::.draw_gig(100L, -0.4, 1e-8, 1e3))),
  all(is.finite(blm:::draw_gig_rcpp_cpp(100L, 5, 1e3, 1e-4)))
)

# Both implementations report progress at 10-percent intervals.
for (sampler_version in c("Rcpp", "R")) {
  progress_amounts <- integer()
  progress_iterations <- integer()
  progress_callback <- function(amount, iteration) {
    progress_amounts <<- c(progress_amounts, amount)
    progress_iterations <<- c(progress_iterations, iteration)
  }
  sampler <- if (sampler_version == "Rcpp") {
    blm:::.blm_gibbs_rcpp
  } else {
    blm:::.blm_gibbs
  }
  invisible(sampler(
    y = y,
    x = x,
    prior_var = c(10, 10),
    residual_shape = 2,
    residual_scale = 1,
    iterations = 100,
    burnin = 20,
    thin = 1,
    seed = 42,
    progress_callback = progress_callback
  ))
  stopifnot(
    isTRUE(all.equal(progress_iterations, seq.int(10L, 100L, by = 10L))),
    isTRUE(all.equal(progress_amounts, rep.int(10L, 10L)))
  )
}

# Invalid designs, priors, and residual specifications are rejected.
stopifnot(
  inherits(try(multiple_blm(y, x[, 1], 10, 1), silent = TRUE), "try-error"),
  inherits(try(multiple_blm(y[-1], x, 10, 1), silent = TRUE), "try-error"),
  inherits(try(multiple_blm(y, x, c(10, 0), 1), silent = TRUE), "try-error"),
  inherits(try(multiple_blm(y, x, c(10, 10, 10), 1), silent = TRUE), "try-error"),
  inherits(try(multiple_blm(y, x, 10), silent = TRUE), "try-error"),
  inherits(
    try(
      multiple_blm(y, x, 10, residual_var = 1, nchains = 2),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(
        y, x, prior_var = 10,
        coefficient_prior = "global_local",
        residual_shape = 2,
        residual_scale = 1
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(
        y, x,
        coefficient_prior = "global_local",
        global_scale = 0,
        residual_shape = 2,
        residual_scale = 1
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(
        y, x,
        coefficient_prior = "global_local",
        local_shape = c(1, 0),
        residual_shape = 2,
        residual_scale = 1
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(y, x, 10, residual_var = 1, standardize = NA),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(
        y, cbind(x, constant = 1), 10, residual_var = 1
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(
        y, x, 10,
        residual_shape = 2,
        residual_scale = 1,
        nchains = 0
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(multiple_blm(y, x, 10, residual_var = 1, version = "C"), silent = TRUE),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(
        y, x, 10,
        residual_var = 1,
        coefficient_prior = "spike_slab"
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(
        y, x, 10,
        coefficient_prior = "spike_slab",
        pi_alpha = 0,
        residual_shape = 2,
        residual_scale = 1
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(y, x, 10, residual_var = 1, verbose = NA),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(
        y, x, 10,
        residual_shape = 2, residual_scale = 1,
        iterations = 10, burnin = 9
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(
      multiple_blm(y, x, 10, residual_var = 1, residual_shape = 2),
      silent = TRUE
    ),
    "try-error"
  )
)
