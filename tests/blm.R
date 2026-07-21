library(BayesLinReg)

x <- cbind(
  first = 1:5,
  second = c(0, 1, 0, 1, 0)
)
y <- 1 + 2 * x[, "first"] - 3 * x[, "second"]

normal_eta <- function(X, var_shape = 2, var_scale = 10,
                       standardize = TRUE) {
  list(
    X = X, model = "Normal", var_shape = var_shape,
    var_scale = var_scale, standardize = standardize
  )
}

# A Normal block samples one shared coefficient variance.
known_fit <- blm(
  y,
  ETA = normal_eta(x, standardize = FALSE),
  residual_var = 1,
  iterations = 1000,
  burnin = 400,
  thin = 2,
  seed = 41
)
known_block <- known_fit$ETA$ETA1
stopifnot(
  identical(names(known_fit$ETA), "ETA1"),
  identical(known_block$model, "Normal"),
  identical(known_block$standardize, FALSE),
  identical(known_block$var_shape, 2),
  identical(known_block$var_scale, 10),
  isTRUE(all.equal(
    known_block$coefficient_var,
    diag(known_block$coefficient_cov)
  )),
  identical(dim(known_block$coefficient_samples), c(300L, 2L)),
  length(known_block$normal_var_samples) == 300L,
  all(known_block$normal_var_samples > 0),
  identical(known_block$normal_var_mean, mean(known_block$normal_var_samples)),
  all(known_fit$residual_var_samples == 1)
)

# Matrix and data-frame inputs use the same fitting engine.
vector_fit <- blm(
  y,
  ETA = list(
    X = x, model = "Normal", var_shape = 2, var_scale = 10,
    standardize = FALSE
  ),
  residual_var = 1, iterations = 1000, burnin = 400, thin = 2, seed = 41
)
data_frame_fit <- blm(
  y,
  ETA = list(
    X = as.data.frame(x), model = "Normal", var_shape = 2, var_scale = 10,
    standardize = FALSE
  ),
  residual_var = 1, iterations = 1000, burnin = 400, thin = 2, seed = 41
)
stopifnot(
  isTRUE(all.equal(known_fit, vector_fit)),
  isTRUE(all.equal(known_fit, data_frame_fit))
)

# A predictor vector and a one-column matrix use the same fitting engine.
simple_y <- 1 + 2 * x[, "first"]
simple_fit <- blm(
  simple_y,
  ETA = list(
    X = x[, "first"], model = "Normal", var_shape = 2, var_scale = 10,
    standardize = FALSE
  ),
  residual_var = 1, iterations = 1000, burnin = 400, thin = 2, seed = 42
)
one_predictor_fit <- blm(
  simple_y,
  ETA = normal_eta(x[, "first", drop = FALSE], standardize = FALSE),
  residual_var = 1, iterations = 1000, burnin = 400, thin = 2, seed = 42
)
stopifnot(
  isTRUE(all.equal(
    unname(one_predictor_fit$ETA$ETA1$coefficient_mean),
    unname(simple_fit$ETA$ETA1$coefficient_mean)
  )),
  isTRUE(all.equal(
    drop(one_predictor_fit$ETA$ETA1$coefficient_cov),
    drop(simple_fit$ETA$ETA1$coefficient_cov)
  )),
  isTRUE(all.equal(one_predictor_fit$intercept_mean, simple_fit$intercept_mean)),
  isTRUE(all.equal(one_predictor_fit$intercept_var, simple_fit$intercept_var))
)

# Standardization is block-specific, and returned coefficients use the
# original predictor scale.
predictor_sd <- apply(x, 2L, stats::sd)
working_x <- sweep(x, 2L, predictor_sd, FUN = "/")
manual_fit <- blm(
  y,
  ETA = normal_eta(working_x, standardize = FALSE),
  residual_var = 1, iterations = 1000, burnin = 400, thin = 2, seed = 43
)
automatic_fit <- blm(
  y, ETA = normal_eta(x), residual_var = 1,
  iterations = 1000, burnin = 400, thin = 2, seed = 43
)
stopifnot(
  isTRUE(all.equal(
    automatic_fit$ETA$ETA1$coefficient_mean,
    manual_fit$ETA$ETA1$coefficient_mean / predictor_sd
  )),
  isTRUE(all.equal(
    automatic_fit$ETA$ETA1$coefficient_cov,
    manual_fit$ETA$ETA1$coefficient_cov / outer(predictor_sd, predictor_sd)
  )),
  isTRUE(all.equal(automatic_fit$intercept_mean, manual_fit$intercept_mean)),
  isTRUE(all.equal(automatic_fit$intercept_var, manual_fit$intercept_var))
)

# Learned residual variance uses Gibbs sampling and remains reproducible.
learned_fit <- blm(
  y,
  ETA = normal_eta(x),
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1500,
  burnin = 500,
  thin = 2,
  seed = 42
)
repeated_fit <- blm(
  y,
  ETA = normal_eta(x),
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1500,
  burnin = 500,
  thin = 2,
  seed = 42
)
stopifnot(
  isTRUE(all.equal(learned_fit, repeated_fit)),
  identical(dim(learned_fit$ETA$ETA1$coefficient_samples), c(500L, 2L)),
  all(learned_fit$residual_var_samples > 0),
  isTRUE(all.equal(
    learned_fit$ETA$ETA1$coefficient_mean,
    colMeans(learned_fit$ETA$ETA1$coefficient_samples)
  ))
)

# The R and Rcpp implementations target the same Normal posterior.
r_fit <- blm(
  y,
  ETA = normal_eta(x),
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1500,
  burnin = 500,
  thin = 2,
  seed = 42,
  version = "R"
)
stopifnot(
  max(abs(
    r_fit$ETA$ETA1$coefficient_mean -
      learned_fit$ETA$ETA1$coefficient_mean
  )) < 0.2,
  abs(r_fit$residual_var_mean - learned_fit$residual_var_mean) < 0.2
)

# Multiple ETA blocks can use distinct models and hyperparameters.
multi_n <- 80
multi_X <- cbind(
  signal = seq(-2, 2, length.out = multi_n),
  nuisance = rep(c(-1, 1), 40),
  noise1 = sin(seq_len(multi_n)),
  noise2 = cos(seq_len(multi_n))
)
multi_y <- 1 + 2.5 * multi_X[, "signal"] +
  0.15 * rep(c(-1, 1), 40)
multi_eta <- list(
  fixed = list(
    X = multi_X[, "signal", drop = FALSE],
    model = "Normal",
    var_shape = 2,
    var_scale = 20
  ),
  selection = list(
    X = multi_X[, c("nuisance", "noise1")],
    model = "SpikeSlab",
    slab_shape = 2,
    slab_scale = 10,
    pi = c(a = 1, b = 2)
  ),
  shrinkage = list(
    X = multi_X[, "noise2", drop = FALSE],
    model = "GlobalLocal",
    local_shape = c(a = 1, b = 0.5),
    global_scale = 0.5
  )
)
multi_rcpp <- blm(
  multi_y, ETA = multi_eta,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1200,
  burnin = 400,
  thin = 2,
  seed = 77,
  version = "Rcpp"
)
multi_r <- blm(
  multi_y, ETA = multi_eta,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1200,
  burnin = 400,
  thin = 2,
  seed = 77,
  version = "R"
)
stopifnot(
  identical(names(multi_rcpp$ETA), c("fixed", "selection", "shrinkage")),
  identical(
    unname(vapply(multi_rcpp$ETA, `[[`, character(1), "model")),
    c("Normal", "SpikeSlab", "GlobalLocal")
  ),
  identical(dim(multi_rcpp$ETA$fixed$coefficient_samples), c(400L, 1L)),
  length(multi_rcpp$ETA$fixed$normal_var_samples) == 400L,
  all(multi_rcpp$ETA$fixed$normal_var_samples > 0),
  identical(dim(multi_rcpp$ETA$selection$inclusion_samples), c(400L, 2L)),
  length(multi_rcpp$ETA$selection$pi_samples) == 400L,
  length(multi_rcpp$ETA$selection$slab_var_samples) == 400L,
  all(multi_rcpp$ETA$selection$slab_var_samples > 0),
  identical(dim(multi_rcpp$ETA$shrinkage$local_var_samples), c(400L, 1L)),
  length(multi_rcpp$ETA$shrinkage$tau_sq_samples) == 400L,
  identical(multi_rcpp$ETA$selection$pi, c(a = 1, b = 2)),
  identical(multi_rcpp$ETA$selection$slab_shape, 2),
  identical(multi_rcpp$ETA$selection$slab_scale, 10),
  identical(multi_rcpp$ETA$shrinkage$global_scale, 0.5),
  abs(multi_rcpp$ETA$fixed$coefficient_mean["signal"] - 2.5) < 0.1,
  max(abs(
    multi_rcpp$ETA$fixed$coefficient_mean -
      multi_r$ETA$fixed$coefficient_mean
  )) < 0.1,
  abs(
    multi_rcpp$ETA$selection$slab_var_mean -
      multi_r$ETA$selection$slab_var_mean
  ) < 5,
  abs(multi_rcpp$residual_var_mean - multi_r$residual_var_mean) < 0.1
)

# Summary-only fits match stored-draw summaries without returning samples.
for (sampler_version in c("Rcpp", "R")) {
  stored_fit <- if (sampler_version == "Rcpp") multi_rcpp else multi_r
  summary_fit <- blm(
    multi_y, ETA = multi_eta,
    residual_shape = 2,
    residual_scale = 1,
    iterations = 1200,
    burnin = 400,
    thin = 2,
    seed = 77,
    version = sampler_version,
    store_samples = FALSE
  )
  variance_only_fit <- blm(
    multi_y, ETA = multi_eta,
    residual_shape = 2,
    residual_scale = 1,
    iterations = 1200,
    burnin = 400,
    thin = 2,
    seed = 77,
    version = sampler_version,
    store_samples = FALSE,
    store_coefficient_cov = FALSE
  )
  stopifnot(
    identical(stored_fit$store_samples, TRUE),
    identical(summary_fit$store_samples, FALSE),
    is.null(summary_fit$intercept_samples),
    is.null(summary_fit$residual_var_samples),
    is.null(summary_fit$chain_id),
    is.null(summary_fit$ETA$fixed$coefficient_samples),
    is.null(summary_fit$ETA$fixed$normal_var_samples),
    is.null(summary_fit$ETA$selection$inclusion_samples),
    is.null(summary_fit$ETA$selection$pi_samples),
    is.null(summary_fit$ETA$selection$slab_var_samples),
    is.null(summary_fit$ETA$shrinkage$local_var_samples),
    is.null(summary_fit$ETA$shrinkage$tau_sq_samples),
    identical(summary_fit$store_coefficient_cov, TRUE),
    identical(variance_only_fit$store_coefficient_cov, FALSE),
    is.null(variance_only_fit$ETA$fixed$coefficient_cov),
    is.null(variance_only_fit$ETA$selection$coefficient_cov),
    is.null(variance_only_fit$ETA$shrinkage$coefficient_cov),
    isTRUE(all.equal(
      stored_fit$ETA$fixed$coefficient_mean,
      summary_fit$ETA$fixed$coefficient_mean,
      tolerance = 1e-10
    )),
    isTRUE(all.equal(
      stored_fit$ETA$selection$coefficient_cov,
      summary_fit$ETA$selection$coefficient_cov,
      tolerance = 1e-10
    )),
    isTRUE(all.equal(
      stored_fit$ETA$selection$coefficient_var,
      diag(stored_fit$ETA$selection$coefficient_cov),
      tolerance = 1e-10
    )),
    isTRUE(all.equal(
      stored_fit$ETA$selection$coefficient_var,
      variance_only_fit$ETA$selection$coefficient_var,
      tolerance = 1e-10
    )),
    isTRUE(all.equal(
      stored_fit$ETA$selection$inclusion_probability,
      summary_fit$ETA$selection$inclusion_probability,
      tolerance = 1e-10
    )),
    isTRUE(all.equal(
      stored_fit$ETA$shrinkage$local_var_mean,
      summary_fit$ETA$shrinkage$local_var_mean,
      tolerance = 1e-10
    )),
    isTRUE(all.equal(
      stored_fit$residual_var_var,
      summary_fit$residual_var_var,
      tolerance = 1e-10
    )),
    object.size(summary_fit) < object.size(stored_fit),
    object.size(variance_only_fit) < object.size(summary_fit)
  )
}

# Covariance storage can be disabled while retaining individual draws.
stored_variance_only_fit <- blm(
  y,
  ETA = normal_eta(x),
  residual_var = 1,
  iterations = 100,
  burnin = 40,
  seed = 78,
  store_coefficient_cov = FALSE
)
stopifnot(
  !is.null(stored_variance_only_fit$ETA$ETA1$coefficient_samples),
  is.null(stored_variance_only_fit$ETA$ETA1$coefficient_cov),
  isTRUE(all.equal(
    stored_variance_only_fit$ETA$ETA1$coefficient_var,
    apply(
      stored_variance_only_fit$ETA$ETA1$coefficient_samples,
      2L,
      stats::var
    )
  ))
)

# GlobalLocal defaults to Strawderman-Berger and can recover the horseshoe.
global_fit <- blm(
  multi_y,
  ETA = list(X = multi_X, model = "GlobalLocal"),
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1000,
  burnin = 400,
  thin = 2,
  seed = 109
)
horseshoe_fit <- blm(
  multi_y,
  ETA = list(
    X = multi_X,
    model = "GlobalLocal",
    local_shape = c(a = 0.5, b = 0.5)
  ),
  residual_var = 0.25,
  iterations = 500,
  burnin = 200,
  seed = 110
)
stopifnot(
  identical(global_fit$ETA$ETA1$local_shape, c(a = 1, b = 0.5)),
  all(global_fit$ETA$ETA1$local_var_samples > 0),
  all(global_fit$ETA$ETA1$tau_sq_samples > 0),
  identical(horseshoe_fit$ETA$ETA1$local_shape, c(a = 0.5, b = 0.5)),
  all(horseshoe_fit$residual_var_samples == 0.25),
  identical(horseshoe_fit$residual_var_var, 0)
)

# A single SpikeSlab block accepts a vector or matrix of predictors.
spike_fit <- blm(
  multi_y,
  ETA = list(
    X = multi_X[, c("signal", "noise1")],
    model = "SpikeSlab",
    slab_shape = 3,
    slab_scale = 4,
    pi = c(a = 1, b = 2)
  ),
  residual_shape = 2,
  residual_scale = 1,
  iterations = 500,
  burnin = 200,
  seed = 111
)
stopifnot(
  identical(spike_fit$ETA$ETA1$model, "SpikeSlab"),
  identical(spike_fit$ETA$ETA1$pi, c(a = 1, b = 2)),
  identical(spike_fit$ETA$ETA1$slab_shape, 3),
  identical(spike_fit$ETA$ETA1$slab_scale, 4),
  all(spike_fit$ETA$ETA1$slab_var_samples > 0),
  identical(
    spike_fit$ETA$ETA1$slab_var_mean,
    mean(spike_fit$ETA$ETA1$slab_var_samples)
  ),
  identical(dim(spike_fit$ETA$ETA1$coefficient_samples), c(300L, 2L))
)

# SpikeSlab supports the same fixed residual-variance interface as other priors.
fixed_spike_rcpp <- blm(
  multi_y,
  ETA = list(
    X = multi_X[, c("signal", "noise1")],
    model = "SpikeSlab",
    slab_shape = 3,
    slab_scale = 4,
    pi = c(a = 1, b = 2)
  ),
  residual_var = 0.25,
  iterations = 500,
  burnin = 200,
  seed = 112,
  version = "Rcpp"
)
fixed_spike_r <- blm(
  multi_y,
  ETA = list(
    X = multi_X[, c("signal", "noise1")],
    model = "SpikeSlab",
    slab_shape = 3,
    slab_scale = 4,
    pi = c(a = 1, b = 2)
  ),
  residual_var = 0.25,
  iterations = 500,
  burnin = 200,
  seed = 112,
  version = "R"
)
stopifnot(
  all(fixed_spike_rcpp$residual_var_samples == 0.25),
  identical(fixed_spike_rcpp$residual_var_mean, 0.25),
  identical(fixed_spike_rcpp$residual_var_var, 0),
  all(fixed_spike_rcpp$ETA$ETA1$slab_var_samples > 0),
  all(fixed_spike_r$residual_var_samples == 0.25),
  max(abs(
    fixed_spike_rcpp$ETA$ETA1$coefficient_mean -
      fixed_spike_r$ETA$ETA1$coefficient_mean
  )) < 0.2,
  abs(
    fixed_spike_rcpp$ETA$ETA1$pi_mean -
      fixed_spike_r$ETA$ETA1$pi_mean
  ) < 0.2
)

# The R and Rcpp GIG entry points match a known moment and challenging cases.
gig_lambda <- 1
gig_chi <- 2
gig_psi <- 3
gig_argument <- sqrt(gig_chi * gig_psi)
gig_expected_mean <- sqrt(gig_chi / gig_psi) *
  besselK(gig_argument, gig_lambda + 1) /
  besselK(gig_argument, gig_lambda)
set.seed(301)
gig_r <- BayesLinReg:::.draw_gig(20000L, gig_lambda, gig_chi, gig_psi)
set.seed(302)
gig_rcpp <- BayesLinReg:::draw_gig_rcpp_cpp(
  20000L, gig_lambda, gig_chi, gig_psi
)
stopifnot(
  abs(mean(gig_r) - gig_expected_mean) / gig_expected_mean < 0.03,
  abs(mean(gig_rcpp) - gig_expected_mean) / gig_expected_mean < 0.03,
  all(is.finite(BayesLinReg:::.draw_gig(100L, -0.4, 1e-8, 1e3))),
  all(is.finite(BayesLinReg:::draw_gig_rcpp_cpp(100L, 5, 1e3, 1e-4)))
)

# Both low-level implementations report progress at 10-percent intervals.
for (sampler_version in c("Rcpp", "R")) {
  progress_amounts <- integer()
  progress_iterations <- integer()
  callback <- function(amount, iteration) {
    progress_amounts <<- c(progress_amounts, amount)
    progress_iterations <<- c(progress_iterations, iteration)
  }
  sampler <- if (sampler_version == "Rcpp") {
    BayesLinReg:::.blm_gibbs_rcpp
  } else {
    BayesLinReg:::.blm_gibbs
  }
  invisible(sampler(
    y = y,
    x = x,
    residual_shape = 2,
    residual_scale = 1,
    iterations = 100,
    burnin = 20,
    thin = 1,
    seed = 42,
    normal_shape = 2,
    normal_scale = 10,
    progress_callback = callback
  ))
  stopifnot(
    isTRUE(all.equal(progress_iterations, seq.int(10L, 100L, by = 10L))),
    isTRUE(all.equal(progress_amounts, rep.int(10L, 10L)))
  )
}

# Chain combination retains block-specific hyperparameter matrices.
mock_chain <- list(
  coefficient_samples = matrix(1:8, nrow = 2),
  intercept_samples = c(1, 2),
  residual_var_samples = c(3, 4),
  normal_var_samples = matrix(c(4, NA, NA, 5, NA, NA), 2, 3, byrow = TRUE),
  inclusion_samples = matrix(1L, 2, 4),
  pi_samples = matrix(c(NA, 0.4, NA, NA, 0.5, NA), 2, 3, byrow = TRUE),
  slab_var_samples = matrix(c(NA, 4, NA, NA, 5, NA), 2, 3, byrow = TRUE),
  local_var_samples = matrix(1, 2, 4),
  tau_sq_samples = matrix(c(NA, NA, 1, NA, NA, 2), 2, 3, byrow = TRUE)
)
mock_combined <- BayesLinReg:::.combine_blm_chains(
  list(mock_chain, mock_chain),
  block_model = 0:2
)
stopifnot(
  identical(mock_combined$chain_id, c(1L, 1L, 2L, 2L)),
  identical(dim(mock_combined$normal_var_samples), c(4L, 3L)),
  identical(dim(mock_combined$pi_samples), c(4L, 3L)),
  identical(dim(mock_combined$slab_var_samples), c(4L, 3L)),
  identical(dim(mock_combined$tau_sq_samples), c(4L, 3L))
)

# Real multisession tests are enabled outside restricted check environments.
parallel_test_flags <- Sys.getenv(c(
  "BAYESLINREG_TEST_FUTURE",
  "BLM_TEST_FUTURE"
))
if (any(parallel_test_flags == "true")) {
  parallel_fit <- blm(
    multi_y, ETA = multi_eta,
    residual_shape = 2,
    residual_scale = 1,
    iterations = 100,
    burnin = 20,
    seed = 123,
    nchains = 2,
    verbose = TRUE
  )
  repeated_parallel_fit <- blm(
    multi_y, ETA = multi_eta,
    residual_shape = 2,
    residual_scale = 1,
    iterations = 100,
    burnin = 20,
    seed = 123,
    nchains = 2,
    verbose = TRUE
  )
  stopifnot(
    isTRUE(all.equal(parallel_fit, repeated_parallel_fit)),
    identical(parallel_fit$chain_id, rep.int(1:2, c(80L, 80L))),
    identical(dim(parallel_fit$ETA$selection$pi_samples), NULL),
    length(parallel_fit$ETA$selection$pi_samples) == 160L,
    length(parallel_fit$ETA$selection$slab_var_samples) == 160L,
    length(parallel_fit$ETA$shrinkage$tau_sq_samples) == 160L
  )
}

# Invalid ETA specifications and model parameters are rejected.
invalid_calls <- list(
  function() blm(y, residual_var = 1),
  function() blm(
    y, ETA = normal_eta(x), X = x, residual_var = 1
  ),
  function() blm(y, ETA = list(), residual_var = 1),
  function() blm(y, ETA = list(X = x), residual_var = 1),
  function() blm(
    y, ETA = list(X = x, model = "normal", var_scale = 10), residual_var = 1
  ),
  function() blm(
    y, ETA = list(X = x, model = "Normal", prior_var = 10), residual_var = 1
  ),
  function() blm(
    y, ETA = list(X = x, model = "Normal", var_shape = 0), residual_var = 1
  ),
  function() blm(
    y, ETA = list(X = x, model = "Normal", var_scale = 0), residual_var = 1
  ),
  function() blm(
    y,
    ETA = list(
      X = cbind(x, constant = 1), model = "Normal", var_scale = 10
    ),
    residual_var = 1
  ),
  function() blm(
    y,
    ETA = list(X = x, model = "SpikeSlab", var = 10),
    residual_var = 1
  ),
  function() blm(
    y,
    ETA = list(
      X = x, model = "SpikeSlab", pi_alpha = 1, pi_beta = 1
    ),
    residual_shape = 2,
    residual_scale = 1
  ),
  function() blm(
    y,
    ETA = list(X = x, model = "SpikeSlab", pi = c(1, 0)),
    residual_shape = 2,
    residual_scale = 1
  ),
  function() blm(
    y,
    ETA = list(X = x, model = "SpikeSlab", slab_shape = 0),
    residual_shape = 2,
    residual_scale = 1
  ),
  function() blm(
    y,
    ETA = list(X = x, model = "SpikeSlab", slab_scale = 0),
    residual_shape = 2,
    residual_scale = 1
  ),
  function() blm(
    y,
    ETA = list(
      X = x, model = "GlobalLocal", local_shape = c(1, 0)
    ),
    residual_shape = 2,
    residual_scale = 1
  ),
  function() blm(
    y, ETA = normal_eta(x), residual_shape = 2, residual_scale = 1,
    nchains = 0
  ),
  function() blm(
    y, ETA = normal_eta(x), residual_var = 1, residual_shape = 2
  ),
  function() blm(
    y, ETA = normal_eta(x), residual_var = 1, store_samples = NA
  ),
  function() blm(
    y, ETA = normal_eta(x), residual_var = 1,
    store_coefficient_cov = NA
  )
)
stopifnot(all(vapply(
  invalid_calls,
  function(call) inherits(try(call(), silent = TRUE), "try-error"),
  logical(1)
)))
