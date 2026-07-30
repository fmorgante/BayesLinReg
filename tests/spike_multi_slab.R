library(BayesLinReg)

set.seed(701)
n <- 80
X <- matrix(rnorm(n * 6), nrow = n)
colnames(X) <- paste0("x", seq_len(ncol(X)))
y <- drop(1 + X %*% c(2, -1, rep(0, 4)) + rnorm(n))

custom_eta <- list(
  X = X,
  model = "SpikeMultiSlab",
  gamma = c(0, 0.05, 0.5, 2),
  alpha = c(4, 2, 1, 1),
  var_shape = 3,
  var_scale = 2
)

for (sampler_version in c("Rcpp", "R")) {
  stored_fit <- blm(
    y, ETA = custom_eta, residual_var = 1,
    iterations = 500, burnin = 200, seed = 702,
    version = sampler_version
  )
  summary_fit <- blm(
    y, ETA = custom_eta, residual_var = 1,
    iterations = 500, burnin = 200, seed = 702,
    version = sampler_version, store_samples = FALSE
  )
  block <- stored_fit$ETA$ETA1
  stopifnot(
    identical(block$model, "SpikeMultiSlab"),
    identical(block$gamma, c(
      spike = 0, slab_1 = 0.05, slab_2 = 0.5, slab_3 = 2
    )),
    identical(block$alpha, c(
      spike = 4, slab_1 = 2, slab_2 = 1, slab_3 = 1
    )),
    identical(block$var_shape, 3),
    identical(block$var_scale, 2),
    identical(dim(block$component_samples), c(300L, 6L)),
    identical(dim(block$pi_samples), c(300L, 4L)),
    length(block$var_samples) == 300L,
    all(block$var_samples > 0),
    all(block$component_samples %in% 1:4),
    isTRUE(all.equal(unname(rowSums(block$component_probability)), rep(1, 6))),
    isTRUE(all.equal(
      block$inclusion_probability,
      1 - block$component_probability[, "spike"]
    )),
    block$inclusion_probability["x1"] > 0.9,
    isTRUE(all.equal(block$coefficient_mean,
                     summary_fit$ETA$ETA1$coefficient_mean)),
    isTRUE(all.equal(block$component_probability,
                     summary_fit$ETA$ETA1$component_probability)),
    isTRUE(all.equal(block$pi_mean, summary_fit$ETA$ETA1$pi_mean)),
    isTRUE(all.equal(block$pi_var, summary_fit$ETA$ETA1$pi_var)),
    isTRUE(all.equal(block$var_mean, summary_fit$ETA$ETA1$var_mean)),
    isTRUE(all.equal(block$var_var, summary_fit$ETA$ETA1$var_var))
  )

  ss_fit <- blm_ss(
    n, crossprod(X), crossprod(X, y),
    ETA = custom_eta[names(custom_eta) != "X"],
    X_means = colMeans(X), y_mean = mean(y), residual_var = 1,
    iterations = 500, burnin = 200, seed = 702,
    version = sampler_version
  )
  stopifnot(isTRUE(all.equal(stored_fit, ss_fit, tolerance = 1e-8)))

  diagnostics <- assess_convergence(stored_fit, plot = FALSE)
  stopifnot(all(c(
    "pi_ETA1_spike", "pi_ETA1_slab_1", "pi_ETA1_slab_2",
    "pi_ETA1_slab_3", "var_ETA1"
  ) %in% names(diagnostics$rhat)))
}

# Different multi-slab blocks may use different numbers of components.
mixed_fit <- blm(
  y,
  ETA = list(
    three = list(
      X = X[, 1:3], model = "SpikeMultiSlab", gamma = c(0, 0.1, 1)
    ),
    five = list(
      X = X[, 4:6], model = "SpikeMultiSlab",
      gamma = c(0, 0.01, 0.1, 1, 10)
    )
  ),
  residual_var = 1, iterations = 100, burnin = 40, seed = 703
)
stopifnot(
  identical(dim(mixed_fit$ETA$three$pi_samples), c(60L, 3L)),
  identical(dim(mixed_fit$ETA$five$pi_samples), c(60L, 5L))
)

# Mixed prior families preserve identical stored and online summaries after
# precomputing block membership and allocating only family-specific state.
all_prior_eta <- list(
  normal = list(X = X[, 1, drop = FALSE], model = "Normal"),
  selection = list(X = X[, 2, drop = FALSE], model = "SpikeSlab"),
  shrinkage = list(X = X[, 3, drop = FALSE], model = "GlobalLocal"),
  mixture = list(X = X[, 4:6, drop = FALSE], model = "SpikeMultiSlab")
)
all_prior_stored <- blm(
  y, ETA = all_prior_eta, residual_var = 1,
  iterations = 100, burnin = 40, seed = 705
)
all_prior_online <- blm(
  y, ETA = all_prior_eta, residual_var = 1,
  iterations = 100, burnin = 40, seed = 705,
  store_samples = FALSE
)
stopifnot(
  all(vapply(names(all_prior_eta), function(block) {
    isTRUE(all.equal(
      all_prior_stored$ETA[[block]]$coefficient_mean,
      all_prior_online$ETA[[block]]$coefficient_mean
    )) && isTRUE(all.equal(
      all_prior_stored$ETA[[block]]$coefficient_cov,
      all_prior_online$ETA[[block]]$coefficient_cov
    ))
  }, logical(1))),
  isTRUE(all.equal(
    all_prior_stored$ETA$normal$normal_var_mean,
    all_prior_online$ETA$normal$normal_var_mean
  )),
  isTRUE(all.equal(
    all_prior_stored$ETA$selection$inclusion_probability,
    all_prior_online$ETA$selection$inclusion_probability
  )),
  isTRUE(all.equal(
    all_prior_stored$ETA$shrinkage$local_var_mean,
    all_prior_online$ETA$shrinkage$local_var_mean
  )),
  isTRUE(all.equal(
    all_prior_stored$ETA$mixture$component_probability,
    all_prior_online$ETA$mixture$component_probability
  ))
)

# Low-level family-specific state uses only the predictors assigned to that
# family, rather than allocating every family across all predictors.
compact_x <- scale(X, center = TRUE, scale = FALSE)
compact_y <- y - mean(y)
compact_arguments <- list(
  y = compact_y,
  x = compact_x,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 30,
  burnin = 20,
  thin = 1,
  seed = 706,
  block_id = c(1L, 2L, 3L, 4L, 4L, 4L),
  block_model = 0:3,
  normal_shape = rep(2, 4),
  normal_scale = rep(1, 4),
  pi_alpha = rep(1, 4),
  pi_beta = rep(1, 4),
  spike_var_shape = rep(2, 4),
  spike_var_scale = rep(1, 4),
  global_scale = rep(1, 4),
  residual_var = 1,
  local_a = rep(1, 4),
  local_b = rep(0.5, 4),
  multi_gamma = rep(list(c(0, 0.1, 1)), 4),
  multi_pi_alpha = rep(list(rep(1, 3)), 4),
  multi_var_shape = rep(2, 4),
  multi_var_scale = rep(1, 4),
  store_coefficient_cov = FALSE,
  effective_n = n,
  fit_intercept = FALSE,
  intercept_x_mean = rep(0, ncol(X)),
  intercept_y_mean = 0,
  center_observations = FALSE
)
for (sampler in list(
  BayesLinReg:::.blm_gibbs_rcpp,
  BayesLinReg:::.blm_gibbs
)) {
  compact_stored <- do.call(
    sampler,
    c(compact_arguments, list(store_samples = TRUE))
  )
  compact_online <- do.call(
    sampler,
    c(compact_arguments, list(store_samples = FALSE))
  )
  stopifnot(
    identical(dim(compact_stored$coefficient_samples), c(10L, 6L)),
    identical(dim(compact_stored$inclusion_samples), c(10L, 1L)),
    identical(dim(compact_stored$local_var_samples), c(10L, 1L)),
    identical(dim(compact_stored$multi_component_samples), c(10L, 3L)),
    length(compact_online$inclusion_sum) == 1L,
    length(compact_online$local_var_sum) == 1L,
    identical(dim(compact_online$multi_component_sum[[4L]]), c(3L, 3L))
  )
}

# SpikeSlab now uses the same variance-prior field names as other models.
renamed_fit <- blm(
  y,
  ETA = list(
    X = X, model = "SpikeSlab", var_shape = 3, var_scale = 4
  ),
  residual_var = 1, iterations = 60, burnin = 20, seed = 704
)
stopifnot(
  identical(renamed_fit$ETA$ETA1$var_shape, 3),
  identical(renamed_fit$ETA$ETA1$var_scale, 4)
)

invalid_eta <- list(
  list(X = X, model = "SpikeMultiSlab", gamma = c(0, 1, 0.1)),
  list(X = X, model = "SpikeMultiSlab", gamma = c(0, 0, 1)),
  list(
    X = X, model = "SpikeMultiSlab", gamma = c(0, 0.1, 1),
    alpha = c(1, 1)
  ),
  list(X = X, model = "SpikeMultiSlab", alpha = c(1, 1, 1, 0)),
  list(X = X, model = "SpikeMultiSlab", pi = rep(1, 4)),
  list(X = X, model = "SpikeMultiSlab", var_shape = 0),
  list(X = X, model = "SpikeMultiSlab", var_scale = 0),
  list(X = X, model = "SpikeSlab", slab_shape = 2),
  list(X = X, model = "SpikeSlab", slab_scale = 2)
)
stopifnot(all(vapply(invalid_eta, function(ETA) {
  inherits(try(blm(y, ETA, residual_var = 1), silent = TRUE), "try-error")
}, logical(1))))
