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
  residual_var = 1
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
vector_prior_fit <- multiple_blm(y, x, c(10, 10), residual_var = 1)
stopifnot(isTRUE(all.equal(known_fit, vector_prior_fit)))

# Data-frame inputs retain their predictor names.
data_frame_fit <- multiple_blm(y, as.data.frame(x), 10, residual_var = 1)
stopifnot(isTRUE(all.equal(known_fit, data_frame_fit)))

# With one predictor, multiple_blm agrees with simple_blm.
simple_y <- 1 + 2 * x[, "first"]
simple_fit <- simple_blm(simple_y, x[, "first"], 10, residual_var = 1)
one_predictor_fit <- multiple_blm(
  simple_y, x[, "first", drop = FALSE], 10, residual_var = 1
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

# Both implementations display an iteration counter when requested.
for (sampler_version in c("Rcpp", "R")) {
  progress_output <- capture.output(
    multiple_blm(
      y, x, 10,
      residual_shape = 2,
      residual_scale = 1,
      iterations = 6,
      burnin = 2,
      seed = 42,
      version = sampler_version,
      verbose = TRUE
    )
  )
  stopifnot(grepl(
    "Gibbs iteration 6/6",
    paste(progress_output, collapse = "\n"),
    fixed = TRUE
  ))
}

# Invalid designs, priors, and residual specifications are rejected.
stopifnot(
  inherits(try(multiple_blm(y, x[, 1], 10, 1), silent = TRUE), "try-error"),
  inherits(try(multiple_blm(y[-1], x, 10, 1), silent = TRUE), "try-error"),
  inherits(try(multiple_blm(y, x, c(10, 0), 1), silent = TRUE), "try-error"),
  inherits(try(multiple_blm(y, x, c(10, 10, 10), 1), silent = TRUE), "try-error"),
  inherits(try(multiple_blm(y, x, 10), silent = TRUE), "try-error"),
  inherits(
    try(multiple_blm(y, x, 10, residual_var = 1, version = "C"), silent = TRUE),
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
