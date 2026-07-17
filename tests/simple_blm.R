library(blm)

fit <- simple_blm(
  y = c(3, 5, 7, 9, 11),
  x = 1:5,
  prior_var = 10,
  residual_var = 1
)

stopifnot(
  identical(
    names(fit),
    c("slope_mean", "slope_var", "intercept_mean", "intercept_var")
  ),
  isTRUE(all.equal(fit$slope_var, 1 / 10.1)),
  isTRUE(all.equal(fit$slope_mean, 20 / 10.1)),
  isTRUE(all.equal(fit$intercept_mean, 7 - 3 * (20 / 10.1))),
  isTRUE(all.equal(fit$intercept_var, 1 / 5 + 3^2 / 10.1))
)

# Centering makes the slope posterior invariant to shifts in x and y.
shifted_fit <- simple_blm(
  y = c(3, 5, 7, 9, 11) + 100,
  x = 1:5 + 50,
  prior_var = 10,
  residual_var = 1
)
stopifnot(
  isTRUE(all.equal(fit$slope_mean, shifted_fit$slope_mean)),
  isTRUE(all.equal(fit$slope_var, shifted_fit$slope_var)),
  isTRUE(all.equal(
    shifted_fit$intercept_mean,
    fit$intercept_mean + 100 - 50 * fit$slope_mean
  ))
)

# With an inverse-gamma prior, the residual variance is learned from the data.
learned_fit <- simple_blm(
  y = c(3, 5, 7, 9, 11),
  x = 1:5,
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
      "slope_mean", "slope_var", "intercept_mean", "intercept_var",
      "residual_var_mean", "residual_var_var", "slope_samples",
      "intercept_samples", "residual_var_samples"
    )
  ),
  length(learned_fit$slope_samples) == 500,
  length(learned_fit$intercept_samples) == 500,
  length(learned_fit$residual_var_samples) == 500,
  all(learned_fit$residual_var_samples > 0),
  identical(learned_fit$slope_mean, mean(learned_fit$slope_samples)),
  identical(learned_fit$slope_var, var(learned_fit$slope_samples)),
  identical(
    learned_fit$residual_var_mean,
    mean(learned_fit$residual_var_samples)
  )
)

learned_shifted_fit <- simple_blm(
  y = c(3, 5, 7, 9, 11) + 100,
  x = 1:5 + 50,
  prior_var = 10,
  residual_shape = 2,
  residual_scale = 1,
  iterations = 1500,
  burnin = 500,
  thin = 2,
  seed = 42
)
stopifnot(
  identical(learned_fit$slope_samples, learned_shifted_fit$slope_samples),
  identical(
    learned_fit$residual_var_samples,
    learned_shifted_fit$residual_var_samples
  ),
  isTRUE(all.equal(
    learned_shifted_fit$intercept_samples,
    learned_fit$intercept_samples + 100 - 50 * learned_fit$slope_samples
  ))
)

# Invalid variances and incompatible inputs are rejected.
stopifnot(
  inherits(try(simple_blm(1:3, 1:2, 1, 1), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 0, 1), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 1, NA_real_), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 1), silent = TRUE), "try-error"),
  inherits(
    try(
      simple_blm(
        1:3, 1:3, 1,
        residual_shape = 2, residual_scale = 1,
        iterations = 10, burnin = 9
      ),
      silent = TRUE
    ),
    "try-error"
  ),
  inherits(
    try(simple_blm(1:3, 1:3, 1, 1, residual_shape = 2), silent = TRUE),
    "try-error"
  )
)
