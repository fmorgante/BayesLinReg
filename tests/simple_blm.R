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
  residual_scale = 1
)
expected_relative_var <- 1 / 10.1
expected_mean <- 20 / 10.1
expected_shape <- 4
expected_intercept_relative_var <- 1 / 5 + 3^2 * expected_relative_var
expected_scale <- 1 + 0.5 * (
  sum((c(-4, -2, 0, 2, 4) - expected_mean * c(-2, -1, 0, 1, 2))^2) +
    expected_mean^2 / 10
)
stopifnot(
  identical(
    names(learned_fit),
    c(
      "slope_mean", "slope_var", "slope_scale", "slope_df",
      "intercept_mean", "intercept_var", "intercept_scale", "intercept_df",
      "residual_var_shape", "residual_var_scale"
    )
  ),
  isTRUE(all.equal(learned_fit$slope_mean, expected_mean)),
  isTRUE(all.equal(learned_fit$intercept_mean, 7 - 3 * expected_mean)),
  isTRUE(all.equal(learned_fit$residual_var_shape, expected_shape)),
  isTRUE(all.equal(learned_fit$residual_var_scale, expected_scale)),
  isTRUE(all.equal(learned_fit$slope_df, 2 * expected_shape)),
  isTRUE(all.equal(learned_fit$intercept_df, 2 * expected_shape)),
  isTRUE(all.equal(
    learned_fit$slope_var,
    expected_scale / (expected_shape - 1) * expected_relative_var
  )),
  isTRUE(all.equal(
    learned_fit$intercept_var,
    expected_scale / (expected_shape - 1) *
      expected_intercept_relative_var
  )),
  isTRUE(all.equal(
    learned_fit$intercept_scale,
    sqrt(expected_scale / expected_shape * expected_intercept_relative_var)
  ))
)

learned_shifted_fit <- simple_blm(
  y = c(3, 5, 7, 9, 11) + 100,
  x = 1:5 + 50,
  prior_var = 10,
  residual_shape = 2,
  residual_scale = 1
)
stopifnot(
  isTRUE(all.equal(
    learned_fit$slope_mean,
    learned_shifted_fit$slope_mean
  )),
  isTRUE(all.equal(
    learned_shifted_fit$intercept_mean,
    learned_fit$intercept_mean + 100 - 50 * learned_fit$slope_mean
  ))
)

# Invalid variances and incompatible inputs are rejected.
stopifnot(
  inherits(try(simple_blm(1:3, 1:2, 1, 1), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 0, 1), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 1, NA_real_), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 1), silent = TRUE), "try-error"),
  inherits(
    try(simple_blm(1:3, 1:3, 1, 1, residual_shape = 2), silent = TRUE),
    "try-error"
  )
)
