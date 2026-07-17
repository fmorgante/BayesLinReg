library(blm)

fit <- simple_blm(
  y = c(3, 5, 7, 9, 11),
  x = 1:5,
  prior_var = 10,
  residual_var = 1
)

stopifnot(
  identical(names(fit), c("posterior_mean", "posterior_var")),
  isTRUE(all.equal(fit$posterior_var, 1 / 10.1)),
  isTRUE(all.equal(fit$posterior_mean, 20 / 10.1))
)

# Centering makes the result invariant to shifts in x and y.
shifted_fit <- simple_blm(
  y = c(3, 5, 7, 9, 11) + 100,
  x = 1:5 + 50,
  prior_var = 10,
  residual_var = 1
)
stopifnot(isTRUE(all.equal(fit, shifted_fit)))

# Invalid variances and incompatible inputs are rejected.
stopifnot(
  inherits(try(simple_blm(1:3, 1:2, 1, 1), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 0, 1), silent = TRUE), "try-error"),
  inherits(try(simple_blm(1:3, 1:3, 1, NA_real_), silent = TRUE), "try-error")
)
