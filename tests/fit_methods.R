library(BayesLinReg)

X <- cbind(first = 1:20, second = rep(c(0, 1), 10))
y <- 1 + 2 * X[, "first"] - X[, "second"]

single_fit <- blm(
  y, ETA = list(X = X, model = "Normal"), residual_var = 1,
  iterations = 100, burnin = 40, seed = 601
)
single_coef <- coef(single_fit)
stopifnot(
  inherits(single_fit, "blm_fit"),
  identical(names(single_coef), c("(Intercept)", "first", "second")),
  isTRUE(all.equal(
    predict(single_fit, X[1:4, , drop = FALSE]),
    unname(drop(single_fit$intercept_mean +
      X[1:4, , drop = FALSE] %*% single_fit$ETA$ETA1$coefficient_mean))
  )),
  isTRUE(all.equal(
    predict(single_fit, X[1, c("second", "first"), drop = FALSE]),
    predict(single_fit, X[1, , drop = FALSE])
  ))
)

summary_fit <- blm(
  y, ETA = list(X = X, model = "Normal"), residual_var = 1,
  iterations = 60, burnin = 20, seed = 603, store_samples = FALSE,
  store_coefficient_cov = FALSE
)
stopifnot(
  length(coef(summary_fit)) == 3L,
  length(predict(summary_fit, X[1:2, , drop = FALSE])) == 2L
)

multi_fit <- blm(
  y,
  ETA = list(
    first_block = list(X = X[, "first", drop = FALSE], model = "Normal"),
    second_block = list(X = X[, "second", drop = FALSE], model = "Normal")
  ),
  residual_var = 1, iterations = 100, burnin = 40, seed = 602
)
stopifnot(
  identical(
    names(coef(multi_fit)),
    c("(Intercept)", "first_block::first", "second_block::second")
  ),
  isTRUE(all.equal(
    predict(multi_fit, list(
      second_block = X[1:3, "second", drop = FALSE],
      first_block = X[1:3, "first", drop = FALSE]
    )),
    unname(multi_fit$intercept_mean +
      X[1:3, "first"] * multi_fit$ETA$first_block$coefficient_mean +
      X[1:3, "second"] * multi_fit$ETA$second_block$coefficient_mean)
  ))
)

constant_X <- cbind(signal = 1:20, constant = 1)
for (standardize in c(TRUE, FALSE)) {
  error <- try(blm(
    y,
    ETA = list(
      X = constant_X, model = "Normal", standardize = standardize
    ),
    residual_var = 1
  ), silent = TRUE)
  stopifnot(
    inherits(error, "try-error"),
    grepl("constant predictor", as.character(error), fixed = TRUE)
  )
}

constant_ss_error <- try(blm_ss(
  n = nrow(constant_X),
  XtX = crossprod(constant_X),
  Xty = drop(crossprod(constant_X, y)),
  X_means = colMeans(constant_X),
  y_mean = mean(y),
  ETA = list(model = "Normal", standardize = FALSE),
  residual_var = 1
), silent = TRUE)
stopifnot(
  inherits(constant_ss_error, "try-error"),
  grepl("constant predictor", as.character(constant_ss_error), fixed = TRUE)
)

invalid_predictions <- list(
  function() predict(single_fit),
  function() predict(single_fit, matrix(1, nrow = 2, ncol = 3)),
  function() predict(single_fit, matrix(NA_real_, nrow = 1, ncol = 2)),
  function() predict(multi_fit, list(first_block = X[, "first"])),
  function() predict(multi_fit, list(
    first_block = X[1:2, "first", drop = FALSE],
    second_block = X[1:3, "second", drop = FALSE]
  ))
)
stopifnot(all(vapply(
  invalid_predictions,
  function(call) inherits(try(call(), silent = TRUE), "try-error"),
  logical(1)
)))
