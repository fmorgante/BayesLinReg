library(BayesLinReg)

blm_public <- BayesLinReg::blm
blm <- function(..., store_samples = TRUE, store_coefficient_cov = TRUE) {
  blm_public(..., store_samples = store_samples,
             store_coefficient_cov = store_coefficient_cov)
}

X <- cbind(first = 1:20, second = rep(c(0, 1), 10))
y <- 1 + 2 * X[, "first"] - X[, "second"]

set.seed(601)
single_fit <- blm(
  y, ETA = list(X = X, model = "Normal"), residual_var = 1,
  iterations = 100, burnin = 40, compute_pve = TRUE
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

# Credible intervals use fitted-mean draws; prediction intervals integrate the
# normal observation model over the retained posterior draws.
new_X <- X[1:4, , drop = FALSE]
credible <- predict(single_fit, new_X, interval = "credible", level = 0.8)
predictive <- predict(single_fit, new_X, interval = "prediction", level = 0.8)
credible_one_at_a_time <- predict(
  single_fit, new_X, interval = "credible", level = 0.8, chunk_size = 1L
)
predictive_two_at_a_time <- predict(
  single_fit, new_X, interval = "prediction", level = 0.8, chunk_size = 2L
)
prediction_draws <- sweep(
  tcrossprod(new_X, single_fit$ETA$ETA1$coefficient_samples),
  2L, single_fit$intercept_samples, FUN = "+"
)
expected_bounds <- t(apply(
  prediction_draws, 1L, stats::quantile,
  probs = c(0.1, 0.9), names = FALSE
))
stopifnot(
  identical(colnames(credible), c("fit", "lwr", "upr")),
  isTRUE(all.equal(credible[, "fit"], predict(single_fit, new_X))),
  isTRUE(all.equal(
    unname(credible[, c("lwr", "upr")]), unname(expected_bounds)
  )),
  identical(credible_one_at_a_time, credible),
  identical(predictive_two_at_a_time, predictive),
  all(predictive[, "lwr"] < predictive[, "upr"]),
  identical(
    predictive,
    predict(single_fit, new_X, interval = "prediction", level = 0.8)
  )
)

stored_summary <- summary(single_fit, probs = c(0.1, 0.5, 0.9))
printed_summary <- NULL
invisible(utils::capture.output(printed_summary <- print(stored_summary)))
stopifnot(
  inherits(stored_summary, "summary.blm_fit"),
  identical(
    colnames(stored_summary$coefficients),
    c("Mean", "SD", "10%", "50%", "90%")
  ),
  isTRUE(all.equal(
    unname(stored_summary$coefficients["first", "Mean"]),
    unname(single_fit$ETA$ETA1$coefficient_mean["first"])
  )),
  !is.null(stored_summary$hyperparameters),
  !is.null(stored_summary$variance_explained),
  identical(printed_summary, stored_summary)
)

compact_summary <- summary(single_fit, max_coefficients = 1)
top_summary <- summary(
  single_fit, coefficients = "top", max_coefficients = 1,
  rank_by = "absolute_mean"
)
all_summary <- summary(single_fit, coefficients = "all",
                       max_coefficients = 0)
stopifnot(
  is.null(compact_summary$coefficients),
  identical(compact_summary$coefficient_summary, "none"),
  identical(compact_summary$total_coefficients, 2L),
  identical(compact_summary$reported_coefficients, 0L),
  identical(rownames(compact_summary$blocks), "ETA1"),
  identical(compact_summary$blocks$Predictors, 2L),
  identical(rownames(compact_summary$intercept), "(Intercept)"),
  identical(rownames(top_summary$coefficients), "first"),
  identical(top_summary$coefficient_summary, "top"),
  nrow(all_summary$coefficients) == 2L
)

set.seed(604)
spike_fit <- blm(
  y, ETA = list(X = X, model = "SpikeSlab"), residual_var = 1,
  iterations = 80, burnin = 30
)
spike_summary <- summary(
  spike_fit, coefficients = "top", max_coefficients = 1,
  rank_by = "inclusion_probability"
)
expected_top <- names(which.max(
  spike_fit$ETA$ETA1$inclusion_probability
))
stopifnot(
  identical(rownames(spike_summary$coefficients), expected_top),
  isTRUE(all.equal(
    spike_summary$blocks["ETA1", "Expected Included"],
    sum(spike_fit$ETA$ETA1$inclusion_probability)
  ))
)

set.seed(603)
summary_fit <- blm(
  y, ETA = list(X = X, model = "Normal"), residual_var = 1,
  iterations = 60, burnin = 20, store_samples = FALSE,
  store_coefficient_cov = FALSE
)
stopifnot(
  length(coef(summary_fit)) == 3L,
  length(predict(summary_fit, X[1:2, , drop = FALSE])) == 2L,
  identical(colnames(summary(summary_fit)$coefficients), c("Mean", "SD")),
  inherits(
    try(
      predict(summary_fit, X[1:2, , drop = FALSE], interval = "credible"),
      silent = TRUE
    ),
    "try-error"
  )
)

set.seed(602)
multi_fit <- blm(
  y,
  ETA = list(
    first_block = list(X = X[, "first", drop = FALSE], model = "Normal"),
    second_block = list(X = X[, "second", drop = FALSE], model = "Normal")
  ),
  residual_var = 1, iterations = 100, burnin = 40
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
  function() predict(single_fit, X[1, , drop = FALSE], level = 1),
  function() predict(single_fit, X[1, , drop = FALSE], chunk_size = 0),
  function() predict(
    single_fit, X[1, , drop = FALSE],
    chunk_size = .Machine$integer.max + 1
  ),
  function() predict(single_fit, X[1, , drop = FALSE], interval = "invalid"),
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

invalid_summaries <- list(
  function() summary(single_fit, probs = 0),
  function() summary(single_fit, probs = c(0.9, 0.1)),
  function() summary(single_fit, probs = c(0.5, 0.5)),
  function() summary(single_fit, max_coefficients = -1),
  function() summary(
    single_fit, coefficients = "top", max_coefficients = 1,
    rank_by = "inclusion_probability"
  )
)
stopifnot(all(vapply(
  invalid_summaries,
  function(call) inherits(try(call(), silent = TRUE), "try-error"),
  logical(1)
)))
