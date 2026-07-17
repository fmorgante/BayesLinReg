#' Bayesian multiple linear regression
#'
#' Computes posterior distributions for the coefficients and intercept in a
#' multiple linear regression. The intercept is integrated out by centering the
#' response and every predictor. The residual variance can be fixed or learned
#' using an inverse-gamma prior and Gibbs sampling.
#'
#' @param y A numeric vector containing the response values.
#' @param X A numeric matrix or data frame with observations in rows and
#'   predictors in columns.
#' @param prior_var A positive numeric scalar giving a common prior variance
#'   for all coefficients, or a positive numeric vector with one variance per
#'   predictor.
#' @param residual_var A positive numeric scalar giving the known residual
#'   variance, or `NULL` to learn it from the data.
#' @param residual_shape A positive numeric scalar giving the shape of the
#'   inverse-gamma prior. Required when `residual_var = NULL`.
#' @param residual_scale A positive numeric scalar giving the scale of the
#'   inverse-gamma prior. Required when `residual_var = NULL`.
#' @param iterations A positive integer giving the total number of Gibbs
#'   iterations when the residual variance is learned.
#' @param burnin A non-negative integer giving the number of initial Gibbs
#'   iterations to discard.
#' @param thin A positive integer giving the interval between retained draws.
#' @param seed `NULL` or an integer used to initialize the random-number
#'   generator.
#'
#' @return A named list containing `coefficient_mean`, `coefficient_cov`,
#'   `intercept_mean`, and `intercept_var`. When the residual variance is
#'   learned, the list additionally contains `residual_var_mean`,
#'   `residual_var_var`, `coefficient_samples`, `intercept_samples`, and
#'   `residual_var_samples`.
#'
#' @details With known residual variance, the coefficients have independent
#'   zero-mean normal priors with variances given by `prior_var`. These priors
#'   are independent of the inverse-gamma residual-variance prior. Posterior
#'   summaries are computed from retained Gibbs draws when the residual
#'   variance is learned.
#' @export
#'
#' @examples
#' X <- cbind(x1 = 1:5, x2 = c(0, 1, 0, 1, 0))
#' y <- 1 + 2 * X[, "x1"] - 3 * X[, "x2"]
#' multiple_blm(y, X, prior_var = 10, residual_var = 1)
#' multiple_blm(
#'   y, X,
#'   prior_var = c(10, 5),
#'   residual_shape = 2,
#'   residual_scale = 1,
#'   seed = 123
#' )
multiple_blm <- function(y, X, prior_var, residual_var = NULL,
                         residual_shape = NULL, residual_scale = NULL,
                         iterations = 4000L, burnin = 1000L, thin = 1L,
                         seed = NULL) {
  if (!is.numeric(y) || !is.atomic(y) || is.object(y) || !is.null(dim(y))) {
    stop("`y` must be a numeric vector.", call. = FALSE)
  }
  if (length(y) < 2L) {
    stop("`y` must contain at least two observations.", call. = FALSE)
  }
  if (anyNA(y) || any(!is.finite(y))) {
    stop("`y` must contain only finite, non-missing values.", call. = FALSE)
  }

  X <- .as_predictor_matrix(X, length(y))
  number_of_predictors <- ncol(X)
  predictor_names <- colnames(X)
  prior_var <- .validate_prior_var(prior_var, number_of_predictors)
  prior_precision <- diag(1 / prior_var, nrow = number_of_predictors)

  x_mean <- colMeans(X)
  y_mean <- mean(y)
  x_centered <- sweep(X, 2L, x_mean, FUN = "-")
  y_centered <- y - y_mean

  if (!is.null(residual_var)) {
    if (!is.null(residual_shape) || !is.null(residual_scale)) {
      stop(
        "Supply either `residual_var` or the inverse-gamma prior, not both.",
        call. = FALSE
      )
    }
    .validate_variance(residual_var, "residual_var")

    posterior_precision <- crossprod(x_centered) / residual_var +
      prior_precision
    coefficient_cov <- chol2inv(chol(posterior_precision))
    coefficient_mean <- drop(
      coefficient_cov %*% crossprod(x_centered, y_centered) / residual_var
    )
    names(coefficient_mean) <- predictor_names
    dimnames(coefficient_cov) <- list(predictor_names, predictor_names)
    intercept_mean <- drop(y_mean - crossprod(x_mean, coefficient_mean))
    intercept_var <- drop(
      residual_var / length(y) +
        crossprod(x_mean, coefficient_cov %*% x_mean)
    )

    return(list(
      coefficient_mean = coefficient_mean,
      coefficient_cov = coefficient_cov,
      intercept_mean = intercept_mean,
      intercept_var = intercept_var
    ))
  }

  if (is.null(residual_shape) || is.null(residual_scale)) {
    stop(
      paste0(
        "`residual_shape` and `residual_scale` are required when ",
        "`residual_var` is NULL."
      ),
      call. = FALSE
    )
  }
  .validate_variance(residual_shape, "residual_shape")
  .validate_variance(residual_scale, "residual_scale")

  samples <- .blm_gibbs(
    y = y,
    x = X,
    prior_var = prior_var,
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    seed = seed
  )

  list(
    coefficient_mean = colMeans(samples$coefficient_samples),
    coefficient_cov = stats::cov(samples$coefficient_samples),
    intercept_mean = mean(samples$intercept_samples),
    intercept_var = stats::var(samples$intercept_samples),
    residual_var_mean = mean(samples$residual_var_samples),
    residual_var_var = stats::var(samples$residual_var_samples),
    coefficient_samples = samples$coefficient_samples,
    intercept_samples = samples$intercept_samples,
    residual_var_samples = samples$residual_var_samples
  )
}
