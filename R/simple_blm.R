#' Bayesian simple linear regression
#'
#' Computes the posterior distribution of the slope in a simple linear
#' regression. The intercept is integrated out by centering both the predictor
#' and response. The residual variance can either be fixed or learned using an
#' inverse-gamma prior and Gibbs sampling.
#'
#' @param y A numeric vector containing the response values.
#' @param x A numeric vector containing the predictor values.
#' @param prior_var A positive numeric scalar giving the prior variance of the
#'   regression coefficient.
#' @param residual_var A positive numeric scalar giving the known residual
#'   variance, or `NULL` to learn it from the data.
#' @param residual_shape A positive numeric scalar giving the shape of the
#'   inverse-gamma prior for the residual variance. Required when
#'   `residual_var = NULL`.
#' @param residual_scale A positive numeric scalar giving the scale of the
#'   inverse-gamma prior for the residual variance. Required when
#'   `residual_var = NULL`. The inverse-gamma density is proportional to
#'   \eqn{v^{-a-1} \exp(-b/v)}, where \eqn{a} is `residual_shape` and \eqn{b} is
#'   `residual_scale`.
#' @param iterations A positive integer giving the total number of Gibbs
#'   iterations when the residual variance is learned.
#' @param burnin A non-negative integer giving the number of initial Gibbs
#'   iterations to discard.
#' @param thin A positive integer giving the interval between retained draws.
#' @param seed `NULL` or an integer used to initialize the random-number
#'   generator.
#'
#' @return A named list containing `slope_mean`, `slope_var`,
#'   `intercept_mean`, and `intercept_var`. When the residual variance is
#'   learned, the list additionally contains `residual_var_mean`,
#'   `residual_var_var`, `slope_samples`, `intercept_samples`, and
#'   `residual_var_samples`.
#'
#' @details The slope and residual variance priors are independent. The slope
#'   has a zero-mean normal prior with variance `prior_var`, and the residual
#'   variance has an inverse-gamma prior. Posterior summaries are computed from
#'   retained Gibbs draws when the residual variance is learned.
#' @export
#'
#' @examples
#' x <- 1:5
#' y <- 1 + 2 * x
#' simple_blm(y, x, prior_var = 10, residual_var = 1)
#' simple_blm(
#'   y, x,
#'   prior_var = 10,
#'   residual_shape = 2,
#'   residual_scale = 1,
#'   seed = 123
#' )
simple_blm <- function(y, x, prior_var, residual_var = NULL,
                       residual_shape = NULL, residual_scale = NULL,
                       iterations = 4000L, burnin = 1000L, thin = 1L,
                       seed = NULL) {
  if (!is.numeric(y) || !is.atomic(y) || is.object(y) || !is.null(dim(y))) {
    stop("`y` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(x) || !is.atomic(x) || is.object(x) || !is.null(dim(x))) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  if (length(y) != length(x)) {
    stop("`y` and `x` must have the same length.", call. = FALSE)
  }
  if (length(y) < 2L) {
    stop("`y` and `x` must contain at least two observations.", call. = FALSE)
  }
  if (anyNA(y) || anyNA(x) || any(!is.finite(y)) || any(!is.finite(x))) {
    stop("`y` and `x` must contain only finite, non-missing values.", call. = FALSE)
  }

  .validate_variance(prior_var, "prior_var")

  x_mean <- mean(x)
  y_mean <- mean(y)
  x_centered <- x - x_mean
  y_centered <- y - y_mean

  if (!is.null(residual_var)) {
    if (!is.null(residual_shape) || !is.null(residual_scale)) {
      stop(
        "Supply either `residual_var` or the inverse-gamma prior, not both.",
        call. = FALSE
      )
    }
    .validate_variance(residual_var, "residual_var")

    posterior_var <- 1 / (sum(x_centered^2) / residual_var + 1 / prior_var)
    posterior_mean <- posterior_var *
      sum(x_centered * y_centered) / residual_var
    intercept_mean <- y_mean - posterior_mean * x_mean
    intercept_var <- residual_var / length(y) + x_mean^2 * posterior_var

    return(list(
      slope_mean = posterior_mean,
      slope_var = posterior_var,
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
    x = matrix(x, ncol = 1L, dimnames = list(NULL, "x")),
    prior_var = prior_var,
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    seed = seed
  )
  slope_samples <- drop(samples$coefficient_samples)

  list(
    slope_mean = mean(slope_samples),
    slope_var = stats::var(slope_samples),
    intercept_mean = mean(samples$intercept_samples),
    intercept_var = stats::var(samples$intercept_samples),
    residual_var_mean = mean(samples$residual_var_samples),
    residual_var_var = stats::var(samples$residual_var_samples),
    slope_samples = slope_samples,
    intercept_samples = samples$intercept_samples,
    residual_var_samples = samples$residual_var_samples
  )
}
