#' Bayesian simple linear regression
#'
#' Computes the posterior distribution of the slope in a simple linear
#' regression. The intercept is integrated out by centering both the predictor
#' and response. The residual variance can either be fixed or learned using a
#' conjugate inverse-gamma prior.
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
#'
#' @return A named list containing `slope_mean`, `slope_var`,
#'   `intercept_mean`, and `intercept_var`. When the residual variance is
#'   learned, the list additionally contains `slope_scale`, `slope_df`,
#'   `intercept_scale`, and `intercept_df` for the marginal Student t
#'   posteriors, plus `residual_var_shape` and `residual_var_scale` for the
#'   inverse-gamma posterior of the residual variance.
#'
#' @details When the residual variance is learned, the conjugate prior is
#'   \eqn{\beta \mid \sigma^2 \sim N(0, \sigma^2 V_0)}, where \eqn{V_0} is
#'   `prior_var`, and \eqn{\sigma^2 \sim IG(a, b)}. Thus, in this case,
#'   `prior_var` is the slope's variance relative to the residual variance.
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
#'   residual_scale = 1
#' )
simple_blm <- function(y, x, prior_var, residual_var = NULL,
                       residual_shape = NULL, residual_scale = NULL) {
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

  # Normal-inverse-gamma update. Centering removes one residual degree of
  # freedom because the intercept has been integrated out.
  posterior_relative_var <- 1 / (sum(x_centered^2) + 1 / prior_var)
  posterior_mean <- posterior_relative_var * sum(x_centered * y_centered)
  posterior_residual_shape <- residual_shape + (length(y) - 1) / 2
  posterior_residual_scale <- residual_scale + 0.5 * (
    sum((y_centered - posterior_mean * x_centered)^2) +
      posterior_mean^2 / prior_var
  )
  posterior_scale <- sqrt(
    posterior_residual_scale / posterior_residual_shape *
      posterior_relative_var
  )
  intercept_mean <- y_mean - posterior_mean * x_mean
  intercept_relative_var <- 1 / length(y) +
    x_mean^2 * posterior_relative_var
  intercept_scale <- sqrt(
    posterior_residual_scale / posterior_residual_shape *
      intercept_relative_var
  )
  posterior_var <- if (posterior_residual_shape > 1) {
    posterior_residual_scale / (posterior_residual_shape - 1) *
      posterior_relative_var
  } else {
    Inf
  }
  intercept_var <- if (posterior_residual_shape > 1) {
    posterior_residual_scale / (posterior_residual_shape - 1) *
      intercept_relative_var
  } else {
    Inf
  }

  list(
    slope_mean = posterior_mean,
    slope_var = posterior_var,
    slope_scale = posterior_scale,
    slope_df = 2 * posterior_residual_shape,
    intercept_mean = intercept_mean,
    intercept_var = intercept_var,
    intercept_scale = intercept_scale,
    intercept_df = 2 * posterior_residual_shape,
    residual_var_shape = posterior_residual_shape,
    residual_var_scale = posterior_residual_scale
  )
}
