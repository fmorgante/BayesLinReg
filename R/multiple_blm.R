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
#' @param prior_var For the normal and spike-and-slab priors, a positive numeric
#'   scalar giving a common prior variance
#'   for all coefficients, or a positive numeric vector with one variance per
#'   predictor. Must be `NULL` for the global-local prior.
#' @param coefficient_prior The coefficient-prior family: `"normal"`,
#'   `"spike_slab"`, or `"global_local"`. The global-local family uses
#'   beta-prime local variances and learns a shared global variance.
#' @param pi_alpha,pi_beta Positive shape parameters for the Beta prior on the
#'   shared inclusion probability \eqn{\pi}. Used with the spike-and-slab prior.
#' @param local_shape A positive numeric vector of length two containing the
#'   beta-prime shape parameters `a` and `b`. The default `c(1, 0.5)` gives the
#'   Strawderman-Berger prior; `c(0.5, 0.5)` gives the horseshoe prior. Used
#'   with the global-local prior.
#' @param global_scale A positive numeric scalar giving the scale of the
#'   half-Cauchy prior on the global standard deviation \eqn{\tau}. The global
#'   variance \eqn{\tau^2} is learned from the data. Used with the global-local
#'   prior.
#' @param standardize A logical scalar. If `TRUE`, the predictors are centered
#'   and scaled before fitting. Coefficients and intercepts are always returned
#'   on the original predictor scale.
#' @param residual_var A positive numeric scalar giving the known residual
#'   variance, or `NULL` to learn it from the data.
#' @param residual_shape A positive numeric scalar giving the shape of the
#'   inverse-gamma prior. Required when `residual_var = NULL`.
#' @param residual_scale A positive numeric scalar giving the scale of the
#'   inverse-gamma prior. Required when `residual_var = NULL`.
#' @param iterations A positive integer giving the total number of Gibbs
#'   iterations when Gibbs sampling is used.
#' @param burnin A non-negative integer giving the number of initial Gibbs
#'   iterations to discard.
#' @param thin A positive integer giving the interval between retained draws.
#' @param seed `NULL` or an integer used to initialize the random-number
#'   generator.
#' @param version The Gibbs-sampler implementation to use: `"Rcpp"` for the
#'   compiled C++ implementation or `"R"` for the reference R implementation.
#' @param verbose A logical scalar. If `TRUE`, display aggregate Gibbs
#'   progress using [progressr::with_progress()], updated every 10 percent per
#'   chain.
#' @param nchains A positive integer giving the number of independent MCMC
#'   chains. Values greater than one run chains in parallel using a temporary
#'   \code{future::multisession} plan.
#'
#' @return A named list containing `coefficient_mean`, `coefficient_cov`,
#'   `intercept_mean`, and `intercept_var`. When the residual variance is
#'   sampled, the list additionally contains `residual_var_mean`,
#'   `residual_var_var`, `coefficient_samples`, `intercept_samples`, and
#'   `residual_var_samples`.
#'   With the spike-and-slab prior, the list also contains
#'   `inclusion_probability`, `pi_mean`, `pi_var`, `inclusion_samples`,
#'   and `pi_samples`.
#'   With the global-local prior, it also contains `local_var_mean`,
#'   `local_var_var`, `tau_sq_mean`, `tau_sq_var`, `local_var_samples`,
#'   `tau_sq_samples`, and `local_shape`. These variance parameters refer to
#'   the internally standardized coefficient scale when `standardize = TRUE`.
#'   When `nchains > 1`, retained draws are combined across chains and the
#'   result additionally contains `nchains` and a corresponding `chain_id`
#'   for each draw.
#'
#' @details With known residual variance, the coefficients have independent
#'   zero-mean normal priors with variances given by `prior_var`. These priors
#'   are independent of the inverse-gamma residual-variance prior. Posterior
#'   summaries are computed from retained Gibbs draws when the residual
#'   variance is learned. During each Gibbs sweep, coefficients are updated one
#'   at a time from their univariate conditional normal distributions.
#'   For the spike-and-slab prior, inclusion indicators are sampled using the
#'   coefficient-marginalized conditional odds, and the shared inclusion
#'   probability is updated from its full conditional Beta distribution.
#'   The global-local hierarchy is
#'   \deqn{\beta_j \mid \tau^2, \psi_j \sim N(0, \tau^2\psi_j),}
#'   \deqn{\psi_j \sim \mathrm{BetaPrime}(a,b), \qquad
#'   \tau \sim C^+(0, \mathrm{global\_scale}).}
#'   The local variances are updated from generalized inverse Gaussian full
#'   conditionals, and \eqn{\tau^2} is updated from its inverse-gamma full
#'   conditional. This coefficient prior is independent of the residual
#'   variance. Global-local fits use Gibbs sampling even when `residual_var` is
#'   known.
#'   When `standardize = TRUE`, `prior_var` and `global_scale` describe the
#'   coefficient scale after predictor standardization. All returned regression
#'   coefficients, their covariance matrix, and intercept quantities are
#'   transformed to the scale of the supplied `X`.
#' @export
#'
#' @examples
#' X <- cbind(x1 = 1:5, x2 = c(0, 1, 0, 1, 0))
#' y <- 1 + 2 * X[, "x1"] - 3 * X[, "x2"]
#' multiple_blm(y, X, prior_var = 10, residual_var = 1)
#' multiple_blm(
#'   y, X,
#'   coefficient_prior = "global_local",
#'   local_shape = c(a = 1, b = 0.5),
#'   residual_shape = 2,
#'   residual_scale = 1,
#'   seed = 123,
#'   version = "Rcpp",
#'   verbose = FALSE,
#'   nchains = 1
#' )
multiple_blm <- function(y, X, prior_var = NULL, residual_var = NULL,
                         coefficient_prior = c(
                           "normal", "spike_slab", "global_local"
                         ),
                         pi_alpha = 1, pi_beta = 1,
                         local_shape = c(a = 1, b = 0.5),
                         global_scale = 1, standardize = TRUE,
                         residual_shape = NULL, residual_scale = NULL,
                         iterations = 4000L, burnin = 1000L, thin = 1L,
                         seed = NULL, version = c("Rcpp", "R"),
                         verbose = FALSE, nchains = 1L) {
  version <- match.arg(version)
  coefficient_prior <- match.arg(coefficient_prior)
  nchains <- .validate_nchains(nchains)
  if (!is.logical(standardize) || length(standardize) != 1L ||
      is.na(standardize)) {
    stop("`standardize` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
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
  sampler_local_shape <- c(a = 1, b = 0.5)
  if (coefficient_prior == "global_local") {
    if (!is.null(prior_var)) {
      stop(
        "`prior_var` must be NULL for the global-local prior.",
        call. = FALSE
      )
    }
    local_shape <- .validate_local_shape(local_shape)
    sampler_local_shape <- local_shape
    .validate_variance(global_scale, "global_scale")
    sampler_prior_var <- rep(1, number_of_predictors)
  } else {
    if (is.null(prior_var)) {
      stop(
        "`prior_var` is required for the normal and spike-and-slab priors.",
        call. = FALSE
      )
    }
    prior_var <- .validate_prior_var(prior_var, number_of_predictors)
    sampler_prior_var <- prior_var
  }
  if (coefficient_prior == "spike_slab") {
    .validate_variance(pi_alpha, "pi_alpha")
    .validate_variance(pi_beta, "pi_beta")
    if (!is.null(residual_var)) {
      stop(
        "`coefficient_prior = \"spike_slab\"` requires learning the residual variance.",
        call. = FALSE
      )
    }
  }

  x_mean <- colMeans(X)
  x_centered <- sweep(X, 2L, x_mean, FUN = "-")
  predictor_scale <- if (standardize) {
    sqrt(colSums(x_centered^2) / (nrow(X) - 1))
  } else {
    rep(1, number_of_predictors)
  }
  if (any(!is.finite(predictor_scale)) || any(predictor_scale <= 0)) {
    stop(
      "`X` cannot contain constant predictors when `standardize = TRUE`.",
      call. = FALSE
    )
  }
  X_sampler <- sweep(X, 2L, predictor_scale, FUN = "/")
  x_sampler_mean <- colMeans(X_sampler)

  y_mean <- mean(y)
  x_centered <- sweep(X_sampler, 2L, x_sampler_mean, FUN = "-")
  y_centered <- y - y_mean

  if (!is.null(residual_var) && coefficient_prior == "normal") {
    if (nchains > 1L) {
      stop(
        "`nchains > 1` is only available when the residual variance is learned.",
        call. = FALSE
      )
    }
    if (!is.null(residual_shape) || !is.null(residual_scale)) {
      stop(
        "Supply either `residual_var` or the inverse-gamma prior, not both.",
        call. = FALSE
      )
    }
    .validate_variance(residual_var, "residual_var")

    prior_precision <- diag(
      1 / sampler_prior_var,
      nrow = number_of_predictors
    )
    posterior_precision <- crossprod(x_centered) / residual_var +
      prior_precision
    sampler_coefficient_cov <- chol2inv(chol(posterior_precision))
    sampler_coefficient_mean <- drop(
      sampler_coefficient_cov %*%
        crossprod(x_centered, y_centered) / residual_var
    )
    coefficient_mean <- sampler_coefficient_mean / predictor_scale
    coefficient_cov <- sampler_coefficient_cov /
      outer(predictor_scale, predictor_scale)
    names(coefficient_mean) <- predictor_names
    dimnames(coefficient_cov) <- list(predictor_names, predictor_names)
    intercept_mean <- drop(
      y_mean - crossprod(x_sampler_mean, sampler_coefficient_mean)
    )
    intercept_var <- drop(
      residual_var / length(y) +
        crossprod(
          x_sampler_mean,
          sampler_coefficient_cov %*% x_sampler_mean
        )
    )

    return(list(
      coefficient_mean = coefficient_mean,
      coefficient_cov = coefficient_cov,
      intercept_mean = intercept_mean,
      intercept_var = intercept_var
    ))
  }

  if (!is.null(residual_var)) {
    if (!is.null(residual_shape) || !is.null(residual_scale)) {
      stop(
        "Supply either `residual_var` or the inverse-gamma prior, not both.",
        call. = FALSE
      )
    }
    .validate_variance(residual_var, "residual_var")
  } else {
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
  }

  sampler_arguments <- list(
    y = y,
    x = X_sampler,
    prior_var = sampler_prior_var,
    residual_shape = if (is.null(residual_shape)) 1 else residual_shape,
    residual_scale = if (is.null(residual_scale)) 1 else residual_scale,
    residual_var = residual_var,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    coefficient_prior = coefficient_prior,
    pi_alpha = pi_alpha,
    pi_beta = pi_beta,
    global_scale = global_scale,
    local_shape = sampler_local_shape
  )

  if (verbose) {
    samples <- progressr::with_progress({
      progress <- progressr::progressor(steps = nchains * iterations)
      .run_blm_chains(
        sampler_arguments = sampler_arguments,
        version = version,
        nchains = nchains,
        seed = seed,
        coefficient_prior = coefficient_prior,
        progressor = progress
      )
    }, enable = TRUE)
  } else {
    samples <- .run_blm_chains(
      sampler_arguments = sampler_arguments,
      version = version,
      nchains = nchains,
      seed = seed,
      coefficient_prior = coefficient_prior
    )
  }

  samples$coefficient_samples <- sweep(
    samples$coefficient_samples,
    2L,
    predictor_scale,
    FUN = "/"
  )
  colnames(samples$coefficient_samples) <- predictor_names

  result <- list(
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
  if (coefficient_prior == "spike_slab") {
    result$inclusion_probability <- colMeans(samples$inclusion_samples)
    result$pi_mean <- mean(samples$pi_samples)
    result$pi_var <- stats::var(samples$pi_samples)
    result$inclusion_samples <- samples$inclusion_samples
    result$pi_samples <- samples$pi_samples
  }
  if (coefficient_prior == "global_local") {
    result$local_var_mean <- colMeans(samples$local_var_samples)
    result$local_var_var <- apply(
      samples$local_var_samples,
      2L,
      stats::var
    )
    result$tau_sq_mean <- mean(samples$tau_sq_samples)
    result$tau_sq_var <- stats::var(samples$tau_sq_samples)
    result$local_var_samples <- samples$local_var_samples
    result$tau_sq_samples <- samples$tau_sq_samples
    result$local_shape <- local_shape
  }
  if (nchains > 1L) {
    result$nchains <- nchains
    result$chain_id <- samples$chain_id
  }
  result
}
