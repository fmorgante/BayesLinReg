#' Assess Gibbs sampler convergence
#'
#' Produces trace plots and computes convergence diagnostics for the retained
#' intercept and residual-variance draws in a sampled fit returned by
#' [multiple_blm()]. The spike-and-slab inclusion probability `pi` is also
#' assessed when present.
#'
#' @param fit A sampled fit returned by [multiple_blm()]. Fits using a known
#'   residual variance do not contain samples and cannot be assessed.
#' @param plot A logical scalar indicating whether to draw trace plots for the
#'   assessed parameters.
#'
#' @return A named list containing `rhat`, a named vector of classical
#'   Gelman-Rubin statistics; `geweke`, a parameter-by-chain matrix of Geweke
#'   z-scores; `effective_sample_size`, a named vector; `nchains`; and
#'   `draws_per_chain`. R-hat is `NA` for single-chain fits and parameters
#'   with zero within-chain variance.
#' @export
#'
#' @examples
#' X <- cbind(x1 = 1:20, x2 = rep(c(0, 1), 10))
#' y <- 1 + 2 * X[, "x1"] - X[, "x2"]
#' fit <- multiple_blm(
#'   y, X,
#'   prior_var = 10,
#'   residual_shape = 2,
#'   residual_scale = 1,
#'   iterations = 100,
#'   burnin = 50,
#'   seed = 123
#' )
#' diagnostics <- assess_convergence(fit, plot = FALSE)
assess_convergence <- function(fit, plot = TRUE) {
  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    stop("`plot` must be TRUE or FALSE.", call. = FALSE)
  }

  chains <- .as_blm_mcmc_list(fit)
  parameter_names <- coda::varnames(chains)
  number_of_chains <- coda::nchain(chains)

  if (plot) {
    coda::traceplot(chains)
  }

  geweke <- vapply(chains, function(chain) {
    z_scores <- suppressWarnings(tryCatch(
      coda::geweke.diag(chain)$z,
      error = function(error) rep(NA_real_, length(parameter_names))
    ))
    stats::setNames(z_scores, parameter_names)
  }, numeric(length(parameter_names)))
  if (is.null(dim(geweke))) {
    geweke <- matrix(
      geweke,
      ncol = 1L,
      dimnames = list(parameter_names, "chain_1")
    )
  } else {
    rownames(geweke) <- parameter_names
    colnames(geweke) <- paste0("chain_", seq_len(number_of_chains))
  }

  effective_sample_size <- suppressWarnings(tryCatch(
    coda::effectiveSize(chains),
    error = function(error) {
      stats::setNames(rep(NA_real_, length(parameter_names)), parameter_names)
    }
  ))

  list(
    rhat = .classical_rhat(chains),
    geweke = geweke,
    effective_sample_size = effective_sample_size,
    nchains = number_of_chains,
    draws_per_chain = coda::niter(chains)
  )
}
