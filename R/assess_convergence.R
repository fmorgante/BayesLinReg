#' Assess Gibbs sampler convergence
#'
#' Produces trace plots and computes convergence diagnostics for the retained
#' non-coefficient draws in a sampled fit returned by [blm()]. These include
#' the intercept and residual variance; the shared coefficient variance for
#' each normal block; the inclusion probability and slab variance for each
#' spike-and-slab block; and the global variance for each global-local block.
#' All per-predictor quantities, including regression coefficients, local
#' variances, and spike-and-slab inclusion indicators, are intentionally
#' excluded.
#'
#' @param fit A fit returned by [blm()] with `store_samples = TRUE`.
#' @param plot A logical scalar indicating whether to draw trace plots for the
#'   assessed parameters.
#'
#' @return A named list containing `rhat`, a named vector of classical
#'   Gelman-Rubin statistics; `geweke`, a parameter-by-chain matrix of Geweke
#'   z-scores; `effective_sample_size`, a named vector; `nchains`; and
#'   `draws_per_chain`. Rows and names identify blocks where applicable.
#'   R-hat is `NA` for single-chain fits and parameters with zero
#'   within-chain variance.
#' @export
#'
#' @examples
#' X <- cbind(x1 = 1:20, x2 = rep(c(0, 1), 10))
#' y <- 1 + 2 * X[, "x1"] - X[, "x2"]
#' fit <- blm(
#'   y,
#'   ETA = list(X = X, model = "Normal", var_shape = 2, var_scale = 10),
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
    .plot_blm_traces(chains)
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

.plot_blm_traces <- function(chains) {
  parameter_names <- coda::varnames(chains)
  titles <- paste("Trace of", parameter_names)
  for (parameter_index in seq_along(parameter_names)) {
    coda::traceplot(
      chains[, parameter_index, drop = FALSE],
      main = titles[parameter_index]
    )
  }
  invisible(titles)
}
