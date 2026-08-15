# MCMC and execution-control validation helpers.

.validate_bounded_integer <- function(
    value, name, minimum = 0, maximum = .Machine$integer.max,
    message = NULL) {
  valid <- is.numeric(value) && length(value) == 1L && !is.na(value) &&
    is.finite(value) && value == floor(value) &&
    isTRUE(value >= minimum && value <= maximum)
  if (!valid) {
    if (is.null(message)) {
      message <- sprintf(
        "`%s` must be an integer between %s and %s.",
        name, format(minimum, scientific = FALSE),
        format(maximum, scientific = FALSE)
      )
    }
    stop(message, call. = FALSE)
  }
  as.integer(value)
}

.validate_mcmc <- function(iterations, burnin, thin) {
  iterations <- .validate_bounded_integer(
    iterations, "iterations", minimum = 2,
    message =
      paste0(
        "`iterations` must be an integer greater than one and no larger ",
        "than `.Machine$integer.max`."
      )
  )
  burnin <- .validate_bounded_integer(
    burnin, "burnin", maximum = iterations - 1L,
    message = "`burnin` must be a non-negative integer smaller than `iterations`."
  )
  thin <- .validate_bounded_integer(
    thin, "thin", minimum = 1,
    message =
      paste0(
        "`thin` must be a positive integer no larger than ",
        "`.Machine$integer.max`."
      )
  )

  retained_iterations <- seq.int(burnin + 1L, iterations, by = thin)
  if (length(retained_iterations) < 2L) {
    stop("The MCMC settings must retain at least two draws.", call. = FALSE)
  }

  retained_iterations
}

.validate_nchains <- function(nchains) {
  .validate_bounded_integer(
    nchains, "nchains", minimum = 1,
    message =
      paste0(
        "`nchains` must be a positive integer no larger than ",
        "`.Machine$integer.max`."
      )
  )
}

.validate_nthreads <- function(nthreads) {
  .validate_bounded_integer(
    nthreads, "nthreads", minimum = 1,
    message =
      paste0(
        "`nthreads` must be a positive integer no larger than ",
        "`.Machine$integer.max`."
      )
  )
}

.resolve_likelihood_df <- function(effective_n, fit_intercept,
                                   likelihood_df = NULL) {
  if (is.null(likelihood_df)) {
    likelihood_df <- effective_n - as.integer(fit_intercept)
  }
  .validate_bounded_integer(
    likelihood_df, "likelihood_df", minimum = 1,
    maximum = min(effective_n, .Machine$integer.max),
    message = "`likelihood_df` must be a positive integer no greater than `n`."
  )
}

.validate_pve_controls <- function(compute_pve, pve_type) {
  if (!is.logical(compute_pve) || length(compute_pve) != 1L ||
      is.na(compute_pve)) {
    stop("`compute_pve` must be TRUE or FALSE.", call. = FALSE)
  }
  list(
    compute_pve = compute_pve,
    pve_type = match.arg(pve_type, c("standalone", "allocated"))
  )
}
