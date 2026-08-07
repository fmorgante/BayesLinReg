# MCMC and execution-control validation helpers.

.validate_mcmc <- function(iterations, burnin, thin) {
  is_whole_number <- function(value) {
    is.numeric(value) && length(value) == 1L && !is.na(value) &&
      is.finite(value) && value == floor(value)
  }

  if (!is_whole_number(iterations) || iterations < 2) {
    stop("`iterations` must be an integer greater than one.", call. = FALSE)
  }
  if (!is_whole_number(burnin) || burnin < 0 || burnin >= iterations) {
    stop(
      "`burnin` must be a non-negative integer smaller than `iterations`.",
      call. = FALSE
    )
  }
  if (!is_whole_number(thin) || thin < 1) {
    stop("`thin` must be a positive integer.", call. = FALSE)
  }

  retained_iterations <- seq.int(burnin + 1L, iterations, by = thin)
  if (length(retained_iterations) < 2L) {
    stop("The MCMC settings must retain at least two draws.", call. = FALSE)
  }

  retained_iterations
}

.validate_nchains <- function(nchains) {
  if (!is.numeric(nchains) || length(nchains) != 1L || is.na(nchains) ||
      !is.finite(nchains) || nchains != floor(nchains) || nchains < 1) {
    stop("`nchains` must be a positive integer.", call. = FALSE)
  }
  as.integer(nchains)
}

.validate_nthreads <- function(nthreads) {
  if (!is.numeric(nthreads) || length(nthreads) != 1L || is.na(nthreads) ||
      !is.finite(nthreads) || nthreads != floor(nthreads) || nthreads < 1) {
    stop("`nthreads` must be a positive integer.", call. = FALSE)
  }
  as.integer(nthreads)
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
