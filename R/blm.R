#' Bayesian linear regression
#'
#' Fits a Bayesian linear regression with one or more predictor blocks using a
#' BGLR-style `ETA` interface. Each block specifies its predictors,
#' coefficient-prior family, standardization, and prior parameters. The
#' intercept is integrated out during coefficient sampling by centering the
#' response and predictors.
#'
#' @param y A finite numeric response vector.
#' @param ETA A predictor specification or named list of predictor
#'   specifications. Each block must contain `X` and `model`. Available models
#'   are `"Normal"`, `"SpikeSlab"`, and `"GlobalLocal"`. See Details.
#' @param residual_var A positive known residual variance, or `NULL` to learn
#'   it using an inverse-gamma prior.
#' @param residual_shape,residual_scale Positive shape and scale parameters for
#'   the inverse-gamma residual-variance prior. Required when
#'   `residual_var = NULL`.
#' @param iterations Total Gibbs iterations when sampling is required.
#' @param burnin Number of initial iterations to discard.
#' @param thin Interval between retained draws.
#' @param seed `NULL` or an integer random-number seed.
#' @param version Gibbs implementation: `"Rcpp"` or `"R"`.
#' @param verbose If `TRUE`, display aggregate progress at 10-percent intervals
#'   per chain using [progressr::with_progress()].
#' @param nchains Number of independent chains. Multiple chains use a temporary
#'   [future::multisession] plan.
#'
#' @return A list containing `ETA`, a named list of block-specific posterior
#'   summaries, plus intercept and residual-variance summaries. Sampled fits
#'   also contain intercept and residual-variance draws. With multiple chains,
#'   `chain_id` identifies the origin of each retained draw.
#'
#' @details `ETA` may be a single-block specification such as
#'   `list(X = X, model = "Normal")`, or a named list of blocks.
#'   A numeric vector supplied as a block's `X` is treated as a one-column
#'   matrix.
#'   Every block accepts `standardize`, which defaults to `TRUE`. Returned
#'   coefficients are always transformed to the original scale of that block's
#'   supplied `X`.
#'
#'   A `"Normal"` block optionally accepts `var_shape = 2` and
#'   `var_scale = 1`. Its coefficients share a variance sampled from an
#'   inverse-gamma prior with the supplied shape and scale.
#'   A `"SpikeSlab"` block optionally accepts `pi = c(a = 1, b = 1)`,
#'   `slab_shape = 2`, and `slab_scale = 1`, and requires the residual variance
#'   to be learned. Its shared slab variance has an inverse-gamma prior with
#'   the supplied shape and scale. A `"GlobalLocal"` block optionally accepts
#'   `local_shape = c(a = 1, b = 0.5)` and `global_scale = 1`. Its hierarchy is
#'   \deqn{\beta_j \mid \tau^2,\psi_j \sim N(0,\tau^2\psi_j),\qquad
#'   \psi_j \sim \mathrm{BetaPrime}(a,b),}
#'   with \eqn{\tau \sim C^+(0,\mathrm{global\_scale})}. Thus the default is
#'   Strawderman-Berger, while `local_shape = c(0.5, 0.5)` gives the horseshoe.
#'   The coefficient priors are independent of the residual variance.
#' @export
#'
#' @examples
#' X <- cbind(x1 = 1:20, x2 = rep(c(0, 1), 10))
#' y <- 1 + 2 * X[, "x1"] - X[, "x2"]
#' blm(
#'   y,
#'   ETA = list(X = X, model = "Normal", var_shape = 2, var_scale = 10),
#'   residual_var = 1
#' )
#' blm(
#'   y,
#'   ETA = list(markers = list(X = X, model = "GlobalLocal")),
#'   residual_shape = 2,
#'   residual_scale = 1,
#'   iterations = 100,
#'   burnin = 50,
#'   seed = 123
#' )
blm <- function(y, ETA, residual_var = NULL,
                residual_shape = NULL, residual_scale = NULL,
                iterations = 4000L, burnin = 1000L, thin = 1L,
                seed = NULL, version = c("Rcpp", "R"),
                verbose = FALSE, nchains = 1L) {
  version <- match.arg(version)
  nchains <- .validate_nchains(nchains)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(y) || !is.atomic(y) || is.object(y) || !is.null(dim(y))) {
    stop("`y` must be a numeric vector.", call. = FALSE)
  }
  if (length(y) < 2L || anyNA(y) || any(!is.finite(y))) {
    stop(
      "`y` must contain at least two finite, non-missing values.",
      call. = FALSE
    )
  }

  blocks <- .normalize_eta(ETA, length(y), residual_var)
  block_sizes <- vapply(blocks, function(block) ncol(block$x), integer(1))
  block_ends <- cumsum(block_sizes)
  block_starts <- block_ends - block_sizes + 1L
  block_indices <- Map(seq.int, block_starts, block_ends)
  block_model <- vapply(blocks, `[[`, integer(1), "model_code")
  block_id <- rep.int(seq_along(blocks), block_sizes)

  x <- do.call(cbind, lapply(seq_along(blocks), function(block_index) {
    block_x <- blocks[[block_index]]$x
    colnames(block_x) <- paste0(
      names(blocks)[block_index], "::", blocks[[block_index]]$predictor_names
    )
    block_x
  }))
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

  normal_shape <- vapply(blocks, `[[`, numeric(1), "normal_shape")
  normal_scale <- vapply(blocks, `[[`, numeric(1), "normal_scale")
  pi_alpha <- vapply(blocks, `[[`, numeric(1), "pi_alpha")
  pi_beta <- vapply(blocks, `[[`, numeric(1), "pi_beta")
  slab_shape <- vapply(blocks, `[[`, numeric(1), "slab_shape")
  slab_scale <- vapply(blocks, `[[`, numeric(1), "slab_scale")
  global_scale <- vapply(blocks, `[[`, numeric(1), "global_scale")
  local_a <- vapply(blocks, function(block) block$local_shape[1L], numeric(1))
  local_b <- vapply(blocks, function(block) block$local_shape[2L], numeric(1))
  sampler_arguments <- list(
    y = y,
    x = x,
    residual_shape = if (is.null(residual_shape)) 1 else residual_shape,
    residual_scale = if (is.null(residual_scale)) 1 else residual_scale,
    residual_var = residual_var,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    block_id = block_id,
    block_model = block_model,
    normal_shape = normal_shape,
    normal_scale = normal_scale,
    pi_alpha = pi_alpha,
    pi_beta = pi_beta,
    slab_shape = slab_shape,
    slab_scale = slab_scale,
    global_scale = global_scale,
    local_a = local_a,
    local_b = local_b
  )
  run_chains <- function(progressor = NULL) {
    .run_blm_chains(
      sampler_arguments = sampler_arguments,
      version = version,
      nchains = nchains,
      seed = seed,
      block_model = block_model,
      progressor = progressor
    )
  }
  samples <- if (verbose) {
    progressr::with_progress({
      progress <- progressr::progressor(steps = nchains * iterations)
      run_chains(progress)
    }, enable = TRUE)
  } else {
    run_chains()
  }

  eta_result <- lapply(seq_along(blocks), function(block_index) {
    block <- blocks[[block_index]]
    indices <- block_indices[[block_index]]
    coefficient_samples <- sweep(
      samples$coefficient_samples[, indices, drop = FALSE],
      2L,
      block$predictor_scale,
      FUN = "/"
    )
    colnames(coefficient_samples) <- block$predictor_names
    result <- list(
      model = block$model,
      standardize = block$standardize,
      coefficient_mean = colMeans(coefficient_samples),
      coefficient_cov = stats::cov(coefficient_samples),
      coefficient_samples = coefficient_samples
    )
    if (block$model == "Normal") {
      normal_var_samples <- samples$normal_var_samples[, block_index]
      result$normal_var_mean <- mean(normal_var_samples)
      result$normal_var_var <- stats::var(normal_var_samples)
      result$normal_var_samples <- normal_var_samples
      result$var_shape <- block$normal_shape
      result$var_scale <- block$normal_scale
    }
    if (block$model == "SpikeSlab") {
      inclusion_samples <- samples$inclusion_samples[, indices, drop = FALSE]
      colnames(inclusion_samples) <- block$predictor_names
      pi_samples <- samples$pi_samples[, block_index]
      slab_var_samples <- samples$slab_var_samples[, block_index]
      result$inclusion_probability <- colMeans(inclusion_samples)
      result$pi_mean <- mean(pi_samples)
      result$pi_var <- stats::var(pi_samples)
      result$slab_var_mean <- mean(slab_var_samples)
      result$slab_var_var <- stats::var(slab_var_samples)
      result$inclusion_samples <- inclusion_samples
      result$pi_samples <- pi_samples
      result$slab_var_samples <- slab_var_samples
      result$pi <- c(a = block$pi_alpha, b = block$pi_beta)
      result$slab_shape <- block$slab_shape
      result$slab_scale <- block$slab_scale
    }
    if (block$model == "GlobalLocal") {
      local_var_samples <- samples$local_var_samples[, indices, drop = FALSE]
      colnames(local_var_samples) <- block$predictor_names
      tau_sq_samples <- samples$tau_sq_samples[, block_index]
      result$local_var_mean <- colMeans(local_var_samples)
      result$local_var_var <- apply(local_var_samples, 2L, stats::var)
      result$tau_sq_mean <- mean(tau_sq_samples)
      result$tau_sq_var <- stats::var(tau_sq_samples)
      result$local_var_samples <- local_var_samples
      result$tau_sq_samples <- tau_sq_samples
      result$local_shape <- block$local_shape
      result$global_scale <- block$global_scale
    }
    result
  })
  names(eta_result) <- names(blocks)

  result <- list(
    ETA = eta_result,
    intercept_mean = mean(samples$intercept_samples),
    intercept_var = stats::var(samples$intercept_samples),
    residual_var_mean = mean(samples$residual_var_samples),
    residual_var_var = stats::var(samples$residual_var_samples),
    intercept_samples = samples$intercept_samples,
    residual_var_samples = samples$residual_var_samples
  )
  if (nchains > 1L) {
    result$nchains <- nchains
    result$chain_id <- samples$chain_id
  }
  result
}
