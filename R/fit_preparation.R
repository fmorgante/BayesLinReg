# Shared fitting-interface preparation helpers.

.prepare_block_layout <- function(blocks, source_indices = NULL,
                                  predictor_scales = NULL) {
  block_sizes <- if (is.null(source_indices)) {
    vapply(blocks, function(block) length(block$predictor_names), integer(1))
  } else {
    lengths(source_indices)
  }
  block_ends <- cumsum(block_sizes)
  block_starts <- block_ends - block_sizes + 1L
  block_indices <- Map(seq.int, block_starts, block_ends)
  source_order <- if (is.null(source_indices)) {
    seq_len(sum(block_sizes))
  } else {
    unlist(source_indices, use.names = FALSE)
  }
  scale_order <- if (is.null(predictor_scales)) {
    rep(1, sum(block_sizes))
  } else {
    unlist(predictor_scales, use.names = FALSE)
  }
  internal_names <- unlist(lapply(seq_along(blocks), function(block_index) {
    paste0(
      names(blocks)[block_index], "::", blocks[[block_index]]$predictor_names
    )
  }))

  list(
    block_sizes = block_sizes,
    block_indices = block_indices,
    block_model = vapply(blocks, `[[`, integer(1), "model_code"),
    block_id = rep.int(seq_along(blocks), block_sizes),
    source_order = source_order,
    scale_order = scale_order,
    internal_names = internal_names
  )
}

.prior_sampler_arguments <- function(blocks) {
  list(
    normal_shape = vapply(blocks, `[[`, numeric(1), "normal_shape"),
    normal_scale = vapply(blocks, `[[`, numeric(1), "normal_scale"),
    pi_alpha = vapply(blocks, `[[`, numeric(1), "pi_alpha"),
    pi_beta = vapply(blocks, `[[`, numeric(1), "pi_beta"),
    spike_var_shape = vapply(blocks, `[[`, numeric(1), "spike_var_shape"),
    spike_var_scale = vapply(blocks, `[[`, numeric(1), "spike_var_scale"),
    global_scale = vapply(blocks, `[[`, numeric(1), "global_scale"),
    fixed_global_var = vapply(
      blocks, `[[`, numeric(1), "fixed_global_var"
    ),
    local_a = vapply(
      blocks, function(block) block$local_shape[1L], numeric(1)
    ),
    local_b = vapply(
      blocks, function(block) block$local_shape[2L], numeric(1)
    ),
    multi_gamma = lapply(blocks, `[[`, "multi_gamma"),
    multi_pi_alpha = lapply(blocks, `[[`, "multi_pi_alpha"),
    multi_var_shape = vapply(blocks, `[[`, numeric(1), "multi_var_shape"),
    multi_var_scale = vapply(blocks, `[[`, numeric(1), "multi_var_scale"),
    fixed_var = vapply(blocks, `[[`, numeric(1), "fixed_var")
  )
}

.prepare_sampler_arguments <- function(
    blocks, layout, y, x, residual_shape, residual_scale, residual_var,
    iterations, burnin, thin, store_samples, store_coefficient_cov,
    compute_pve, pve_type, effective_n, likelihood_df, fit_intercept,
    intercept_x_mean, intercept_y_mean, center_observations = FALSE,
    residual_sse_offset = 0) {
  c(
    list(
      y = y,
      x = x,
      residual_shape = if (is.null(residual_shape)) 1 else residual_shape,
      residual_scale = if (is.null(residual_scale)) 1 else residual_scale,
      residual_var = residual_var,
      iterations = iterations,
      burnin = burnin,
      thin = thin,
      block_id = layout$block_id,
      block_model = layout$block_model
    ),
    .prior_sampler_arguments(blocks),
    list(
      store_samples = store_samples,
      store_coefficient_cov = store_coefficient_cov,
      compute_pve = compute_pve,
      pve_type = pve_type,
      effective_n = effective_n,
      likelihood_df = likelihood_df,
      fit_intercept = fit_intercept,
      intercept_x_mean = intercept_x_mean,
      intercept_y_mean = intercept_y_mean,
      center_observations = center_observations,
      residual_sse_offset = residual_sse_offset
    )
  )
}

.run_prepared_sampler <- function(sampler_arguments, version, nchains,
                                  block_model, verbose, iterations) {
  run_chains <- function(progressor = NULL) {
    .run_blm_chains(
      sampler_arguments, version, nchains, block_model, progressor
    )
  }
  if (verbose) {
    progressr::with_progress({
      progress <- progressr::progressor(steps = nchains * iterations)
      run_chains(progress)
    }, enable = TRUE)
  } else {
    run_chains()
  }
}
