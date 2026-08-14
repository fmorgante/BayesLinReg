.assemble_blm_result <- function(blocks, block_indices, samples, nchains,
                                 store_samples, store_coefficient_cov,
                                 fit_intercept = TRUE,
                                 compute_pve = FALSE,
                                 pve_type = "standalone",
                                 residual_shape = NULL,
                                 residual_scale = NULL,
                                 residual_scale_calibrated = FALSE,
                                 expected_pve_total = NULL,
                                 reference_response_var = NULL,
                                 reference_residual_var = NULL,
                                 likelihood_df = NULL,
                                 sampler_block_id = NULL) {
  block_models <- vapply(blocks, `[[`, character(1), "model")
  model_names <- c(
    "Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab", "Fixed"
  )
  family_indices <- if (is.null(sampler_block_id)) {
    lapply(model_names, function(model) {
      unlist(block_indices[block_models == model], use.names = FALSE)
    })
  } else {
    lapply(model_names, function(model) {
      which(block_models[sampler_block_id] == model)
    })
  }
  eta_result <- lapply(seq_along(blocks), function(block_index) {
    block <- blocks[[block_index]]
    indices <- block_indices[[block_index]]
    block_local_order <- match(indices, sort(indices))
    if (store_samples) {
      coefficient_samples <- sweep(
        samples$coefficient_samples[, indices, drop = FALSE],
        2L, block$predictor_scale, FUN = "/"
      )
      colnames(coefficient_samples) <- block$predictor_names
      coefficient_mean <- colMeans(coefficient_samples)
      coefficient_var <- apply(coefficient_samples, 2L, stats::var)
      if (store_coefficient_cov) {
        coefficient_cov <- stats::cov(coefficient_samples)
      }
    } else {
      coefficient_mean <- samples$coefficient_mean[indices] /
        block$predictor_scale
      coefficient_var <- .variance_from_m2(
        samples$coefficient_m2[indices], samples$number_of_draws
      ) / block$predictor_scale^2
      if (store_coefficient_cov) {
        coefficient_cov <- .covariance_from_m2(
          samples$coefficient_cov_m2[[block_index]][
            block_local_order, block_local_order, drop = FALSE
          ],
          samples$number_of_draws
        ) / outer(block$predictor_scale, block$predictor_scale)
      }
      names(coefficient_mean) <- block$predictor_names
      names(coefficient_var) <- block$predictor_names
      if (store_coefficient_cov) {
        dimnames(coefficient_cov) <- list(
          block$predictor_names, block$predictor_names
        )
      }
    }
    result <- list(
      model = block$model,
      standardize = block$standardize,
      coefficient_mean = coefficient_mean,
      coefficient_var = coefficient_var
    )
    if (store_coefficient_cov) result$coefficient_cov <- coefficient_cov
    if (store_samples) result$coefficient_samples <- coefficient_samples
    if (compute_pve) {
      if (store_samples) {
        pve_samples <- samples$block_pve_samples[, block_index]
        result$pve_mean <- mean(pve_samples)
        result$pve_var <- stats::var(pve_samples)
        result$pve_samples <- pve_samples
      } else {
        result$pve_mean <- samples$block_pve_mean[block_index]
        result$pve_var <- .variance_from_m2(
          samples$block_pve_m2[block_index],
          samples$number_of_draws
        )
      }
    }
    if (block$model == "Normal") {
      if (store_samples) {
        normal_var_samples <- samples$normal_var_samples[, block_index]
        result$normal_var_mean <- mean(normal_var_samples)
        result$normal_var_var <- stats::var(normal_var_samples)
        result$normal_var_samples <- normal_var_samples
      } else {
        result$normal_var_mean <- samples$normal_var_mean[block_index]
        result$normal_var_var <- .variance_from_m2(
          samples$normal_var_m2[block_index],
          samples$number_of_draws
        )
      }
      if (!is.na(block$fixed_var)) {
        result$normal_var_mean <- block$fixed_var
        result$normal_var_var <- 0
      }
      if (is.na(block$fixed_var)) {
        result$var_shape <- block$normal_shape
        result$var_scale <- block$normal_scale
        result$var_scale_calibrated <- block$normal_scale_calibrated
        if (block$normal_scale_calibrated) {
          result$expected_pve <- block$expected_pve
          result$reference_response_var <- block$reference_response_var
        }
      } else {
        result$var <- block$fixed_var
      }
    }
    if (block$model == "SpikeSlab") {
      local_indices <- match(indices, family_indices[[2L]])
      if (store_samples) {
        inclusion_samples <-
          samples$inclusion_samples[, local_indices, drop = FALSE]
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
      } else {
        result$inclusion_probability <- samples$inclusion_sum[local_indices] /
          samples$number_of_draws
        names(result$inclusion_probability) <- block$predictor_names
        result$pi_mean <- samples$pi_mean[block_index]
        result$pi_var <- .variance_from_m2(
          samples$pi_m2[block_index],
          samples$number_of_draws
        )
        result$slab_var_mean <- samples$slab_var_mean[block_index]
        result$slab_var_var <- .variance_from_m2(
          samples$slab_var_m2[block_index], samples$number_of_draws
        )
      }
      if (!is.na(block$fixed_var)) {
        result$slab_var_mean <- block$fixed_var
        result$slab_var_var <- 0
      }
      result$pi <- c(a = block$pi_alpha, b = block$pi_beta)
      if (is.na(block$fixed_var)) {
        result$var_shape <- block$spike_var_shape
        result$var_scale <- block$spike_var_scale
        result$var_scale_calibrated <- block$spike_var_scale_calibrated
        if (block$spike_var_scale_calibrated) {
          result$expected_pve <- block$expected_pve
          result$reference_response_var <- block$reference_response_var
        }
      } else {
        result$var <- block$fixed_var
      }
    }
    if (block$model == "GlobalLocal") {
      local_indices <- match(indices, family_indices[[3L]])
      if (store_samples) {
        local_var_samples <-
          samples$local_var_samples[, local_indices, drop = FALSE]
        colnames(local_var_samples) <- block$predictor_names
        tau_sq_samples <- samples$tau_sq_samples[, block_index]
        result$local_var_mean <- colMeans(local_var_samples)
        result$local_var_var <- apply(local_var_samples, 2L, stats::var)
        result$tau_sq_mean <- mean(tau_sq_samples)
        result$tau_sq_var <- stats::var(tau_sq_samples)
        result$local_var_samples <- local_var_samples
        result$tau_sq_samples <- tau_sq_samples
      } else {
        result$local_var_mean <- samples$local_var_mean[local_indices]
        result$local_var_var <- .variance_from_m2(
          samples$local_var_m2[local_indices], samples$number_of_draws
        )
        names(result$local_var_mean) <- block$predictor_names
        names(result$local_var_var) <- block$predictor_names
        result$tau_sq_mean <- samples$tau_sq_mean[block_index]
        result$tau_sq_var <- .variance_from_m2(
          samples$tau_sq_m2[block_index], samples$number_of_draws
        )
      }
      if (!is.na(block$fixed_global_var)) {
        result$tau_sq_mean <- block$fixed_global_var
        result$tau_sq_var <- 0
      }
      result$local_shape <- block$local_shape
      if (is.na(block$fixed_global_var)) {
        result$global_scale <- block$global_scale
        result$global_scale_calibrated <- block$global_scale_calibrated
        result$global_scale_calibration <- block$global_scale_calibration
        if (block$global_scale_calibrated) {
          result$expected_nonzero <- block$expected_nonzero
          result$reference_residual_var <- block$reference_residual_var
        }
        if (identical(block$global_scale_calibration, "expected_pve")) {
          result$expected_pve <- block$expected_pve
          result$reference_response_var <- block$reference_response_var
        }
      } else {
        result$global_var <- block$fixed_global_var
      }
    }
    if (block$model == "SpikeMultiSlab") {
      local_indices <- match(indices, family_indices[[4L]])
      component_names <- c(
        "spike", paste0("slab_", seq_len(length(block$multi_gamma) - 1L))
      )
      if (store_samples) {
        component_samples <-
          samples$multi_component_samples[, local_indices, drop = FALSE]
        colnames(component_samples) <- block$predictor_names
        component_probability <- matrix(
          vapply(
            seq_along(block$multi_gamma),
            function(component) colMeans(component_samples == component),
            numeric(length(indices))
          ),
          nrow = length(indices),
          ncol = length(block$multi_gamma)
        )
        dimnames(component_probability) <- list(
          block$predictor_names, component_names
        )
        multi_pi_samples <- samples$multi_pi_samples[[block_index]]
        colnames(multi_pi_samples) <- component_names
        multi_var_samples <- samples$multi_var_samples[, block_index]
        result$component_samples <- component_samples
        result$pi_samples <- multi_pi_samples
        result$var_samples <- multi_var_samples
        result$pi_mean <- colMeans(multi_pi_samples)
        result$pi_var <- apply(multi_pi_samples, 2L, stats::var)
        result$var_mean <- mean(multi_var_samples)
        result$var_var <- stats::var(multi_var_samples)
      } else {
        component_probability <-
          samples$multi_component_sum[[block_index]][
            block_local_order, , drop = FALSE
          ] /
            samples$number_of_draws
        dimnames(component_probability) <- list(
          block$predictor_names, component_names
        )
        result$pi_mean <- samples$multi_pi_mean[[block_index]]
        result$pi_var <- .variance_from_m2(
          samples$multi_pi_m2[[block_index]], samples$number_of_draws
        )
        names(result$pi_mean) <- names(result$pi_var) <- component_names
        result$var_mean <- samples$multi_var_mean[block_index]
        result$var_var <- .variance_from_m2(
          samples$multi_var_m2[block_index], samples$number_of_draws
        )
      }
      if (!is.na(block$fixed_var)) {
        result$var_mean <- block$fixed_var
        result$var_var <- 0
      }
      result$component_probability <- component_probability
      result$inclusion_probability <- 1 - component_probability[, "spike"]
      result$gamma <- stats::setNames(block$multi_gamma, component_names)
      result$alpha <- stats::setNames(block$multi_pi_alpha, component_names)
      if (is.na(block$fixed_var)) {
        result$var_shape <- block$multi_var_shape
        result$var_scale <- block$multi_var_scale
        result$var_scale_calibrated <- block$multi_var_scale_calibrated
        if (block$multi_var_scale_calibrated) {
          result$expected_pve <- block$expected_pve
          result$reference_response_var <- block$reference_response_var
        }
      } else {
        result$var <- block$fixed_var
      }
    }
    result
  })
  names(eta_result) <- names(blocks)

  if (store_samples) {
    residual_var_mean <- mean(samples$residual_var_samples)
    residual_var_var <- stats::var(samples$residual_var_samples)
  } else {
    residual_var_mean <- samples$residual_var_mean
    residual_var_var <- .variance_from_m2(
      samples$residual_var_m2, samples$number_of_draws
    )
  }
  result <- list(ETA = eta_result)
  if (fit_intercept) {
    if (store_samples) {
      result$intercept_mean <- mean(samples$intercept_samples)
      result$intercept_var <- stats::var(samples$intercept_samples)
      result$intercept_samples <- samples$intercept_samples
    } else {
      result$intercept_mean <- samples$intercept_mean
      result$intercept_var <- .variance_from_m2(
        samples$intercept_m2, samples$number_of_draws
      )
    }
  }
  result$residual_var_mean <- residual_var_mean
  result$residual_var_var <- residual_var_var
  if (compute_pve) {
    if (store_samples) {
      result$total_pve_mean <- mean(samples$total_pve_samples)
      result$total_pve_var <- stats::var(samples$total_pve_samples)
      result$cross_block_pve_mean <- mean(samples$cross_block_pve_samples)
      result$cross_block_pve_var <- stats::var(samples$cross_block_pve_samples)
      result$total_pve_samples <- samples$total_pve_samples
      result$cross_block_pve_samples <- samples$cross_block_pve_samples
    } else {
      result$total_pve_mean <- samples$total_pve_mean
      result$total_pve_var <- .variance_from_m2(
        samples$total_pve_m2,
        samples$number_of_draws
      )
      result$cross_block_pve_mean <- samples$cross_block_pve_mean
      result$cross_block_pve_var <- .variance_from_m2(
        samples$cross_block_pve_m2,
        samples$number_of_draws
      )
    }
    result$pve_type <- pve_type
  }
  if (!is.null(residual_shape)) {
    result$residual_shape <- residual_shape
    result$residual_scale <- residual_scale
    result$residual_scale_calibrated <- residual_scale_calibrated
    if (residual_scale_calibrated) {
      result$expected_pve_total <- expected_pve_total
      result$reference_response_var <- reference_response_var
      result$reference_residual_var <- reference_residual_var
    }
  }
  result$store_samples <- store_samples
  result$store_coefficient_cov <- store_coefficient_cov
  if (!is.null(likelihood_df)) result$likelihood_df <- likelihood_df
  if (store_samples) {
    result$residual_var_samples <- samples$residual_var_samples
  }
  if (nchains > 1L) {
    result$nchains <- nchains
    if (store_samples) result$chain_id <- samples$chain_id
  }
  structure(result, class = "blm_fit")
}
