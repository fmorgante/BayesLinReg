.assemble_blm_result <- function(blocks, block_indices, samples, nchains,
                                 store_samples, store_coefficient_cov,
                                 fit_intercept = TRUE,
                                 residual_shape = NULL,
                                 residual_scale = NULL,
                                 residual_scale_calibrated = FALSE,
                                 expected_pve_total = NULL,
                                 reference_response_var = NULL,
                                 reference_residual_var = NULL) {
  eta_result <- lapply(seq_along(blocks), function(block_index) {
    block <- blocks[[block_index]]
    indices <- block_indices[[block_index]]
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
      coefficient_mean <- samples$coefficient_sum[indices] /
        samples$number_of_draws / block$predictor_scale
      coefficient_var <- vapply(indices, function(index) {
        .variance_from_sums(
          samples$coefficient_sum[index],
          samples$coefficient_sum_sq[index],
          samples$number_of_draws
        )
      }, numeric(1)) / block$predictor_scale^2
      if (store_coefficient_cov) {
        coefficient_cov <- .covariance_from_sums(
          samples$coefficient_sum[indices],
          samples$coefficient_crossprod[indices, indices, drop = FALSE],
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
    if (block$model == "Normal") {
      if (store_samples) {
        normal_var_samples <- samples$normal_var_samples[, block_index]
        result$normal_var_mean <- mean(normal_var_samples)
        result$normal_var_var <- stats::var(normal_var_samples)
        result$normal_var_samples <- normal_var_samples
      } else {
        result$normal_var_mean <- samples$normal_var_sum[block_index] /
          samples$number_of_draws
        result$normal_var_var <- .variance_from_sums(
          samples$normal_var_sum[block_index],
          samples$normal_var_sum_sq[block_index],
          samples$number_of_draws
        )
      }
      result$var_shape <- block$normal_shape
      result$var_scale <- block$normal_scale
      result$var_scale_calibrated <- block$normal_scale_calibrated
      if (block$normal_scale_calibrated) {
        result$expected_pve <- block$expected_pve
        result$reference_response_var <- block$reference_response_var
      }
    }
    if (block$model == "SpikeSlab") {
      if (store_samples) {
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
      } else {
        result$inclusion_probability <- samples$inclusion_sum[indices] /
          samples$number_of_draws
        names(result$inclusion_probability) <- block$predictor_names
        result$pi_mean <- samples$pi_sum[block_index] / samples$number_of_draws
        result$pi_var <- .variance_from_sums(
          samples$pi_sum[block_index], samples$pi_sum_sq[block_index],
          samples$number_of_draws
        )
        result$slab_var_mean <- samples$slab_var_sum[block_index] /
          samples$number_of_draws
        result$slab_var_var <- .variance_from_sums(
          samples$slab_var_sum[block_index],
          samples$slab_var_sum_sq[block_index], samples$number_of_draws
        )
      }
      result$pi <- c(a = block$pi_alpha, b = block$pi_beta)
      result$var_shape <- block$spike_var_shape
      result$var_scale <- block$spike_var_scale
      result$var_scale_calibrated <- block$spike_var_scale_calibrated
      if (block$spike_var_scale_calibrated) {
        result$expected_pve <- block$expected_pve
        result$reference_response_var <- block$reference_response_var
      }
    }
    if (block$model == "GlobalLocal") {
      if (store_samples) {
        local_var_samples <- samples$local_var_samples[, indices, drop = FALSE]
        colnames(local_var_samples) <- block$predictor_names
        tau_sq_samples <- samples$tau_sq_samples[, block_index]
        result$local_var_mean <- colMeans(local_var_samples)
        result$local_var_var <- apply(local_var_samples, 2L, stats::var)
        result$tau_sq_mean <- mean(tau_sq_samples)
        result$tau_sq_var <- stats::var(tau_sq_samples)
        result$local_var_samples <- local_var_samples
        result$tau_sq_samples <- tau_sq_samples
      } else {
        result$local_var_mean <- samples$local_var_sum[indices] /
          samples$number_of_draws
        result$local_var_var <- vapply(indices, function(index) {
          .variance_from_sums(
            samples$local_var_sum[index], samples$local_var_sum_sq[index],
            samples$number_of_draws
          )
        }, numeric(1))
        names(result$local_var_mean) <- block$predictor_names
        names(result$local_var_var) <- block$predictor_names
        result$tau_sq_mean <- samples$tau_sq_sum[block_index] /
          samples$number_of_draws
        result$tau_sq_var <- .variance_from_sums(
          samples$tau_sq_sum[block_index],
          samples$tau_sq_sum_sq[block_index], samples$number_of_draws
        )
      }
      result$local_shape <- block$local_shape
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
    }
    if (block$model == "SpikeMultiSlab") {
      component_names <- c(
        "spike", paste0("slab_", seq_len(length(block$multi_gamma) - 1L))
      )
      if (store_samples) {
        component_samples <-
          samples$multi_component_samples[, indices, drop = FALSE]
        colnames(component_samples) <- block$predictor_names
        component_probability <- vapply(
          seq_along(block$multi_gamma),
          function(component) colMeans(component_samples == component),
          numeric(length(indices))
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
          samples$multi_component_sum[[block_index]] /
            samples$number_of_draws
        dimnames(component_probability) <- list(
          block$predictor_names, component_names
        )
        result$pi_mean <- samples$multi_pi_sum[[block_index]] /
          samples$number_of_draws
        result$pi_var <- vapply(seq_along(block$multi_gamma), function(index) {
          .variance_from_sums(
            samples$multi_pi_sum[[block_index]][index],
            samples$multi_pi_sum_sq[[block_index]][index],
            samples$number_of_draws
          )
        }, numeric(1))
        names(result$pi_mean) <- names(result$pi_var) <- component_names
        result$var_mean <- samples$multi_var_sum[block_index] /
          samples$number_of_draws
        result$var_var <- .variance_from_sums(
          samples$multi_var_sum[block_index],
          samples$multi_var_sum_sq[block_index], samples$number_of_draws
        )
      }
      result$component_probability <- component_probability
      result$inclusion_probability <- 1 - component_probability[, "spike"]
      result$gamma <- stats::setNames(block$multi_gamma, component_names)
      result$alpha <- stats::setNames(block$multi_pi_alpha, component_names)
      result$var_shape <- block$multi_var_shape
      result$var_scale <- block$multi_var_scale
      result$var_scale_calibrated <- block$multi_var_scale_calibrated
      if (block$multi_var_scale_calibrated) {
        result$expected_pve <- block$expected_pve
        result$reference_response_var <- block$reference_response_var
      }
    }
    result
  })
  names(eta_result) <- names(blocks)

  if (store_samples) {
    residual_var_mean <- mean(samples$residual_var_samples)
    residual_var_var <- stats::var(samples$residual_var_samples)
  } else {
    residual_var_mean <- samples$residual_var_sum / samples$number_of_draws
    residual_var_var <- .variance_from_sums(
      samples$residual_var_sum, samples$residual_var_sum_sq,
      samples$number_of_draws
    )
  }
  result <- list(ETA = eta_result)
  if (fit_intercept) {
    if (store_samples) {
      result$intercept_mean <- mean(samples$intercept_samples)
      result$intercept_var <- stats::var(samples$intercept_samples)
      result$intercept_samples <- samples$intercept_samples
    } else {
      result$intercept_mean <- samples$intercept_sum / samples$number_of_draws
      result$intercept_var <- .variance_from_sums(
        samples$intercept_sum, samples$intercept_sum_sq,
        samples$number_of_draws
      )
    }
  }
  result$residual_var_mean <- residual_var_mean
  result$residual_var_var <- residual_var_var
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
  if (store_samples) {
    result$residual_var_samples <- samples$residual_var_samples
  }
  if (nchains > 1L) {
    result$nchains <- nchains
    if (store_samples) result$chain_id <- samples$chain_id
  }
  structure(result, class = "blm_fit")
}
