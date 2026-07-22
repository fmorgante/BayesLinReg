.assemble_blm_result <- function(blocks, block_indices, samples, nchains,
                                 store_samples, store_coefficient_cov,
                                 fit_intercept = TRUE) {
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
      result$slab_shape <- block$slab_shape
      result$slab_scale <- block$slab_scale
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
