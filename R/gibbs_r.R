# Reference R Gibbs sampler and its distribution helpers.

.draw_gig <- function(n, lambda, chi, psi) {
  values <- c(n = n, lambda = lambda, chi = chi, psi = psi)
  if (!is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
      n != floor(n) || n < 1 || chi <= 0 || psi <= 0) {
    stop(
      paste0(
        "GIG parameters require a positive integer `n`, finite `lambda`, ",
        "and positive finite `chi` and `psi`."
      ),
      call. = FALSE
    )
  }
  GIGrvg::rgig(n = n, lambda = lambda, chi = chi, psi = psi)
}

.guard_reconstructed_sse <- function(sse, yty, linear, quadratic,
                                     iteration) {
  scale <- max(1, abs(yty), 2 * abs(linear), abs(quadratic))
  tolerance <- sqrt(.Machine$double.eps) * scale
  if (!is.finite(sse)) {
    stop(
      sprintf(
        "The reconstructed residual SSE is non-finite at Gibbs iteration %d.",
        iteration
      ),
      call. = FALSE
    )
  }
  if (sse < -tolerance) {
    stop(
      sprintf(
        paste0(
          "The reconstructed residual SSE is materially negative at Gibbs ",
          "iteration %d (SSE %.6g; tolerance %.6g). The supplied ",
          "sufficient statistics are incompatible."
        ),
        iteration, sse, tolerance
      ),
      call. = FALSE
    )
  }
  max(0, sse)
}

.blm_gibbs <- function(y, x, residual_shape, residual_scale,
                       iterations, burnin, thin,
                       progress_callback = NULL,
                       block_id = NULL, block_model = 0L,
                       normal_shape = 2, normal_scale = 1,
                       pi_alpha = 1, pi_beta = 1,
                       spike_var_shape = 2, spike_var_scale = 1,
                       global_scale = 1, fixed_global_var = NULL,
                       residual_var = NULL,
                       local_a = 1, local_b = 0.5,
                       multi_gamma = list(c(0, 0.01, 0.1, 1)),
                       multi_pi_alpha = list(rep(1, 4)),
                       multi_var_shape = 2, multi_var_scale = 1,
                       fixed_var = NULL,
                       store_samples = TRUE,
                       store_coefficient_cov = TRUE,
                       effective_n = NULL, fit_intercept = TRUE,
                       intercept_x_mean = NULL, intercept_y_mean = NULL,
                       XtX = NULL, XtX_center = NULL, Xty = NULL, yty = NULL,
                       center_observations = TRUE,
                       residual_sse_offset = 0, compute_pve = FALSE,
                       pve_type = c("standalone", "allocated")) {
  pve_controls <- .validate_pve_controls(compute_pve, pve_type)
  compute_pve <- pve_controls$compute_pve
  pve_type <- pve_controls$pve_type
  retained_iterations <- .validate_mcmc(iterations, burnin, thin)
  use_sufficient_statistics <- !is.null(XtX)
  number_of_predictors <- if (use_sufficient_statistics) {
    ncol(XtX)
  } else {
    ncol(x)
  }
  predictor_names <- if (use_sufficient_statistics) {
    colnames(XtX)
  } else {
    colnames(x)
  }
  if (is.null(block_id)) {
    block_id <- rep.int(1L, number_of_predictors)
  }
  number_of_blocks <- length(block_model)
  if (is.null(fixed_var)) fixed_var <- rep(NA_real_, number_of_blocks)
  if (is.null(fixed_global_var)) {
    fixed_global_var <- rep(NA_real_, number_of_blocks)
  }
  block_predictors <- lapply(seq_len(number_of_blocks), function(block) {
    which(block_id == block)
  })
  block_XtX <- if (compute_pve && use_sufficient_statistics) {
    lapply(block_predictors, function(predictors) {
      XtX[predictors, predictors, drop = FALSE]
    })
  } else {
    NULL
  }
  model_predictors <- lapply(0:4, function(model) {
    unlist(block_predictors[block_model == model], use.names = FALSE)
  })
  model_local_index <- integer(number_of_predictors)
  for (model in 0:4) {
    predictors <- model_predictors[[model + 1L]]
    model_local_index[predictors] <- seq_along(predictors)
  }

  if (use_sufficient_statistics) {
    x_squared <- diag(XtX)
    corrected_rhs <- as.numeric(Xty)
    residual_sse <- yty
  } else {
    x_mean <- if (center_observations) colMeans(x) else rep(0, ncol(x))
    y_mean <- if (center_observations) mean(y) else 0
    if (is.null(effective_n)) effective_n <- length(y)
    if (is.null(intercept_x_mean)) intercept_x_mean <- x_mean
    if (is.null(intercept_y_mean)) intercept_y_mean <- y_mean
    x_centered <- if (center_observations) {
      sweep(x, 2L, x_mean, FUN = "-")
    } else {
      x
    }
    y_centered <- if (center_observations) y - y_mean else y
    x_squared <- colSums(x_centered^2)
  }

  number_of_draws <- length(retained_iterations)
  if (store_samples) {
    coefficient_samples <- matrix(
      NA_real_,
      nrow = number_of_draws,
      ncol = number_of_predictors,
      dimnames = list(NULL, predictor_names)
    )
    intercept_samples <- numeric(number_of_draws)
    residual_var_samples <- numeric(number_of_draws)
  } else {
    coefficient_sum <- numeric(number_of_predictors)
    coefficient_sum_sq <- numeric(number_of_predictors)
    if (store_coefficient_cov) {
      coefficient_crossprod <- lapply(block_predictors, function(predictors) {
        matrix(0, length(predictors), length(predictors))
      })
    }
    intercept_sum <- intercept_sum_sq <- 0
    residual_var_sum <- residual_var_sum_sq <- 0
  }
  if (compute_pve) {
    if (store_samples) {
      block_pve_samples <- matrix(
        NA_real_, number_of_draws, number_of_blocks
      )
      total_pve_samples <- cross_block_pve_samples <-
        numeric(number_of_draws)
    } else {
      block_pve_sum <- block_pve_sum_sq <- numeric(number_of_blocks)
      total_pve_sum <- total_pve_sum_sq <- 0
      cross_block_pve_sum <- cross_block_pve_sum_sq <- 0
    }
  }
  has_normal <- any(block_model == 0L)
  has_spike_slab <- any(block_model == 1L)
  has_global_local <- any(block_model == 2L)
  has_spike_multi_slab <- any(block_model == 3L)
  if (has_normal) {
    if (store_samples) {
      normal_var_samples <- matrix(NA_real_, number_of_draws, number_of_blocks)
    } else {
      normal_var_sum <- normal_var_sum_sq <- numeric(number_of_blocks)
    }
    normal_var <- ifelse(
      is.na(fixed_var), normal_scale / (normal_shape + 1), fixed_var
    )
  }
  if (has_spike_slab) {
    if (store_samples) {
      inclusion_samples <- matrix(
        NA_integer_,
        nrow = number_of_draws,
        ncol = length(model_predictors[[2L]]),
        dimnames = list(NULL, predictor_names[model_predictors[[2L]]])
      )
      pi_samples <- matrix(NA_real_, number_of_draws, number_of_blocks)
      slab_var_samples <- matrix(NA_real_, number_of_draws, number_of_blocks)
    } else {
      inclusion_sum <- numeric(length(model_predictors[[2L]]))
      pi_sum <- pi_sum_sq <- numeric(number_of_blocks)
      slab_var_sum <- slab_var_sum_sq <- numeric(number_of_blocks)
    }
    inclusion <- rep.int(1L, length(model_predictors[[2L]]))
    pi <- pi_alpha / (pi_alpha + pi_beta)
    slab_var <- ifelse(
      is.na(fixed_var), spike_var_scale / (spike_var_shape + 1), fixed_var
    )
  }
  if (has_global_local) {
    if (store_samples) {
      local_var_samples <- matrix(
        NA_real_,
        nrow = number_of_draws,
        ncol = length(model_predictors[[3L]]),
        dimnames = list(NULL, predictor_names[model_predictors[[3L]]])
      )
      tau_sq_samples <- matrix(NA_real_, number_of_draws, number_of_blocks)
    } else {
      local_var_sum <- local_var_sum_sq <-
        numeric(length(model_predictors[[3L]]))
      tau_sq_sum <- tau_sq_sum_sq <- numeric(number_of_blocks)
    }
    local_var <- rep(1, length(model_predictors[[3L]]))
    local_aux <- rep(1, length(model_predictors[[3L]]))
    tau_sq <- ifelse(
      is.na(fixed_global_var), global_scale^2, fixed_global_var
    )
    global_aux <- rep(1, number_of_blocks)
  }
  if (has_spike_multi_slab) {
    multi_component <- rep.int(1L, length(model_predictors[[4L]]))
    multi_pi <- lapply(seq_len(number_of_blocks), function(block) {
      multi_pi_alpha[[block]] / sum(multi_pi_alpha[[block]])
    })
    multi_var <- ifelse(
      is.na(fixed_var), multi_var_scale / (multi_var_shape + 1), fixed_var
    )
    if (store_samples) {
      multi_component_samples <- matrix(
        NA_integer_, number_of_draws, length(model_predictors[[4L]]),
        dimnames = list(NULL, predictor_names[model_predictors[[4L]]])
      )
      multi_pi_samples <- lapply(seq_len(number_of_blocks), function(block) {
        if (block_model[block] == 3L) {
          matrix(NA_real_, number_of_draws, length(multi_gamma[[block]]))
        } else {
          NULL
        }
      })
      multi_var_samples <- matrix(NA_real_, number_of_draws, number_of_blocks)
    } else {
      multi_component_sum <- lapply(seq_len(number_of_blocks), function(block) {
        if (block_model[block] == 3L) {
          matrix(
            0, sum(block_id == block), length(multi_gamma[[block]])
          )
        } else {
          NULL
        }
      })
      multi_pi_sum <- lapply(multi_pi, function(value) numeric(length(value)))
      multi_pi_sum_sq <- lapply(multi_pi, function(value) numeric(length(value)))
      multi_var_sum <- multi_var_sum_sq <- numeric(number_of_blocks)
    }
  }

  coefficient <- numeric(number_of_predictors)
  if (!use_sufficient_statistics) residuals <- y_centered
  learn_residual_var <- is.null(residual_var)
  if (learn_residual_var) {
    residual_var <- residual_scale / (residual_shape + 1)
    residual_posterior_shape <- residual_shape +
      (effective_n - as.integer(fit_intercept)) / 2
  }
  retained_index <- 1L
  progress_thresholds <- if (!is.null(progress_callback)) {
    unique(pmax(
      1L,
      as.integer((iterations * seq_len(10L) + 9) %/% 10)
    ))
  } else {
    integer(0)
  }
  progress_index <- 1L
  last_reported_iteration <- 0L
  residual_refresh_interval <- 100L

  for (iteration in seq_len(iterations)) {
    for (predictor in seq_len(number_of_predictors)) {
      block <- block_id[predictor]
      model <- block_model[block]
      local_index <- model_local_index[predictor]
      old_coefficient <- coefficient[predictor]
      if (use_sufficient_statistics) {
        partial_rhs <- corrected_rhs[predictor] +
          x_squared[predictor] * old_coefficient
        conditional_numerator <- partial_rhs / residual_var
      } else {
        partial_residuals <- residuals +
          x_centered[, predictor] * old_coefficient
        conditional_numerator <-
          sum(x_centered[, predictor] * partial_residuals) / residual_var
      }
      if (model == 3L) {
        gamma <- multi_gamma[[block]]
        log_weights <- log(pmax(multi_pi[[block]], .Machine$double.xmin))
        conditional_vars <- conditional_means <- rep(0, length(gamma))
        for (component in seq.int(2L, length(gamma))) {
          prior_var <- gamma[component] * multi_var[block]
          conditional_vars[component] <- 1 / (
            x_squared[predictor] / residual_var + 1 / prior_var
          )
          conditional_means[component] <-
            conditional_vars[component] * conditional_numerator
          log_weights[component] <- log_weights[component] +
            0.5 * log(conditional_vars[component] / prior_var) +
            conditional_means[component]^2 /
              (2 * conditional_vars[component])
        }
        probabilities <- exp(log_weights - max(log_weights))
        probabilities <- probabilities / sum(probabilities)
        multi_component[local_index] <- sample.int(
          length(gamma), 1L, prob = probabilities
        )
        component <- multi_component[local_index]
        coefficient[predictor] <- if (component == 1L) {
          0
        } else {
          stats::rnorm(
            1L, conditional_means[component],
            sqrt(conditional_vars[component])
          )
        }
      } else {
        prior_precision <- if (model == 4L) {
          0
        } else if (model == 2L) {
          1 / tau_sq[block] / local_var[local_index]
        } else if (model == 1L) {
          1 / slab_var[block]
        } else {
          1 / normal_var[block]
        }
        conditional_var <- 1 / (
          x_squared[predictor] / residual_var + prior_precision
        )
        conditional_mean <- conditional_var * conditional_numerator
        if (model == 1L) {
          bounded_pi <- min(
            max(pi[block], .Machine$double.eps),
            1 - .Machine$double.eps
          )
          log_inclusion_odds <- stats::qlogis(bounded_pi) +
            0.5 * log(conditional_var / slab_var[block]) +
            conditional_mean^2 / (2 * conditional_var)
          inclusion[local_index] <- stats::rbinom(
            1L,
            size = 1L,
            prob = stats::plogis(log_inclusion_odds)
          )
        }
        if (model != 1L || inclusion[local_index] == 1L) {
          coefficient[predictor] <- stats::rnorm(
            1L,
            mean = conditional_mean,
            sd = sqrt(conditional_var)
          )
        } else {
          coefficient[predictor] <- 0
        }
      }
      if (use_sufficient_statistics) {
        coefficient_change <- coefficient[predictor] - old_coefficient
        if (coefficient_change != 0) {
          corrected_rhs <- corrected_rhs -
            XtX[, predictor] * coefficient_change
          if (learn_residual_var) {
            residual_sse <- residual_sse -
              2 * coefficient_change * partial_rhs +
              (coefficient[predictor]^2 - old_coefficient^2) *
                x_squared[predictor]
          }
        }
      } else {
        residuals <- partial_residuals -
          x_centered[, predictor] * coefficient[predictor]
      }
    }

    # Reconstruct accumulated state periodically to limit floating-point drift.
    if (iteration %% residual_refresh_interval == 0L) {
      if (use_sufficient_statistics) {
        fitted_crossproducts <- drop(XtX %*% coefficient)
        corrected_rhs <- Xty - fitted_crossproducts
        if (learn_residual_var) {
          linear <- sum(coefficient * Xty)
          quadratic <- sum(coefficient * fitted_crossproducts)
          residual_sse <- .guard_reconstructed_sse(
            yty - 2 * linear + quadratic,
            yty, linear, quadratic, iteration
          )
        }
      } else {
        residuals <- y_centered - drop(x_centered %*% coefficient)
      }
    }

    if (has_normal) {
      for (block in which(block_model == 0L)) {
        if (!is.na(fixed_var[block])) next
        predictors <- block_predictors[[block]]
        normal_var[block] <- 1 / stats::rgamma(
          1L,
          shape = normal_shape[block] + length(predictors) / 2,
          rate = normal_scale[block] + sum(coefficient[predictors]^2) / 2
        )
      }
    }

    if (has_spike_slab) {
      for (block in which(block_model == 1L)) {
        predictors <- block_predictors[[block]]
        local_indices <- model_local_index[predictors]
        number_included <- sum(inclusion[local_indices])
        pi[block] <- stats::rbeta(
          1L,
          shape1 = pi_alpha[block] + number_included,
          shape2 = pi_beta[block] + length(predictors) - number_included
        )
        if (is.na(fixed_var[block])) {
          included_predictors <- predictors[inclusion[local_indices] == 1L]
          slab_var[block] <- 1 / stats::rgamma(
            1L,
            shape = spike_var_shape[block] + length(included_predictors) / 2,
            rate = spike_var_scale[block] +
              sum(coefficient[included_predictors]^2) / 2
          )
        }
      }
    }

    if (has_spike_multi_slab) {
      for (block in which(block_model == 3L)) {
        predictors <- block_predictors[[block]]
        components <- multi_component[model_local_index[predictors]]
        counts <- tabulate(components, nbins = length(multi_gamma[[block]]))
        gamma_draws <- stats::rgamma(
          length(counts), shape = multi_pi_alpha[[block]] + counts
        )
        multi_pi[[block]] <- gamma_draws / sum(gamma_draws)
        if (is.na(fixed_var[block])) {
          nonzero <- components > 1L
          scaled_sum_of_squares <- if (any(nonzero)) {
            sum(
              coefficient[predictors[nonzero]]^2 /
                multi_gamma[[block]][components[nonzero]]
            )
          } else {
            0
          }
          multi_var[block] <- 1 / stats::rgamma(
            1L,
            shape = multi_var_shape[block] + sum(nonzero) / 2,
            rate = multi_var_scale[block] + scaled_sum_of_squares / 2
          )
        }
      }
    }

    if (has_global_local) {
      for (block in which(block_model == 2L)) {
        predictors <- block_predictors[[block]]
        local_indices <- model_local_index[predictors]
        gig_chi <- pmax(
          coefficient[predictors]^2 / tau_sq[block],
          .Machine$double.xmin
        )
        local_var[local_indices] <- vapply(
          seq_along(predictors),
          function(index) {
            .draw_gig(
              n = 1L,
              lambda = local_a[block] - 0.5,
              chi = gig_chi[index],
              psi = 2 * local_aux[local_indices[index]]
            )
          },
          numeric(1)
        )
        local_aux[local_indices] <- stats::rgamma(
          length(predictors),
          shape = local_a[block] + local_b[block],
          rate = 1 + local_var[local_indices]
        )
        if (is.na(fixed_global_var[block])) {
          tau_sq[block] <- 1 / stats::rgamma(
            1L,
            shape = (length(predictors) + 1) / 2,
            rate = 1 / global_aux[block] +
              sum(coefficient[predictors]^2 / local_var[local_indices]) / 2
          )
          global_aux[block] <- 1 / stats::rgamma(
            1L,
            shape = 1,
            rate = 1 / global_scale[block]^2 + 1 / tau_sq[block]
          )
        }
      }
    }

    if (learn_residual_var) {
      sum_squared_residuals <- if (use_sufficient_statistics) {
        preliminary_tolerance <- sqrt(.Machine$double.eps) * max(1, abs(yty))
        if (residual_sse < -preliminary_tolerance) {
          fitted_crossproducts <- drop(XtX %*% coefficient)
          corrected_rhs <- Xty - fitted_crossproducts
          linear <- sum(coefficient * Xty)
          quadratic <- sum(coefficient * fitted_crossproducts)
          residual_sse <- .guard_reconstructed_sse(
            yty - 2 * linear + quadratic,
            yty, linear, quadratic, iteration
          )
        }
        max(0, residual_sse)
      } else {
        residual_sse_offset + sum(residuals^2)
      }
      residual_posterior_scale <- residual_scale +
        0.5 * sum_squared_residuals
      residual_var <- 1 / stats::rgamma(
        1L,
        shape = residual_posterior_shape,
        rate = residual_posterior_scale
      )
    }

    if (retained_index <= number_of_draws &&
        iteration == retained_iterations[retained_index]) {
      if (compute_pve) {
        if (use_sufficient_statistics) {
          fitted_crossproducts <- Xty - corrected_rhs
          total_sum_squares <- max(
            0, sum(coefficient * fitted_crossproducts)
          )
          standalone_sum_squares <- vapply(
            seq_along(block_predictors),
            function(block) {
              predictors <- block_predictors[[block]]
              block_coefficient <- coefficient[predictors]
              max(0, sum(
                block_coefficient *
                  drop(block_XtX[[block]] %*% block_coefficient)
              ))
            },
            numeric(1)
          )
          allocated_sum_squares <- vapply(
            block_predictors,
            function(predictors) {
              sum(coefficient[predictors] *
                fitted_crossproducts[predictors])
            },
            numeric(1)
          )
        } else {
          total_fitted <- y_centered - residuals
          total_sum_squares <- max(0, sum(total_fitted^2))
          standalone_sum_squares <- allocated_sum_squares <-
            numeric(number_of_blocks)
          for (block in seq_len(number_of_blocks)) {
            predictors <- block_predictors[[block]]
            block_fitted <- drop(
              x_centered[, predictors, drop = FALSE] %*%
                coefficient[predictors]
            )
            standalone_sum_squares[block] <- max(0, sum(block_fitted^2))
            allocated_sum_squares[block] <- sum(block_fitted * total_fitted)
          }
        }
        variance_df <- effective_n - as.integer(fit_intercept)
        total_signal_variance <- total_sum_squares / variance_df
        pve_denominator <- total_signal_variance + residual_var
        block_pve <- if (pve_type == "standalone") {
          standalone_sum_squares / variance_df / pve_denominator
        } else {
          allocated_sum_squares / variance_df / pve_denominator
        }
        total_pve <- total_signal_variance / pve_denominator
        cross_block_pve <-
          (total_sum_squares - sum(standalone_sum_squares)) /
            variance_df / pve_denominator
      }
      intercept_draw <- if (fit_intercept) {
        stats::rnorm(
          1L,
          mean = intercept_y_mean - sum(intercept_x_mean * coefficient),
          sd = sqrt(residual_var / effective_n)
        )
      } else {
        0
      }
      if (store_samples) {
        coefficient_samples[retained_index, ] <- coefficient
        intercept_samples[retained_index] <- intercept_draw
        residual_var_samples[retained_index] <- residual_var
        if (compute_pve) {
          block_pve_samples[retained_index, ] <- block_pve
          total_pve_samples[retained_index] <- total_pve
          cross_block_pve_samples[retained_index] <- cross_block_pve
        }
        if (has_normal) {
          normal_var_samples[retained_index, ] <- normal_var
        }
        if (has_spike_slab) {
          inclusion_samples[retained_index, ] <- inclusion
          pi_samples[retained_index, ] <- pi
          slab_var_samples[retained_index, ] <- slab_var
        }
        if (has_global_local) {
          local_var_samples[retained_index, ] <- local_var
          tau_sq_samples[retained_index, ] <- tau_sq
        }
        if (has_spike_multi_slab) {
          multi_component_samples[retained_index, ] <- multi_component
          for (block in which(block_model == 3L)) {
            multi_pi_samples[[block]][retained_index, ] <- multi_pi[[block]]
            multi_var_samples[retained_index, block] <- multi_var[block]
          }
        }
      } else {
        coefficient_sum <- coefficient_sum + coefficient
        coefficient_sum_sq <- coefficient_sum_sq + coefficient^2
        if (store_coefficient_cov) {
          for (block in seq_len(number_of_blocks)) {
            predictors <- block_predictors[[block]]
            coefficient_crossprod[[block]] <-
              coefficient_crossprod[[block]] +
                tcrossprod(coefficient[predictors])
          }
        }
        intercept_sum <- intercept_sum + intercept_draw
        intercept_sum_sq <- intercept_sum_sq + intercept_draw^2
        residual_var_sum <- residual_var_sum + residual_var
        residual_var_sum_sq <- residual_var_sum_sq + residual_var^2
        if (compute_pve) {
          block_pve_sum <- block_pve_sum + block_pve
          block_pve_sum_sq <- block_pve_sum_sq + block_pve^2
          total_pve_sum <- total_pve_sum + total_pve
          total_pve_sum_sq <- total_pve_sum_sq + total_pve^2
          cross_block_pve_sum <- cross_block_pve_sum + cross_block_pve
          cross_block_pve_sum_sq <-
            cross_block_pve_sum_sq + cross_block_pve^2
        }
        if (has_normal) {
          normal_var_sum <- normal_var_sum + normal_var
          normal_var_sum_sq <- normal_var_sum_sq + normal_var^2
        }
        if (has_spike_slab) {
          inclusion_sum <- inclusion_sum + inclusion
          pi_sum <- pi_sum + pi
          pi_sum_sq <- pi_sum_sq + pi^2
          slab_var_sum <- slab_var_sum + slab_var
          slab_var_sum_sq <- slab_var_sum_sq + slab_var^2
        }
        if (has_global_local) {
          local_var_sum <- local_var_sum + local_var
          local_var_sum_sq <- local_var_sum_sq + local_var^2
          tau_sq_sum <- tau_sq_sum + tau_sq
          tau_sq_sum_sq <- tau_sq_sum_sq + tau_sq^2
        }
        if (has_spike_multi_slab) {
          for (block in which(block_model == 3L)) {
            predictors <- block_predictors[[block]]
            components <-
              multi_component[model_local_index[predictors]]
            for (component in seq_along(multi_gamma[[block]])) {
              selected <- components == component
              multi_component_sum[[block]][selected, component] <-
                multi_component_sum[[block]][selected, component] + 1
            }
            multi_pi_sum[[block]] <-
              multi_pi_sum[[block]] + multi_pi[[block]]
            multi_pi_sum_sq[[block]] <-
              multi_pi_sum_sq[[block]] + multi_pi[[block]]^2
            multi_var_sum[block] <- multi_var_sum[block] + multi_var[block]
            multi_var_sum_sq[block] <-
              multi_var_sum_sq[block] + multi_var[block]^2
          }
        }
      }
      retained_index <- retained_index + 1L
    }

    if (progress_index <= length(progress_thresholds) &&
        iteration >= progress_thresholds[progress_index]) {
      progress_callback(
        amount = iteration - last_reported_iteration,
        iteration = iteration
      )
      last_reported_iteration <- iteration
      progress_index <- progress_index + 1L
    }
  }

  samples <- if (store_samples) {
    list(
      coefficient_samples = coefficient_samples,
      intercept_samples = intercept_samples,
      residual_var_samples = residual_var_samples
    )
  } else {
    list(
      number_of_draws = number_of_draws,
      coefficient_sum = coefficient_sum,
      coefficient_sum_sq = coefficient_sum_sq,
      intercept_sum = intercept_sum,
      intercept_sum_sq = intercept_sum_sq,
      residual_var_sum = residual_var_sum,
      residual_var_sum_sq = residual_var_sum_sq
    )
  }
  if (!store_samples && store_coefficient_cov) {
    samples$coefficient_crossprod <- coefficient_crossprod
  }
  if (compute_pve) {
    if (store_samples) {
      samples$block_pve_samples <- block_pve_samples
      samples$total_pve_samples <- total_pve_samples
      samples$cross_block_pve_samples <- cross_block_pve_samples
    } else {
      samples$block_pve_sum <- block_pve_sum
      samples$block_pve_sum_sq <- block_pve_sum_sq
      samples$total_pve_sum <- total_pve_sum
      samples$total_pve_sum_sq <- total_pve_sum_sq
      samples$cross_block_pve_sum <- cross_block_pve_sum
      samples$cross_block_pve_sum_sq <- cross_block_pve_sum_sq
    }
  }
  if (has_normal) {
    if (store_samples) {
      samples$normal_var_samples <- normal_var_samples
    } else {
      samples$normal_var_sum <- normal_var_sum
      samples$normal_var_sum_sq <- normal_var_sum_sq
    }
  }
  if (has_spike_slab) {
    if (store_samples) {
      samples$inclusion_samples <- inclusion_samples
      samples$pi_samples <- pi_samples
      samples$slab_var_samples <- slab_var_samples
    } else {
      samples$inclusion_sum <- inclusion_sum
      samples$pi_sum <- pi_sum
      samples$pi_sum_sq <- pi_sum_sq
      samples$slab_var_sum <- slab_var_sum
      samples$slab_var_sum_sq <- slab_var_sum_sq
    }
  }
  if (has_global_local) {
    if (store_samples) {
      samples$local_var_samples <- local_var_samples
      samples$tau_sq_samples <- tau_sq_samples
    } else {
      samples$local_var_sum <- local_var_sum
      samples$local_var_sum_sq <- local_var_sum_sq
      samples$tau_sq_sum <- tau_sq_sum
      samples$tau_sq_sum_sq <- tau_sq_sum_sq
    }
  }
  if (has_spike_multi_slab) {
    if (store_samples) {
      samples$multi_component_samples <- multi_component_samples
      samples$multi_pi_samples <- multi_pi_samples
      samples$multi_var_samples <- multi_var_samples
    } else {
      samples$multi_component_sum <- multi_component_sum
      samples$multi_pi_sum <- multi_pi_sum
      samples$multi_pi_sum_sq <- multi_pi_sum_sq
      samples$multi_var_sum <- multi_var_sum
      samples$multi_var_sum_sq <- multi_var_sum_sq
    }
  }
  samples
}
