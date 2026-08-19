# Posterior draw conversion and convergence calculations.

.fit_sample_matrix <- function(fit) {
  if (is.list(fit) && identical(fit$store_samples, FALSE)) {
    stop(
      paste0(
        "Convergence diagnostics require individual posterior draws; ",
        "refit with `store_samples = TRUE`."
      ),
      call. = FALSE
    )
  }
  required <- c("ETA", "residual_var_samples")
  if (!is.list(fit) || !all(required %in% names(fit))) {
    stop(
      "`fit` must be a sampled fit returned by `blm()`.",
      call. = FALSE
    )
  }

  if (!is.list(fit$ETA) || length(fit$ETA) < 1L ||
      is.null(fit$ETA[[1L]]$coefficient_samples)) {
    stop(
      "`fit` must contain posterior samples from `blm()`.",
      call. = FALSE
    )
  }
  number_of_draws <- nrow(as.matrix(fit$ETA[[1L]]$coefficient_samples))
  has_intercept <- !is.null(fit$intercept_samples)
  if ((has_intercept && length(fit$intercept_samples) != number_of_draws) ||
      length(fit$residual_var_samples) != number_of_draws) {
    stop("`fit` contains sample components with incompatible lengths.",
         call. = FALSE)
  }

  sample_matrix <- matrix(numeric(), nrow = number_of_draws, ncol = 0L)
  if (has_intercept) {
    sample_matrix <- cbind(sample_matrix, intercept = fit$intercept_samples)
  }
  if (!is.null(fit$residual_shape)) {
    sample_matrix <- cbind(
      sample_matrix, residual_var = fit$residual_var_samples
    )
  }

  for (block_name in names(fit$ETA)) {
    block <- fit$ETA[[block_name]]
    coefficient_samples <- as.matrix(block$coefficient_samples)
    if (is.null(block$coefficient_samples) ||
        nrow(coefficient_samples) != number_of_draws) {
      stop("`fit` contains incompatible ETA samples.", call. = FALSE)
    }
    if (!is.null(block$pve_samples)) {
      if (length(block$pve_samples) != number_of_draws) {
        stop("`fit` contains incompatible block-PVE samples.",
             call. = FALSE)
      }
      sample_matrix <- cbind(sample_matrix, block$pve_samples)
      colnames(sample_matrix)[ncol(sample_matrix)] <- paste0(
        "pve_", block_name
      )
    }
    if (identical(block$model, "Normal") && is.null(block$var)) {
      if (is.null(block$normal_var_samples) ||
          length(block$normal_var_samples) != number_of_draws) {
        stop("`fit` contains incompatible normal-variance samples.",
             call. = FALSE)
      }
      sample_matrix <- cbind(sample_matrix, block$normal_var_samples)
      colnames(sample_matrix)[ncol(sample_matrix)] <- paste0(
        "normal_var_", block_name
      )
    }

    if (identical(block$model, "SpikeSlab")) {
      if (is.null(block$pi_samples) ||
          length(block$pi_samples) != number_of_draws) {
        stop("`fit` contains incompatible pi samples.", call. = FALSE)
      }
      sample_matrix <- cbind(sample_matrix, block$pi_samples)
      colnames(sample_matrix)[ncol(sample_matrix)] <- paste0("pi_", block_name)

      if (is.null(block$var)) {
        if (is.null(block$slab_var_samples) ||
            length(block$slab_var_samples) != number_of_draws) {
          stop("`fit` contains incompatible slab-variance samples.",
               call. = FALSE)
        }
        sample_matrix <- cbind(sample_matrix, block$slab_var_samples)
        colnames(sample_matrix)[ncol(sample_matrix)] <- paste0(
          "slab_var_", block_name
        )
      }
    }

    if (identical(block$model, "GlobalLocal") &&
        is.null(block$global_var)) {
      if (length(block$tau_sq_samples) != number_of_draws) {
        stop("`fit` contains incompatible global-variance samples.",
             call. = FALSE)
      }
      sample_matrix <- cbind(sample_matrix, block$tau_sq_samples)
      colnames(sample_matrix)[ncol(sample_matrix)] <- paste0(
        "tau_sq_", block_name
      )
    }
    if (identical(block$model, "SpikeMultiSlab")) {
      pi_samples <- as.matrix(block$pi_samples)
      if (nrow(pi_samples) != number_of_draws ||
          ncol(pi_samples) != length(block$gamma)) {
        stop("`fit` contains incompatible multi-slab pi samples.",
             call. = FALSE)
      }
      sample_matrix <- cbind(sample_matrix, pi_samples)
      component_names <- names(block$gamma)
      new_columns <- seq.int(
        ncol(sample_matrix) - ncol(pi_samples) + 1L, ncol(sample_matrix)
      )
      colnames(sample_matrix)[new_columns] <- paste0(
        "pi_", block_name, "_", component_names
      )
      if (is.null(block$var)) {
        if (is.null(block$sampled_slab_var_samples) ||
            length(block$sampled_slab_var_samples) != number_of_draws) {
          stop("`fit` contains incompatible multi-slab variance samples.",
               call. = FALSE)
        }
        sample_matrix <- cbind(sample_matrix, block$sampled_slab_var_samples)
        colnames(sample_matrix)[ncol(sample_matrix)] <- paste0(
          "var_", block_name
        )
      }
    }
  }
  if (!is.null(fit$total_pve_samples)) {
    if (length(fit$total_pve_samples) != number_of_draws ||
        length(fit$cross_block_pve_samples) != number_of_draws) {
      stop("`fit` contains incompatible total-PVE samples.", call. = FALSE)
    }
    sample_matrix <- cbind(
      sample_matrix,
      total_pve = fit$total_pve_samples,
      cross_block_pve = fit$cross_block_pve_samples
    )
  }
  if (ncol(sample_matrix) == 0L) {
    stop(
      paste0(
        "The fit has no sampled non-coefficient parameters or posterior ",
        "PVE quantities to assess."
      ),
      call. = FALSE
    )
  }
  sample_matrix
}

.as_blm_mcmc_list <- function(fit) {
  sample_matrix <- .fit_sample_matrix(fit)
  number_of_draws <- nrow(sample_matrix)
  chain_id <- if (is.null(fit$chain_id)) {
    rep.int(1L, number_of_draws)
  } else {
    fit$chain_id
  }
  if (length(chain_id) != number_of_draws || anyNA(chain_id) ||
      any(chain_id < 1) || any(chain_id != floor(chain_id))) {
    stop("`fit$chain_id` is invalid.", call. = FALSE)
  }

  split_indices <- split(seq_len(number_of_draws), chain_id)
  chain_lengths <- vapply(split_indices, length, integer(1))
  if (length(unique(chain_lengths)) != 1L) {
    stop("All chains must contain the same number of retained draws.",
         call. = FALSE)
  }
  if (chain_lengths[1] < 20L) {
    stop("At least 20 retained draws per chain are required.",
         call. = FALSE)
  }

  coda::mcmc.list(lapply(
    split_indices,
    function(indices) coda::mcmc(sample_matrix[indices, , drop = FALSE])
  ))
}

.classical_rhat <- function(chains) {
  parameter_names <- coda::varnames(chains)
  number_of_chains <- coda::nchain(chains)
  if (number_of_chains < 2L) {
    return(stats::setNames(rep(NA_real_, length(parameter_names)),
                           parameter_names))
  }

  chain_matrices <- lapply(chains, as.matrix)
  draws_per_chain <- nrow(chain_matrices[[1]])
  rhat <- vapply(seq_along(parameter_names), function(parameter) {
    chain_means <- vapply(
      chain_matrices,
      function(chain) mean(chain[, parameter]),
      numeric(1)
    )
    within_variance <- mean(vapply(
      chain_matrices,
      function(chain) stats::var(chain[, parameter]),
      numeric(1)
    ))
    if (!is.finite(within_variance) || within_variance <= 0) {
      return(NA_real_)
    }
    between_variance <- draws_per_chain * stats::var(chain_means)
    pooled_variance <- (draws_per_chain - 1) / draws_per_chain *
      within_variance + between_variance / draws_per_chain
    sqrt(pooled_variance / within_variance)
  }, numeric(1))
  stats::setNames(rhat, parameter_names)
}
