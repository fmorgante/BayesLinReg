# Internal validation helpers.

.validate_variance <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value <= 0) {
    stop(
      sprintf("`%s` must be a positive, finite numeric scalar.", name),
      call. = FALSE
    )
  }
}

.variance_from_m2 <- function(m2, number_of_draws) {
  pmax(0, m2 / (number_of_draws - 1))
}

.covariance_from_m2 <- function(m2, number_of_draws) {
  m2 / (number_of_draws - 1)
}

.validate_local_shape <- function(local_shape) {
  if (!is.numeric(local_shape) || !is.atomic(local_shape) ||
      is.object(local_shape) || !is.null(dim(local_shape)) ||
      length(local_shape) != 2L || anyNA(local_shape) ||
      any(!is.finite(local_shape)) || any(local_shape <= 0)) {
    stop(
      "`local_shape` must contain two positive, finite numeric values.",
      call. = FALSE
    )
  }
  stats::setNames(as.numeric(local_shape), c("a", "b"))
}

.validate_pi <- function(pi) {
  if (!is.numeric(pi) || !is.atomic(pi) || is.object(pi) ||
      !is.null(dim(pi)) || length(pi) != 2L || anyNA(pi) ||
      any(!is.finite(pi)) || any(pi <= 0)) {
    stop(
      "`pi` must contain two positive, finite numeric values.",
      call. = FALSE
    )
  }
  stats::setNames(as.numeric(pi), c("a", "b"))
}

.validate_multi_alpha <- function(alpha, number_of_components) {
  if (!is.numeric(alpha) || !is.atomic(alpha) || is.object(alpha) ||
      !is.null(dim(alpha)) || length(alpha) != number_of_components ||
      anyNA(alpha) || any(!is.finite(alpha)) || any(alpha <= 0)) {
    stop(
      sprintf(
        "`alpha` must contain %d positive, finite Dirichlet concentrations.",
        number_of_components
      ),
      call. = FALSE
    )
  }
  as.numeric(alpha)
}

.validate_gamma <- function(gamma) {
  if (!is.numeric(gamma) || !is.atomic(gamma) || is.object(gamma) ||
      !is.null(dim(gamma)) || length(gamma) < 2L || anyNA(gamma) ||
      any(!is.finite(gamma)) || gamma[1L] != 0 ||
      any(gamma[-1L] <= 0) || is.unsorted(gamma, strictly = TRUE)) {
    stop(
      paste0(
        "`gamma` must start with zero and continue with strictly increasing, ",
        "positive, finite variance multipliers."
      ),
      call. = FALSE
    )
  }
  as.numeric(gamma)
}

.validate_expected_pve <- function(value, block_name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value <= 0 || value >= 1) {
    stop(
      sprintf(
        "ETA block `%s`: `expected_pve` must be between zero and one.",
        block_name
      ),
      call. = FALSE
    )
  }
  as.numeric(value)
}

.normalize_eta <- function(ETA, number_of_observations, residual_var,
                           calibration_n = number_of_observations) {
  if (!is.list(ETA) || length(ETA) < 1L) {
    stop("`ETA` must be a non-empty list.", call. = FALSE)
  }
  if (all(c("X", "model") %in% names(ETA))) {
    ETA <- list(ETA1 = ETA)
  }
  if (!all(vapply(ETA, is.list, logical(1)))) {
    stop(
      paste0(
        "`ETA` must be a predictor specification or a list of predictor ",
        "specifications."
      ),
      call. = FALSE
    )
  }

  block_names <- names(ETA)
  if (is.null(block_names)) {
    block_names <- paste0("ETA", seq_along(ETA))
  } else {
    missing_name <- is.na(block_names) | block_names == ""
    block_names[missing_name] <- paste0("ETA", which(missing_name))
    block_names <- make.unique(block_names)
  }

  blocks <- lapply(seq_along(ETA), function(block_index) {
    specification <- ETA[[block_index]]
    block_name <- block_names[block_index]
    if (!all(c("X", "model") %in% names(specification))) {
      stop(
        sprintf("ETA block `%s` must contain `X` and `model`.", block_name),
        call. = FALSE
      )
    }
    if (!is.character(specification$model) ||
        length(specification$model) != 1L || is.na(specification$model) ||
        !specification$model %in%
          c("Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab",
            "Fixed")) {
      stop(
        sprintf(
          paste0(
            "ETA block `%s` has an invalid `model`; use `Normal`, ",
            "`SpikeSlab`, `GlobalLocal`, `SpikeMultiSlab`, or `Fixed`."
          ),
          block_name
        ),
        call. = FALSE
      )
    }
    model <- specification$model
    allowed <- switch(
      model,
      Normal = c(
        "X", "model", "standardize", "var", "var_shape", "var_scale",
        "expected_pve"
      ),
      SpikeSlab = c(
        "X", "model", "standardize", "var", "var_shape", "var_scale", "pi",
        "expected_pve"
      ),
      GlobalLocal = c(
        "X", "model", "standardize", "local_shape", "global_var",
        "global_scale", "expected_nonzero", "reference_residual_var",
        "expected_pve"
      ),
      SpikeMultiSlab = c(
        "X", "model", "standardize", "gamma", "alpha", "var", "var_shape",
        "var_scale", "expected_pve"
      ),
      Fixed = c("X", "model", "standardize")
    )
    unknown <- setdiff(names(specification), allowed)
    if (length(unknown) > 0L) {
      stop(
        sprintf(
          "ETA block `%s` contains unsupported field(s): %s.",
          block_name,
          paste(unknown, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    X <- .as_predictor_matrix(specification$X, number_of_observations)
    standardize <- specification$standardize
    if (is.null(standardize)) {
      standardize <- TRUE
    }
    if (!is.logical(standardize) || length(standardize) != 1L ||
        is.na(standardize)) {
      stop(
        sprintf("ETA block `%s`: `standardize` must be TRUE or FALSE.",
                block_name),
        call. = FALSE
      )
    }
    predictor_mean <- colMeans(X)
    x_centered <- sweep(X, 2L, predictor_mean, FUN = "-")
    constant_predictors <- vapply(
      seq_len(ncol(X)),
      function(index) all(X[, index] == X[1L, index]),
      logical(1)
    )
    if (any(constant_predictors)) {
      stop(
        sprintf(
          "ETA block `%s` contains constant predictor(s): %s.",
          block_name,
          paste(colnames(X)[constant_predictors], collapse = ", ")
        ),
        call. = FALSE
      )
    }
    predictor_scale <- if (standardize) {
      sqrt(colSums(x_centered^2) / (nrow(X) - 1))
    } else {
      rep(1, ncol(X))
    }
    if (any(!is.finite(predictor_scale)) || any(predictor_scale <= 0)) {
      stop(
        sprintf(
          paste0(
            "ETA block `%s` contains a predictor with nonpositive variance."
          ),
          block_name
        ),
        call. = FALSE
      )
    }

    normal_shape <- 2
    normal_scale <- 1
    normal_scale_calibrated <- FALSE
    pi_alpha <- pi_beta <- 1
    spike_var_shape <- 2
    spike_var_scale <- 1
    spike_var_scale_calibrated <- FALSE
    local_shape <- c(a = 1, b = 0.5)
    global_scale <- 1
    expected_nonzero <- NULL
    reference_residual_var <- NULL
    global_scale_calibrated <- FALSE
    global_scale_calibration <- "none"
    fixed_global_var <- NA_real_
    multi_gamma <- c(0, 0.01, 0.1, 1)
    multi_pi_alpha <- rep(1, length(multi_gamma))
    multi_var_shape <- 2
    multi_var_scale <- 1
    multi_var_scale_calibrated <- FALSE
    fixed_var <- NA_real_
    expected_pve <- if (is.null(specification$expected_pve)) {
      NULL
    } else {
      .validate_expected_pve(specification$expected_pve, block_name)
    }
    if (model == "Normal") {
      if (!is.null(specification[["var"]])) {
        conflicting <- intersect(
          c("var_shape", "var_scale", "expected_pve"), names(specification)
        )
        if (length(conflicting)) {
          stop(sprintf(
            "ETA block `%s`: `var` cannot be combined with %s.",
            block_name, paste(sprintf("`%s`", conflicting), collapse = ", ")
          ), call. = FALSE)
        }
        .validate_variance(specification[["var"]], "var")
        fixed_var <- as.numeric(specification[["var"]])
      }
      if (!is.null(expected_pve) && !is.null(specification$var_scale)) {
        stop(
          sprintf(
            "ETA block `%s`: supply either `var_scale` or `expected_pve`, not both.",
            block_name
          ),
          call. = FALSE
        )
      }
      normal_shape <- if (is.null(specification$var_shape)) {
        2
      } else {
        specification$var_shape
      }
      normal_scale <- if (is.null(specification$var_scale)) {
        1
      } else {
        specification$var_scale
      }
      .validate_variance(normal_shape, "var_shape")
      if (!is.null(expected_pve) && normal_shape <= 1) {
        stop(
          sprintf(
            "ETA block `%s`: `var_shape` must exceed one with `expected_pve`.",
            block_name
          ),
          call. = FALSE
        )
      }
      if (is.null(expected_pve)) .validate_variance(normal_scale, "var_scale")
    }
    if (model == "SpikeSlab") {
      if (!is.null(specification[["var"]])) {
        conflicting <- intersect(
          c("var_shape", "var_scale", "expected_pve"), names(specification)
        )
        if (length(conflicting)) {
          stop(sprintf(
            "ETA block `%s`: `var` cannot be combined with %s.",
            block_name, paste(sprintf("`%s`", conflicting), collapse = ", ")
          ), call. = FALSE)
        }
        .validate_variance(specification[["var"]], "var")
        fixed_var <- as.numeric(specification[["var"]])
      }
      if (!is.null(expected_pve) && !is.null(specification$var_scale)) {
        stop(
          sprintf(
            "ETA block `%s`: supply either `var_scale` or `expected_pve`, not both.",
            block_name
          ),
          call. = FALSE
        )
      }
      pi <- if (is.null(specification$pi)) {
        c(a = 1, b = 1)
      } else {
        .validate_pi(specification$pi)
      }
      pi_alpha <- unname(pi["a"])
      pi_beta <- unname(pi["b"])
      spike_var_shape <- if (is.null(specification$var_shape)) {
        2
      } else {
        specification$var_shape
      }
      spike_var_scale <- if (is.null(specification$var_scale)) {
        1
      } else {
        specification$var_scale
      }
      .validate_variance(spike_var_shape, "var_shape")
      if (!is.null(expected_pve) && spike_var_shape <= 1) {
        stop(
          sprintf(
            "ETA block `%s`: `var_shape` must exceed one with `expected_pve`.",
            block_name
          ),
          call. = FALSE
        )
      }
      if (is.null(expected_pve)) {
        .validate_variance(spike_var_scale, "var_scale")
      }
    }
    if (model == "GlobalLocal") {
      local_shape <- if (is.null(specification$local_shape)) {
        c(a = 1, b = 0.5)
      } else {
        .validate_local_shape(specification$local_shape)
      }
      has_global_scale <- !is.null(specification$global_scale)
      has_expected_nonzero <- !is.null(specification$expected_nonzero)
      has_reference_residual_var <-
        !is.null(specification$reference_residual_var)
      if (!is.null(specification[["global_var"]])) {
        conflicting <- intersect(
          c(
            "global_scale", "expected_nonzero", "reference_residual_var",
            "expected_pve"
          ),
          names(specification)
        )
        if (length(conflicting)) {
          stop(sprintf(
            "ETA block `%s`: `global_var` cannot be combined with %s.",
            block_name, paste(sprintf("`%s`", conflicting), collapse = ", ")
          ), call. = FALSE)
        }
        .validate_variance(specification[["global_var"]], "global_var")
        fixed_global_var <- as.numeric(specification[["global_var"]])
      }
      if (has_global_scale &&
          (has_expected_nonzero || has_reference_residual_var ||
           !is.null(expected_pve))) {
        stop(
          sprintf(
            paste0(
              "ETA block `%s`: supply either `global_scale` or ",
              "an expected-sparsity calibration, not both."
            ),
            block_name
          ),
          call. = FALSE
        )
      }
      if (has_reference_residual_var && !is.null(expected_pve)) {
        stop(
          sprintf(
            paste0(
              "ETA block `%s`: supply either `reference_residual_var` or ",
              "`expected_pve`, not both."
            ),
            block_name
          ),
          call. = FALSE
        )
      }
      has_scale_reference <- has_reference_residual_var ||
        !is.null(expected_pve)
      if (xor(has_expected_nonzero, has_scale_reference)) {
        stop(
          sprintf(
            paste0(
              "ETA block `%s`: `expected_nonzero` and ",
              "either `reference_residual_var` or `expected_pve` must be ",
              "supplied together."
            ),
            block_name
          ),
          call. = FALSE
        )
      }
      if (has_expected_nonzero) {
        expected_nonzero <- specification$expected_nonzero
        if (!is.numeric(expected_nonzero) ||
            length(expected_nonzero) != 1L || is.na(expected_nonzero) ||
            !is.finite(expected_nonzero) || expected_nonzero <= 0 ||
            expected_nonzero >= ncol(X)) {
          stop(
            sprintf(
              paste0(
                "ETA block `%s`: `expected_nonzero` must be a positive ",
                "finite number smaller than the block's %d predictors."
              ),
              block_name, ncol(X)
            ),
            call. = FALSE
          )
        }
        if (has_reference_residual_var) {
          reference_residual_var <- specification$reference_residual_var
          .validate_variance(
            reference_residual_var, "reference_residual_var"
          )
          global_scale <- expected_nonzero /
            (ncol(X) - expected_nonzero) *
            sqrt(reference_residual_var / calibration_n)
          global_scale_calibrated <- TRUE
          global_scale_calibration <- "reference_residual_var"
        }
      } else {
        global_scale <- if (has_global_scale) specification$global_scale else 1
        .validate_variance(global_scale, "global_scale")
      }
    }
    if (model == "SpikeMultiSlab") {
      if (!is.null(specification[["var"]])) {
        conflicting <- intersect(
          c("var_shape", "var_scale", "expected_pve"), names(specification)
        )
        if (length(conflicting)) {
          stop(sprintf(
            "ETA block `%s`: `var` cannot be combined with %s.",
            block_name, paste(sprintf("`%s`", conflicting), collapse = ", ")
          ), call. = FALSE)
        }
        .validate_variance(specification[["var"]], "var")
        fixed_var <- as.numeric(specification[["var"]])
      }
      if (!is.null(expected_pve) && !is.null(specification$var_scale)) {
        stop(
          sprintf(
            "ETA block `%s`: supply either `var_scale` or `expected_pve`, not both.",
            block_name
          ),
          call. = FALSE
        )
      }
      multi_gamma <- if (is.null(specification$gamma)) {
        c(0, 0.01, 0.1, 1)
      } else {
        .validate_gamma(specification$gamma)
      }
      multi_pi_alpha <- if (is.null(specification$alpha)) {
        rep(1, length(multi_gamma))
      } else {
        .validate_multi_alpha(specification$alpha, length(multi_gamma))
      }
      multi_var_shape <- if (is.null(specification$var_shape)) {
        2
      } else {
        specification$var_shape
      }
      multi_var_scale <- if (is.null(specification$var_scale)) {
        1
      } else {
        specification$var_scale
      }
      .validate_variance(multi_var_shape, "var_shape")
      if (!is.null(expected_pve) && multi_var_shape <= 1) {
        stop(
          sprintf(
            "ETA block `%s`: `var_shape` must exceed one with `expected_pve`.",
            block_name
          ),
          call. = FALSE
        )
      }
      if (is.null(expected_pve)) {
        .validate_variance(multi_var_scale, "var_scale")
      }
    }

    list(
      name = block_name,
      model = model,
      model_code = match(
        model,
        c("Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab", "Fixed")
      ) - 1L,
      x = sweep(x_centered, 2L, predictor_scale, FUN = "/"),
      predictor_names = colnames(X),
      predictor_scale = predictor_scale,
      predictor_mean = predictor_mean / predictor_scale,
      standardize = standardize,
      normal_shape = normal_shape,
      normal_scale = normal_scale,
      normal_scale_calibrated = normal_scale_calibrated,
      pi_alpha = pi_alpha,
      pi_beta = pi_beta,
      spike_var_shape = spike_var_shape,
      spike_var_scale = spike_var_scale,
      spike_var_scale_calibrated = spike_var_scale_calibrated,
      local_shape = local_shape,
      global_scale = global_scale,
      expected_nonzero = expected_nonzero,
      reference_residual_var = reference_residual_var,
      global_scale_calibrated = global_scale_calibrated,
      global_scale_calibration = global_scale_calibration,
      fixed_global_var = fixed_global_var,
      multi_gamma = multi_gamma,
      multi_pi_alpha = multi_pi_alpha,
      multi_var_shape = multi_var_shape,
      multi_var_scale = multi_var_scale,
      multi_var_scale_calibrated = multi_var_scale_calibrated,
      fixed_var = fixed_var,
      expected_pve = expected_pve,
      reference_response_var = NULL
    )
  })
  names(blocks) <- block_names
  blocks
}

.calibrate_eta_priors <- function(blocks, predictor_variance_sums,
                                  reference_response_var, n) {
  has_expected_pve <- vapply(
    blocks, function(block) !is.null(block$expected_pve), logical(1)
  )
  if (!any(has_expected_pve)) return(blocks)

  .validate_variance(reference_response_var, "reference_response_var")
  expected_pve_total <- sum(vapply(
    blocks[has_expected_pve], `[[`, numeric(1), "expected_pve"
  ))
  if (expected_pve_total >= 1) {
    stop(
      "The sum of supplied `expected_pve` values must be smaller than one.",
      call. = FALSE
    )
  }

  global_local_from_pve <- has_expected_pve & vapply(
    blocks, function(block) block$model == "GlobalLocal", logical(1)
  )
  penalized <- vapply(
    blocks, function(block) block$model != "Fixed", logical(1)
  )
  if (any(global_local_from_pve) && !all(has_expected_pve[penalized])) {
    stop(
      paste0(
        "When a GlobalLocal block uses `expected_pve`, every penalized ETA ",
        "block must supply `expected_pve` so the common reference residual ",
        "variance can be determined."
      ),
      call. = FALSE
    )
  }
  reference_residual_var <- if (any(global_local_from_pve)) {
    (1 - expected_pve_total) * reference_response_var
  } else {
    NULL
  }

  for (block_index in which(has_expected_pve)) {
    block <- blocks[[block_index]]
    predictor_variance_sum <- unname(predictor_variance_sums[block_index])
    if (!is.finite(predictor_variance_sum) ||
        predictor_variance_sum <= 0) {
      stop(
        sprintf(
          "ETA block `%s` has an invalid predictor-variance sum.",
          names(blocks)[block_index]
        ),
        call. = FALSE
      )
    }
    expected_signal_var <-
      block$expected_pve * reference_response_var
    if (block$model == "Normal") {
      prior_mean_var <- expected_signal_var / predictor_variance_sum
      block$normal_scale <- (block$normal_shape - 1) * prior_mean_var
      block$normal_scale_calibrated <- TRUE
    } else if (block$model == "SpikeSlab") {
      expected_inclusion <-
        block$pi_alpha / (block$pi_alpha + block$pi_beta)
      prior_mean_var <- expected_signal_var /
        (predictor_variance_sum * expected_inclusion)
      block$spike_var_scale <-
        (block$spike_var_shape - 1) * prior_mean_var
      block$spike_var_scale_calibrated <- TRUE
    } else if (block$model == "SpikeMultiSlab") {
      expected_component_probability <-
        block$multi_pi_alpha / sum(block$multi_pi_alpha)
      expected_gamma <- sum(
        expected_component_probability * block$multi_gamma
      )
      prior_mean_var <- expected_signal_var /
        (predictor_variance_sum * expected_gamma)
      block$multi_var_scale <-
        (block$multi_var_shape - 1) * prior_mean_var
      block$multi_var_scale_calibrated <- TRUE
    } else {
      block$reference_residual_var <- reference_residual_var
      block$global_scale <- block$expected_nonzero /
        (length(block$predictor_names) - block$expected_nonzero) *
        sqrt(reference_residual_var / n)
      block$global_scale_calibrated <- TRUE
      block$global_scale_calibration <- "expected_pve"
    }
    block$reference_response_var <- reference_response_var
    blocks[[block_index]] <- block
  }
  blocks
}

.fixed_predictor_indices <- function(blocks, block_indices) {
  fixed <- vapply(blocks, function(block) block$model == "Fixed", logical(1))
  unlist(block_indices[fixed], use.names = FALSE)
}

.validate_fixed_design <- function(x, fixed_indices, predictor_names) {
  if (length(fixed_indices) == 0L) return(invisible(NULL))
  fixed_x <- x[, fixed_indices, drop = FALSE]
  fixed_rank <- qr(fixed_x)$rank
  if (fixed_rank < ncol(fixed_x)) {
    stop(
      sprintf(
        paste0(
          "The Fixed predictors are jointly rank-deficient (rank %d for %d ",
          "predictors): %s."
        ),
        fixed_rank, ncol(fixed_x), paste(predictor_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.validate_fixed_gram <- function(gram, fixed_indices, predictor_names) {
  if (length(fixed_indices) == 0L) return(invisible(NULL))
  fixed_gram <- as.matrix(gram[fixed_indices, fixed_indices, drop = FALSE])
  fixed_gram <- (fixed_gram + t(fixed_gram)) / 2
  eigenvalues <- eigen(fixed_gram, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(eigenvalues)))
  fixed_rank <- sum(eigenvalues > tolerance)
  if (fixed_rank < length(fixed_indices)) {
    stop(
      sprintf(
        paste0(
          "The Fixed predictors are jointly rank-deficient (rank %d for %d ",
          "predictors): %s."
        ),
        fixed_rank, length(fixed_indices), paste(predictor_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(NULL)
}

.prepare_residual_prior <- function(residual_var, residual_shape,
                                    residual_scale, blocks,
                                    reference_response_var) {
  if (!is.null(residual_var)) {
    if (!is.null(residual_shape) || !is.null(residual_scale)) {
      stop(
        "Supply either `residual_var` or the inverse-gamma prior, not both.",
        call. = FALSE
      )
    }
    .validate_variance(residual_var, "residual_var")
    return(list(
      residual_var = residual_var,
      residual_shape = NULL,
      residual_scale = NULL,
      residual_scale_calibrated = FALSE,
      expected_pve_total = NULL,
      reference_residual_var = NULL
    ))
  }

  if (is.null(residual_shape)) {
    stop(
      "`residual_shape` is required when `residual_var` is NULL.",
      call. = FALSE
    )
  }
  .validate_variance(residual_shape, "residual_shape")

  residual_scale_calibrated <- is.null(residual_scale)
  expected_pve_total <- NULL
  reference_residual_var <- NULL
  if (residual_scale_calibrated) {
    has_expected_pve <- vapply(
      blocks, function(block) !is.null(block$expected_pve), logical(1)
    )
    if (!all(has_expected_pve)) {
      stop(
        paste0(
          "`residual_scale` is required unless every ETA block supplies ",
          "`expected_pve`."
        ),
        call. = FALSE
      )
    }
    if (residual_shape <= 1) {
      stop(
        paste0(
          "`residual_shape` must exceed one when `residual_scale` is ",
          "calibrated from `expected_pve`."
        ),
        call. = FALSE
      )
    }
    .validate_variance(reference_response_var, "reference_response_var")
    expected_pve_total <- sum(vapply(
      blocks, `[[`, numeric(1), "expected_pve"
    ))
    if (expected_pve_total >= 1) {
      stop(
        "The sum of supplied `expected_pve` values must be smaller than one.",
        call. = FALSE
      )
    }
    reference_residual_var <-
      (1 - expected_pve_total) * reference_response_var
    residual_scale <-
      (residual_shape - 1) * reference_residual_var
  }
  .validate_variance(residual_scale, "residual_scale")

  list(
    residual_var = NULL,
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    residual_scale_calibrated = residual_scale_calibrated,
    expected_pve_total = expected_pve_total,
    reference_residual_var = reference_residual_var
  )
}

.as_predictor_matrix <- function(x, number_of_observations) {
  if (is.numeric(x) && is.atomic(x) && !is.object(x) && is.null(dim(x))) {
    x <- matrix(as.numeric(x), ncol = 1L, dimnames = list(NULL, "x"))
  } else if (is.data.frame(x)) {
    if (ncol(x) < 1L || !all(vapply(x, is.numeric, logical(1)))) {
      stop("`X` must contain at least one numeric predictor.", call. = FALSE)
    }
    x <- as.matrix(x)
  } else if (!is.matrix(x) || !is.numeric(x)) {
    stop("`X` must be a numeric matrix or data frame.", call. = FALSE)
  }

  if (nrow(x) != number_of_observations) {
    stop("`y` and `X` must have the same number of observations.",
         call. = FALSE)
  }
  if (ncol(x) < 1L) {
    stop("`X` must contain at least one predictor.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`X` must contain only finite, non-missing values.", call. = FALSE)
  }

  predictor_names <- colnames(x)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("x", seq_len(ncol(x)))
  } else {
    missing_name <- is.na(predictor_names) | predictor_names == ""
    predictor_names[missing_name] <- paste0("x", which(missing_name))
    predictor_names <- make.unique(predictor_names)
  }

  storage.mode(x) <- "double"
  colnames(x) <- predictor_names
  x
}
