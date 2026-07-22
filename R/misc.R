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

.variance_from_sums <- function(sum, sum_sq, number_of_draws) {
  variance <- (sum_sq - sum^2 / number_of_draws) / (number_of_draws - 1)
  max(0, variance)
}

.covariance_from_sums <- function(sum, crossprod, number_of_draws) {
  (crossprod - tcrossprod(sum) / number_of_draws) / (number_of_draws - 1)
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

.normalize_eta <- function(ETA, number_of_observations, residual_var) {
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
          c("Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab")) {
      stop(
        sprintf(
          paste0(
            "ETA block `%s` has an invalid `model`; use `Normal`, ",
            "`SpikeSlab`, `GlobalLocal`, or `SpikeMultiSlab`."
          ),
          block_name
        ),
        call. = FALSE
      )
    }
    model <- specification$model
    allowed <- switch(
      model,
      Normal = c("X", "model", "standardize", "var_shape", "var_scale"),
      SpikeSlab = c(
        "X", "model", "standardize", "var_shape", "var_scale", "pi"
      ),
      GlobalLocal = c(
        "X", "model", "standardize", "local_shape", "global_scale"
      ),
      SpikeMultiSlab = c(
        "X", "model", "standardize", "gamma", "alpha", "var_shape",
        "var_scale"
      )
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
    x_centered <- sweep(X, 2L, colMeans(X), FUN = "-")
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
    pi_alpha <- pi_beta <- 1
    spike_var_shape <- 2
    spike_var_scale <- 1
    local_shape <- c(a = 1, b = 0.5)
    global_scale <- 1
    multi_gamma <- c(0, 0.01, 0.1, 1)
    multi_pi_alpha <- rep(1, length(multi_gamma))
    multi_var_shape <- 2
    multi_var_scale <- 1
    if (model == "Normal") {
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
      .validate_variance(normal_scale, "var_scale")
    }
    if (model == "SpikeSlab") {
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
      .validate_variance(spike_var_scale, "var_scale")
    }
    if (model == "GlobalLocal") {
      local_shape <- if (is.null(specification$local_shape)) {
        c(a = 1, b = 0.5)
      } else {
        .validate_local_shape(specification$local_shape)
      }
      global_scale <- if (is.null(specification$global_scale)) 1 else
        specification$global_scale
      .validate_variance(global_scale, "global_scale")
    }
    if (model == "SpikeMultiSlab") {
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
      .validate_variance(multi_var_scale, "var_scale")
    }

    list(
      name = block_name,
      model = model,
      model_code = match(
        model, c("Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab")
      ) - 1L,
      X = X,
      x = sweep(X, 2L, predictor_scale, FUN = "/"),
      predictor_names = colnames(X),
      predictor_scale = predictor_scale,
      standardize = standardize,
      normal_shape = normal_shape,
      normal_scale = normal_scale,
      pi_alpha = pi_alpha,
      pi_beta = pi_beta,
      spike_var_shape = spike_var_shape,
      spike_var_scale = spike_var_scale,
      local_shape = local_shape,
      global_scale = global_scale,
      multi_gamma = multi_gamma,
      multi_pi_alpha = multi_pi_alpha,
      multi_var_shape = multi_var_shape,
      multi_var_scale = multi_var_scale
    )
  })
  names(blocks) <- block_names
  blocks
}

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

.validate_mcmc <- function(iterations, burnin, thin, seed) {
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

  if (!is.null(seed)) {
    if (!is_whole_number(seed)) {
      stop("`seed` must be NULL or a finite integer.", call. = FALSE)
    }
    set.seed(seed)
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

.blm_gibbs <- function(y, x, residual_shape, residual_scale,
                       iterations, burnin, thin, seed,
                       progress_callback = NULL,
                       block_id = NULL, block_model = 0L,
                       normal_shape = 2, normal_scale = 1,
                       pi_alpha = 1, pi_beta = 1,
                       spike_var_shape = 2, spike_var_scale = 1,
                       global_scale = 1, residual_var = NULL,
                       local_a = 1, local_b = 0.5,
                       multi_gamma = list(c(0, 0.01, 0.1, 1)),
                       multi_pi_alpha = list(rep(1, 4)),
                       multi_var_shape = 2, multi_var_scale = 1,
                       store_samples = TRUE,
                       store_coefficient_cov = TRUE,
                       effective_n = NULL, fit_intercept = TRUE,
                       intercept_x_mean = NULL, intercept_y_mean = NULL,
                       XtX = NULL, XtX_center = NULL, Xty = NULL, yty = NULL) {
  retained_iterations <- .validate_mcmc(iterations, burnin, thin, seed)
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

  if (use_sufficient_statistics) {
    x_squared <- diag(XtX)
    corrected_rhs <- as.numeric(Xty)
    residual_sse <- yty
  } else {
    x_mean <- colMeans(x)
    y_mean <- mean(y)
    if (is.null(effective_n)) effective_n <- length(y)
    if (is.null(intercept_x_mean)) intercept_x_mean <- x_mean
    if (is.null(intercept_y_mean)) intercept_y_mean <- y_mean
    x_centered <- sweep(x, 2L, x_mean, FUN = "-")
    y_centered <- y - y_mean
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
      coefficient_crossprod <- matrix(
        0, number_of_predictors, number_of_predictors
      )
    }
    intercept_sum <- intercept_sum_sq <- 0
    residual_var_sum <- residual_var_sum_sq <- 0
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
    normal_var <- normal_scale / (normal_shape + 1)
  }
  if (has_spike_slab) {
    if (store_samples) {
      inclusion_samples <- matrix(
        NA_integer_,
        nrow = number_of_draws,
        ncol = number_of_predictors,
        dimnames = list(NULL, predictor_names)
      )
      pi_samples <- matrix(NA_real_, number_of_draws, number_of_blocks)
      slab_var_samples <- matrix(NA_real_, number_of_draws, number_of_blocks)
    } else {
      inclusion_sum <- numeric(number_of_predictors)
      pi_sum <- pi_sum_sq <- numeric(number_of_blocks)
      slab_var_sum <- slab_var_sum_sq <- numeric(number_of_blocks)
    }
    inclusion <- rep.int(1L, number_of_predictors)
    pi <- pi_alpha / (pi_alpha + pi_beta)
    slab_var <- spike_var_scale / (spike_var_shape + 1)
  }
  if (has_global_local) {
    if (store_samples) {
      local_var_samples <- matrix(
        NA_real_,
        nrow = number_of_draws,
        ncol = number_of_predictors,
        dimnames = list(NULL, predictor_names)
      )
      tau_sq_samples <- matrix(NA_real_, number_of_draws, number_of_blocks)
    } else {
      local_var_sum <- local_var_sum_sq <- numeric(number_of_predictors)
      tau_sq_sum <- tau_sq_sum_sq <- numeric(number_of_blocks)
    }
    local_var <- rep(1, number_of_predictors)
    local_aux <- rep(1, number_of_predictors)
    tau_sq <- global_scale^2
    global_aux <- rep(1, number_of_blocks)
  }
  if (has_spike_multi_slab) {
    multi_component <- rep.int(1L, number_of_predictors)
    multi_pi <- lapply(seq_len(number_of_blocks), function(block) {
      multi_pi_alpha[[block]] / sum(multi_pi_alpha[[block]])
    })
    multi_var <- multi_var_scale / (multi_var_shape + 1)
    if (store_samples) {
      multi_component_samples <- matrix(
        NA_integer_, number_of_draws, number_of_predictors,
        dimnames = list(NULL, predictor_names)
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

  for (iteration in seq_len(iterations)) {
    for (predictor in seq_len(number_of_predictors)) {
      block <- block_id[predictor]
      model <- block_model[block]
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
        multi_component[predictor] <- sample.int(
          length(gamma), 1L, prob = probabilities
        )
        component <- multi_component[predictor]
        coefficient[predictor] <- if (component == 1L) {
          0
        } else {
          stats::rnorm(
            1L, conditional_means[component],
            sqrt(conditional_vars[component])
          )
        }
      } else {
        prior_precision <- if (model == 2L) {
          1 / tau_sq[block] / local_var[predictor]
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
          inclusion[predictor] <- stats::rbinom(
            1L,
            size = 1L,
            prob = stats::plogis(log_inclusion_odds)
          )
        }
        if (model != 1L || inclusion[predictor] == 1L) {
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

    if (use_sufficient_statistics && iteration %% 100L == 0L) {
      fitted_crossproducts <- drop(XtX %*% coefficient)
      corrected_rhs <- Xty - fitted_crossproducts
      if (learn_residual_var) {
        residual_sse <- yty - 2 * sum(coefficient * Xty) +
          sum(coefficient * fitted_crossproducts)
      }
    }

    if (has_normal) {
      for (block in which(block_model == 0L)) {
        predictors <- which(block_id == block)
        normal_var[block] <- 1 / stats::rgamma(
          1L,
          shape = normal_shape[block] + length(predictors) / 2,
          rate = normal_scale[block] + sum(coefficient[predictors]^2) / 2
        )
      }
    }

    if (has_spike_slab) {
      for (block in which(block_model == 1L)) {
        predictors <- which(block_id == block)
        number_included <- sum(inclusion[predictors])
        pi[block] <- stats::rbeta(
          1L,
          shape1 = pi_alpha[block] + number_included,
          shape2 = pi_beta[block] + length(predictors) - number_included
        )
        included_predictors <- predictors[inclusion[predictors] == 1L]
        slab_var[block] <- 1 / stats::rgamma(
          1L,
          shape = spike_var_shape[block] + length(included_predictors) / 2,
          rate = spike_var_scale[block] +
            sum(coefficient[included_predictors]^2) / 2
        )
      }
    }

    if (has_spike_multi_slab) {
      for (block in which(block_model == 3L)) {
        predictors <- which(block_id == block)
        components <- multi_component[predictors]
        counts <- tabulate(components, nbins = length(multi_gamma[[block]]))
        gamma_draws <- stats::rgamma(
          length(counts), shape = multi_pi_alpha[[block]] + counts
        )
        multi_pi[[block]] <- gamma_draws / sum(gamma_draws)
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

    if (has_global_local) {
      for (block in which(block_model == 2L)) {
        predictors <- which(block_id == block)
        gig_chi <- pmax(
          coefficient[predictors]^2 / tau_sq[block],
          .Machine$double.xmin
        )
        local_var[predictors] <- vapply(
          seq_along(predictors),
          function(index) {
            predictor <- predictors[index]
            .draw_gig(
              n = 1L,
              lambda = local_a[block] - 0.5,
              chi = gig_chi[index],
              psi = 2 * local_aux[predictor]
            )
          },
          numeric(1)
        )
        local_aux[predictors] <- stats::rgamma(
          length(predictors),
          shape = local_a[block] + local_b[block],
          rate = 1 + local_var[predictors]
        )
        tau_sq[block] <- 1 / stats::rgamma(
          1L,
          shape = (length(predictors) + 1) / 2,
          rate = 1 / global_aux[block] +
            sum(coefficient[predictors]^2 / local_var[predictors]) / 2
        )
        global_aux[block] <- 1 / stats::rgamma(
          1L,
          shape = 1,
          rate = 1 / global_scale[block]^2 + 1 / tau_sq[block]
        )
      }
    }

    if (learn_residual_var) {
      sum_squared_residuals <- if (use_sufficient_statistics) {
        max(0, residual_sse)
      } else {
        sum(residuals^2)
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
          coefficient_crossprod <- coefficient_crossprod +
            tcrossprod(coefficient)
        }
        intercept_sum <- intercept_sum + intercept_draw
        intercept_sum_sq <- intercept_sum_sq + intercept_draw^2
        residual_var_sum <- residual_var_sum + residual_var
        residual_var_sum_sq <- residual_var_sum_sq + residual_var^2
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
            predictors <- which(block_id == block)
            for (component in seq_along(multi_gamma[[block]])) {
              selected <- multi_component[predictors] == component
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

.blm_gibbs_rcpp <- function(y, x, residual_shape, residual_scale,
                            iterations, burnin, thin, seed,
                            progress_callback = NULL,
                            block_id = NULL, block_model = 0L,
                            normal_shape = 2, normal_scale = 1,
                            pi_alpha = 1, pi_beta = 1,
                            spike_var_shape = 2, spike_var_scale = 1,
                            global_scale = 1, residual_var = NULL,
                            local_a = 1, local_b = 0.5,
                            multi_gamma = list(c(0, 0.01, 0.1, 1)),
                            multi_pi_alpha = list(rep(1, 4)),
                            multi_var_shape = 2, multi_var_scale = 1,
                            store_samples = TRUE,
                            store_coefficient_cov = TRUE,
                            effective_n = NULL, fit_intercept = TRUE,
                            intercept_x_mean = NULL,
                            intercept_y_mean = NULL,
                            XtX = NULL, XtX_center = NULL, Xty = NULL,
                            yty = NULL) {
  .validate_mcmc(iterations, burnin, thin, seed)
  use_sufficient_statistics <- !is.null(XtX)
  if (is.null(block_id)) {
    block_id <- rep.int(
      1L, if (use_sufficient_statistics) ncol(XtX) else ncol(x)
    )
  }
  if (is.null(progress_callback)) {
    progress_callback <- function(amount, iteration) invisible(NULL)
  }
  if (is.null(effective_n)) effective_n <- length(y)
  if (is.null(intercept_x_mean)) intercept_x_mean <- colMeans(x)
  if (is.null(intercept_y_mean)) intercept_y_mean <- mean(y)
  common_arguments <- list(
    y = y,
    X = x,
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    progress_callback = progress_callback,
    block_id = block_id,
    block_model = block_model,
    normal_shape = normal_shape,
    normal_scale = normal_scale,
    pi_alpha = pi_alpha,
    pi_beta = pi_beta,
    spike_var_shape = spike_var_shape,
    spike_var_scale = spike_var_scale,
    global_scale = global_scale,
    local_a = local_a,
    local_b = local_b,
    multi_gamma_list = multi_gamma,
    multi_pi_alpha_list = multi_pi_alpha,
    multi_var_shape = multi_var_shape,
    multi_var_scale = multi_var_scale,
    learn_residual_var = is.null(residual_var),
    fixed_residual_var = if (is.null(residual_var)) 1 else residual_var,
    store_samples = store_samples,
    store_coefficient_cov = store_coefficient_cov,
    effective_n = effective_n,
    fit_intercept = fit_intercept,
    intercept_x_mean = intercept_x_mean,
    intercept_y_mean = intercept_y_mean,
    use_sufficient_statistics = use_sufficient_statistics,
    summary_Xty = if (use_sufficient_statistics) Xty else numeric(),
    summary_yty = if (use_sufficient_statistics && !is.null(yty)) yty else 0
  )
  samples <- if (inherits(XtX, "dgCMatrix")) {
    do.call(
      blm_gibbs_sparse_rcpp_cpp,
      c(
        common_arguments[c(
          "residual_shape", "residual_scale", "iterations", "burnin", "thin",
          "progress_callback", "block_id", "block_model", "normal_shape",
          "normal_scale", "pi_alpha", "pi_beta", "spike_var_shape",
          "spike_var_scale", "global_scale", "local_a", "local_b",
          "multi_gamma_list", "multi_pi_alpha_list", "multi_var_shape",
          "multi_var_scale", "learn_residual_var", "fixed_residual_var",
          "store_samples", "store_coefficient_cov", "effective_n",
          "fit_intercept", "intercept_x_mean", "intercept_y_mean"
        )],
        list(
          summary_XtX = XtX,
          summary_center = XtX_center,
          summary_Xty = Xty,
          summary_yty = if (is.null(yty)) 0 else yty
        )
      )
    )
  } else {
    common_arguments$summary_XtX <- if (use_sufficient_statistics) {
      XtX
    } else {
      matrix(0, 0L, 0L)
    }
    do.call(blm_gibbs_rcpp_cpp, common_arguments)
  }
  if (store_samples) {
    predictor_names <- if (use_sufficient_statistics) {
      colnames(XtX)
    } else {
      colnames(x)
    }
    colnames(samples$coefficient_samples) <- predictor_names
    if (!any(block_model == 0L)) {
      samples$normal_var_samples <- NULL
    }
    if (any(block_model == 1L)) {
      colnames(samples$inclusion_samples) <- predictor_names
    } else {
      samples$inclusion_samples <- NULL
      samples$pi_samples <- NULL
      samples$slab_var_samples <- NULL
    }
    if (any(block_model == 2L)) {
      colnames(samples$local_var_samples) <- predictor_names
    } else {
      samples$local_var_samples <- NULL
      samples$tau_sq_samples <- NULL
    }
    if (any(block_model == 3L)) {
      colnames(samples$multi_component_samples) <- predictor_names
    } else {
      samples$multi_component_samples <- NULL
      samples$multi_pi_samples <- NULL
      samples$multi_var_samples <- NULL
    }
  }
  samples
}

.chain_progress_callback <- function(progressor, chain, nchains) {
  if (is.null(progressor)) {
    return(NULL)
  }
  force(progressor)
  force(chain)
  force(nchains)
  function(amount, iteration) {
    progressor(
      amount = amount,
      message = sprintf(
        "Chain %d/%d: iteration %d",
        chain,
        nchains,
        iteration
      )
    )
  }
}

.run_blm_chains <- function(sampler_arguments, version, nchains, seed,
                            block_model, progressor = NULL) {
  if (nchains == 1L) {
    sampler <- if (version == "Rcpp") .blm_gibbs_rcpp else .blm_gibbs
    return(do.call(
      sampler,
      c(
        sampler_arguments,
        list(
          seed = seed,
          progress_callback = .chain_progress_callback(
            progressor,
            chain = 1L,
            nchains = 1L
          )
        )
      )
    ))
  }

  chain_seeds <- if (is.null(seed)) {
    rep(list(NULL), nchains)
  } else {
    as.list((abs(as.double(seed)) + seq_len(nchains) - 1) %% 2147483647)
  }
  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)
  future::plan(future::multisession, workers = nchains)

  chain_futures <- lapply(seq_len(nchains), function(chain) {
    chain_seed <- chain_seeds[[chain]]
    chain_progress <- .chain_progress_callback(progressor, chain, nchains)
    future::future({
      namespace <- asNamespace("BayesLinReg")
      chain_sampler <- if (version == "Rcpp") {
        get(".blm_gibbs_rcpp", envir = namespace)
      } else {
        get(".blm_gibbs", envir = namespace)
      }
      do.call(
        chain_sampler,
        c(
          sampler_arguments,
          list(
            seed = chain_seed,
            progress_callback = chain_progress
          )
        )
      )
    }, seed = TRUE)
  })
  chain_samples <- lapply(chain_futures, future::value)
  .combine_blm_chains(
    chain_samples,
    block_model = block_model,
    store_samples = sampler_arguments$store_samples,
    store_coefficient_cov = sampler_arguments$store_coefficient_cov
  )
}

.combine_blm_chains <- function(chain_samples, block_model = 0L,
                                store_samples = TRUE,
                                store_coefficient_cov = TRUE) {
  if (!store_samples) {
    summary_names <- c(
      "number_of_draws", "coefficient_sum", "coefficient_sum_sq",
      "intercept_sum", "intercept_sum_sq", "residual_var_sum",
      "residual_var_sum_sq"
    )
    if (store_coefficient_cov) {
      summary_names <- c(summary_names, "coefficient_crossprod")
    }
    if (any(block_model == 0L)) {
      summary_names <- c(summary_names, "normal_var_sum", "normal_var_sum_sq")
    }
    if (any(block_model == 1L)) {
      summary_names <- c(
        summary_names, "inclusion_sum", "pi_sum", "pi_sum_sq",
        "slab_var_sum", "slab_var_sum_sq"
      )
    }
    if (any(block_model == 2L)) {
      summary_names <- c(
        summary_names, "local_var_sum", "local_var_sum_sq",
        "tau_sq_sum", "tau_sq_sum_sq"
      )
    }
    if (any(block_model == 3L)) {
      summary_names <- c(
        summary_names, "multi_var_sum", "multi_var_sum_sq"
      )
    }
    combined <- stats::setNames(lapply(summary_names, function(name) {
      Reduce(`+`, lapply(chain_samples, `[[`, name))
    }), summary_names)
    if (any(block_model == 3L)) {
      list_names <- c(
        "multi_component_sum", "multi_pi_sum", "multi_pi_sum_sq"
      )
      for (name in list_names) {
        combined[[name]] <- lapply(seq_along(block_model), function(block) {
          values <- lapply(chain_samples, function(samples) {
            samples[[name]][[block]]
          })
          if (all(vapply(values, is.null, logical(1)))) NULL else Reduce(`+`, values)
        })
      }
    }
    return(combined)
  }
  number_of_draws <- vapply(
    chain_samples,
    function(samples) nrow(samples$coefficient_samples),
    integer(1)
  )
  combined <- list(
    coefficient_samples = do.call(
      rbind,
      lapply(chain_samples, `[[`, "coefficient_samples")
    ),
    intercept_samples = unlist(
      lapply(chain_samples, `[[`, "intercept_samples"),
      use.names = FALSE
    ),
    residual_var_samples = unlist(
      lapply(chain_samples, `[[`, "residual_var_samples"),
      use.names = FALSE
    ),
    chain_id = rep.int(seq_along(chain_samples), number_of_draws)
  )
  if (any(block_model == 0L)) {
    combined$normal_var_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "normal_var_samples")
    )
  }
  if (any(block_model == 1L)) {
    combined$inclusion_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "inclusion_samples")
    )
    combined$pi_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "pi_samples")
    )
    combined$slab_var_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "slab_var_samples")
    )
  }
  if (any(block_model == 2L)) {
    combined$local_var_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "local_var_samples")
    )
    combined$tau_sq_samples <- do.call(
      rbind,
      lapply(chain_samples, `[[`, "tau_sq_samples")
    )
  }
  if (any(block_model == 3L)) {
    combined$multi_component_samples <- do.call(
      rbind, lapply(chain_samples, `[[`, "multi_component_samples")
    )
    combined$multi_pi_samples <- lapply(
      seq_along(block_model),
      function(block) {
        if (block_model[block] != 3L) return(NULL)
        do.call(rbind, lapply(chain_samples, function(samples) {
          samples$multi_pi_samples[[block]]
        }))
      }
    )
    combined$multi_var_samples <- do.call(
      rbind, lapply(chain_samples, `[[`, "multi_var_samples")
    )
  }
  combined
}

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

  sample_matrix <- if (has_intercept) {
    cbind(
      intercept = fit$intercept_samples,
      residual_var = fit$residual_var_samples
    )
  } else {
    cbind(residual_var = fit$residual_var_samples)
  }

  for (block_name in names(fit$ETA)) {
    block <- fit$ETA[[block_name]]
    coefficient_samples <- as.matrix(block$coefficient_samples)
    if (is.null(block$coefficient_samples) ||
        nrow(coefficient_samples) != number_of_draws) {
      stop("`fit` contains incompatible ETA samples.", call. = FALSE)
    }
    if (identical(block$model, "Normal")) {
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

    if (identical(block$model, "GlobalLocal")) {
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
      if (is.null(block$var_samples) ||
          length(block$var_samples) != number_of_draws) {
        stop("`fit` contains incompatible multi-slab variance samples.",
             call. = FALSE)
      }
      sample_matrix <- cbind(sample_matrix, block$var_samples)
      colnames(sample_matrix)[ncol(sample_matrix)] <- paste0(
        "var_", block_name
      )
    }
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
