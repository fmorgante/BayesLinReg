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

.validate_prior_var <- function(prior_var, number_of_predictors) {
  if (length(prior_var) == 1L) {
    .validate_variance(prior_var, "prior_var")
    return(rep(prior_var, number_of_predictors))
  }

  if (!is.numeric(prior_var) || !is.atomic(prior_var) ||
      is.object(prior_var) || !is.null(dim(prior_var)) ||
      length(prior_var) != number_of_predictors ||
      anyNA(prior_var) || any(!is.finite(prior_var)) ||
      any(prior_var <= 0)) {
    stop(
      paste0(
        "`prior_var` must be a positive, finite numeric scalar or have ",
        "one value per predictor."
      ),
      call. = FALSE
    )
  }

  prior_var
}

.as_predictor_matrix <- function(x, number_of_observations) {
  if (is.data.frame(x)) {
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

.blm_gibbs <- function(y, x, prior_var, residual_shape, residual_scale,
                       iterations, burnin, thin, seed, verbose = FALSE,
                       coefficient_prior = "normal",
                       pi_alpha = 1, pi_beta = 1) {
  retained_iterations <- .validate_mcmc(iterations, burnin, thin, seed)
  number_of_predictors <- ncol(x)
  predictor_names <- colnames(x)

  x_mean <- colMeans(x)
  y_mean <- mean(y)
  x_centered <- sweep(x, 2L, x_mean, FUN = "-")
  y_centered <- y - y_mean
  x_squared <- colSums(x_centered^2)

  number_of_draws <- length(retained_iterations)
  coefficient_samples <- matrix(
    NA_real_,
    nrow = number_of_draws,
    ncol = number_of_predictors,
    dimnames = list(NULL, predictor_names)
  )
  intercept_samples <- numeric(number_of_draws)
  residual_var_samples <- numeric(number_of_draws)
  use_spike_slab <- coefficient_prior == "spike_slab"
  if (use_spike_slab) {
    inclusion_samples <- matrix(
      NA_integer_,
      nrow = number_of_draws,
      ncol = number_of_predictors,
      dimnames = list(NULL, predictor_names)
    )
    pi_samples <- numeric(number_of_draws)
    inclusion <- rep.int(1L, number_of_predictors)
    pi <- pi_alpha / (pi_alpha + pi_beta)
  }

  coefficient <- numeric(number_of_predictors)
  residuals <- y_centered
  residual_var <- residual_scale / (residual_shape + 1)
  residual_posterior_shape <- residual_shape + (length(y) - 1) / 2
  retained_index <- 1L

  for (iteration in seq_len(iterations)) {
    if (verbose) {
      cat(sprintf("\rGibbs iteration %d/%d", iteration, iterations))
      utils::flush.console()
    }

    for (predictor in seq_len(number_of_predictors)) {
      partial_residuals <- residuals +
        x_centered[, predictor] * coefficient[predictor]
      conditional_var <- 1 / (
        x_squared[predictor] / residual_var + 1 / prior_var[predictor]
      )
      conditional_mean <- conditional_var *
        sum(x_centered[, predictor] * partial_residuals) / residual_var
      if (use_spike_slab) {
        bounded_pi <- min(
          max(pi, .Machine$double.eps),
          1 - .Machine$double.eps
        )
        log_inclusion_odds <- stats::qlogis(bounded_pi) +
          0.5 * log(conditional_var / prior_var[predictor]) +
          conditional_mean^2 / (2 * conditional_var)
        inclusion[predictor] <- stats::rbinom(
          1L,
          size = 1L,
          prob = stats::plogis(log_inclusion_odds)
        )
      }
      if (!use_spike_slab || inclusion[predictor] == 1L) {
        coefficient[predictor] <- stats::rnorm(
          1L,
          mean = conditional_mean,
          sd = sqrt(conditional_var)
        )
      } else {
        coefficient[predictor] <- 0
      }
      residuals <- partial_residuals -
        x_centered[, predictor] * coefficient[predictor]
    }

    if (use_spike_slab) {
      number_included <- sum(inclusion)
      pi <- stats::rbeta(
        1L,
        shape1 = pi_alpha + number_included,
        shape2 = pi_beta + number_of_predictors - number_included
      )
    }

    residual_posterior_scale <- residual_scale +
      0.5 * sum(residuals^2)
    residual_var <- 1 / stats::rgamma(
      1L,
      shape = residual_posterior_shape,
      rate = residual_posterior_scale
    )

    if (retained_index <= number_of_draws &&
        iteration == retained_iterations[retained_index]) {
      coefficient_samples[retained_index, ] <- coefficient
      intercept_samples[retained_index] <- stats::rnorm(
        1L,
        mean = y_mean - sum(x_mean * coefficient),
        sd = sqrt(residual_var / length(y))
      )
      residual_var_samples[retained_index] <- residual_var
      if (use_spike_slab) {
        inclusion_samples[retained_index, ] <- inclusion
        pi_samples[retained_index] <- pi
      }
      retained_index <- retained_index + 1L
    }
  }

  if (verbose) {
    cat("\n")
  }

  samples <- list(
    coefficient_samples = coefficient_samples,
    intercept_samples = intercept_samples,
    residual_var_samples = residual_var_samples
  )
  if (use_spike_slab) {
    samples$inclusion_samples <- inclusion_samples
    samples$pi_samples <- pi_samples
  }
  samples
}

.blm_gibbs_rcpp <- function(y, x, prior_var, residual_shape, residual_scale,
                            iterations, burnin, thin, seed, verbose = FALSE,
                            coefficient_prior = "normal",
                            pi_alpha = 1, pi_beta = 1) {
  .validate_mcmc(iterations, burnin, thin, seed)
  samples <- blm_gibbs_rcpp_cpp(
    y = y,
    X = x,
    prior_var = prior_var,
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    verbose = verbose,
    use_spike_slab = coefficient_prior == "spike_slab",
    pi_alpha = pi_alpha,
    pi_beta = pi_beta
  )
  colnames(samples$coefficient_samples) <- colnames(x)
  if (coefficient_prior == "spike_slab") {
    colnames(samples$inclusion_samples) <- colnames(x)
  } else {
    samples$inclusion_samples <- NULL
    samples$pi_samples <- NULL
  }
  samples
}
