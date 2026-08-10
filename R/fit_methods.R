#' Extract coefficients from a BayesLinReg fit
#'
#' Returns posterior means of the intercept and regression coefficients.
#'
#' @param object A fitted model returned by [blm()], [blm_ss()],
#'   [blm_ss_eigen()], or [blm_gwas()].
#' @param ... Additional arguments. Currently unused.
#'
#' @return A named numeric vector. For multi-block fits, coefficient names use
#'   the form `block::predictor`. The intercept, when fitted, is named
#'   `(Intercept)`.
#' @export
coef.blm_fit <- function(object, ...) {
  if (length(list(...)) > 0L) {
    stop("Additional arguments are not supported.", call. = FALSE)
  }
  coefficients <- lapply(object$ETA, `[[`, "coefficient_mean")
  if (length(coefficients) > 1L) {
    coefficients <- Map(function(values, block_name) {
      names(values) <- paste0(block_name, "::", names(values))
      values
    }, coefficients, names(coefficients))
  }
  coefficients <- do.call(c, unname(coefficients))
  if (!is.null(object$intercept_mean)) {
    coefficients <- c(`(Intercept)` = object$intercept_mean, coefficients)
  }
  coefficients
}

#' Predict from a BayesLinReg fit
#'
#' Computes posterior-mean predictions for new predictor values.
#'
#' @param object A fitted model returned by [blm()], [blm_ss()],
#'   [blm_ss_eigen()], or [blm_gwas()].
#' @param newdata For a single-block fit, a numeric vector, matrix, or data
#'   frame containing that block's predictors. A vector represents one
#'   observation unless the fitted block has one predictor, in which case it
#'   represents multiple observations. For a multi-block fit, a named list of
#'   matrices or data frames, one for each `ETA` block. Named columns are
#'   reordered to match the fitted predictors.
#' @param interval Interval to return. `"none"` returns only posterior-mean
#'   predictions. `"credible"` returns an equal-tailed posterior interval for
#'   the conditional mean. `"prediction"` returns an equal-tailed posterior
#'   predictive interval that also includes residual variation. Intervals
#'   require a fit created with `store_samples = TRUE`.
#' @param level Coverage probability for credible or prediction intervals.
#' @param chunk_size Positive integer giving the maximum number of new
#'   observations whose retained fitted-value draws are materialized at once.
#'   Chunking bounds peak memory for interval prediction and does not affect
#'   point predictions.
#' @param ... Additional arguments. Currently unused.
#'
#' @return With `interval = "none"`, a numeric vector of posterior-mean
#'   predictions. Otherwise, a matrix with columns `fit`, `lwr`, and `upr`.
#'
#' @details Prediction intervals use deterministic quantiles of the
#'   equal-weight normal mixture induced by the retained coefficient and
#'   residual-variance draws; they do not add a second layer of Monte Carlo
#'   simulation.
#' @export
predict.blm_fit <- function(object, newdata,
                            interval = c("none", "credible", "prediction"),
                            level = 0.95, chunk_size = 1000L, ...) {
  if (missing(newdata)) {
    stop("`newdata` is required because fitted predictor data are not stored.",
         call. = FALSE)
  }
  if (length(list(...)) > 0L) {
    stop("Additional arguments are not supported.", call. = FALSE)
  }
  interval <- match.arg(interval)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must be a number strictly between zero and one.",
         call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
      is.na(chunk_size) || !is.finite(chunk_size) ||
      chunk_size != floor(chunk_size) || chunk_size < 1 ||
      chunk_size > .Machine$integer.max) {
    stop("`chunk_size` must be a positive integer.", call. = FALSE)
  }
  chunk_size <- as.integer(chunk_size)

  block_names <- names(object$ETA)
  if (length(object$ETA) == 1L &&
      (is.numeric(newdata) || is.matrix(newdata) || is.data.frame(newdata))) {
    newdata <- stats::setNames(list(newdata), block_names)
  } else {
    if (!is.list(newdata) || is.data.frame(newdata)) {
      stop("Multi-block predictions require a named list of predictor inputs.",
           call. = FALSE)
    }
    if (is.null(names(newdata)) || anyDuplicated(names(newdata)) ||
        !setequal(names(newdata), block_names)) {
      stop("`newdata` names must match the fitted ETA block names.",
           call. = FALSE)
    }
    newdata <- newdata[block_names]
  }

  matrices <- Map(function(values, block) {
    .prediction_matrix(values, names(block$coefficient_mean))
  }, newdata, object$ETA)
  row_counts <- vapply(matrices, nrow, integer(1))
  if (length(unique(row_counts)) != 1L) {
    stop("All `newdata` blocks must have the same number of rows.",
         call. = FALSE)
  }

  prediction <- if (is.null(object$intercept_mean)) {
    numeric(row_counts[1L])
  } else {
    rep(object$intercept_mean, row_counts[1L])
  }
  for (block_index in seq_along(matrices)) {
    prediction <- prediction + drop(
      matrices[[block_index]] %*%
        object$ETA[[block_index]]$coefficient_mean
    )
  }
  prediction <- unname(prediction)
  if (interval == "none") return(prediction)
  if (!isTRUE(object$store_samples)) {
    stop(
      paste0(
        "Prediction intervals require retained posterior draws; refit with ",
        "`store_samples = TRUE`."
      ),
      call. = FALSE
    )
  }

  number_of_draws <- nrow(object$ETA[[1L]]$coefficient_samples)
  tail_probability <- (1 - level) / 2
  probabilities <- c(tail_probability, 1 - tail_probability)
  number_of_observations <- row_counts[1L]
  bounds <- matrix(NA_real_, nrow = number_of_observations, ncol = 2L)
  residual_sd <- if (interval == "prediction") {
    sqrt(object$residual_var_samples)
  } else {
    NULL
  }
  for (start in seq.int(1L, number_of_observations, by = chunk_size)) {
    indices <- seq.int(
      start, min(number_of_observations, start + chunk_size - 1L)
    )
    prediction_draws <- matrix(
      0, nrow = length(indices), ncol = number_of_draws
    )
    if (!is.null(object$intercept_samples)) {
      prediction_draws <- sweep(
        prediction_draws, 2L, object$intercept_samples, FUN = "+"
      )
    }
    for (block_index in seq_along(matrices)) {
      prediction_draws <- prediction_draws + tcrossprod(
        matrices[[block_index]][indices, , drop = FALSE],
        object$ETA[[block_index]]$coefficient_samples
      )
    }
    bounds[indices, ] <- if (interval == "credible") {
      t(apply(prediction_draws, 1L, stats::quantile,
              probs = probabilities, names = FALSE))
    } else {
      t(apply(prediction_draws, 1L, function(draws) {
        vapply(
          probabilities,
          .normal_mixture_quantile,
          numeric(1),
          means = draws,
          standard_deviations = residual_sd
        )
      }))
    }
  }
  result <- cbind(fit = prediction, lwr = bounds[, 1L], upr = bounds[, 2L])
  rownames(result) <- rownames(matrices[[1L]])
  result
}

.normal_mixture_quantile <- function(probability, means,
                                     standard_deviations) {
  component_quantiles <- means +
    standard_deviations * stats::qnorm(probability)
  lower <- min(component_quantiles)
  upper <- max(component_quantiles)
  if (lower == upper) return(lower)
  stats::uniroot(
    function(value) {
      mean(stats::pnorm(
        value, mean = means, sd = standard_deviations
      )) - probability
    },
    interval = c(lower, upper),
    extendInt = "yes"
  )$root
}

#' Summarize a BayesLinReg fit
#'
#' Organizes scalable block-level, residual-variance, hyperparameter, and
#' variance-explained summaries. Coefficient tables are adaptive and opt-in
#' for high-dimensional fits. Empirical posterior quantiles are computed only
#' for reported coefficients when the fit retained individual draws.
#'
#' @param object A fitted model returned by [blm()], [blm_ss()],
#'   [blm_ss_eigen()], or [blm_gwas()].
#' @param probs Numeric vector of posterior probabilities used for empirical
#'   quantiles. Quantiles are unavailable when `store_samples = FALSE`.
#' @param coefficients Coefficient-table policy. `"auto"` reports all
#'   coefficients only when their number does not exceed `max_coefficients`;
#'   otherwise it omits the table. `"none"` always omits it, `"top"` reports
#'   at most `max_coefficients` ranked predictors, and `"all"` reports every
#'   predictor.
#' @param max_coefficients Nonnegative integer controlling the `"auto"` and
#'   `"top"` modes.
#' @param rank_by Ranking used with `coefficients = "top"`.
#'   `"absolute_mean"` ranks by absolute posterior coefficient mean.
#'   `"inclusion_probability"` is available when at least one block has a
#'   discrete spike component; coefficients without inclusion probabilities
#'   sort after those with probabilities. Absolute-mean ranking uses
#'   coefficients transformed back to the original predictor scale. It is
#'   therefore not invariant to predictor rescaling and is most meaningful
#'   when predictors use comparable units. Internal standardization does not
#'   alter this behavior. For spike models with differently scaled predictors,
#'   inclusion-probability ranking is generally preferable.
#' @param ... Additional arguments. Currently unused.
#'
#' @return An object of class `summary.blm_fit` containing block-level and
#'   residual-variance summaries, an intercept summary when fitted, the
#'   requested coefficient table, and, when available, hyperparameter and
#'   variance-explained tables. The print method returns its input invisibly.
#' @export
summary.blm_fit <- function(
    object, probs = c(0.025, 0.5, 0.975),
    coefficients = c("auto", "none", "top", "all"),
    max_coefficients = 20L,
    rank_by = c("absolute_mean", "inclusion_probability"), ...) {
  if (length(list(...)) > 0L) {
    stop("Additional arguments are not supported.", call. = FALSE)
  }
  probs <- .validate_summary_probs(probs)
  coefficients <- match.arg(coefficients)
  rank_by <- match.arg(rank_by)
  if (!is.numeric(max_coefficients) || length(max_coefficients) != 1L ||
      is.na(max_coefficients) || !is.finite(max_coefficients) ||
      max_coefficients != floor(max_coefficients) || max_coefficients < 0) {
    stop("`max_coefficients` must be a nonnegative integer.", call. = FALSE)
  }
  max_coefficients <- as.integer(max_coefficients)
  stored <- isTRUE(object$store_samples)

  intercept <- if (is.null(object$intercept_mean)) {
    NULL
  } else {
    intercept_row <- .posterior_summary_row(
      object$intercept_mean, object$intercept_var,
      if (stored) object$intercept_samples else NULL, probs
    )
    matrix(
      intercept_row, nrow = 1L,
      dimnames = list("(Intercept)", names(intercept_row))
    )
  }

  metadata <- .coefficient_summary_metadata(object)
  selection <- .select_summary_coefficients(
    metadata, coefficients, max_coefficients, rank_by
  )
  coefficient_table <- .summarize_selected_coefficients(
    object, metadata, selection$indices, probs
  )

  residual_row <- .posterior_summary_row(
    object$residual_var_mean, object$residual_var_var,
    if (stored) object$residual_var_samples else NULL, probs
  )
  residual_variance <- matrix(
    residual_row,
    nrow = 1L,
    dimnames = list("Residual variance", names(residual_row))
  )
  blocks <- .summarize_blocks(object)
  hyperparameters <- .summarize_hyperparameters(object, probs)
  variance_explained <- .summarize_pve(object, probs)

  result <- list(
    blocks = blocks,
    intercept = intercept,
    coefficients = coefficient_table,
    residual_variance = residual_variance,
    hyperparameters = hyperparameters,
    variance_explained = variance_explained,
    coefficient_summary = selection$mode,
    coefficient_rank_by = if (selection$mode == "top") rank_by else NULL,
    total_coefficients = nrow(metadata),
    reported_coefficients = length(selection$indices),
    store_samples = stored,
    number_of_draws = if (stored) {
      length(object$residual_var_samples)
    } else {
      NA_integer_
    },
    nchains = if (is.null(object$nchains)) 1L else object$nchains,
    probs = probs
  )
  class(result) <- "summary.blm_fit"
  result
}

.coefficient_summary_metadata <- function(object) {
  multiple_blocks <- length(object$ETA) > 1L
  rows <- lapply(seq_along(object$ETA), function(block_index) {
    block_name <- names(object$ETA)[block_index]
    block <- object$ETA[[block_index]]
    predictor_names <- names(block$coefficient_mean)
    inclusion <- if (is.null(block$inclusion_probability)) {
      rep(NA_real_, length(predictor_names))
    } else {
      unname(block$inclusion_probability[predictor_names])
    }
    data.frame(
      block_index = block_index,
      predictor = predictor_names,
      row_name = if (multiple_blocks) {
        paste0(block_name, "::", predictor_names)
      } else {
        predictor_names
      },
      mean = unname(block$coefficient_mean[predictor_names]),
      variance = unname(block$coefficient_var[predictor_names]),
      inclusion = inclusion,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.select_summary_coefficients <- function(metadata, mode, maximum, rank_by) {
  resolved_mode <- if (mode == "auto") {
    if (nrow(metadata) <= maximum) "all" else "none"
  } else {
    mode
  }
  indices <- if (resolved_mode == "all") {
    seq_len(nrow(metadata))
  } else if (resolved_mode == "none" || maximum == 0L) {
    integer()
  } else {
    score <- if (rank_by == "absolute_mean") {
      abs(metadata$mean)
    } else {
      if (all(is.na(metadata$inclusion))) {
        stop(
          paste0(
            "`rank_by = \"inclusion_probability\"` requires a SpikeSlab ",
            "or SpikeMultiSlab block."
          ),
          call. = FALSE
        )
      }
      metadata$inclusion
    }
    score[is.na(score)] <- -Inf
    ordered <- order(score, decreasing = TRUE)
    ordered[seq_len(min(maximum, length(ordered)))]
  }
  list(mode = resolved_mode, indices = indices)
}

.summarize_selected_coefficients <- function(object, metadata, indices,
                                             probs) {
  if (length(indices) == 0L) return(NULL)
  stored <- isTRUE(object$store_samples)
  rows <- lapply(indices, function(index) {
    item <- metadata[index, ]
    block <- object$ETA[[item$block_index]]
    .posterior_summary_row(
      item$mean, item$variance,
      if (stored) block$coefficient_samples[, item$predictor] else NULL,
      probs
    )
  })
  names(rows) <- metadata$row_name[indices]
  result <- do.call(rbind, rows)
  inclusion <- metadata$inclusion[indices]
  if (any(!is.na(inclusion))) {
    result <- cbind(result, `Inclusion Probability` = inclusion)
  }
  result
}

.summarize_blocks <- function(object) {
  result <- data.frame(
    Model = vapply(object$ETA, `[[`, character(1), "model"),
    Standardized = vapply(object$ETA, `[[`, logical(1), "standardize"),
    Predictors = vapply(
      object$ETA, function(block) length(block$coefficient_mean), integer(1)
    ),
    check.names = FALSE
  )
  has_inclusion <- vapply(
    object$ETA,
    function(block) !is.null(block$inclusion_probability),
    logical(1)
  )
  if (any(has_inclusion)) {
    result$`Expected Included` <- vapply(object$ETA, function(block) {
      if (is.null(block$inclusion_probability)) {
        NA_real_
      } else {
        sum(block$inclusion_probability)
      }
    }, numeric(1))
    result$`Mean Inclusion Probability` <- vapply(object$ETA, function(block) {
      if (is.null(block$inclusion_probability)) {
        NA_real_
      } else {
        mean(block$inclusion_probability)
      }
    }, numeric(1))
  }
  if (!is.null(object$total_pve_mean)) {
    result$`PVE Mean` <- vapply(object$ETA, `[[`, numeric(1), "pve_mean")
    result$`PVE SD` <- sqrt(vapply(
      object$ETA, `[[`, numeric(1), "pve_var"
    ))
  }
  result
}

.validate_summary_probs <- function(probs) {
  if (!is.numeric(probs) || !is.atomic(probs) || is.object(probs) ||
      !is.null(dim(probs)) || length(probs) < 1L || anyNA(probs) ||
      any(!is.finite(probs)) || any(probs <= 0 | probs >= 1) ||
      is.unsorted(probs, strictly = TRUE)) {
    stop(
      "`probs` must contain strictly increasing numbers between zero and one.",
      call. = FALSE
    )
  }
  as.numeric(probs)
}

.posterior_quantile_names <- function(probs) {
  names(stats::quantile(0, probs = probs))
}

.posterior_summary_row <- function(mean, variance, samples, probs) {
  result <- c(Mean = unname(mean), SD = sqrt(max(0, unname(variance))))
  if (!is.null(samples)) {
    quantiles <- stats::quantile(samples, probs = probs, names = FALSE)
    names(quantiles) <- .posterior_quantile_names(probs)
    result <- c(result, quantiles)
  }
  result
}

.summarize_hyperparameters <- function(object, probs) {
  rows <- list()
  stored <- isTRUE(object$store_samples)
  add_row <- function(name, mean, variance, samples = NULL) {
    rows[[name]] <<- .posterior_summary_row(
      mean, variance, if (stored) samples else NULL, probs
    )
  }
  for (block_name in names(object$ETA)) {
    block <- object$ETA[[block_name]]
    if (block$model == "Normal") {
      add_row(paste0(block_name, "::coefficient_var"),
              block$normal_var_mean, block$normal_var_var,
              block$normal_var_samples)
    } else if (block$model == "SpikeSlab") {
      add_row(paste0(block_name, "::pi"), block$pi_mean, block$pi_var,
              block$pi_samples)
      add_row(paste0(block_name, "::slab_var"),
              block$slab_var_mean, block$slab_var_var,
              block$slab_var_samples)
    } else if (block$model == "GlobalLocal") {
      add_row(paste0(block_name, "::global_var"),
              block$tau_sq_mean, block$tau_sq_var,
              block$tau_sq_samples)
    } else if (block$model == "SpikeMultiSlab") {
      add_row(paste0(block_name, "::shared_var"),
              block$var_mean, block$var_var, block$var_samples)
      for (component in names(block$pi_mean)) {
        add_row(
          paste0(block_name, "::pi::", component),
          block$pi_mean[component], block$pi_var[component],
          if (stored) block$pi_samples[, component] else NULL
        )
      }
    }
  }
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

.summarize_pve <- function(object, probs) {
  if (is.null(object$total_pve_mean)) return(NULL)
  stored <- isTRUE(object$store_samples)
  rows <- lapply(names(object$ETA), function(block_name) {
    block <- object$ETA[[block_name]]
    .posterior_summary_row(
      block$pve_mean, block$pve_var,
      if (stored) block$pve_samples else NULL, probs
    )
  })
  names(rows) <- names(object$ETA)
  rows[["Total"]] <- .posterior_summary_row(
    object$total_pve_mean, object$total_pve_var,
    if (stored) object$total_pve_samples else NULL, probs
  )
  rows[["Cross-block"]] <- .posterior_summary_row(
    object$cross_block_pve_mean, object$cross_block_pve_var,
    if (stored) object$cross_block_pve_samples else NULL, probs
  )
  do.call(rbind, rows)
}

#' @param x An object returned by `summary.blm_fit()`.
#' @rdname summary.blm_fit
#' @export
print.summary.blm_fit <- function(x, ...) {
  if (length(list(...)) > 0L) {
    stop("Additional arguments are not supported.", call. = FALSE)
  }
  cat("Bayesian linear regression posterior summary\n")
  if (x$store_samples) {
    cat(sprintf("Retained draws: %d across %d chain(s)\n\n",
                x$number_of_draws, x$nchains))
  } else {
    cat("Online posterior moments (individual draws were not retained)\n\n")
  }
  cat("Blocks:\n")
  print(x$blocks)
  if (!is.null(x$intercept)) {
    cat("\nIntercept:\n")
    print(x$intercept)
  }
  if (!is.null(x$coefficients)) {
    label <- if (x$coefficient_summary == "top") {
      sprintf(
        "Coefficients (top %d by %s):\n",
        x$reported_coefficients,
        gsub("_", " ", x$coefficient_rank_by, fixed = TRUE)
      )
    } else {
      "Coefficients:\n"
    }
    cat("\n", label, sep = "")
    print(x$coefficients)
  } else {
    cat(sprintf(
      paste0(
        "\nIndividual coefficient table omitted (%d predictors). ",
        "Use `coefficients = \"top\"` or `coefficients = \"all\"` ",
        "to request one.\n"
      ),
      x$total_coefficients
    ))
  }
  cat("\nResidual variance:\n")
  print(x$residual_variance)
  if (!is.null(x$hyperparameters)) {
    cat("\nHyperparameters:\n")
    print(x$hyperparameters)
  }
  if (!is.null(x$variance_explained)) {
    cat("\nVariance explained:\n")
    print(x$variance_explained)
  }
  invisible(x)
}

.prediction_matrix <- function(values, predictor_names) {
  number_of_predictors <- length(predictor_names)
  if (is.numeric(values) && is.atomic(values) && !is.object(values) &&
      is.null(dim(values))) {
    values <- if (number_of_predictors == 1L) {
      matrix(values, ncol = 1L)
    } else {
      matrix(values, nrow = 1L)
    }
  } else if (is.data.frame(values)) {
    if (!all(vapply(values, is.numeric, logical(1)))) {
      stop("Prediction data frames must contain only numeric columns.",
           call. = FALSE)
    }
    values <- as.matrix(values)
  } else if (!is.matrix(values) || !is.numeric(values)) {
    stop("Prediction inputs must be numeric vectors, matrices, or data frames.",
         call. = FALSE)
  }
  if (ncol(values) != number_of_predictors || nrow(values) < 1L) {
    stop("Prediction inputs must match the fitted number of predictors.",
         call. = FALSE)
  }
  if (anyNA(values) || any(!is.finite(values))) {
    stop("Prediction inputs must contain only finite, non-missing values.",
         call. = FALSE)
  }
  supplied_names <- colnames(values)
  if (!is.null(supplied_names)) {
    if (anyNA(supplied_names) || any(supplied_names == "") ||
        anyDuplicated(supplied_names) ||
        !setequal(supplied_names, predictor_names)) {
      stop("Prediction column names must match the fitted predictors.",
           call. = FALSE)
    }
    values <- values[, predictor_names, drop = FALSE]
  }
  storage.mode(values) <- "double"
  values
}
