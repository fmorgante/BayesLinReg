#' Extract coefficients from a BayesLinReg fit
#'
#' Returns posterior means of the intercept and regression coefficients.
#'
#' @param object A fitted model returned by [blm()] or [blm_ss()].
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
#' @param object A fitted model returned by [blm()] or [blm_ss()].
#' @param newdata For a single-block fit, a numeric vector, matrix, or data
#'   frame containing that block's predictors. A vector represents one
#'   observation unless the fitted block has one predictor, in which case it
#'   represents multiple observations. For a multi-block fit, a named list of
#'   matrices or data frames, one for each `ETA` block. Named columns are
#'   reordered to match the fitted predictors.
#' @param ... Additional arguments. Currently unused.
#'
#' @return A numeric vector of posterior-mean predictions.
#' @export
predict.blm_fit <- function(object, newdata, ...) {
  if (missing(newdata)) {
    stop("`newdata` is required because fitted predictor data are not stored.",
         call. = FALSE)
  }
  if (length(list(...)) > 0L) {
    stop("Additional arguments are not supported.", call. = FALSE)
  }

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
  unname(prediction)
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
