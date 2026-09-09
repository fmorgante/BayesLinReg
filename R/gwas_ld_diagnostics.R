#' Diagnose variant-level consistency between GWAS statistics and LD
#'
#' Performs a read-only, blockwise conditional z-score diagnostic after the
#' same position and allele harmonization used by [match_gwas_ld()]. No variant
#' is removed and no allele is changed based on a statistical diagnostic.
#'
#' @param gwas A GWAS summary-statistics object accepted by [blm_gwas()].
#' @param ld A native `blm_ld` object returned by [as_blm_ld()]. Diagnostics
#'   should normally be run before conversion to `blm_ld_eigen`.
#' @param window_variants Maximum number of central variants assigned to one
#'   diagnostic window.
#' @param overlap_variants Number of neighboring variants added on either side
#'   of each central window for prediction.
#' @param ld_shrink Numeric scalar in `[0, 1)`. The diagnostic covariance is
#'   `(1 - ld_shrink) * R + ld_shrink * I`. A positive value stabilizes the
#'   conditional solves without modifying `ld`.
#' @param p_threshold Conditional-discrepancy p-value threshold used to flag
#'   variants.
#' @param min_tagging Minimum conditional tagging strength required before a
#'   variant can be flagged statistically.
#' @param conditional_variance_floor Positive numerical floor for conditional
#'   variances.
#' @param n_variation_tolerance Nonnegative relative range allowed within an LD
#'   block before its sample sizes are marked dissimilar.
#' @param nthreads Number of threads used to process independent diagnostic
#'   windows. When the package uses external BLAS, values greater than one
#'   require the applicable BLAS thread setting to be explicitly equal to one.
#' @param n_partitions Number of independent random balanced partitions evaluated
#'   for every window. With more than one partition, the reported diagnostic
#'   values come from the successfully assessed partition having the largest
#'   statistic, and `p_value` is its Bonferroni-adjusted conditional p-value.
#' @param store_variant_report Whether to retain only statistically flagged
#'   variants or the complete per-variant diagnostic table.
#'
#' @return An object of class `blm_gwas_ld_diagnostics` containing
#'   `block_report`, `variant_failures`, optionally `variant_report`, the
#'   harmonization counts, controls, and the score-correlation assumption.
#'   Per-variant output includes the raw minimum-partition p-value, its adjusted
#'   `p_value`, and the numbers of assessed and requested partitions. The block
#'   report separates covariance-factorization failures, non-finite solves, and
#'   invalid conditional calculations. With repeated partitions,
#'   `failed_windows` counts failed partition-window evaluations and can exceed
#'   `windows`.
#'
#' @details Let `z` be oriented to the LD dosage allele. Within each overlapping
#'   window, variants are randomly divided into two balanced groups. Each group
#'   is predicted from the other. For target variant `j` and predictor set `t`,
#'   the diagnostic calculates
#'   `predicted_z = R[j,t] solve(R[t,t], z[t])` and divides the squared residual
#'   by `1 - R[j,t] solve(R[t,t], R[t,j])`. The supplied `ld_shrink` is applied
#'   to all these covariance terms. Set the R seed with [set.seed()] to make the
#'   random partitions reproducible. Repeated partitions provide more than one
#'   opportunity to identify a discrepancy caused by the random predictor split.
#'   Their minimum raw conditional p-value is multiplied by `n_partitions`, so
#'   `p_value` controls the familywise error rate across the attempted splits by
#'   the Bonferroni inequality without assuming independence.
#'
#'   This is a scalable, conservative DENTIST-style first-pass diagnostic, not
#'   a reimplementation of the complete iterative DENTIST procedure. Windows
#'   are bounded by variant count and placed in bounded cross-block parallel
#'   work queues. Linear solves are batched so a global dense LD matrix is never
#'   created. Optimized BLAS can accelerate the Cholesky and triangular-solve
#'   operations.
#'
#'   Conditional p-values assume that correlations among GWAS z-scores equal
#'   the supplied reference LD. Marginal sample sizes do not identify pairwise
#'   sample overlap. Blocks exceeding `n_variation_tolerance` are marked in
#'   `block_report`, but their diagnostics are still returned as sensitivity
#'   measures and should not be used for automatic filtering.
#'
#'   The external-BLAS guard checks the documented vendor environment settings
#'   present in the R process. Vendor-specific calls that change BLAS threads
#'   after R starts cannot be detected portably and should not be used together
#'   with `nthreads > 1`.
#'
#'   `possible_allele_flip` is only statistical evidence: alleles should be
#'   changed automatically only when variant metadata establish the change.
#'   To filter explicitly, remove selected IDs from `gwas` and then call
#'   [match_gwas_ld()] to construct the corresponding native LD subset.
#'
#' @examples
#' R <- matrix(c(1, 0.95, 0.95, 1), 2)
#' variants <- data.frame(
#'   CHR = 1, ID = c("rs1", "rs2"), POS = 1:2,
#'   A1 = c("A", "A"), A0 = c("C", "C")
#' )
#' ld <- as_blm_ld(R, variants)
#' gwas <- transform(
#'   variants, N = 1000, BETA = c(0.3, -0.3), SE = 0.05
#' )
#' set.seed(1)
#' qc <- diagnose_gwas_ld(gwas, ld, window_variants = 2)
#' qc$variant_failures
#' @export
diagnose_gwas_ld <- function(
    gwas, ld, window_variants = 2000L, overlap_variants = 200L,
    ld_shrink = 0.05, p_threshold = 5e-8, min_tagging = 0.1,
    conditional_variance_floor = sqrt(.Machine$double.eps),
    n_variation_tolerance = 0.1, nthreads = 1L,
    n_partitions = 1L,
    store_variant_report = c("failures", "all")) {
  if (inherits(ld, "blm_ld_eigen")) {
    stop(
      paste0(
        "`ld` must be a native `blm_ld` object. Run diagnostics before ",
        "conversion with `as_blm_ld_eigen()`."
      ),
      call. = FALSE
    )
  }
  .validate_blm_ld_object(ld)
  gwas <- .validate_blm_gwas(gwas)
  store_variant_report <- match.arg(store_variant_report)
  nthreads <- .validate_nthreads(nthreads)
  .validate_diagnostic_blas_threads(nthreads)
  n_partitions <- .validate_qc_count(
    n_partitions, "n_partitions", minimum = 1L
  )
  window_variants <- .validate_qc_count(
    window_variants, "window_variants", minimum = 2L
  )
  overlap_variants <- .validate_qc_count(
    overlap_variants, "overlap_variants", minimum = 0L
  )
  .validate_qc_probability(ld_shrink, "ld_shrink", include_zero = TRUE)
  .validate_qc_probability(p_threshold, "p_threshold")
  .validate_qc_probability(min_tagging, "min_tagging", include_zero = TRUE)
  if (!is.numeric(conditional_variance_floor) ||
      length(conditional_variance_floor) != 1L ||
      is.na(conditional_variance_floor) ||
      !is.finite(conditional_variance_floor) ||
      conditional_variance_floor <= 0 || conditional_variance_floor >= 1) {
    stop(
      "`conditional_variance_floor` must be finite and in (0, 1).",
      call. = FALSE
    )
  }
  if (!is.numeric(n_variation_tolerance) ||
      length(n_variation_tolerance) != 1L ||
      is.na(n_variation_tolerance) || !is.finite(n_variation_tolerance) ||
      n_variation_tolerance < 0) {
    stop(
      "`n_variation_tolerance` must be nonnegative and finite.",
      call. = FALSE
    )
  }

  input_gwas_variants <- nrow(gwas)
  input_ld_variants <- nrow(ld$variants)
  harmonized <- .harmonize_gwas_ld(gwas, ld)
  matched_gwas <- harmonized$gwas
  matched_ld <- harmonized$ld
  orientation <- harmonized$orientation
  aligned_z <- matched_gwas$BETA / matched_gwas$SE * orientation
  input_z <- matched_gwas$BETA / matched_gwas$SE

  block_sizes <- vapply(matched_ld$blocks, `[[`, integer(1), "size")
  window_counts <- as.integer(ceiling(block_sizes / window_variants))
  block_ends <- cumsum(block_sizes)
  block_starts <- block_ends - block_sizes + 1L
  batches <- .diagnostic_block_batches(
    block_sizes, window_counts, nthreads
  )

  retain_all_variants <- identical(store_variant_report, "all")
  block_results <- if (retain_all_variants) {
    vector("list", length(matched_ld$blocks))
  } else {
    NULL
  }
  failure_results <- list()
  failure_index <- 0L
  variant_template <- NULL
  block_report <- vector("list", length(matched_ld$blocks))
  assessed_total <- 0
  conditional_outlier_total <- 0
  possible_flip_total <- 0
  for (batch in batches) {
    batch_global <- seq.int(
      block_starts[[batch[[1L]]]], block_ends[[batch[[length(batch)]]]]
    )
    diagnosed_batch <- .diagnose_gwas_ld_batch(
      matched_ld$blocks[batch], aligned_z[batch_global], window_variants,
      overlap_variants, ld_shrink, conditional_variance_floor, nthreads,
      n_partitions
    )
    batch_offset <- 0L
    for (local_block_index in seq_along(batch)) {
      block_index <- batch[[local_block_index]]
      block <- matched_ld$blocks[[block_index]]
      global <- seq.int(block_starts[[block_index]], block_ends[[block_index]])
      local <- seq.int(batch_offset + 1L, batch_offset + block$size)
      block_n <- matched_gwas$N[global]
      n_range <- range(block_n)
      n_ratio <- n_range[[2L]] / n_range[[1L]]
      similar_n <- n_ratio <= 1 + n_variation_tolerance
      diagnosed <- lapply(
        diagnosed_batch[c(
          "predicted_z", "conditional_variance", "tagging", "conditional_z",
          "statistic", "minimum_partition_p_value", "p_value",
          "flip_log_likelihood_ratio", "predictors_used",
          "partitions_assessed"
        )],
        function(value) value[local]
      )
      failures_by_type <- diagnosed_batch$failure_counts[local_block_index, ]
      result <- data.frame(
        CHR = matched_gwas$CHR[global],
        ID = matched_gwas$ID[global],
        POS = matched_gwas$POS[global],
        gwas_A1 = matched_gwas$A1[global],
        gwas_A0 = matched_gwas$A0[global],
        ld_A1 = matched_ld$variants$A1[global],
        ld_A0 = matched_ld$variants$A0[global],
        effect_orientation = orientation[global],
        block = block$name,
        parent = block$parent,
        N = block_n,
        input_z = input_z[global],
        ld_aligned_z = aligned_z[global],
        predicted_z = diagnosed$predicted_z,
        conditional_variance = diagnosed$conditional_variance,
        tagging = diagnosed$tagging,
        conditional_z = diagnosed$conditional_z,
        statistic = diagnosed$statistic,
        minimum_partition_p_value = diagnosed$minimum_partition_p_value,
        p_value = diagnosed$p_value,
        flip_log_likelihood_ratio = diagnosed$flip_log_likelihood_ratio,
        predictors_used = diagnosed$predictors_used,
        partitions_assessed = diagnosed$partitions_assessed,
        partitions_requested = n_partitions,
        stringsAsFactors = FALSE
      )
      assessed <- is.finite(result$statistic)
      sufficiently_tagged <- assessed & result$tagging >= min_tagging
      result$conditional_outlier <- sufficiently_tagged &
        result$p_value < p_threshold
      result$possible_allele_flip <- sufficiently_tagged &
        result$flip_log_likelihood_ratio > 2 &
        abs(result$ld_aligned_z) > 2
      result$similar_sample_sizes <- similar_n
      result$status <- ifelse(
        !assessed, "not_assessed",
        ifelse(
          result$tagging < min_tagging, "low_tagging",
          ifelse(
            result$possible_allele_flip, "possible_allele_flip",
            ifelse(result$conditional_outlier, "conditional_outlier", "ok")
          )
        )
      )
      if (is.null(variant_template)) variant_template <- result[FALSE, ]
      if (retain_all_variants) block_results[[block_index]] <- result
      failure <- result[
        result$conditional_outlier | result$possible_allele_flip,
        , drop = FALSE
      ]
      if (nrow(failure)) {
        failure_index <- failure_index + 1L
        failure_results[[failure_index]] <- failure
      }
      assessed_total <- assessed_total + sum(assessed)
      conditional_outlier_total <- conditional_outlier_total +
        sum(result$conditional_outlier)
      possible_flip_total <- possible_flip_total +
        sum(result$possible_allele_flip)
      block_report[[block_index]] <- data.frame(
        block = block$name,
        parent = block$parent,
        variants = block$size,
        windows = window_counts[[block_index]],
        partition_evaluations =
          as.double(window_counts[[block_index]]) * n_partitions,
        failed_windows = failures_by_type[["failed"]],
        factorization_failures = failures_by_type[["factorization"]],
        solve_failures = failures_by_type[["solve"]],
        invalid_conditional_windows = failures_by_type[["conditional"]],
        assessed = sum(assessed),
        low_tagging = sum(assessed & result$tagging < min_tagging),
        conditional_outliers = sum(result$conditional_outlier),
        possible_allele_flips = sum(result$possible_allele_flip),
        minimum_p_value = .finite_summary(result$p_value, min),
        median_abs_conditional_z = .finite_summary(
          abs(result$conditional_z), stats::median
        ),
        minimum_n = n_range[[1L]],
        maximum_n = n_range[[2L]],
        n_ratio = n_ratio,
        similar_sample_sizes = similar_n,
        stringsAsFactors = FALSE
      )
      batch_offset <- batch_offset + block$size
    }
  }
  failures <- if (length(failure_results)) {
    do.call(rbind, failure_results)
  } else {
    variant_template
  }
  rownames(failures) <- NULL
  answer <- list(
    summary = c(
      input_gwas_variants = input_gwas_variants,
      input_ld_variants = input_ld_variants,
      retained_variants = nrow(matched_gwas),
      assessed_variants = assessed_total,
      flagged_variants = nrow(failures),
      conditional_outliers = conditional_outlier_total,
      possible_allele_flips = possible_flip_total,
      dissimilar_n_blocks = sum(!vapply(
        block_report, `[[`, logical(1), "similar_sample_sizes"
      ))
    ),
    block_report = do.call(rbind, block_report),
    variant_failures = failures,
    harmonization_report = harmonized$counts,
    controls = list(
      window_variants = window_variants,
      overlap_variants = overlap_variants,
      ld_shrink = ld_shrink,
      p_threshold = p_threshold,
      min_tagging = min_tagging,
      conditional_variance_floor = conditional_variance_floor,
      n_variation_tolerance = n_variation_tolerance,
      nthreads = nthreads,
      n_partitions = n_partitions,
      partition_p_value_adjustment = "bonferroni",
      flip_log_likelihood_ratio_threshold = 2,
      flip_absolute_z_threshold = 2
    ),
    score_correlation_assumption =
      "gwas_score_correlation_equals_reference_ld",
    filtering_applied = FALSE
  )
  rownames(answer$block_report) <- NULL
  factorization_failure_count <- sum(
    answer$block_report$factorization_failures
  )
  if (factorization_failure_count > 0L) {
    warning(sprintf(
      paste0(
        "%d diagnostic partition-window evaluation(s) could not be ",
        "factorized. Increase ",
        "`ld_shrink` or regularize the native LD object."
      ),
      factorization_failure_count
    ), call. = FALSE)
  }
  solve_failure_count <- sum(answer$block_report$solve_failures)
  if (solve_failure_count > 0L) {
    warning(sprintf(
      paste0(
        "%d diagnostic partition-window evaluation(s) returned non-finite ",
        "solve values."
      ),
      solve_failure_count
    ), call. = FALSE)
  }
  conditional_failure_count <- sum(
    answer$block_report$invalid_conditional_windows
  )
  if (conditional_failure_count > 0L) {
    warning(sprintf(
      paste0(
        "%d diagnostic partition-window evaluation(s) produced invalid ",
        "conditional tagging or fitted values. The supplied LD may be ",
        "jointly incompatible; increase `ld_shrink` or regularize the native ",
        "LD object."
      ),
      conditional_failure_count
    ), call. = FALSE)
  }
  if (answer$summary[["dissimilar_n_blocks"]] > 0) {
    warning(sprintf(
      paste0(
        "%d LD block(s) exceed `n_variation_tolerance`; their conditional ",
        "diagnostics should be treated as sensitivity measures because ",
        "pairwise GWAS sample overlap is unknown."
      ),
      answer$summary[["dissimilar_n_blocks"]]
    ), call. = FALSE)
  }
  if (retain_all_variants) {
    variant_report <- do.call(rbind, block_results)
    rownames(variant_report) <- NULL
    answer$variant_report <- variant_report
  }
  structure(answer, class = "blm_gwas_ld_diagnostics")
}

#' @export
print.blm_gwas_ld_diagnostics <- function(x, ...) {
  cat("GWAS-LD consistency diagnostics\n")
  cat(sprintf("  Retained variants: %s\n", x$summary[["retained_variants"]]))
  cat(sprintf("  Assessed variants: %s\n", x$summary[["assessed_variants"]]))
  cat(sprintf("  Flagged variants: %s\n", x$summary[["flagged_variants"]]))
  cat(sprintf(
    "  Blocks with dissimilar N: %s\n",
    x$summary[["dissimilar_n_blocks"]]
  ))
  invisible(x)
}

.diagnostic_block_batches <- function(
    block_sizes, window_counts, nthreads, max_variants = 250000) {
  batches <- list()
  first <- 1L
  variants <- 0
  windows <- 0
  minimum_windows <- max(1, 2 * as.double(nthreads))
  for (index in seq_along(block_sizes)) {
    variants <- variants + block_sizes[[index]]
    windows <- windows + window_counts[[index]]
    final <- index == length(block_sizes)
    if (final || (variants >= max_variants && windows >= minimum_windows)) {
      batches[[length(batches) + 1L]] <- seq.int(first, index)
      first <- index + 1L
      variants <- 0
      windows <- 0
    }
  }
  batches
}

.diagnose_gwas_ld_batch <- function(
    blocks, z, window_variants, overlap_variants, ld_shrink,
    conditional_variance_floor, nthreads, n_partitions) {
  sizes <- vapply(blocks, `[[`, integer(1), "size")
  block_offsets <- cumsum(sizes) - sizes
  output_names <- c(
    "predicted_z", "conditional_variance", "tagging", "conditional_z",
    "statistic", "flip_log_likelihood_ratio"
  )
  best <- stats::setNames(
    replicate(length(output_names), rep(NA_real_, length(z)), simplify = FALSE),
    output_names
  )
  best$predictors_used <- integer(length(z))
  best_statistic <- rep(-Inf, length(z))
  partitions_assessed <- integer(length(z))
  failure_counts <- matrix(
    0L, nrow = length(blocks), ncol = 4L,
    dimnames = list(NULL, c("failed", "factorization", "solve", "conditional"))
  )

  for (partition in seq_len(n_partitions)) {
    plans <- lapply(sizes, function(size) {
      .prepare_gwas_ld_windows(size, window_variants, overlap_variants)
    })
    window_counts <- vapply(plans, function(plan) {
      length(plan$core_start)
    }, integer(1))
    window_block <- rep.int(seq_along(plans) - 1L, window_counts)
    expanded_lengths <- unlist(lapply(plans, function(plan) {
      plan$expanded_end - plan$expanded_start + 1L
    }), use.names = FALSE)
    group_offset <- c(0, cumsum(expanded_lengths))
    if (group_offset[[length(group_offset)]] > .Machine$integer.max) {
      stop(
        "A diagnostic batch exceeds the native integer indexing limit.",
        call. = FALSE
      )
    }
    diagnosed <- diagnose_gwas_ld_cpp(
      blocks = blocks,
      block_offset = as.integer(block_offsets),
      z = z,
      window_block = as.integer(window_block),
      core_start = as.integer(unlist(lapply(
        plans, `[[`, "core_start"
      ), use.names = FALSE)),
      core_end = as.integer(unlist(lapply(
        plans, `[[`, "core_end"
      ), use.names = FALSE)),
      expanded_start = as.integer(unlist(lapply(
        plans, `[[`, "expanded_start"
      ), use.names = FALSE)),
      expanded_end = as.integer(unlist(lapply(
        plans, `[[`, "expanded_end"
      ), use.names = FALSE)),
      group_offset = as.integer(group_offset),
      group = as.integer(unlist(lapply(
        plans, `[[`, "group"
      ), use.names = FALSE)),
      ld_shrink = ld_shrink,
      conditional_variance_floor = conditional_variance_floor,
      nthreads = nthreads
    )

    window_ends <- cumsum(window_counts)
    window_starts <- window_ends - window_counts + 1L
    for (block_index in seq_along(blocks)) {
      tasks <- seq.int(window_starts[[block_index]], window_ends[[block_index]])
      status <- diagnosed$window_status[tasks]
      failure_counts[block_index, "failed"] <-
        failure_counts[block_index, "failed"] + sum(status != 0L)
      failure_counts[block_index, "factorization"] <-
        failure_counts[block_index, "factorization"] +
        sum(bitwAnd(status, 1L) != 0L)
      failure_counts[block_index, "solve"] <-
        failure_counts[block_index, "solve"] +
        sum(bitwAnd(status, 2L) != 0L)
      failure_counts[block_index, "conditional"] <-
        failure_counts[block_index, "conditional"] +
        sum(bitwAnd(status, 4L) != 0L)
    }

    assessed <- is.finite(diagnosed$statistic)
    partitions_assessed <- partitions_assessed + assessed
    replace <- assessed & diagnosed$statistic > best_statistic
    if (any(replace)) {
      best_statistic[replace] <- diagnosed$statistic[replace]
      for (name in output_names) {
        best[[name]][replace] <- diagnosed[[name]][replace]
      }
      best$predictors_used[replace] <- diagnosed$predictors_used[replace]
    }
  }

  best$minimum_partition_p_value <- stats::pchisq(
    best$statistic, df = 1, lower.tail = FALSE
  )
  best$p_value <- pmin(
    1, n_partitions * best$minimum_partition_p_value
  )
  best$partitions_assessed <- partitions_assessed
  best$failure_counts <- failure_counts
  best
}

.prepare_gwas_ld_windows <- function(
    size, window_variants, overlap_variants) {
  starts <- seq.int(1L, size, by = window_variants)
  core_end <- pmin(
    size, as.double(starts) + window_variants - 1
  )
  expanded_start <- pmax(1, as.double(starts) - overlap_variants)
  expanded_end <- pmin(size, core_end + overlap_variants)
  group_lengths <- expanded_end - expanded_start + 1
  group_offset <- c(0, cumsum(group_lengths))
  group <- integer(group_offset[[length(group_offset)]])
  for (window_index in seq_along(starts)) {
    length_window <- group_lengths[[window_index]]
    permutation <- sample.int(length_window)
    window_group <- integer(length_window)
    window_group[permutation] <- rep_len(1:2, length_window)
    destination <- seq.int(
      group_offset[[window_index]] + 1,
      group_offset[[window_index + 1L]]
    )
    group[destination] <- window_group
  }
  list(
    core_start = as.integer(starts - 1L),
    core_end = as.integer(core_end - 1),
    expanded_start = as.integer(expanded_start - 1),
    expanded_end = as.integer(expanded_end - 1),
    group = group
  )
}

.validate_diagnostic_blas_threads <- function(nthreads) {
  if (nthreads <= 1L) return(invisible(NULL))
  information <- blm_build_info()
  if (!isTRUE(information$eigen_blas)) return(invisible(NULL))
  requested <- .requested_blas_threads(information$blas, Sys.getenv())
  if (!identical(requested, 1L)) {
    detail <- if (is.na(requested)) {
      "the external BLAS thread count could not be verified"
    } else {
      sprintf("the external BLAS requests %d threads", requested)
    }
    stop(sprintf(
      paste0(
        "`nthreads > 1` cannot be combined with multithreaded or ",
        "unverified external BLAS (%s). Configure the BLAS to use exactly ",
        "one thread before starting R, or use `nthreads = 1`."
      ),
      detail
    ), call. = FALSE)
  }
  invisible(NULL)
}

.requested_blas_threads <- function(blas, environment) {
  blas <- tolower(paste(blas, collapse = " "))
  get_setting <- function(name) {
    value <- unname(environment[name])
    if (!length(value) || is.na(value)) return("")
    value[[1L]]
  }
  read_setting <- function(names) {
    for (name in names) {
      value <- get_setting(name)
      if (!nzchar(value)) next
      match <- regexec("^[[:space:]]*([0-9]+)", value)
      pieces <- regmatches(value, match)[[1L]]
      if (length(pieces) == 2L) {
        count <- suppressWarnings(as.integer(pieces[[2L]]))
        if (!is.na(count) && count >= 1L) return(count)
      }
      return(NA_integer_)
    }
    NA_integer_
  }
  if (grepl("mkl|oneapi", blas)) {
    domain <- get_setting("MKL_DOMAIN_NUM_THREADS")
    if (nzchar(domain)) {
      for (scope in c("BLAS", "ALL")) {
        match <- regexec(
          paste0("MKL_", scope, "[[:space:]]*=[[:space:]]*([0-9]+)"),
          toupper(domain)
        )
        pieces <- regmatches(toupper(domain), match)[[1L]]
        if (length(pieces) == 2L) {
          count <- suppressWarnings(as.integer(pieces[[2L]]))
          if (!is.na(count) && count >= 1L) return(count)
        }
      }
    }
    return(read_setting(c("MKL_NUM_THREADS", "OMP_NUM_THREADS")))
  }
  if (grepl("openblas", blas)) {
    values <- vapply(
      c("OPENBLAS_NUM_THREADS", "GOTO_NUM_THREADS", "OMP_NUM_THREADS"),
      function(name) read_setting(name), integer(1)
    )
    values <- values[!is.na(values)]
    if (!length(values)) return(NA_integer_)
    if (any(values > 1L)) return(max(values))
    return(1L)
  }
  if (grepl("accelerate|veclib", blas)) {
    return(read_setting("VECLIB_MAXIMUM_THREADS"))
  }
  if (grepl("blis", blas)) {
    return(read_setting(c("BLIS_NUM_THREADS", "OMP_NUM_THREADS")))
  }
  if (grepl("librblas|reference blas", blas)) return(1L)
  read_setting("OMP_NUM_THREADS")
}

.validate_qc_count <- function(value, name, minimum) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop(sprintf("`%s` must be an integer of at least %d.", name, minimum),
         call. = FALSE)
  }
  as.integer(value)
}

.validate_qc_probability <- function(value, name, include_zero = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value)) {
    interval <- if (include_zero) "[0, 1)" else "(0, 1)"
    stop(sprintf("`%s` must be finite and in %s.", name, interval),
         call. = FALSE)
  }
  lower_valid <- if (include_zero) value >= 0 else value > 0
  if (!lower_valid || value >= 1) {
    interval <- if (include_zero) "[0, 1)" else "(0, 1)"
    stop(sprintf("`%s` must be finite and in %s.", name, interval),
         call. = FALSE)
  }
  invisible(value)
}

.finite_summary <- function(value, function_) {
  value <- value[is.finite(value)]
  if (!length(value)) return(NA_real_)
  function_(value)
}
