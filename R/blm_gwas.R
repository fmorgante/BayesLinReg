#' Bayesian linear regression from GWAS summary statistics and LD
#'
#' Fits the same coefficient-prior models as [blm_ss()] from additive
#' quantitative-trait GWAS results and a reusable LD object created by
#' [as_blm_ld()] or its truncated-eigen analogue [as_blm_ld_eigen()]. The
#' working Gram matrix is applied without materializing a global dense matrix.
#'
#' @param gwas A data frame with columns `CHR`, `ID`, `POS`, `A1`, `A0`, `N`,
#'   `BETA`, and `SE`. `BETA` is the additive marginal effect per `A1` allele.
#'   `N` may vary among retained variants; noninteger effective sample sizes
#'   are allowed.
#' @param ld A `blm_ld` object returned by [as_blm_ld()] or a
#'   `blm_ld_eigen` object returned by [as_blm_ld_eigen()].
#' @param ETA Prior specifications in the same format as [blm_ss()]. Character
#'   `indices` refer to variant IDs.
#' @param residual_var,residual_shape,residual_scale Residual-variance controls
#'   with the same meaning and defaults as in [blm_ss()].
#' @param reference_response_var Reference response variance used for original-
#'   scale reconstruction and expected-PVE calibration. It has the same prior-
#'   calibration meaning as in [blm_ss()].
#' @param iterations,burnin,thin MCMC controls; see [blm_ss()].
#' @param verbose,nchains,nthreads Progress, chain, and within-chain block
#'   parallelism controls; see [blm_ss()].
#' @param store_samples,store_coefficient_cov Posterior-storage controls; see
#'   [blm_ss()].
#' @param check_psd Whether to perform the optional blockwise PSD and joint
#'   compatibility validation; see [blm_ss()].
#' @param compute_pve,pve_type Posterior PVE controls; see [blm_ss()].
#' @param scale Working scale. `"auto"` uses the original response scale when
#'   `reference_response_var` is supplied and otherwise standardizes the
#'   response and predictors. `"original"` requires
#'   `reference_response_var`.
#' @param residual_df_gwas Positive residual degrees of freedom used by the
#'   marginal GWAS regressions. It may be a scalar or one value per input GWAS
#'   row. The default is `N - 2` elementwise.
#' @param ld_shrink Numeric scalar in `[0, 1)`. Off-diagonal LD correlations
#'   are multiplied by `1 - ld_shrink` while the unit diagonal is retained.
#'   This is applied by the native LD operator without copying or modifying
#'   `ld`. Positive values can stabilize analyses using external LD. It is not
#'   available for `blm_ld_eigen` input because the eigen object already
#'   represents an explicit LD approximation.
#'
#' @return An object of class `blm_fit`. Coefficients are oriented to the input
#'   GWAS `A1` alleles. No intercept is fitted because centered GWAS summary
#'   statistics do not identify the phenotype mean. Nevertheless, their
#'   centered likelihood has `likelihood_df = N - 1`; this value is recorded
#'   in the result and is separate from `residual_df_gwas`.
#'   `gwas_sample_size_summary` describes retained `N`, and `gwas_likelihood`
#'   identifies the likelihood construction. With heterogeneous `N`, the fit
#'   also records `gwas_reference_n` and `gwas_sample_overlap_assumption`, and
#'   it has no single `likelihood_df`. When `ld` was created by
#'   [regularize_blm_ld()], `ld_regularization_report` preserves that object's
#'   original regularization audit report unchanged, while
#'   `ld_regularization_block_map` maps fitted post-harmonization blocks to the
#'   source blocks in that report. `ld_harmonization` reports
#'   retained and flipped variants together with separate counts for GWAS-only,
#'   LD-only, location-mismatched, allele-mismatched, and ambiguous variants.
#'   Its `excluded` element counts excluded table entries: unmatched entries
#'   contribute one and matched-but-incompatible variant pairs contribute two.
#'   Eigen LD fits additionally report `ld_eigen_rank`, block-specific ranks
#'   and retained trace fractions, aggregate `ld_prop_var`, and whether the
#'   representation is approximate. Trace fractions use the unit-diagonal
#'   trace under the strict eigenvalue policy and the available positive trace
#'   under `negative_eigenvalues = "discard"`.
#'
#' @details Variants are matched by `ID` and checked against chromosome,
#'   position, and alleles. Reversed alleles are handled by changing effect
#'   orientation. With native LD, unmatched, incompatible, and unresolved
#'   strand-ambiguous variants are excluded with a warning. Excluded IDs are
#'   also removed from character-indexed `ETA` blocks, with a block-specific
#'   warning; a block that becomes empty is rejected. If harmonization excludes
#'   an LD variant, numeric `ETA` indices are rejected because their intended
#'   position is ambiguous; use character variant IDs instead. Eigen LD instead
#'   requires complete compatible coverage of its LD variant panel. Gibbs
#'   coordinates retain LD order, while `ETA` blocks remain independent of LD
#'   blocks.
#'
#'   A `blm_ld_eigen` object uses the pure truncated representation
#'   \eqn{R_q=U_q\Lambda_qU_q'} stored in that object. The omitted eigenspace
#'   is treated as having zero LD variance; no diagonal correction is added.
#'   Eigen LD requires complete coverage of its variant panel. GWAS rows may be
#'   reordered and may contain additional variants, but every eigen-LD variant
#'   must have compatible GWAS position and allele metadata. Match a native
#'   `blm_ld` object with [match_gwas_ld()] before eigen decomposition when
#'   variants need to be removed.
#'   Predictor standardization and expected-PVE prior calibration continue to
#'   use the original unit-diagonal LD scale, so changing `prop_var` does not
#'   silently redefine prior scales.
#'
#'   With `scale = "standardized"`, the working statistics use
#'   `XtX = (N - 1) R` and response variance one. With `scale = "original"`,
#'   predictor cross-product diagonals are reconstructed from `BETA`, `SE`,
#'   `residual_df_gwas`, and `reference_response_var`. Reference-panel LD and
#'   GWAS results not obtained by common-sample ordinary least squares define
#'   approximate working sufficient statistics.
#'   The centered working likelihood uses `N - 1` degrees of freedom for the
#'   residual inverse-gamma update and PVE normalization even though no
#'   intercept can be returned from the summary statistics.
#'
#'   When LD was estimated outside the GWAS sample, fixing `residual_var` is
#'   recommended. Learning it requires the reconstructed `XtX`, `Xty`, and
#'   `yty` to be mutually compatible; reference-panel mismatch can violate
#'   that requirement. For standardized quantitative traits,
#'   `residual_var = 1` is a conservative robust choice. `ld_shrink` addresses
#'   LD regularization but does not make approximate statistics exact.
#'
#'   When `N` differs among variants (or contains noninteger effective sample
#'   sizes), `blm_gwas()` uses a heterogeneous-sample-size RSS approximation.
#'   Each marginal statistic receives its own likelihood scaling, while the
#'   score-correlation matrix is assumed to equal the supplied LD matrix. This
#'   assumes essentially common sample overlap; marginal `N` values alone do
#'   not identify pairwise sample overlap. The residual variance must be fixed
#'   at one and `check_psd = TRUE` is unavailable in this mode. PVE continues
#'   to use predictor covariance rather than likelihood information, and a
#'   separate PVE scaling is applied internally. The median `N`, rounded to an
#'   integer, is used only as a computational and GlobalLocal prior-calibration
#'   reference.
#' @export
blm_gwas <- function(
    gwas, ld, ETA, residual_var = NULL,
    residual_shape = NULL, residual_scale = NULL,
    reference_response_var = NULL,
    scale = c("auto", "standardized", "original"),
    residual_df_gwas = NULL,
    ld_shrink = 0,
    iterations = 4000L, burnin = 1000L, thin = 1L,
    verbose = FALSE, nchains = 1L, nthreads = 1L,
    store_samples = FALSE, store_coefficient_cov = FALSE,
    check_psd = FALSE, compute_pve = FALSE,
    pve_type = c("standalone", "allocated")) {
  scale <- match.arg(scale)
  if (!is.numeric(ld_shrink) || length(ld_shrink) != 1L ||
      is.na(ld_shrink) || !is.finite(ld_shrink) ||
      ld_shrink < 0 || ld_shrink >= 1) {
    stop("`ld_shrink` must be a finite numeric scalar in [0, 1).",
         call. = FALSE)
  }
  ld_shrink <- as.numeric(ld_shrink)
  pve_controls <- .validate_pve_controls(compute_pve, pve_type)
  compute_pve <- pve_controls$compute_pve
  pve_type <- pve_controls$pve_type
  eigen_ld <- inherits(ld, "blm_ld_eigen")
  if (eigen_ld) {
    .validate_blm_ld_eigen_object(ld)
    if (ld_shrink != 0) {
      stop("`ld_shrink` must be zero for `blm_ld_eigen` input.",
           call. = FALSE)
    }
  } else {
    .validate_blm_ld_object(ld)
  }
  input_ld_regularization_report <- ld$regularization_report
  controls <- list(
    verbose = verbose,
    store_samples = store_samples,
    store_coefficient_cov = store_coefficient_cov,
    check_psd = check_psd
  )
  for (name in names(controls)) {
    value <- controls[[name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
    }
  }
  nchains <- .validate_nchains(nchains)
  nthreads <- .validate_nthreads(nthreads)
  if (nthreads > 1L && nchains != 1L) {
    stop("`nthreads > 1` requires `nchains = 1`.", call. = FALSE)
  }
  if (!is.null(reference_response_var)) {
    .validate_variance(reference_response_var, "reference_response_var")
  }
  if (scale == "auto") {
    scale <- if (is.null(reference_response_var)) "standardized" else "original"
  }
  if (scale == "original" && is.null(reference_response_var)) {
    stop(
      "`reference_response_var` is required for `scale = \"original\"`.",
      call. = FALSE
    )
  }
  if (scale == "standardized" && !is.null(reference_response_var) &&
      reference_response_var != 1) {
    stop(
      paste0(
        "`reference_response_var` must be NULL or one for ",
        "`scale = \"standardized\"`."
      ),
      call. = FALSE
    )
  }

  gwas <- .validate_blm_gwas(gwas)
  input_gwas_ids <- gwas$ID
  residual_df_by_id <- NULL
  if (!is.null(residual_df_gwas) && length(residual_df_gwas) != 1L) {
    if (length(residual_df_gwas) != nrow(gwas)) {
      stop(
        paste0(
          "`residual_df_gwas` must be a scalar or have one value per input ",
          "GWAS row."
        ),
        call. = FALSE
      )
    }
    residual_df_by_id <- stats::setNames(residual_df_gwas, input_gwas_ids)
  }
  harmonized <- .harmonize_gwas_ld(
    gwas, ld, require_complete_ld = eigen_ld
  )
  gwas <- harmonized$gwas
  if (!is.null(residual_df_by_id)) {
    residual_df_gwas <- unname(residual_df_by_id[gwas$ID])
  }
  orientation <- harmonized$orientation
  ETA <- .harmonize_gwas_eta(
    ETA, gwas$ID, input_gwas_ids, ld$variants$ID
  )
  ld <- harmonized$ld
  n_by_variant <- gwas$N
  n_range <- range(n_by_variant)
  heterogeneous_n <- n_range[[1L]] != n_range[[2L]] ||
    n_range[[1L]] != floor(n_range[[1L]])
  sample_size_quantiles <- if (heterogeneous_n) {
    unname(stats::quantile(n_by_variant, c(0.25, 0.5, 0.75)))
  } else {
    rep(n_range[[1L]], 3L)
  }
  reference_n <- if (heterogeneous_n) {
    as.integer(round(sample_size_quantiles[[2L]]))
  } else {
    as.integer(n_range[[1L]])
  }
  likelihood_df <- .resolve_likelihood_df(
    reference_n, FALSE, reference_n - 1L
  )
  if (heterogeneous_n) {
    if (is.null(residual_var) || !is.numeric(residual_var) ||
        length(residual_var) != 1L ||
        is.na(residual_var) || !is.finite(residual_var) || residual_var != 1) {
      stop(
        paste0(
          "Heterogeneous `gwas$N` requires `residual_var = 1`; residual ",
          "variance learning is not identified by this RSS approximation."
        ),
        call. = FALSE
      )
    }
    if (check_psd) {
      stop(
        "`check_psd = TRUE` is unavailable with heterogeneous `gwas$N`.",
        call. = FALSE
      )
    }
    warning(
      paste0(
        "Heterogeneous `gwas$N` uses an RSS approximation that assumes ",
        "score correlations equal the supplied reference LD; pairwise ",
        "sample overlap is not modeled."
      ),
      call. = FALSE
    )
  }
  if (is.null(residual_df_gwas)) {
    residual_df_gwas <- if (heterogeneous_n) {
      n_by_variant - 2
    } else {
      reference_n - 2
    }
  }
  if (!is.numeric(residual_df_gwas) ||
      !length(residual_df_gwas) %in% c(1L, length(n_by_variant)) ||
      anyNA(residual_df_gwas) || any(!is.finite(residual_df_gwas)) ||
      any(residual_df_gwas <= 0) ||
      any(residual_df_gwas > n_by_variant - 1)) {
    stop(
      paste0(
        "`residual_df_gwas` must be positive, with every value no greater ",
        "than `N - 1`."
      ),
      call. = FALSE
    )
  }
  if (length(residual_df_gwas) == 1L && heterogeneous_n) {
    residual_df_gwas <- rep(residual_df_gwas, length(n_by_variant))
  }

  working_reference_var <- if (scale == "standardized") {
    1
  } else {
    reference_response_var
  }
  components <- .gwas_working_components(
    gwas$BETA * orientation, gwas$SE, n_by_variant, residual_df_gwas,
    working_reference_var, scale, reference_n
  )
  predictor_names <- gwas$ID
  names(components$Xty) <- predictor_names
  normalized <- .normalize_ss_eta(
    ETA, predictor_names, residual_var, reference_n
  )
  blocks <- normalized$blocks
  source_indices <- normalized$source_indices
  has_expected_pve <- any(vapply(
    blocks, function(block) !is.null(block$expected_pve), logical(1)
  ))
  predictor_scales <- lapply(seq_along(blocks), function(block_index) {
    indices <- source_indices[[block_index]]
    if (!blocks[[block_index]]$standardize) return(rep(1, length(indices)))
    sqrt(components$diagonal[indices] / (n_by_variant[indices] - 1))
  })
  for (block_index in seq_along(blocks)) {
    blocks[[block_index]]$predictor_scale <- predictor_scales[[block_index]]
  }
  if (has_expected_pve) {
    predictor_variance_sums <- vapply(seq_along(blocks), function(block_index) {
      indices <- source_indices[[block_index]]
      sum(
        components$diagonal[indices] /
          predictor_scales[[block_index]]^2 /
          (n_by_variant[indices] - 1)
      )
    }, numeric(1))
    blocks <- .calibrate_eta_priors(
      blocks, predictor_variance_sums, components$reference_response_var,
      reference_n
    )
  }
  learn_residual_var <- is.null(residual_var)
  residual_prior <- .prepare_residual_prior(
    residual_var, residual_shape, residual_scale, blocks,
    components$reference_response_var
  )
  residual_var <- residual_prior$residual_var
  residual_shape <- residual_prior$residual_shape
  residual_scale <- residual_prior$residual_scale

  layout <- .prepare_block_layout(blocks, source_indices, predictor_scales)
  block_model <- layout$block_model
  p <- length(predictor_names)
  sampler_block_id <- integer(p)
  sampler_predictor_scale <- numeric(p)
  sampler_internal_names <- character(p)
  for (block_index in seq_along(blocks)) {
    indices <- source_indices[[block_index]]
    sampler_block_id[indices] <- block_index
    sampler_predictor_scale[indices] <- predictor_scales[[block_index]]
    sampler_internal_names[indices] <- paste0(
      names(blocks)[[block_index]], "::", blocks[[block_index]]$predictor_names
    )
  }
  sampler_layout <- layout
  sampler_layout$block_id <- sampler_block_id
  sampler_layout$internal_names <- sampler_internal_names
  source_scale <- sqrt(components$diagonal)
  sampler_scale <- source_scale / sampler_predictor_scale
  pve_scale <- if (heterogeneous_n) {
    sqrt(
      likelihood_df * components$diagonal / (n_by_variant - 1)
    ) / sampler_predictor_scale
  } else {
    sampler_scale
  }
  working_Xty <- components$Xty / sampler_predictor_scale
  names(working_Xty) <- sampler_internal_names
  ld_block_ends <- cumsum(vapply(ld$blocks, `[[`, integer(1), "size"))
  ld_block_starts <- ld_block_ends -
    vapply(ld$blocks, `[[`, integer(1), "size") + 1L
  ld_indices <- Map(seq.int, ld_block_starts, ld_block_ends)

  fixed_blocks <- vapply(
    blocks, function(block) block$model == "Fixed", logical(1)
  )
  fixed_source <- unlist(source_indices[fixed_blocks], use.names = FALSE)
  if (eigen_ld) {
    transformed_X_blocks <- vector("list", length(ld$blocks))
    separate_pve_factor <- heterogeneous_n && compute_pve
    pve_X_blocks <- if (separate_pve_factor) {
      vector("list", length(ld$blocks))
    } else {
      NULL
    }
    transformed_y_blocks <- vector("list", length(ld$blocks))
    projected_Xty <- numeric(p)
    approximate_diagonal <- numeric(p)
    for (block_index in seq_along(ld$blocks)) {
      block <- ld$blocks[[block_index]]
      indices <- ld_indices[[block_index]]
      block_scale <- sampler_scale[indices]
      scale_tolerance <- 100 * .Machine$double.eps *
        max(1, max(abs(block_scale)))
      common_scale <- max(abs(block_scale - block_scale[[1L]])) <=
        scale_tolerance
      factor_scale <- if (common_scale) {
        rep(block_scale[[1L]], length(block_scale))
      } else {
        block_scale
      }
      factor <- build_scaled_eigen_factor_cpp(
        block$eigenvectors, block$eigenvalues,
        1 / factor_scale, integer()
      )
      if (separate_pve_factor) {
        pve_X_blocks[[block_index]] <- build_scaled_eigen_factor_cpp(
          block$eigenvectors, block$eigenvalues,
          1 / pve_scale[indices], integer()
        )
      }
      prepared <- if (common_scale) {
        prepare_eigen_statistics_cpp(
          block$eigenvectors,
          block$eigenvalues * block_scale[[1L]]^2,
          working_Xty[indices]
        )
      } else {
        prepare_eigen_factor_statistics_cpp(factor, working_Xty[indices])
      }
      transformed_X_blocks[[block_index]] <- factor
      transformed_y_blocks[[block_index]] <- prepared$transformed_response
      projected_Xty[indices] <- prepared$projected_crossproduct
      approximate_diagonal[indices] <- prepared$diagonal
      if (block$complete_eigenspace &&
          block$discarded_negative_eigenvalues == 0L) {
        projection_error <- working_Xty[indices] -
          prepared$projected_crossproduct
        projection_tolerance <- sqrt(.Machine$double.eps) *
          max(1, sqrt(sum(working_Xty[indices]^2)))
        if (sqrt(sum(projection_error^2)) > projection_tolerance) {
          stop(sprintf(
            paste0(
              "The GWAS cross-products have a component outside the exact ",
              "eigenspace of LD block `%s`."
            ),
            block$name
          ), call. = FALSE)
        }
      }
    }
    if (!separate_pve_factor) pve_X_blocks <- transformed_X_blocks
    rm(factor, prepared)
    diagonal_tolerance <- 100 * .Machine$double.eps *
      pmax(1, approximate_diagonal)
    constant_predictors <- approximate_diagonal <= diagonal_tolerance
    if (any(constant_predictors)) {
      stop(sprintf(
        "The eigen LD representation contains constant predictor(s): %s.",
        paste(predictor_names[constant_predictors], collapse = ", ")
      ), call. = FALSE)
    }
    transformed_y_norm <- sum(vapply(
      transformed_y_blocks, function(value) sum(value^2), numeric(1)
    ))
    compatibility_tolerance <- sqrt(.Machine$double.eps) *
      max(1, components$yty)
    if ((learn_residual_var || check_psd) &&
        transformed_y_norm > components$yty + compatibility_tolerance) {
      stop(
        paste0(
          "The GWAS cross-products and eigen LD representation are ",
          "incompatible with the reconstructed response sum of squares."
        ),
        call. = FALSE
      )
    }
    if (length(fixed_source)) {
      fixed_design <- do.call(rbind, lapply(seq_along(ld$blocks), function(b) {
        answer <- matrix(
          0, nrow = nrow(transformed_X_blocks[[b]]),
          ncol = length(fixed_source)
        )
        selected <- match(ld_indices[[b]], fixed_source, nomatch = 0L)
        keep <- selected > 0L
        answer[, selected[keep]] <-
          transformed_X_blocks[[b]][, keep, drop = FALSE]
        answer
      }))
      .validate_fixed_design(
        fixed_design, seq_along(fixed_source),
        sampler_internal_names[fixed_source]
      )
    }
  } else if (length(fixed_source)) {
    fixed_R <- .materialize_blm_ld(ld, fixed_source, ld_shrink)
    fixed_scale <- sampler_scale[fixed_source]
    fixed_gram <- fixed_R * tcrossprod(fixed_scale)
    .validate_fixed_gram(
      fixed_gram, seq_along(fixed_source),
      sampler_internal_names[fixed_source]
    )
  }
  if (check_psd && !eigen_ld) {
    .validate_block_working_crossproducts(
      ld$blocks, ld_indices, working_Xty, components$yty,
      transform = function(block, block_indices) {
        matrix <- .materialize_ld_block(block)
        if (ld_shrink > 0) {
          matrix <- (1 - ld_shrink) * matrix
          diag(matrix) <- 1
        }
        matrix * tcrossprod(sampler_scale[block_indices])
      }
    )
  }

  sampler_arguments <- .prepare_sampler_arguments(
    blocks = blocks,
    layout = sampler_layout,
    y = numeric(),
    x = matrix(numeric(), nrow = 0L, ncol = length(working_Xty)),
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    residual_var = residual_var,
    iterations = iterations,
    burnin = burnin,
    thin = thin,
    store_samples = store_samples,
    store_coefficient_cov = store_coefficient_cov,
    compute_pve = compute_pve,
    pve_type = pve_type,
    effective_n = reference_n,
    likelihood_df = likelihood_df,
    fit_intercept = FALSE,
    intercept_x_mean = numeric(length(working_Xty)),
    intercept_y_mean = 0
  )
  if (eigen_ld) {
    sampler_arguments$eigen_X <- transformed_X_blocks
    sampler_arguments$eigen_pve_X <- pve_X_blocks
    sampler_arguments$eigen_y <- transformed_y_blocks
    sampler_arguments$eigen_indices <- ld_indices
    sampler_arguments$Xty <- projected_Xty
  } else {
    sampler_arguments$ld_blocks <- lapply(ld$blocks, function(block) {
      block[c("type", "size", "data", "indptr", "row_index")]
    })
    sampler_arguments$ld_indices <- ld_indices
    sampler_arguments$ld_scale <- sampler_scale
    sampler_arguments$ld_pve_scale <- pve_scale
    sampler_arguments$ld_shrink <- ld_shrink
    sampler_arguments$Xty <- working_Xty
  }
  sampler_arguments$yty <- components$yty
  sampler_arguments$nthreads <- nthreads
  samples <- .run_prepared_sampler(
    sampler_arguments, "Rcpp", nchains, block_model, verbose, iterations
  )

  result <- .assemble_blm_result(
    blocks, source_indices, samples, nchains, store_samples,
    store_coefficient_cov, FALSE,
    compute_pve = compute_pve, pve_type = pve_type,
    residual_shape = residual_shape,
    residual_scale = residual_scale,
    residual_scale_calibrated = residual_prior$residual_scale_calibrated,
    expected_pve_total = residual_prior$expected_pve_total,
    reference_response_var = components$reference_response_var,
    reference_residual_var = residual_prior$reference_residual_var,
    likelihood_df = if (heterogeneous_n) NULL else likelihood_df,
    sampler_block_id = sampler_block_id
  )
  result <- .orient_gwas_coefficients(
    result, orientation, source_indices, store_samples,
    store_coefficient_cov
  )
  result$gwas_variants <- gwas[c("CHR", "ID", "POS", "A1", "A0", "N")]
  result$gwas_scale <- scale
  result$residual_df_gwas <- residual_df_gwas
  result$reference_response_var <- components$reference_response_var
  result$gwas_likelihood <- if (heterogeneous_n) {
    "rss_heterogeneous_n"
  } else {
    "reconstructed_sufficient_statistics"
  }
  result$gwas_sample_size_summary <- c(
    minimum = n_range[[1L]],
    first_quartile = sample_size_quantiles[[1L]],
    median = sample_size_quantiles[[2L]],
    mean = mean(n_by_variant),
    third_quartile = sample_size_quantiles[[3L]],
    maximum = n_range[[2L]]
  )
  if (heterogeneous_n) {
    result$gwas_reference_n <- reference_n
    result$gwas_sample_overlap_assumption <-
      "score_correlation_equals_reference_ld"
  }
  result$ld_block_table <- ld$block_table
  result$ld_cross_block_assumption <- ld$cross_block_assumption
  result$ld_representation <- if (eigen_ld) "truncated_eigen" else "explicit"
  if (eigen_ld) {
    block_rank <- vapply(ld$blocks, `[[`, integer(1), "rank")
    block_prop_var <- vapply(ld$blocks, `[[`, numeric(1), "prop_var")
    names(block_rank) <- names(ld$blocks)
    names(block_prop_var) <- names(ld$blocks)
    result$ld_eigen_rank <- sum(block_rank)
    result$ld_eigen_rank_by_block <- block_rank
    result$ld_prop_var <- sum(vapply(
      ld$blocks, `[[`, numeric(1), "retained_trace"
    )) / sum(vapply(ld$blocks, `[[`, numeric(1), "trace_basis"))
    result$ld_prop_var_by_block <- block_prop_var
    result$ld_approximate <- any(
      !vapply(ld$blocks, `[[`, logical(1), "complete_eigenspace") |
      vapply(
        ld$blocks, `[[`, integer(1), "discarded_negative_eigenvalues"
      ) > 0L
    )
    result$ld_eigen_requested_prop_var <- ld$requested_prop_var
  }
  if (!is.null(input_ld_regularization_report)) {
    result$ld_regularization_report <- input_ld_regularization_report
    result$ld_regularization_block_map <- .ld_regularization_block_map(
      input_ld_regularization_report, ld$regularization_report
    )
  }
  result$ld_shrink <- ld_shrink
  result$ld_harmonization <- harmonized$counts
  result$nthreads <- nthreads
  result
}

#' Match GWAS summary statistics to native LD
#'
#' Harmonizes GWAS summary statistics with a native [blm_ld][as_blm_ld()]
#' object and keeps only variants with compatible identifiers, positions, and
#' alleles in both inputs. This is intended to establish the final variant
#' panel before calling [as_blm_ld_eigen()].
#'
#' @param gwas A GWAS summary-statistics object accepted by [blm_gwas()].
#' @param ld A native `blm_ld` object returned by [as_blm_ld()]. Eigen LD input
#'   is rejected because filtering must occur before eigen decomposition.
#'
#' @return A list containing `gwas` and `ld` restricted to the same variants in
#'   LD order, `orientation` giving the GWAS-to-LD effect orientation,
#'   `retained_ids`, `excluded_gwas_ids`, `excluded_ld_ids`, and a named count
#'   vector `report`. GWAS alleles and effect estimates retain their input
#'   orientation; [blm_gwas()] handles any required coefficient reorientation.
#'
#' @details Matching uses `ID`, then checks chromosome, position, and alleles.
#' Reversed and complementary alleles are retained, while unresolved
#' strand-ambiguous variants are excluded. Subsetting preserves the native LD
#' representation, detects newly created exact contiguous sub-blocks, and
#' carries forward regularization provenance.
#'
#' @examples
#' R <- matrix(c(1, 0.2, 0.2, 1), 2)
#' variants <- data.frame(
#'   CHR = 1, ID = c("rs1", "rs2"), POS = 1:2,
#'   A1 = c("A", "C"), A0 = c("C", "T")
#' )
#' ld <- as_blm_ld(R, variants)
#' gwas <- transform(variants, N = 1000, BETA = c(0.1, -0.1), SE = 0.05)
#' matched <- match_gwas_ld(gwas, ld)
#' eigen_ld <- as_blm_ld_eigen(matched$ld)
#' @export
match_gwas_ld <- function(gwas, ld) {
  if (inherits(ld, "blm_ld_eigen")) {
    stop(
      "`ld` must be a native `blm_ld` object, not `blm_ld_eigen`.",
      call. = FALSE
    )
  }
  .validate_blm_ld_object(ld)
  gwas <- .validate_blm_gwas(gwas)
  matched <- .harmonize_gwas_ld(gwas, ld)
  retained_ids <- matched$gwas$ID
  list(
    gwas = matched$gwas,
    ld = matched$ld,
    orientation = stats::setNames(matched$orientation, retained_ids),
    retained_ids = retained_ids,
    excluded_gwas_ids = gwas$ID[!gwas$ID %in% retained_ids],
    excluded_ld_ids = ld$variants$ID[
      !ld$variants$ID %in% retained_ids
    ],
    report = matched$counts
  )
}

.validate_blm_gwas <- function(gwas) {
  if (!is.data.frame(gwas)) {
    stop("`gwas` must be a data frame or data-frame-like R object.",
         call. = FALSE)
  }
  required <- c("CHR", "ID", "POS", "A1", "A0", "N", "BETA", "SE")
  missing <- setdiff(required, names(gwas))
  if (length(missing)) {
    stop(sprintf("`gwas` is missing required column(s): %s.",
                 paste(missing, collapse = ", ")), call. = FALSE)
  }
  result <- as.data.frame(gwas, stringsAsFactors = FALSE)
  for (name in c("CHR", "ID", "A1", "A0")) {
    result[[name]] <- as.character(result[[name]])
    if (anyNA(result[[name]]) || any(result[[name]] == "")) {
      stop(sprintf("`gwas$%s` must not contain missing or empty values.", name),
           call. = FALSE)
    }
  }
  if (anyDuplicated(result$ID)) {
    stop("`gwas$ID` must contain unique variant identifiers.", call. = FALSE)
  }
  if (!is.numeric(result$POS) || anyNA(result$POS) ||
      any(!is.finite(result$POS)) || any(result$POS != floor(result$POS)) ||
      any(result$POS < 1)) {
    stop("`gwas$POS` must contain valid positive integers.", call. = FALSE)
  }
  result$POS <- as.numeric(result$POS)
  if (!is.numeric(result$N) || anyNA(result$N) ||
      any(!is.finite(result$N)) || any(result$N < 3) ||
      any(result$N > .Machine$integer.max)) {
    stop(
      paste0(
        "`gwas$N` must contain positive finite sample sizes of at least ",
        "three and no greater than `.Machine$integer.max`."
      ),
      call. = FALSE
    )
  }
  result$N <- as.numeric(result$N)
  if (!is.numeric(result$BETA) || anyNA(result$BETA) ||
      any(!is.finite(result$BETA))) {
    stop("`gwas$BETA` must contain finite numeric values.", call. = FALSE)
  }
  if (!is.numeric(result$SE) || anyNA(result$SE) ||
      any(!is.finite(result$SE)) || any(result$SE <= 0)) {
    stop("`gwas$SE` must contain positive finite numeric values.",
         call. = FALSE)
  }
  result$A1 <- toupper(result$A1)
  result$A0 <- toupper(result$A0)
  if (any(result$A1 == result$A0)) {
    stop("`gwas$A1` and `gwas$A0` must differ.", call. = FALSE)
  }
  result
}

.complement_allele <- function(allele) {
  unname(c(A = "T", T = "A", C = "G", G = "C")[allele])
}

.harmonize_gwas_eta <- function(ETA, retained_ids, gwas_ids, ld_ids) {
  if (!is.list(ETA)) return(ETA)
  single_block <- "model" %in% names(ETA)
  specifications <- if (single_block) list(ETA = ETA) else ETA
  if (!length(specifications) ||
      !all(vapply(specifications, is.list, logical(1)))) {
    return(ETA)
  }
  block_names <- names(specifications)
  if (is.null(block_names)) {
    block_names <- paste0("ETA", seq_along(specifications))
  } else {
    missing_name <- is.na(block_names) | block_names == ""
    block_names[missing_name] <- paste0("ETA", which(missing_name))
    block_names <- make.unique(block_names)
  }
  removed <- integer(length(specifications))
  ld_variants_excluded <- length(retained_ids) < length(ld_ids)
  for (block_index in seq_along(specifications)) {
    indices <- specifications[[block_index]]$indices
    if (ld_variants_excluded && is.numeric(indices)) {
      stop(sprintf(
        paste0(
          "ETA block `%s` uses numeric `indices`, which cannot be safely ",
          "remapped after GWAS-LD harmonization excludes LD variants; use ",
          "character variant IDs."
        ),
        block_names[[block_index]]
      ), call. = FALSE)
    }
    if (!is.character(indices)) next
    known <- indices %in% gwas_ids | indices %in% ld_ids
    unknown <- unique(indices[!known])
    if (length(unknown)) {
      stop(sprintf(
        "Character `indices` in ETA block `%s` contain unknown variant ID(s): %s.",
        block_names[[block_index]], paste(unknown, collapse = ", ")
      ), call. = FALSE)
    }
    keep <- indices %in% retained_ids
    removed[[block_index]] <- sum(!keep)
    if (!any(keep)) {
      stop(sprintf(
        "ETA block `%s` has no predictors remaining after GWAS-LD harmonization.",
        block_names[[block_index]]
      ), call. = FALSE)
    }
    specifications[[block_index]]$indices <- indices[keep]
  }
  changed <- which(removed > 0L)
  if (length(changed)) {
    detail <- paste0(block_names[changed], " (", removed[changed], ")")
    warning(sprintf(
      paste0(
        "GWAS-LD harmonization removed excluded character-indexed ",
        "predictors from ETA block(s): %s."
      ),
      paste(detail, collapse = ", ")
    ), call. = FALSE)
  }
  if (single_block) specifications[[1L]] else specifications
}

.harmonize_gwas_ld <- function(
    gwas, ld, require_complete_ld = FALSE) {
  match_index <- match(ld$variants$ID, gwas$ID)
  id_match <- !is.na(match_index)
  candidate <- which(id_match)
  if (!length(candidate)) {
    if (require_complete_ld) {
      stop(
        paste0(
          "The GWAS statistics do not cover the eigen-LD variant panel. ",
          "Match the GWAS and native `blm_ld` object with ",
          "`match_gwas_ld()` before calling `as_blm_ld_eigen()`, or ",
          "impute the missing summary statistics."
        ),
        call. = FALSE
      )
    }
    stop("No GWAS variants match the LD object by `ID`.", call. = FALSE)
  }
  rows <- match_index[candidate]
  ld_variants <- ld$variants[candidate, , drop = FALSE]
  input <- gwas[rows, , drop = FALSE]
  location_match <- as.character(input$CHR) == as.character(ld_variants$CHR) &
    input$POS == ld_variants$POS
  palindromic <- paste0(input$A1, input$A0) %in%
    c("AT", "TA", "CG", "GC")
  direct <- input$A1 == ld_variants$A1 & input$A0 == ld_variants$A0
  reversed <- input$A1 == ld_variants$A0 & input$A0 == ld_variants$A1
  complement_a1 <- .complement_allele(input$A1)
  complement_a0 <- .complement_allele(input$A0)
  complement_direct <- !is.na(complement_a1) & !is.na(complement_a0) &
    complement_a1 == ld_variants$A1 & complement_a0 == ld_variants$A0
  complement_reversed <- !is.na(complement_a1) & !is.na(complement_a0) &
    complement_a1 == ld_variants$A0 & complement_a0 == ld_variants$A1
  allele_compatible <- direct | reversed | complement_direct |
    complement_reversed
  compatible <- location_match & !palindromic & allele_compatible
  retained_ld <- candidate[compatible]
  retained_gwas <- input[compatible, , drop = FALSE]
  orientation <- ifelse(
    direct[compatible] | complement_direct[compatible], 1, -1
  )
  gwas_unmatched <- sum(!gwas$ID %in% ld$variants$ID)
  ld_unmatched <- sum(!id_match)
  location_mismatch <- sum(!location_match)
  ambiguous <- sum(location_match & palindromic)
  allele_mismatch <- sum(location_match & !palindromic & !allele_compatible)
  if (require_complete_ld && length(retained_ld) != nrow(ld$variants)) {
    stop(sprintf(
      paste0(
        "The GWAS statistics do not cover the complete eigen-LD variant ",
        "panel (LD-only %d; location mismatch %d; allele mismatch %d; ",
        "ambiguous %d). Match the GWAS and native `blm_ld` object with ",
        "`match_gwas_ld()` before calling `as_blm_ld_eigen()`, or impute ",
        "the missing summary statistics."
      ),
      ld_unmatched, location_mismatch, allele_mismatch, ambiguous
    ), call. = FALSE)
  }
  if (!length(retained_ld)) {
    stop("No GWAS variants remain after position and allele harmonization.",
         call. = FALSE)
  }
  subset_ld <- if (identical(retained_ld, seq_len(nrow(ld$variants)))) {
    ld
  } else {
    .subset_blm_ld(ld, retained_ld)
  }
  excluded <- gwas_unmatched + ld_unmatched +
    2L * (location_mismatch + ambiguous + allele_mismatch)
  if (excluded > 0L) {
    warning(sprintf(
      paste0(
        "GWAS-LD harmonization retained %d variants and excluded %d entries ",
        "(GWAS-only %d; LD-only %d; location mismatch %d; allele mismatch ",
        "%d; ambiguous %d)."
      ),
      length(retained_ld), excluded, gwas_unmatched, ld_unmatched,
      location_mismatch, allele_mismatch, ambiguous
    ), call. = FALSE)
  }
  list(
    gwas = retained_gwas,
    ld = subset_ld,
    orientation = orientation,
    counts = c(
      retained = length(retained_ld),
      flipped = sum(orientation < 0),
      excluded = excluded,
      gwas_only = gwas_unmatched,
      ld_only = ld_unmatched,
      location_mismatch = location_mismatch,
      allele_mismatch = allele_mismatch,
      ambiguous = ambiguous
    )
  )
}

.subset_blm_ld <- function(ld, selected) {
  selected_flag <- logical(nrow(ld$variants))
  selected_flag[selected] <- TRUE
  new_blocks <- list()
  source_blocks <- character()
  complete_blocks <- logical()
  reserved_names <- vapply(ld$blocks, `[[`, character(1), "name")
  offset <- 0L
  parent_counts <- integer()
  for (block in ld$blocks) {
    global <- seq.int(offset + 1L, offset + block$size)
    keep <- which(selected_flag[global])
    offset <- offset + block$size
    if (!length(keep)) next
    if (length(keep) == block$size) {
      output_index <- length(new_blocks) + 1L
      new_blocks[[output_index]] <- block
      names(new_blocks)[output_index] <- block$name
      source_blocks[block$name] <- block$name
      complete_blocks[block$name] <- TRUE
      next
    }
    local_map <- integer(block$size)
    local_map[keep] <- seq_along(keep)
    triplets <- .ld_block_triplets(block)
    edge_keep <- local_map[triplets$column] > 0L &
      local_map[triplets$row] > 0L
    rows <- local_map[triplets$row[edge_keep]]
    columns <- local_map[triplets$column[edge_keep]]
    values <- triplets$value[edge_keep]
    size <- length(keep)
    ranges <- .exact_contiguous_ld_triplet_blocks(
      size, rows, columns, values
    )
    for (range in ranges) {
      parent <- block$parent
      repeat {
        count <- if (parent %in% names(parent_counts)) {
          parent_counts[[parent]] + 1L
        } else {
          1L
        }
        parent_counts[parent] <- count
        child_name <- paste0(parent, ".", count)
        if (!child_name %in% c(reserved_names, names(new_blocks))) break
      }
      first <- range[[1L]]
      last <- range[[length(range)]]
      in_range <- rows >= first & rows <= last &
        columns >= first & columns <= last
      child <- .compress_ld_triplets(
        length(range), rows[in_range] - first + 1L,
        columns[in_range] - first + 1L, values[in_range], parent, child_name
      )
      new_blocks[[child$name]] <- child
      source_blocks[child$name] <- block$name
      complete_blocks[child$name] <- FALSE
    }
  }
  if (!length(new_blocks)) stop("LD subset is empty.", call. = FALSE)
  variants <- ld$variants[selected_flag, , drop = FALSE]
  rownames(variants) <- NULL
  block_table <- .ld_block_table(new_blocks)
  result <- structure(list(
    blocks = new_blocks,
    variants = variants,
    parents = unique(block_table$parent),
    block_table = block_table,
    format_version = .blm_ld_format_version,
    cross_block_assumption = if (length(new_blocks) > 1L) "zero" else NULL
  ), class = "blm_ld")
  if (!is.null(ld$regularization_report)) {
    input_report <- ld$regularization_report
    report_rows <- match(
      unname(source_blocks[names(new_blocks)]), input_report$block
    )
    if (anyNA(report_rows)) {
      stop("`ld` regularization metadata are inconsistent.", call. = FALSE)
    }
    report <- input_report[report_rows, , drop = FALSE]
    if (!"source_block" %in% names(report)) {
      report$source_block <- report$block
    }
    report$block <- names(new_blocks)
    report$parent <- vapply(new_blocks, `[[`, character(1), "parent")
    report$predictors <- vapply(new_blocks, `[[`, integer(1), "size")
    incomplete <- !unname(complete_blocks[names(new_blocks)])
    numerical_fields <- intersect(
      c(
        "minimum_eigenvalue_before", "minimum_eigenvalue_after",
        "positive_definite_after"
      ),
      names(report)
    )
    report[incomplete, numerical_fields] <- NA
    rownames(report) <- NULL
    result$regularization_report <- report
  }
  result
}

.ld_regularization_block_map <- function(source_report, fitted_report) {
  fitted_block <- as.character(fitted_report$block)
  source_row <- match(fitted_block, source_report$block)
  unmatched <- is.na(source_row)
  if (any(unmatched)) {
    source_provenance <- if ("source_block" %in% names(fitted_report)) {
      as.character(fitted_report$source_block)
    } else {
      fitted_block
    }
    source_row[unmatched] <- match(
      source_provenance[unmatched], source_report$block
    )
  }
  if (anyNA(source_row)) {
    stop("`ld` regularization metadata are inconsistent.", call. = FALSE)
  }
  source_block <- as.character(source_report$block[source_row])
  source_predictors <- as.integer(source_report$predictors[source_row])
  fitted_predictors <- as.integer(fitted_report$predictors)
  data.frame(
    fitted_block = fitted_block,
    source_block = source_block,
    fitted_predictors = fitted_predictors,
    source_predictors = source_predictors,
    subset = fitted_block != source_block |
      fitted_predictors != source_predictors,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.orient_gwas_coefficients <- function(result, orientation, source_indices,
                                      store_samples,
                                      store_coefficient_cov) {
  for (block_index in seq_along(result$ETA)) {
    signs <- orientation[source_indices[[block_index]]]
    block <- result$ETA[[block_index]]
    block$coefficient_mean <- block$coefficient_mean * signs
    if (store_samples) {
      block$coefficient_samples <- sweep(
        block$coefficient_samples, 2L, signs, FUN = "*"
      )
    }
    if (store_coefficient_cov) {
      block$coefficient_cov <- block$coefficient_cov * tcrossprod(signs)
    }
    result$ETA[[block_index]] <- block
  }
  result
}
