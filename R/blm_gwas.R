#' Bayesian linear regression from GWAS summary statistics and LD
#'
#' Fits the same coefficient-prior models as [blm_ss()] from additive
#' quantitative-trait GWAS results and a reusable LD object created by
#' [as_blm_ld()]. The working Gram matrix is applied in LD-native form and is
#' not materialized.
#'
#' @param gwas A data frame with columns `CHR`, `ID`, `POS`, `A1`, `A0`, `N`,
#'   `BETA`, and `SE`. `BETA` is the additive marginal effect per `A1` allele.
#'   All retained variants must have a common sample size `N`.
#' @param ld A `blm_ld` object returned by the current [as_blm_ld()].
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
#' @param check_psd Whether to perform the optional full PSD and joint
#'   compatibility validation; see [blm_ss()].
#' @param compute_pve,pve_type Posterior PVE controls; see [blm_ss()].
#' @param scale Working scale. `"auto"` uses the original response scale when
#'   `reference_response_var` is supplied and otherwise standardizes the
#'   response and predictors. `"original"` requires
#'   `reference_response_var`.
#' @param residual_df_gwas Positive residual degrees of freedom used by the
#'   marginal GWAS regressions. The default is `N - 2`.
#' @param ld_shrink Numeric scalar in `[0, 1)`. Off-diagonal LD correlations
#'   are multiplied by `1 - ld_shrink` while the unit diagonal is retained.
#'   This is applied by the native LD operator without copying or modifying
#'   `ld`. Positive values can stabilize analyses using external LD.
#'
#' @return An object of class `blm_fit`. Coefficients are oriented to the input
#'   GWAS `A1` alleles. No intercept is fitted because centered GWAS summary
#'   statistics do not identify the phenotype mean. Nevertheless, their
#'   centered likelihood has `likelihood_df = N - 1`; this value is recorded
#'   in the result and is separate from `residual_df_gwas`. When `ld` was created by
#'   [regularize_blm_ld()], `ld_regularization_report` preserves that object's
#'   original regularization audit report unchanged, while
#'   `ld_regularization_block_map` maps fitted post-harmonization blocks to the
#'   source blocks in that report. `ld_harmonization` reports
#'   retained and flipped variants together with separate counts for GWAS-only,
#'   LD-only, location-mismatched, allele-mismatched, and ambiguous variants.
#'   Its `excluded` element counts excluded table entries: unmatched entries
#'   contribute one and matched-but-incompatible variant pairs contribute two.
#'
#' @details Variants are matched by `ID` and checked against chromosome,
#'   position, and alleles. Reversed alleles are handled by changing effect
#'   orientation. Unmatched, incompatible, and unresolved strand-ambiguous
#'   variants are excluded with a warning. Excluded IDs are also removed from
#'   character-indexed `ETA` blocks, with a block-specific warning; a block
#'   that becomes empty is rejected. If harmonization excludes an LD variant,
#'   numeric `ETA` indices are rejected because their intended position is
#'   ambiguous; use character variant IDs instead. Gibbs coordinates retain
#'   LD order, while `ETA` blocks remain independent of LD blocks.
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
  .validate_blm_ld_object(ld)
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
  harmonized <- .harmonize_gwas_ld(gwas, ld)
  gwas <- harmonized$gwas
  orientation <- harmonized$orientation
  ETA <- .harmonize_gwas_eta(
    ETA, gwas$ID, input_gwas_ids, ld$variants$ID
  )
  ld <- harmonized$ld
  n_values <- unique(gwas$N)
  if (length(n_values) != 1L) {
    stop("Retained GWAS variants must have one common value of `N`.",
         call. = FALSE)
  }
  n <- n_values[[1L]]
  likelihood_df <- .resolve_likelihood_df(n, FALSE, n - 1)
  if (is.null(residual_df_gwas)) residual_df_gwas <- n - 2
  if (!is.numeric(residual_df_gwas) || length(residual_df_gwas) != 1L ||
      is.na(residual_df_gwas) || !is.finite(residual_df_gwas) ||
      residual_df_gwas <= 0 || residual_df_gwas > n - 1) {
    stop(
      "`residual_df_gwas` must be positive and no greater than `N - 1`.",
      call. = FALSE
    )
  }

  working_reference_var <- if (scale == "standardized") {
    1
  } else {
    reference_response_var
  }
  components <- .gwas_working_components(
    gwas$BETA * orientation, gwas$SE, n, residual_df_gwas,
    working_reference_var, scale
  )
  predictor_names <- gwas$ID
  names(components$Xty) <- predictor_names
  normalized <- .normalize_ss_eta(ETA, predictor_names, residual_var, n)
  blocks <- normalized$blocks
  source_indices <- normalized$source_indices
  has_expected_pve <- any(vapply(
    blocks, function(block) !is.null(block$expected_pve), logical(1)
  ))
  predictor_scales <- lapply(seq_along(blocks), function(block_index) {
    indices <- source_indices[[block_index]]
    if (!blocks[[block_index]]$standardize) return(rep(1, length(indices)))
    sqrt(components$diagonal[indices] / (n - 1))
  })
  for (block_index in seq_along(blocks)) {
    blocks[[block_index]]$predictor_scale <- predictor_scales[[block_index]]
  }
  if (has_expected_pve) {
    predictor_variance_sums <- vapply(seq_along(blocks), function(block_index) {
      indices <- source_indices[[block_index]]
      sum(
        components$diagonal[indices] /
          predictor_scales[[block_index]]^2 / (n - 1)
      )
    }, numeric(1))
    blocks <- .calibrate_eta_priors(
      blocks, predictor_variance_sums, components$reference_response_var, n
    )
  }
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
  working_Xty <- components$Xty / sampler_predictor_scale
  names(working_Xty) <- sampler_internal_names

  fixed_blocks <- vapply(
    blocks, function(block) block$model == "Fixed", logical(1)
  )
  fixed_source <- unlist(source_indices[fixed_blocks], use.names = FALSE)
  if (length(fixed_source)) {
    fixed_R <- .materialize_blm_ld(ld, fixed_source, ld_shrink)
    fixed_scale <- sampler_scale[fixed_source]
    fixed_gram <- fixed_R * tcrossprod(fixed_scale)
    .validate_fixed_gram(
      fixed_gram, seq_along(fixed_source),
      sampler_internal_names[fixed_source]
    )
  }
  if (check_psd) {
    source_R <- .materialize_blm_ld(ld, ld_shrink = ld_shrink)
    validation_XtX <- source_R * tcrossprod(sampler_scale)
    .validate_working_crossproducts(
      validation_XtX, working_Xty, components$yty
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
    effective_n = n,
    likelihood_df = likelihood_df,
    fit_intercept = FALSE,
    intercept_x_mean = numeric(length(working_Xty)),
    intercept_y_mean = 0
  )
  sampler_arguments$ld_blocks <- lapply(ld$blocks, function(block) {
    block[c("type", "size", "data", "indptr", "row_index")]
  })
  block_ends <- cumsum(vapply(ld$blocks, `[[`, integer(1), "size"))
  block_starts <- block_ends - vapply(ld$blocks, `[[`, integer(1), "size") + 1L
  sampler_arguments$ld_indices <- Map(seq.int, block_starts, block_ends)
  sampler_arguments$ld_scale <- sampler_scale
  sampler_arguments$ld_shrink <- ld_shrink
  sampler_arguments$Xty <- working_Xty
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
    likelihood_df = likelihood_df,
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
  result$ld_block_table <- ld$block_table
  result$ld_cross_block_assumption <- ld$cross_block_assumption
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
  for (name in c("POS", "N")) {
    value <- result[[name]]
    if (!is.numeric(value) || anyNA(value) || any(!is.finite(value)) ||
        any(value != floor(value)) || any(value < if (name == "N") 3 else 1)) {
      stop(sprintf("`gwas$%s` must contain valid positive integers.", name),
           call. = FALSE)
    }
    result[[name]] <- as.numeric(value)
  }
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

.harmonize_gwas_ld <- function(gwas, ld) {
  match_index <- match(ld$variants$ID, gwas$ID)
  id_match <- !is.na(match_index)
  candidate <- which(id_match)
  if (!length(candidate)) {
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
  if (!length(retained_ld)) {
    stop("No GWAS variants remain after position and allele harmonization.",
         call. = FALSE)
  }
  retained_gwas <- input[compatible, , drop = FALSE]
  orientation <- ifelse(
    direct[compatible] | complement_direct[compatible], 1, -1
  )
  subset_ld <- if (identical(retained_ld, seq_len(nrow(ld$variants)))) {
    ld
  } else {
    .subset_blm_ld(ld, retained_ld)
  }
  gwas_unmatched <- sum(!gwas$ID %in% ld$variants$ID)
  ld_unmatched <- sum(!id_match)
  location_mismatch <- sum(!location_match)
  ambiguous <- sum(location_match & palindromic)
  allele_mismatch <- sum(location_match & !palindromic & !allele_compatible)
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
  source_block <- if ("source_block" %in% names(fitted_report)) {
    fitted_report$source_block
  } else {
    fitted_report$block
  }
  source_row <- match(source_block, source_report$block)
  if (anyNA(source_row)) {
    stop("`ld` regularization metadata are inconsistent.", call. = FALSE)
  }
  source_predictors <- as.integer(source_report$predictors[source_row])
  fitted_predictors <- as.integer(fitted_report$predictors)
  fitted_block <- as.character(fitted_report$block)
  data.frame(
    fitted_block = fitted_block,
    source_block = as.character(source_block),
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
