library(BayesLinReg)

# The constructor detects exact contiguous sub-blocks and reconstructs the
# supplied correlations without changing them.
ids <- paste0("rs", seq_len(6L))
R1 <- matrix(c(
  1, 0.2, 0,
  0.2, 1, 0.1,
  0, 0.1, 1
), 3L)
R2 <- matrix(c(
  1, -0.3, 0,
  -0.3, 1, 0.15,
  0, 0.15, 1
), 3L)
dimnames(R1) <- list(ids[1:3], ids[1:3])
dimnames(R2) <- list(ids[4:6], ids[4:6])
variants1 <- data.frame(
  CHR = 1, ID = ids[1:3], POS = 1:3,
  A1 = c("A", "A", "C"), A0 = c("C", "G", "T")
)
variants2 <- data.frame(
  CHR = 2, ID = ids[4:6], POS = 4:6,
  A1 = c("A", "C", "A"), A0 = c("C", "A", "G")
)
ld <- as_blm_ld(
  list(chr1 = R1, chr2 = R2),
  list(chr1 = variants1, chr2 = variants2)
)
stopifnot(
  inherits(ld, "blm_ld"),
  identical(ld$format_version, 1L),
  nrow(ld$block_table) == 2L,
  identical(ld$parents, c("chr1", "chr2")),
  identical(ld$block_table$variant_start, c(1L, 4L)),
  identical(ld$block_table$variant_end, c(3L, 6L)),
  all(vapply(ld$blocks, function(block) {
    !"variants" %in% names(block)
  }, logical(1))),
  isTRUE(all.equal(
    BayesLinReg:::.materialize_blm_ld(ld),
    as.matrix(Matrix::bdiag(R1, R2)),
    check.attributes = FALSE
  ))
)

# LD diagnostics are read-only. Exact repair is limited to manageable blocks,
# while shrinkage preserves the compressed representation.
indefinite_R <- matrix(c(
  1, 0.9, 0.9,
  0.9, 1, -0.9,
  0.9, -0.9, 1
), 3L, dimnames = list(ids[1:3], ids[1:3]))
indefinite_ld <- as_blm_ld(indefinite_R, variants1)
indefinite_diagnostics <- diagnose_blm_ld(indefinite_ld)
stopifnot(
  identical(indefinite_diagnostics$status, "indefinite"),
  indefinite_diagnostics$minimum_eigenvalue < 0,
  indefinite_diagnostics$minimum_ld_shrink > 0,
  identical(indefinite_ld$blocks$LD$data,
            as_blm_ld(indefinite_R, variants1)$blocks$LD$data)
)
unassessed_diagnostics <- diagnose_blm_ld(
  indefinite_ld, max_block_size = 2L
)
oversized_control_error <- try(
  diagnose_blm_ld(
    indefinite_ld, max_block_size = .Machine$integer.max + 1
  ),
  silent = TRUE
)
stopifnot(
  identical(unassessed_diagnostics$status, "not_assessed"),
  is.na(unassessed_diagnostics$minimum_eigenvalue),
  inherits(oversized_control_error, "try-error")
)
eigen_repaired_ld <- regularize_blm_ld(
  indefinite_ld, method = "eigen", eigen_floor = 1e-6
)
eigen_repaired_R <- BayesLinReg:::.materialize_blm_ld(eigen_repaired_ld)
stopifnot(
  min(eigen(eigen_repaired_R, symmetric = TRUE, only.values = TRUE)$values) >=
    1e-6,
  all(diag(eigen_repaired_R) == 1),
  identical(eigen_repaired_ld$regularization_report$method, "eigen"),
  identical(eigen_repaired_ld$regularization_report$source_block, "LD"),
  eigen_repaired_ld$regularization_report$floor_shrink > 0,
  isTRUE(eigen_repaired_ld$regularization_report$positive_definite_after)
)
corrupt_report_ld <- eigen_repaired_ld
corrupt_report_ld$regularization_report$block <- "wrong"
corrupt_report_error <- try(
  BayesLinReg:::.validate_blm_ld_object(corrupt_report_ld), silent = TRUE
)
corrupt_report_value_ld <- eigen_repaired_ld
corrupt_report_value_ld$regularization_report$shrink <- 0.1
corrupt_report_value_error <- try(
  BayesLinReg:::.validate_blm_ld_object(corrupt_report_value_ld), silent = TRUE
)
stopifnot(
  inherits(corrupt_report_error, "try-error"),
  grepl("regularization metadata", corrupt_report_error),
  inherits(corrupt_report_value_error, "try-error"),
  grepl("regularization metadata", corrupt_report_value_error)
)
legacy_report_ld <- eigen_repaired_ld
legacy_report_ld$regularization_report$source_block <- NULL
legacy_report_ld$regularization_report$floor_shrink <- NULL
legacy_report_ld$regularization_report$predictors <- as.numeric(
  legacy_report_ld$regularization_report$predictors
)
stopifnot(inherits(
  BayesLinReg:::.validate_blm_ld_object(legacy_report_ld), "blm_ld"
))
upgraded_report_ld <- regularize_blm_ld(
  legacy_report_ld, method = "shrink", shrink = 0.01
)
stopifnot(
  all(c("source_block", "floor_shrink") %in%
      names(upgraded_report_ld$regularization_report)),
  identical(upgraded_report_ld$regularization_report$source_block, "LD"),
  identical(upgraded_report_ld$regularization_report$floor_shrink, 0)
)
complete_block_map <- BayesLinReg:::.ld_regularization_block_map(
  eigen_repaired_ld$regularization_report,
  eigen_repaired_ld$regularization_report
)
legacy_block_map <- BayesLinReg:::.ld_regularization_block_map(
  legacy_report_ld$regularization_report,
  legacy_report_ld$regularization_report
)
stopifnot(
  identical(complete_block_map$fitted_block, "LD"),
  identical(complete_block_map$source_block, "LD"),
  !complete_block_map$subset,
  identical(legacy_block_map, complete_block_map)
)

# Automatic repair uses the requested final floor, including for matrices
# that are already positive definite but do not meet that floor.
near_singular_R <- matrix(c(1, 0.999, 0.999, 1), 2L)
dimnames(near_singular_R) <- list(ids[1:2], ids[1:2])
near_singular_ld <- as_blm_ld(near_singular_R, variants1[1:2, ])
auto_floor <- 0.01
auto_repaired_ld <- regularize_blm_ld(
  near_singular_ld, method = "auto", eigen_floor = auto_floor
)
auto_repaired_R <- BayesLinReg:::.materialize_blm_ld(auto_repaired_ld)
auto_unchanged_ld <- regularize_blm_ld(
  as_blm_ld(R1, variants1), method = "auto", eigen_floor = auto_floor
)
stopifnot(
  identical(auto_repaired_ld$regularization_report$method, "eigen"),
  min(eigen(auto_repaired_R, symmetric = TRUE, only.values = TRUE)$values) >=
    auto_floor,
  auto_repaired_ld$regularization_report$floor_shrink > 0,
  identical(auto_unchanged_ld$regularization_report$method, "none"),
  identical(auto_unchanged_ld$regularization_report$floor_shrink, 0)
)
shrunk_ld <- regularize_blm_ld(
  indefinite_ld, method = "shrink", shrink = 0.6
)
stopifnot(
  isTRUE(all.equal(
    BayesLinReg:::.materialize_blm_ld(shrunk_ld),
    0.4 * indefinite_R + 0.6 * diag(3),
    check.attributes = FALSE
  )),
  identical(shrunk_ld$regularization_report$shrink, 0.6),
  is.na(shrunk_ld$regularization_report$positive_definite_after)
)
oversized_repair <- try(
  regularize_blm_ld(
    indefinite_ld, method = "eigen", max_block_size = 2L
  ),
  silent = TRUE
)
stopifnot(
  inherits(oversized_repair, "try-error"),
  grepl("Dense eigen repair exceeds", oversized_repair)
)

# Dimnames with array attributes compare by their character values.
array_named_R1 <- R1
dimnames(array_named_R1) <- list(ids[1:3], array(ids[1:3]))
array_named_ld <- as_blm_ld(array_named_R1, variants1)
stopifnot(isTRUE(all.equal(
  BayesLinReg:::.materialize_blm_ld(array_named_ld),
  R1,
  check.attributes = FALSE
)))

# Small diagonal deviations are normalized before correlation bounds are
# checked, and stale or corrupted serialized LD objects fail in R.
rounded_R1 <- R1
diag(rounded_R1) <- c(1 + 5e-7, 1 - 5e-7, 1)
rounded_ld <- as_blm_ld(rounded_R1, variants1)
stopifnot(all(diag(BayesLinReg:::.materialize_blm_ld(rounded_ld)) == 1))
stale_ld <- ld
stale_ld$format_version <- NULL
corrupt_ld <- ld
corrupt_ld$blocks[[1L]]$indptr[1L] <- 1L
validation_gwas <- transform(
  rbind(variants1, variants2), N = 100L, BETA = 0, SE = 0.1
)
validation_eta <- list(model = "Normal")
stale_error <- try(
  blm_gwas(
    validation_gwas, stale_ld, validation_eta, residual_var = 1,
    iterations = 10L, burnin = 5L
  ),
  silent = TRUE
)
corrupt_error <- try(
  blm_gwas(
    validation_gwas, corrupt_ld, validation_eta, residual_var = 1,
    iterations = 10L, burnin = 5L
  ),
  silent = TRUE
)
stopifnot(
  inherits(stale_error, "try-error"),
  grepl("unsupported internal format", stale_error),
  inherits(corrupt_error, "try-error"),
  grepl("invalid compressed block", corrupt_error)
)

# Fully compatible inputs reuse the original compressed LD representation.
unchanged_harmonization <- BayesLinReg:::.harmonize_gwas_ld(
  BayesLinReg:::.validate_blm_gwas(
    transform(rbind(variants1, variants2),
              N = 100L, BETA = 0, SE = 0.1)
  ),
  ld
)
stopifnot(identical(unchanged_harmonization$ld, ld))

internally_blocked <- matrix(0, 5L, 5L)
internally_blocked[1:3, 1:3] <- R1
internally_blocked[4:5, 4:5] <- R2[1:2, 1:2]
internal_ids <- paste0("i", seq_len(5L))
dimnames(internally_blocked) <- list(internal_ids, internal_ids)
internal_variants <- data.frame(
  CHR = 1, ID = internal_ids, POS = seq_len(5L),
  A1 = c("A", "A", "C", "A", "C"),
  A0 = c("C", "G", "T", "C", "A")
)
internal_ld <- as_blm_ld(internally_blocked, internal_variants)
stopifnot(
  nrow(internal_ld$block_table) == 2L,
  identical(internal_ld$block_table$predictors, c(3L, 2L)),
  identical(internal_ld$block_table$parent, c("LD", "LD"))
)

# Generated computational names remain unique when parent names overlap them.
collision_ids <- paste0("c", seq_len(6L))
collision_first <- internally_blocked
dimnames(collision_first) <- list(collision_ids[1:5], collision_ids[1:5])
collision_second <- matrix(
  1, 1L, 1L, dimnames = list(collision_ids[6L], collision_ids[6L])
)
collision_variants <- data.frame(
  CHR = 1, ID = collision_ids, POS = seq_len(6L), A1 = "A", A0 = "C"
)
collision_ld <- as_blm_ld(
  list(a = collision_first, a.1 = collision_second),
  list(a = collision_variants[1:5, ], a.1 = collision_variants[6L, ])
)
stopifnot(
  identical(names(collision_ld$blocks), c("a.1", "a.2", "a.1.1")),
  !anyDuplicated(collision_ld$block_table$block),
  identical(
    unname(vapply(collision_ld$blocks, `[[`, character(1), "name")),
    names(collision_ld$blocks)
  )
)

# The LD-native kernel reproduces materialized sufficient-statistics sampling.
n <- 100L
beta <- stats::setNames(c(0.1, -0.2, 0.05, 0.3, 0, -0.1), ids)
se <- stats::setNames(rep(0.08, 6L), ids)
gwas <- transform(rbind(variants1, variants2), N = n, BETA = beta, SE = se)
ETA <- list(
  normal = list(
    indices = ids[1:3], model = "Normal"
  ),
  mixture = list(
    indices = ids[4:6], model = "SpikeMultiSlab"
  )
)
ss <- compute_ss_from_gwas(
  beta, se, list(chr1 = R1, chr2 = R2), n, response_var = 2
)
common <- list(
  ETA = ETA,
  residual_var = 1,
  iterations = 100L,
  burnin = 30L,
  store_samples = TRUE,
  store_coefficient_cov = TRUE,
  compute_pve = TRUE
)
set.seed(1101)
gwas_fit <- do.call(
  blm_gwas,
  c(list(gwas = gwas, ld = ld, reference_response_var = 2), common)
)
set.seed(1101)
ss_fit <- suppressWarnings(do.call(
  blm_ss,
  c(list(
    n = n, XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty,
    likelihood_df = n - 1L
  ), common)
))
stopifnot(
  isTRUE(all.equal(gwas_fit$ETA, ss_fit$ETA, tolerance = 1e-10)),
  isTRUE(all.equal(
    gwas_fit$residual_var_mean, ss_fit$residual_var_mean, tolerance = 1e-10
  )),
  isTRUE(all.equal(
    gwas_fit$total_pve_mean, ss_fit$total_pve_mean, tolerance = 1e-10
  )),
  identical(gwas_fit$gwas_scale, "original"),
  identical(gwas_fit$residual_df_gwas, n - 2),
  identical(gwas_fit$likelihood_df, n - 1L),
  identical(ss_fit$likelihood_df, n - 1L)
)

# Runtime shrinkage is mathematically identical to supplying a pre-shrunk LD
# matrix but leaves the original compressed values untouched.
runtime_shrink <- 0.25
explicit_ld <- as_blm_ld(
  list(
    chr1 = (1 - runtime_shrink) * R1 + runtime_shrink * diag(3),
    chr2 = (1 - runtime_shrink) * R2 + runtime_shrink * diag(3)
  ),
  list(chr1 = variants1, chr2 = variants2)
)
original_ld_data <- lapply(ld$blocks, `[[`, "data")
set.seed(1102)
runtime_shrunk_fit <- do.call(
  blm_gwas,
  c(list(
    gwas = gwas, ld = ld, reference_response_var = 2,
    ld_shrink = runtime_shrink
  ), common)
)
set.seed(1102)
explicit_shrunk_fit <- do.call(
  blm_gwas,
  c(list(gwas = gwas, ld = explicit_ld, reference_response_var = 2), common)
)
stopifnot(
  isTRUE(all.equal(
    runtime_shrunk_fit$ETA, explicit_shrunk_fit$ETA, tolerance = 1e-12
  )),
  isTRUE(all.equal(
    runtime_shrunk_fit$total_pve_samples,
    explicit_shrunk_fit$total_pve_samples,
    tolerance = 1e-12
  )),
  identical(runtime_shrunk_fit$ld_shrink, runtime_shrink),
  identical(explicit_shrunk_fit$ld_shrink, 0),
  identical(lapply(ld$blocks, `[[`, "data"), original_ld_data)
)

# Predictor order inside ETA blocks affects only output order. The sampler
# retains LD order, so the underlying draws and hyperparameter updates match.
reordered_ETA <- ETA
reordered_ETA$normal$indices <- rev(reordered_ETA$normal$indices)
reordered_ETA$mixture$indices <- rev(reordered_ETA$mixture$indices)
set.seed(1101)
reordered_fit <- do.call(
  blm_gwas,
  c(list(
    gwas = gwas, ld = ld, ETA = reordered_ETA,
    reference_response_var = 2
  ), common[names(common) != "ETA"])
)
for (block_name in names(gwas_fit$ETA)) {
  reference_block <- gwas_fit$ETA[[block_name]]
  reordered_block <- reordered_fit$ETA[[block_name]]
  requested <- names(reordered_block$coefficient_mean)
  stopifnot(
    identical(
      reference_block$coefficient_samples[, requested, drop = FALSE],
      reordered_block$coefficient_samples
    ),
    identical(
      reference_block$coefficient_cov[requested, requested, drop = FALSE],
      reordered_block$coefficient_cov
    )
  )
}
stopifnot(
  identical(gwas_fit$residual_var_samples, reordered_fit$residual_var_samples),
  identical(gwas_fit$total_pve_samples, reordered_fit$total_pve_samples),
  identical(
    gwas_fit$ETA$mixture$component_samples[
      , names(reordered_fit$ETA$mixture$coefficient_mean), drop = FALSE
    ],
    reordered_fit$ETA$mixture$component_samples
  )
)

# Prior blocks may cross LD blocks without changing the LD scan order. Verify
# the retained-draw PVE values directly from the materialized Gram matrix.
crossed_eta <- list(
  normal = list(
    indices = c("rs5", "rs1", "rs4"), model = "Normal"
  ),
  mixture = list(
    indices = c("rs6", "rs3", "rs2"), model = "SpikeMultiSlab"
  )
)
set.seed(1109)
crossed_fit <- blm_gwas(
  gwas, ld, crossed_eta,
  residual_var = 1, reference_response_var = 2,
  iterations = 70L, burnin = 20L,
  store_samples = TRUE, compute_pve = TRUE
)
crossed_samples <- do.call(
  cbind, lapply(crossed_fit$ETA, `[[`, "coefficient_samples")
)
crossed_samples <- crossed_samples[, ids, drop = FALSE]
materialized_XtX <- as.matrix(Matrix::bdiag(ss$XtX))
signal_variance <- rowSums(
  (crossed_samples %*% materialized_XtX) * crossed_samples
) / (n - 1L)
expected_total_pve <- signal_variance / (signal_variance + 1)
expected_standalone <- vapply(names(crossed_eta), function(block_name) {
  block_ids <- crossed_eta[[block_name]]$indices
  block_samples <- crossed_samples[, block_ids, drop = FALSE]
  block_indices <- match(block_ids, ids)
  block_signal <- rowSums(
    (block_samples %*%
       materialized_XtX[block_indices, block_indices, drop = FALSE]) *
      block_samples
  ) / (n - 1L)
  block_signal / (signal_variance + 1)
}, numeric(nrow(crossed_samples)))
stopifnot(
  isTRUE(all.equal(
    crossed_fit$total_pve_samples, expected_total_pve, tolerance = 1e-10
  )),
  isTRUE(all.equal(
    vapply(crossed_fit$ETA, `[[`, numeric(length(expected_total_pve)),
           "pve_samples"),
    expected_standalone,
    tolerance = 1e-10,
    check.attributes = FALSE
  ))
)

# Harmonization exclusions are removed from character-indexed ETA blocks.
incompatible_gwas <- gwas
incompatible_gwas$POS[incompatible_gwas$ID == "rs2"] <- 200L
regularized_ld <- regularize_blm_ld(ld, method = "eigen")
numeric_eta_error <- try(
  suppressWarnings(blm_gwas(
    incompatible_gwas, ld,
    list(
      first = list(indices = 1:3, model = "Normal"),
      second = list(indices = 4:6, model = "Normal")
    ),
    residual_var = 1, iterations = 10L, burnin = 5L
  )),
  silent = TRUE
)
stopifnot(
  inherits(numeric_eta_error, "try-error"),
  grepl("use character variant IDs", numeric_eta_error)
)
partial_harmonization <- suppressWarnings(
  BayesLinReg:::.harmonize_gwas_ld(
    BayesLinReg:::.validate_blm_gwas(incompatible_gwas), regularized_ld
  )
)
partial_report <- partial_harmonization$ld$regularization_report
legacy_regularized_ld <- regularized_ld
legacy_regularized_ld$regularization_report$source_block <- NULL
legacy_regularized_ld$regularization_report$floor_shrink <- NULL
legacy_partial_harmonization <- suppressWarnings(
  BayesLinReg:::.harmonize_gwas_ld(
    BayesLinReg:::.validate_blm_gwas(incompatible_gwas),
    legacy_regularized_ld
  )
)
stopifnot(
  identical(
    partial_harmonization$ld$blocks$chr2, regularized_ld$blocks$chr2
  ),
  nrow(partial_report) == nrow(partial_harmonization$ld$block_table),
  all(partial_report$source_block %in% c("chr1", "chr2")),
  all(is.na(partial_report$minimum_eigenvalue_after[
    partial_report$source_block == "chr1"
  ])),
  all(!is.na(partial_report$minimum_eigenvalue_after[
    partial_report$source_block == "chr2"
  ])),
  inherits(
    BayesLinReg:::.validate_blm_ld_object(legacy_partial_harmonization$ld),
    "blm_ld"
  ),
  "source_block" %in%
    names(legacy_partial_harmonization$ld$regularization_report),
  !"floor_shrink" %in%
    names(legacy_partial_harmonization$ld$regularization_report)
)
categorized_gwas <- gwas[gwas$ID != "rs6", ]
categorized_gwas$POS[categorized_gwas$ID == "rs2"] <- 200L
categorized_gwas$A1[categorized_gwas$ID == "rs3"] <- "A"
categorized_gwas$A0[categorized_gwas$ID == "rs3"] <- "C"
categorized_gwas <- rbind(
  categorized_gwas,
  transform(gwas[1L, ], ID = "gwas_only", POS = 999L)
)
categorized_harmonization <- suppressWarnings(
  BayesLinReg:::.harmonize_gwas_ld(
    BayesLinReg:::.validate_blm_gwas(categorized_gwas), ld
  )
)
stopifnot(identical(
  categorized_harmonization$counts,
  c(
    retained = 3L, flipped = 0L, excluded = 6L, gwas_only = 1L,
    ld_only = 1L, location_mismatch = 1L, allele_mismatch = 1L,
    ambiguous = 0L
  )
))
harmonization_warnings <- character()
set.seed(1104)
excluded_fit <- withCallingHandlers(
  blm_gwas(
    incompatible_gwas, regularized_ld, ETA,
    residual_var = 1, iterations = 30L, burnin = 10L
  ),
  warning = function(condition) {
    harmonization_warnings <<- c(
      harmonization_warnings, conditionMessage(condition)
    )
    invokeRestart("muffleWarning")
  }
)
stopifnot(
  identical(excluded_fit$gwas_variants$ID, ids[c(1L, 3L:6L)]),
  identical(
    excluded_fit$ld_regularization_report,
    regularized_ld$regularization_report
  ),
  identical(
    names(excluded_fit$ld_regularization_block_map),
    c(
      "fitted_block", "source_block", "fitted_predictors",
      "source_predictors", "subset"
    )
  ),
  identical(
    excluded_fit$ld_regularization_block_map$fitted_block,
    excluded_fit$ld_block_table$block
  ),
  all(
    excluded_fit$ld_regularization_block_map$fitted_predictors ==
      excluded_fit$ld_block_table$predictors
  ),
  all(
    excluded_fit$ld_regularization_block_map$source_predictors ==
      regularized_ld$regularization_report$predictors[match(
        excluded_fit$ld_regularization_block_map$source_block,
        regularized_ld$regularization_report$block
      )]
  ),
  all(excluded_fit$ld_regularization_block_map$subset[
    excluded_fit$ld_regularization_block_map$source_block == "chr1"
  ]),
  all(!excluded_fit$ld_regularization_block_map$subset[
    excluded_fit$ld_regularization_block_map$source_block == "chr2"
  ]),
  identical(
    excluded_fit$ld_harmonization,
    c(
      retained = 5L, flipped = 0L, excluded = 2L, gwas_only = 0L,
      ld_only = 0L, location_mismatch = 1L, allele_mismatch = 0L,
      ambiguous = 0L
    )
  ),
  any(grepl("normal \\(1\\)", harmonization_warnings)),
  any(grepl("location mismatch 1", harmonization_warnings)),
  !"rs2" %in% names(coef(excluded_fit))
)

empty_eta <- list(
  removed = list(indices = "rs2", model = "Normal"),
  retained = list(indices = ids[c(1L, 3L:6L)], model = "Normal")
)
empty_eta_error <- try(
  suppressWarnings(blm_gwas(
    incompatible_gwas, ld, empty_eta,
    residual_var = 1, iterations = 10L, burnin = 5L
  )),
  silent = TRUE
)
stopifnot(
  inherits(empty_eta_error, "try-error"),
  grepl("block `removed` has no predictors remaining", empty_eta_error)
)

# Allele reversal changes the returned coefficient orientation but not the
# internally fitted model.
flipped <- gwas
flipped$A1[1] <- gwas$A0[1]
flipped$A0[1] <- gwas$A1[1]
flipped$BETA[1] <- -gwas$BETA[1]
single_eta <- list(model = "Normal", standardize = FALSE)
set.seed(1102)
original_fit <- blm_gwas(
  gwas, ld, single_eta, residual_var = 1,
  iterations = 60L, burnin = 20L
)
set.seed(1102)
stored_original_fit <- blm_gwas(
  gwas, ld, single_eta, residual_var = 1,
  iterations = 60L, burnin = 20L,
  store_samples = TRUE, store_coefficient_cov = TRUE
)
set.seed(1102)
flipped_fit <- blm_gwas(
  flipped, ld, single_eta, residual_var = 1,
  iterations = 60L, burnin = 20L
)
expected_flipped <- coef(original_fit)
expected_flipped[1] <- -expected_flipped[1]
stopifnot(
  isTRUE(all.equal(
    unname(coef(flipped_fit)), unname(expected_flipped), tolerance = 1e-12
  )),
  flipped_fit$ld_harmonization[["flipped"]] == 1,
  identical(original_fit$store_samples, FALSE),
  identical(original_fit$store_coefficient_cov, FALSE),
  is.null(original_fit$ETA$ETA1$coefficient_samples),
  is.null(original_fit$ETA$ETA1$coefficient_cov),
  isTRUE(all.equal(
    coef(original_fit), coef(stored_original_fit), tolerance = 1e-12
  ))
)

# Complement-strand alleles retain the same effect orientation.
complemented <- gwas
complemented$A1[1] <- "T"
complemented$A0[1] <- "G"
set.seed(1102)
complemented_fit <- blm_gwas(
  complemented, ld, single_eta, residual_var = 1,
  iterations = 60L, burnin = 20L
)
stopifnot(
  complemented_fit$ld_harmonization[["flipped"]] == 0,
  isTRUE(all.equal(
    coef(complemented_fit), coef(original_fit), tolerance = 1e-12
  ))
)

# Strand-ambiguous variants are excluded even when their alleles match LD.
ambiguous_ids <- c("amb1", "amb2")
ambiguous_R <- diag(2L)
dimnames(ambiguous_R) <- list(ambiguous_ids, ambiguous_ids)
ambiguous_variants <- data.frame(
  CHR = 1, ID = ambiguous_ids, POS = 1:2,
  A1 = c("A", "C"), A0 = c("T", "A")
)
ambiguous_ld <- as_blm_ld(ambiguous_R, ambiguous_variants)
ambiguous_gwas <- transform(
  ambiguous_variants, N = n, BETA = c(0.1, 0.2), SE = 0.08
)
ambiguous_warnings <- character()
ambiguous_fit <- withCallingHandlers(
  blm_gwas(
    ambiguous_gwas, ambiguous_ld, single_eta,
    residual_var = 1, iterations = 20L, burnin = 5L
  ),
  warning = function(condition) {
    ambiguous_warnings <<- c(ambiguous_warnings, conditionMessage(condition))
    invokeRestart("muffleWarning")
  }
)
stopifnot(
  identical(ambiguous_fit$gwas_variants$ID, "amb2"),
  ambiguous_fit$ld_harmonization[["ambiguous"]] == 1,
  ambiguous_fit$ld_harmonization[["excluded"]] == 2,
  any(grepl("ambiguous", ambiguous_warnings))
)

# Sparse input and within-chain block parallelism are reproducible across
# repeated runs and thread counts while learning the residual variance.
sparse_ld <- as_blm_ld(
  list(
    chr1 = Matrix::forceSymmetric(Matrix::Matrix(R1, sparse = TRUE), "L"),
    chr2 = Matrix::forceSymmetric(Matrix::Matrix(R2, sparse = TRUE), "L")
  ),
  list(chr1 = variants1, chr2 = variants2)
)
set.seed(1103)
parallel_fit <- blm_gwas(
  gwas, sparse_ld, ETA,
  residual_shape = 2, residual_scale = 1,
  iterations = 60L, burnin = 20L,
  store_samples = TRUE,
  nthreads = 2L,
  check_psd = TRUE
)
set.seed(1103)
parallel_repeat <- blm_gwas(
  gwas, sparse_ld, ETA,
  residual_shape = 2, residual_scale = 1,
  iterations = 60L, burnin = 20L,
  store_samples = TRUE,
  nthreads = 2L,
  check_psd = TRUE
)
set.seed(1103)
parallel_three <- blm_gwas(
  gwas, sparse_ld, ETA,
  residual_shape = 2, residual_scale = 1,
  iterations = 60L, burnin = 20L,
  store_samples = TRUE,
  nthreads = 3L,
  check_psd = TRUE
)
stopifnot(
  identical(parallel_fit, parallel_repeat),
  identical(parallel_fit$ETA, parallel_three$ETA),
  identical(parallel_fit$residual_var_samples,
            parallel_three$residual_var_samples),
  identical(parallel_fit$nthreads, 2L),
  identical(parallel_three$nthreads, 3L),
  parallel_fit$residual_var_mean > 0
)

# Every prior family, expected-PVE calibration, fixed effects, both PVE
# definitions, and multiple chains agree with materialized sufficient
# statistics when both backends use the same Gibbs scan order.
all_ids <- paste0("all", seq_len(10L))
all_R <- outer(seq_len(10L), seq_len(10L), function(i, j) 0.25^abs(i - j))
dimnames(all_R) <- list(all_ids, all_ids)
all_variants <- data.frame(
  CHR = 4, ID = all_ids, POS = seq_len(10L),
  A1 = rep(c("A", "C"), 5L), A0 = rep(c("C", "A"), 5L)
)
all_ld <- as_blm_ld(all_R, all_variants)
all_beta <- stats::setNames(seq(-0.12, 0.15, length.out = 10L), all_ids)
all_se <- stats::setNames(rep(0.07, 10L), all_ids)
all_gwas <- transform(
  all_variants, N = n, BETA = all_beta, SE = all_se
)
all_eta <- list(
  fixed = list(indices = all_ids[1:2], model = "Fixed"),
  normal = list(
    indices = all_ids[3:4], model = "Normal", expected_pve = 0.08
  ),
  spike = list(
    indices = all_ids[5:6], model = "SpikeSlab", expected_pve = 0.08
  ),
  global = list(
    indices = all_ids[7:8], model = "GlobalLocal",
    expected_nonzero = 1, expected_pve = 0.08
  ),
  multi = list(
    indices = all_ids[9:10], model = "SpikeMultiSlab",
    expected_pve = 0.08
  )
)
all_ss <- compute_ss_from_gwas(
  all_beta, all_se, all_R, n, response_var = 2
)
for (pve_definition in c("standalone", "allocated")) {
  all_common <- list(
    ETA = all_eta, residual_shape = 2, residual_scale = 1,
    iterations = 70L, burnin = 20L,
    store_samples = TRUE, store_coefficient_cov = TRUE,
    compute_pve = TRUE, pve_type = pve_definition
  )
  set.seed(1106)
  all_gwas_fit <- do.call(
    blm_gwas,
    c(list(
      gwas = all_gwas, ld = all_ld, reference_response_var = 2
    ), all_common)
  )
  set.seed(1106)
  all_ss_fit <- suppressWarnings(do.call(
    blm_ss,
    c(list(
      n = n, XtX = all_ss$XtX, Xty = all_ss$Xty, yty = all_ss$yty,
      reference_response_var = 2, likelihood_df = n - 1L
    ), all_common)
  ))
  stopifnot(
    isTRUE(all.equal(all_gwas_fit$ETA, all_ss_fit$ETA, tolerance = 1e-10)),
    isTRUE(all.equal(
      all_gwas_fit$residual_var_samples, all_ss_fit$residual_var_samples,
      tolerance = 1e-10
    )),
    isTRUE(all.equal(
      all_gwas_fit$total_pve_samples, all_ss_fit$total_pve_samples,
      tolerance = 1e-10
    )),
    isTRUE(all.equal(
      all_gwas_fit$cross_block_pve_samples,
      all_ss_fit$cross_block_pve_samples, tolerance = 1e-10
    ))
  )
}

all_eta_reordered <- lapply(all_eta, function(specification) {
  specification$indices <- rev(specification$indices)
  specification
})
set.seed(1108)
all_online_fit <- blm_gwas(
  all_gwas, all_ld, all_eta,
  residual_var = 1, reference_response_var = 2,
  iterations = 60L, burnin = 20L,
  store_coefficient_cov = TRUE, compute_pve = TRUE
)
set.seed(1108)
all_online_reordered <- blm_gwas(
  all_gwas, all_ld, all_eta_reordered,
  residual_var = 1, reference_response_var = 2,
  iterations = 60L, burnin = 20L,
  store_coefficient_cov = TRUE, compute_pve = TRUE
)
for (block_name in names(all_online_fit$ETA)) {
  reference_block <- all_online_fit$ETA[[block_name]]
  reordered_block <- all_online_reordered$ETA[[block_name]]
  requested <- names(reordered_block$coefficient_mean)
  stopifnot(
    identical(reference_block$coefficient_mean[requested],
              reordered_block$coefficient_mean),
    identical(
      reference_block$coefficient_cov[requested, requested, drop = FALSE],
      reordered_block$coefficient_cov
    )
  )
}
stopifnot(
  identical(
    all_online_fit$ETA$spike$inclusion_probability[
      names(all_online_reordered$ETA$spike$inclusion_probability)
    ],
    all_online_reordered$ETA$spike$inclusion_probability
  ),
  identical(
    all_online_fit$ETA$global$local_var_mean[
      names(all_online_reordered$ETA$global$local_var_mean)
    ],
    all_online_reordered$ETA$global$local_var_mean
  ),
  identical(
    all_online_fit$ETA$multi$component_probability[
      rownames(all_online_reordered$ETA$multi$component_probability),
      , drop = FALSE
    ],
    all_online_reordered$ETA$multi$component_probability
  )
)

# Real multisession tests are opt-in outside environments that permit workers.
parallel_test_flags <- Sys.getenv(c(
  "BAYESLINREG_TEST_FUTURE",
  "BLM_TEST_FUTURE"
))
if (any(parallel_test_flags == "true")) {
  chain_common <- list(
    ETA = all_eta, residual_var = 1,
    iterations = 40L, burnin = 10L,
    store_samples = TRUE, nchains = 2L
  )
  set.seed(1107)
  chain_gwas_fit <- do.call(
    blm_gwas,
    c(list(
      gwas = all_gwas, ld = all_ld, reference_response_var = 2
    ), chain_common)
  )
  set.seed(1107)
  chain_ss_fit <- suppressWarnings(do.call(
    blm_ss,
    c(list(
      n = n, XtX = all_ss$XtX, Xty = all_ss$Xty, yty = all_ss$yty,
      reference_response_var = 2, likelihood_df = n - 1L
    ),
      chain_common)
  ))
  stopifnot(
    isTRUE(all.equal(chain_gwas_fit$ETA, chain_ss_fit$ETA,
                     tolerance = 1e-10)),
    isTRUE(all.equal(
      chain_gwas_fit$residual_var_samples,
      chain_ss_fit$residual_var_samples, tolerance = 1e-10
    )),
    identical(chain_gwas_fit$chain_id, chain_ss_fit$chain_id)
  )
}

# Irregular sparse lower triangles use indexed storage and match materialized
# sufficient-statistics sampling, including posterior PVE summaries.
irregular_ids <- paste0("ir", seq_len(6L))
irregular_R <- diag(6L)
irregular_R[6L, 1L] <- irregular_R[1L, 6L] <- 0.15
irregular_R[4L, 1L] <- irregular_R[1L, 4L] <- 0.05
irregular_R[5L, 2L] <- irregular_R[2L, 5L] <- -0.1
irregular_R[4L, 3L] <- irregular_R[3L, 4L] <- 0.2
dimnames(irregular_R) <- list(irregular_ids, irregular_ids)
irregular_variants <- data.frame(
  CHR = 3, ID = irregular_ids, POS = seq_len(6L),
  A1 = c("A", "C", "A", "G", "C", "T"),
  A0 = c("C", "A", "G", "A", "T", "C")
)
irregular_ld <- as_blm_ld(
  Matrix::forceSymmetric(Matrix::Matrix(irregular_R, sparse = TRUE), "L"),
  irregular_variants
)
indexed_columns <- rep.int(
  seq_len(irregular_ld$blocks[[1L]]$size),
  diff(irregular_ld$blocks[[1L]]$indptr)
)
duplicate_column <- which(tabulate(indexed_columns) > 1L)[1L]
duplicate_positions <- which(indexed_columns == duplicate_column)[1:2]
duplicate_index_ld <- irregular_ld
duplicate_index_ld$blocks[[1L]]$row_index[duplicate_positions[2L]] <-
  duplicate_index_ld$blocks[[1L]]$row_index[duplicate_positions[1L]]
unsorted_index_ld <- irregular_ld
unsorted_index_ld$blocks[[1L]]$row_index[duplicate_positions] <-
  rev(unsorted_index_ld$blocks[[1L]]$row_index[duplicate_positions])
above_diagonal_ld <- irregular_ld
above_diagonal_ld$blocks[[1L]]$row_index[1L] <-
  indexed_columns[[1L]] - 1L
indexed_validation_errors <- lapply(
  list(duplicate_index_ld, unsorted_index_ld, above_diagonal_ld),
  function(object) try(
    BayesLinReg:::.validate_blm_ld_object(object), silent = TRUE
  )
)
stopifnot(all(vapply(indexed_validation_errors, function(error) {
  inherits(error, "try-error") && grepl("strictly lower triangular", error)
}, logical(1))))
irregular_beta <- stats::setNames(
  c(0.1, -0.05, 0.2, 0, -0.1, 0.15), irregular_ids
)
irregular_se <- stats::setNames(rep(0.09, 6L), irregular_ids)
irregular_gwas <- transform(
  irregular_variants, N = n, BETA = irregular_beta, SE = irregular_se
)
irregular_ss <- compute_ss_from_gwas(
  irregular_beta, irregular_se, irregular_R, n
)
irregular_arguments <- list(
  ETA = list(model = "Normal"), residual_var = 1,
  iterations = 60L, burnin = 20L,
  store_samples = TRUE, compute_pve = TRUE
)
set.seed(1105)
irregular_fit <- do.call(
  blm_gwas,
  c(list(gwas = irregular_gwas, ld = irregular_ld), irregular_arguments)
)
set.seed(1105)
irregular_ss_fit <- suppressWarnings(do.call(
  blm_ss,
  c(list(
    n = n, XtX = irregular_ss$XtX,
    Xty = irregular_ss$Xty, yty = irregular_ss$yty,
    likelihood_df = n - 1L
  ), irregular_arguments)
))
stopifnot(
  identical(irregular_ld$block_table$storage, "indexed_triangular"),
  isTRUE(all.equal(
    irregular_fit$ETA, irregular_ss_fit$ETA, tolerance = 1e-10
  )),
  isTRUE(all.equal(
    irregular_fit$total_pve_mean,
    irregular_ss_fit$total_pve_mean,
    tolerance = 1e-10
  ))
)

# Input validation distinguishes GWAS residual degrees of freedom from the
# fitted regression residual-variance controls.
bad_n <- gwas
bad_n$N[1] <- n - 1L
stopifnot(
  inherits(try(
    blm_gwas(
      bad_n, ld, single_eta, residual_var = 1,
      iterations = 10L, burnin = 5L
    ),
    silent = TRUE
  ), "try-error"),
  inherits(try(
    blm_gwas(
      gwas, ld, single_eta, residual_var = 1,
      residual_df_gwas = n, iterations = 10L, burnin = 5L
    ),
    silent = TRUE
  ), "try-error")
)
