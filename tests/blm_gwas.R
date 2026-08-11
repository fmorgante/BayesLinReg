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

# The LD-native kernel reproduces materialized sufficient-statistics sampling,
# including prior blocks that cross LD blocks.
n <- 100L
beta <- stats::setNames(c(0.1, -0.2, 0.05, 0.3, 0, -0.1), ids)
se <- stats::setNames(rep(0.08, 6L), ids)
gwas <- transform(rbind(variants1, variants2), N = n, BETA = beta, SE = se)
ETA <- list(
  normal = list(
    indices = c("rs5", "rs1", "rs4"), model = "Normal"
  ),
  mixture = list(
    indices = c("rs6", "rs3", "rs2"), model = "SpikeMultiSlab"
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
  c(list(n = n, XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty), common)
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
  identical(gwas_fit$residual_df_gwas, n - 2)
)

# Harmonization exclusions are removed from character-indexed ETA blocks.
incompatible_gwas <- gwas
incompatible_gwas$POS[incompatible_gwas$ID == "rs2"] <- 200L
partial_harmonization <- suppressWarnings(
  BayesLinReg:::.harmonize_gwas_ld(
    BayesLinReg:::.validate_blm_gwas(incompatible_gwas), ld
  )
)
stopifnot(identical(partial_harmonization$ld$blocks$chr2, ld$blocks$chr2))
harmonization_warnings <- character()
set.seed(1104)
excluded_fit <- withCallingHandlers(
  blm_gwas(
    incompatible_gwas, ld, ETA,
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
  any(grepl("mixture \\(1\\)", harmonization_warnings)),
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

# Irregular sparse lower triangles use indexed storage and match materialized
# sufficient-statistics sampling, including posterior PVE summaries.
irregular_ids <- paste0("ir", seq_len(6L))
irregular_R <- diag(6L)
irregular_R[6L, 1L] <- irregular_R[1L, 6L] <- 0.15
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
    Xty = irregular_ss$Xty, yty = irregular_ss$yty
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
