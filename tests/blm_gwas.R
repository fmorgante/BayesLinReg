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
  isTRUE(all.equal(
    BayesLinReg:::.materialize_blm_ld(ld),
    as.matrix(Matrix::bdiag(R1, R2)),
    check.attributes = FALSE
  ))
)

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
  flipped_fit$ld_harmonization[["flipped"]] == 1
)

# Sparse input and within-chain block parallelism are supported.
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
  residual_var = 1,
  iterations = 40L,
  burnin = 10L,
  nthreads = 2L,
  check_psd = TRUE
)
stopifnot(inherits(parallel_fit, "blm_fit"), parallel_fit$nthreads == 2L)

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
