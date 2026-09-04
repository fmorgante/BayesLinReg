library(BayesLinReg)

reconstruct_eigen_ld <- function(object) {
  matrices <- lapply(object$blocks, function(block) {
    tcrossprod(
      sweep(block$eigenvectors, 2L, sqrt(block$eigenvalues), `*`)
    )
  })
  as.matrix(Matrix::bdiag(matrices))
}

ids <- paste0("ev", seq_len(7L))
R1 <- matrix(c(
  1, 0.3, 0.1,
  0.3, 1, 0.2,
  0.1, 0.2, 1
), 3L)
R2 <- matrix(c(
  1, 0.2, 0.1, 0.05,
  0.2, 1, 0.25, 0.1,
  0.1, 0.25, 1, 0.3,
  0.05, 0.1, 0.3, 1
), 4L)
dimnames(R1) <- list(ids[1:3], ids[1:3])
dimnames(R2) <- list(ids[4:7], ids[4:7])
variants1 <- data.frame(
  CHR = 1, ID = ids[1:3], POS = 1:3,
  A1 = c("A", "A", "C"), A0 = c("C", "G", "T")
)
variants2 <- data.frame(
  CHR = 2, ID = ids[4:7], POS = 4:7,
  A1 = c("A", "C", "A", "G"), A0 = c("C", "A", "G", "A")
)
explicit <- as_blm_ld(
  list(chr1 = R1, chr2 = R2),
  list(chr1 = variants1, chr2 = variants2)
)
eigen_ld <- as_blm_ld_eigen(explicit, prop_var = 1)
stopifnot(
  inherits(eigen_ld, "blm_ld_eigen"),
  identical(eigen_ld$format_version, 3L),
  all(vapply(
    eigen_ld$blocks, `[[`, logical(1), "source_eigenspace_complete"
  )),
  identical(eigen_ld$parents, c("chr1", "chr2")),
  identical(eigen_ld$block_table$rank, c(3L, 4L)),
  isTRUE(all.equal(
    reconstruct_eigen_ld(eigen_ld),
    as.matrix(Matrix::bdiag(R1, R2)), tolerance = 1e-12,
    check.attributes = FALSE
  )),
  identical(as_blm_ld_eigen(eigen_ld), eigen_ld),
  inherits(
    BayesLinReg:::.validate_blm_ld_eigen_object(eigen_ld),
    "blm_ld_eigen"
  )
)

# Existing eigenpairs and correlation matrices produce the same object-level
# representation. Independently named batches combine in supplied order.
decomposition1 <- eigen(R1, symmetric = TRUE)
decomposition2 <- eigen(R2, symmetric = TRUE)
from_pairs <- as_blm_ld_eigen(
  variants = list(chr1 = variants1, chr2 = variants2),
  eigenvectors = list(
    chr1 = decomposition1$vectors, chr2 = decomposition2$vectors
  ),
  eigenvalues = list(chr1 = decomposition1$values,
                     chr2 = decomposition2$values),
  prop_var = 1, check_eigenvectors = TRUE
)
batch1 <- as_blm_ld_eigen(
  list(chr1 = R1), list(chr1 = variants1), prop_var = 1
)
batch2 <- as_blm_ld_eigen(
  list(chr2 = R2), list(chr2 = variants2), prop_var = 1
)
combined <- combine_blm_ld_eigen(batch1, batch2)
stopifnot(
  isTRUE(all.equal(
    reconstruct_eigen_ld(from_pairs), reconstruct_eigen_ld(eigen_ld),
    tolerance = 1e-12
  )),
  all(vapply(
    from_pairs$blocks, `[[`, logical(1), "source_eigenspace_complete"
  )),
  identical(combined$variants, eigen_ld$variants),
  identical(names(combined$blocks), names(eigen_ld$blocks)),
  identical(combine_blm_ld_eigen(list(batch1, batch2)), combined),
  identical(combine_blm_ld_eigen(batch1), batch1)
)

regularized_batch1 <- regularize_blm_ld(
  as_blm_ld(list(chr1r = R1), list(chr1r = variants1)),
  method = "shrink", shrink = 0.01
)
regularized_eigen <- as_blm_ld_eigen(regularized_batch1, prop_var = 1)
mixed_eigen <- combine_blm_ld_eigen(regularized_eigen, batch2)
stopifnot(
  identical(
    regularized_eigen$regularization_report,
    regularized_batch1$regularization_report
  ),
  identical(mixed_eigen$regularization_report$method, c("shrink", "none")),
  inherits(
    BayesLinReg:::.validate_blm_ld_eigen_object(mixed_eigen),
    "blm_ld_eigen"
  )
)

# Pure truncation retains the leading eigenspace and records the realized
# trace fraction without adding a diagonal correction.
equicorrelation <- matrix(0.8, 6L, 6L)
diag(equicorrelation) <- 1
equicorrelation_ids <- paste0("eq", seq_len(6L))
dimnames(equicorrelation) <- list(equicorrelation_ids, equicorrelation_ids)
equicorrelation_variants <- data.frame(
  CHR = 3, ID = equicorrelation_ids, POS = seq_len(6L),
  A1 = rep("A", 6L), A0 = rep("C", 6L)
)
truncated <- as_blm_ld_eigen(
  list(chr3 = equicorrelation),
  list(chr3 = equicorrelation_variants), prop_var = 0.8
)
stopifnot(
  identical(truncated$blocks[[1L]]$rank, 1L),
  truncated$blocks[[1L]]$source_eigenspace_complete,
  !truncated$blocks[[1L]]$complete_eigenspace,
  truncated$blocks[[1L]]$prop_var >= 0.8,
  any(abs(diag(reconstruct_eigen_ld(truncated)) - 1) > 1e-6)
)

# Indefinite LD remains an error by default. The explicit discard policy uses
# the retained positive eigenspace directly, with prop_var relative to the
# positive trace and without a unit-diagonal correction.
indefinite <- matrix(c(
  1, 0.9, 0.9,
  0.9, 1, -0.9,
  0.9, -0.9, 1
), 3L)
indefinite_ids <- paste0("ind", seq_len(3L))
dimnames(indefinite) <- list(indefinite_ids, indefinite_ids)
indefinite_variants <- data.frame(
  CHR = 4, ID = indefinite_ids, POS = seq_len(3L),
  A1 = rep("A", 3L), A0 = rep("C", 3L)
)
indefinite_error <- try(as_blm_ld_eigen(
  indefinite, indefinite_variants, prop_var = 1
), silent = TRUE)
discarded <- as_blm_ld_eigen(
  indefinite, indefinite_variants, prop_var = 1,
  negative_eigenvalues = "discard"
)
indefinite_decomposition <- eigen(indefinite, symmetric = TRUE)
positive <- indefinite_decomposition$values > 0
positive_part <- tcrossprod(sweep(
  indefinite_decomposition$vectors[, positive, drop = FALSE], 2L,
  sqrt(indefinite_decomposition$values[positive]), `*`
))
stopifnot(
  inherits(indefinite_error, "try-error"),
  grepl("materially negative eigenvalue", indefinite_error),
  identical(discarded$blocks[[1L]]$rank, 2L),
  identical(
    discarded$blocks[[1L]]$negative_eigenvalues, "discard"
  ),
  identical(discarded$blocks[[1L]]$discarded_negative_eigenvalues, 1L),
  isTRUE(all.equal(
    discarded$blocks[[1L]]$minimum_source_eigenvalue, -0.8,
    tolerance = 1e-12
  )),
  discarded$blocks[[1L]]$complete_eigenspace,
  discarded$blocks[[1L]]$source_eigenspace_complete,
  isTRUE(all.equal(discarded$blocks[[1L]]$prop_var, 1)),
  isTRUE(all.equal(
    reconstruct_eigen_ld(discarded), positive_part, tolerance = 1e-12,
    check.attributes = FALSE
  )),
  any(abs(diag(reconstruct_eigen_ld(discarded)) - 1) > 1e-6),
  inherits(
    BayesLinReg:::.validate_blm_ld_eigen_object(discarded),
    "blm_ld_eigen"
  )
)

discard_gwas <- transform(
  indefinite_variants, N = 100, BETA = c(0.1, -0.05, 0.2), SE = 0.1
)
discard_fit <- blm_gwas(
  discard_gwas, discarded, list(model = "Normal"), residual_var = 1,
  iterations = 10L, burnin = 5L
)
stopifnot(
  discard_fit$ld_approximate,
  isTRUE(all.equal(discard_fit$ld_prop_var, 1))
)

discarded_half <- as_blm_ld_eigen(
  variants = indefinite_variants,
  eigenvectors = indefinite_decomposition$vectors,
  eigenvalues = indefinite_decomposition$values,
  prop_var = 0.5, check_eigenvectors = TRUE,
  negative_eigenvalues = "discard"
)
stopifnot(
  identical(discarded_half$blocks[[1L]]$rank, 1L),
  isTRUE(all.equal(discarded_half$blocks[[1L]]$prop_var, 0.5)),
  !discarded_half$blocks[[1L]]$complete_eigenspace
)

# Already-truncated supplied eigenpairs are not marked complete merely because
# every supplied component is retained. The approximate likelihood projects
# GWAS cross-products rather than applying the exact-eigenspace guard.
truncated_vectors <- qr.Q(qr(cbind(c(1, -1, 0), c(1, 1, -2))))
truncated_pair_variants <- data.frame(
  CHR = 5, ID = paste0("pre", seq_len(3L)), POS = seq_len(3L),
  A1 = c("A", "A", "C"), A0 = c("C", "G", "T")
)
truncated_pairs <- as_blm_ld_eigen(
  variants = truncated_pair_variants,
  eigenvectors = truncated_vectors,
  eigenvalues = c(1.495, 1.495), prop_var = 0.995,
  check_eigenvectors = TRUE
)
truncated_pairs_discard <- as_blm_ld_eigen(
  variants = truncated_pair_variants,
  eigenvectors = truncated_vectors,
  eigenvalues = c(1.495, 1.495), prop_var = 1,
  check_eigenvectors = TRUE, negative_eigenvalues = "discard"
)
truncated_pair_gwas <- transform(
  truncated_pair_variants, N = 200, BETA = c(0.1, -0.05, 0.2), SE = 0.08
)
truncated_pair_fit <- blm_gwas(
  truncated_pair_gwas, truncated_pairs, list(model = "Normal"),
  residual_var = 1, iterations = 10L, burnin = 5L
)
stopifnot(
  identical(truncated_pairs$blocks[[1L]]$rank, 2L),
  !truncated_pairs$blocks[[1L]]$source_eigenspace_complete,
  !truncated_pairs$blocks[[1L]]$complete_eigenspace,
  !truncated_pairs_discard$blocks[[1L]]$source_eigenspace_complete,
  !truncated_pairs_discard$blocks[[1L]]$complete_eigenspace,
  truncated_pair_fit$ld_approximate
)

# The native finite scan rejects non-finite eigenvector storage without an
# equally sized logical temporary.
nonfinite_eigen <- eigen_ld
nonfinite_eigen$blocks[[1L]]$eigenvectors[1L, 1L] <- Inf
nonfinite_error <- try(
  BayesLinReg:::.validate_blm_ld_eigen_object(nonfinite_eigen), silent = TRUE
)
stopifnot(
  inherits(nonfinite_error, "try-error"),
  grepl("invalid eigen block", nonfinite_error)
)
old_format_eigen <- eigen_ld
old_format_eigen$format_version <- 2L
stopifnot(inherits(
  try(BayesLinReg:::.validate_blm_ld_eigen_object(old_format_eigen),
      silent = TRUE),
  "try-error"
))

# Exact eigen LD and explicit LD give the same Markov chain up to floating-
# point roundoff, including PVE calculations.
n <- 200
gwas <- transform(
  rbind(variants1, variants2), N = n,
  BETA = c(0.1, -0.05, 0.2, 0.03, -0.1, 0.15, 0.04),
  SE = c(0.05, 0.07, 0.06, 0.08, 0.06, 0.05, 0.09)
)
common <- list(
  gwas = gwas, ETA = list(model = "Normal"), residual_var = 1,
  iterations = 70L, burnin = 20L, store_samples = TRUE,
  store_coefficient_cov = TRUE, compute_pve = TRUE
)
set.seed(301)
explicit_fit <- do.call(blm_gwas, c(list(ld = explicit), common))
set.seed(301)
eigen_fit <- do.call(blm_gwas, c(list(ld = eigen_ld), common))
stopifnot(
  isTRUE(all.equal(explicit_fit$ETA, eigen_fit$ETA, tolerance = 1e-12)),
  isTRUE(all.equal(
    explicit_fit$total_pve_samples, eigen_fit$total_pve_samples,
    tolerance = 1e-12
  )),
  identical(eigen_fit$ld_representation, "truncated_eigen"),
  identical(eigen_fit$ld_eigen_rank, 7L),
  isTRUE(all.equal(eigen_fit$ld_prop_var, 1, tolerance = 1e-12)),
  !eigen_fit$ld_approximate
)

# Eigen LD permits reordered rows and irrelevant GWAS-only variants while
# retaining the complete LD panel.
extra_gwas <- rbind(
  gwas[rev(seq_len(nrow(gwas))), ],
  transform(gwas[1L, ], ID = "gwas_extra", POS = 999L)
)
extra_fit <- suppressWarnings(blm_gwas(
  extra_gwas, eigen_ld, list(model = "Normal"), residual_var = 1,
  iterations = 20L, burnin = 10L
))
stopifnot(identical(extra_fit$gwas_variants$ID, eigen_ld$variants$ID))

# Native-LD and eigen blocks use the same block-indexed RNG streams and
# parallel coordinate/local-variance workers for every prior family. Check
# backend parity as well as invariance to the requested thread count. The two
# LD blocks ensure that the parallel dispatch is actually exercised.
parallel_prior_models <- c(
  "Normal", "SpikeSlab", "GlobalLocal", "SpikeMultiSlab", "Fixed"
)
parallel_common <- common
parallel_common$residual_var <- NULL
parallel_common$residual_shape <- 2
parallel_common$residual_scale <- 1
parallel_common$iterations <- 55L
parallel_common$burnin <- 20L
parallel_common$pve_type <- "allocated"
compare_parallel_chain <- function(first, second, tolerance = 1e-10) {
  stopifnot(
    isTRUE(all.equal(first$ETA, second$ETA, tolerance = tolerance)),
    isTRUE(all.equal(
      first$residual_var_samples, second$residual_var_samples,
      tolerance = tolerance
    )),
    isTRUE(all.equal(
      first$total_pve_samples, second$total_pve_samples,
      tolerance = tolerance
    )),
    isTRUE(all.equal(
      first$cross_block_pve_samples, second$cross_block_pve_samples,
      tolerance = tolerance
    ))
  )
}
for (model_index in seq_along(parallel_prior_models)) {
  parallel_common$ETA <- list(model = parallel_prior_models[[model_index]])
  fit_backend <- function(ld_input, threads) {
    set.seed(310 + model_index)
    do.call(
      blm_gwas,
      c(list(ld = ld_input, nthreads = threads), parallel_common)
    )
  }
  explicit_two <- fit_backend(explicit, 2L)
  explicit_three <- fit_backend(explicit, 3L)
  eigen_two <- fit_backend(eigen_ld, 2L)
  eigen_three <- fit_backend(eigen_ld, 3L)
  compare_parallel_chain(explicit_two, explicit_three)
  compare_parallel_chain(eigen_two, eigen_three)
  compare_parallel_chain(explicit_two, eigen_two)
}

# Predictor-specific original-scale factors use the general low-rank
# projection and remain equivalent for a full eigen representation.
original_common <- common
original_common$ETA <- list(model = "Normal", standardize = FALSE)
original_common$scale <- "original"
original_common$reference_response_var <- 2
original_common$compute_pve <- FALSE
set.seed(302)
explicit_original <- do.call(
  blm_gwas, c(list(ld = explicit), original_common)
)
set.seed(302)
eigen_original <- do.call(
  blm_gwas, c(list(ld = eigen_ld), original_common)
)
stopifnot(isTRUE(all.equal(
  explicit_original$ETA, eigen_original$ETA, tolerance = 1e-12
)))

# Eigen LD cannot be subset during fitting. Match the native LD and GWAS first,
# then decompose the final shared panel.
subset_gwas <- gwas[-2L, ]
subset_error <- try(suppressWarnings(blm_gwas(
  subset_gwas, eigen_ld, list(model = "Normal"), residual_var = 1,
  iterations = 20L, burnin = 10L
)), silent = TRUE)
matched_subset <- suppressWarnings(match_gwas_ld(subset_gwas, explicit))
subset_object <- as_blm_ld_eigen(matched_subset$ld, prop_var = 1)
subset_fit <- blm_gwas(
  matched_subset$gwas, subset_object, list(model = "Normal"),
  residual_var = 1, iterations = 20L, burnin = 10L
)
stopifnot(
  inherits(subset_error, "try-error"),
  grepl("complete eigen-LD variant panel", subset_error),
  identical(matched_subset$retained_ids, ids[c(1L, 3:7)]),
  identical(subset_fit$gwas_variants$ID, ids[c(1L, 3:7)]),
  identical(subset_fit$ld_eigen_rank, 6L),
  isTRUE(all.equal(
    reconstruct_eigen_ld(subset_object),
    as.matrix(Matrix::bdiag(R1[c(1L, 3L), c(1L, 3L)], R2)),
    tolerance = 1e-12, check.attributes = FALSE
  ))
)

# Invalid combinations and eigenpairs fail before sampling.
bad_values <- decomposition1$values
bad_values[[length(bad_values)]] <- -1
negative_error <- try(as_blm_ld_eigen(
  variants = variants1, eigenvectors = decomposition1$vectors,
  eigenvalues = bad_values, prop_var = 1
), silent = TRUE)
shrink_error <- try(blm_gwas(
  gwas, eigen_ld, list(model = "Normal"), residual_var = 1,
  ld_shrink = 0.01, iterations = 10L, burnin = 5L
), silent = TRUE)
duplicate_parent_error <- try(
  combine_blm_ld_eigen(batch1, batch1), silent = TRUE
)
retruncate_error <- try(
  as_blm_ld_eigen(truncated, prop_var = 0.9), silent = TRUE
)
policy_change_error <- try(
  as_blm_ld_eigen(eigen_ld, negative_eigenvalues = "discard"),
  silent = TRUE
)
singular_R <- matrix(1, 3L, 3L)
dimnames(singular_R) <- list(ids[1:3], ids[1:3])
singular_ld <- as_blm_ld_eigen(
  list(singular = singular_R), list(singular = variants1), prop_var = 1
)
incompatible_singular_gwas <- transform(
  variants1, N = n, BETA = c(0.1, -0.1, 0.2), SE = 0.05
)
projection_error <- try(blm_gwas(
  incompatible_singular_gwas, singular_ld, list(model = "Normal"),
  residual_var = 1, iterations = 10L, burnin = 5L
), silent = TRUE)
stopifnot(
  inherits(negative_error, "try-error"),
  grepl("materially negative eigenvalue", negative_error),
  inherits(shrink_error, "try-error"),
  grepl("must be zero", shrink_error),
  inherits(duplicate_parent_error, "try-error"),
  inherits(retruncate_error, "try-error"),
  grepl("cannot change", retruncate_error),
  inherits(policy_change_error, "try-error"),
  grepl("cannot change", policy_change_error),
  inherits(projection_error, "try-error"),
  grepl("outside the exact eigenspace", projection_error)
)
