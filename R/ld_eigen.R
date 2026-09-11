#' Construct a truncated-eigen LD object
#'
#' Constructs the blockwise low-rank LD representation accepted by
#' [blm_gwas()]. Input may be a correlation matrix, a named list of correlation
#' matrices, an existing [blm_ld][as_blm_ld()] object, or existing LD
#' eigenvectors accompanied by `eigenvalues`.
#'
#' @param R A correlation matrix, a named list of correlation matrices, a
#'   `blm_ld` object, or a `blm_ld_eigen` object. It must be `NULL` when
#'   `eigenvectors` is supplied.
#' @param variants A variant metadata data frame, or a named list of data
#'   frames, with the same required columns and ordering as in [as_blm_ld()].
#'   It is required for matrix, list, and supplied-eigenpair input and must be
#'   `NULL` when `R` is an LD object.
#' @param eigenvectors `NULL`, a numeric matrix, or a named list of numeric
#'   matrices containing precomputed LD eigenvectors.
#' @param eigenvalues `NULL`, a numeric vector, or a named list of numeric
#'   vectors corresponding to `eigenvectors`.
#' @param prop_var Number in `(0, 1]`. The smallest leading set of positive
#'   eigenvalues reaching this proportion of the unit-diagonal LD trace is
#'   retained. If supplied eigenpairs do not reach the requested proportion,
#'   all supplied positive eigenpairs are retained and a warning is produced.
#'   With `negative_eigenvalues = "discard"`, the denominator is instead the
#'   sum of the available positive eigenvalues.
#'   Existing `blm_ld_eigen` input is returned unchanged; if `prop_var` is
#'   supplied explicitly, it must match the object's requested value.
#' @param check_eigenvectors Whether to check mutual orthonormality of supplied
#'   eigenvectors. Eigenvectors calculated internally are not rechecked.
#' @param negative_eigenvalues How materially negative eigenvalues are handled.
#'   `"error"` preserves strict validation. `"discard"` removes them together
#'   with numerically nonpositive components and applies `prop_var` relative to
#'   the remaining positive eigenvalue sum. This is a pure positive-eigenspace
#'   approximation: it does not restore the unit diagonal or add a diagonal
#'   correction.
#'
#' @return A `blm_ld_eigen` object containing pure truncated eigen
#'   representations of block-diagonal LD.
#'
#' @details Correlation-matrix input is first converted with [as_blm_ld()], so
#' exact contiguous sub-blocks are detected before eigendecomposition. Blocks
#' are processed sequentially. By default, materially negative eigenvalues are
#' rejected; `negative_eigenvalues = "discard"` instead provides SBayesRC-style
#' positive-eigenspace truncation. Tiny nonpositive eigenvalues at the numerical
#' tolerance are always discarded.
#'
#' Existing eigenpairs must describe correlation matrices, not scaled
#' cross-products. Their rows must follow `variants`, their columns must be
#' eigenvectors, and their eigenvalues may be supplied in any order. Supplied
#' eigenpairs with fewer columns than rows are conservatively recorded as an
#' incomplete source eigenspace unless, under the strict unit-diagonal policy,
#' their eigenvalues account for the complete correlation trace within
#' numerical tolerance. This prevents already-truncated input from being
#' treated as an exact representation merely because every supplied component
#' was retained.
#'
#' @examples
#' R <- matrix(c(1, 0.4, 0.4, 1), 2)
#' variants <- data.frame(
#'   CHR = 1, ID = c("rs1", "rs2"), POS = 1:2,
#'   A1 = c("A", "C"), A0 = c("C", "T")
#' )
#' ld_eigen <- as_blm_ld_eigen(R, variants, prop_var = 0.8)
#' @export
as_blm_ld_eigen <- function(
    R = NULL, variants = NULL, eigenvectors = NULL, eigenvalues = NULL,
    prop_var = 0.995,
    check_eigenvectors = FALSE,
    negative_eigenvalues = c("error", "discard")) {
  prop_var_missing <- missing(prop_var)
  negative_eigenvalues_missing <- missing(negative_eigenvalues)
  if (!is.numeric(prop_var) || length(prop_var) != 1L || is.na(prop_var) ||
      !is.finite(prop_var) || prop_var <= 0 || prop_var > 1) {
    stop("`prop_var` must be a finite number in (0, 1].", call. = FALSE)
  }
  if (!is.logical(check_eigenvectors) ||
      length(check_eigenvectors) != 1L || is.na(check_eigenvectors)) {
    stop("`check_eigenvectors` must be TRUE or FALSE.", call. = FALSE)
  }
  prop_var <- as.numeric(prop_var)
  negative_eigenvalues <- match.arg(negative_eigenvalues)

  supplied_pairs <- !is.null(eigenvectors) || !is.null(eigenvalues)
  if (supplied_pairs) {
    if (!is.null(R)) {
      stop("`R` must be NULL when eigenpairs are supplied.", call. = FALSE)
    }
    if (is.null(eigenvectors) || is.null(eigenvalues)) {
      stop("Supply both `eigenvectors` and `eigenvalues`.", call. = FALSE)
    }
    return(.as_blm_ld_eigen_from_pairs(
      eigenvectors, eigenvalues, variants, prop_var, check_eigenvectors,
      negative_eigenvalues
    ))
  }
  if (is.null(R)) {
    stop("Supply `R` or precomputed `eigenvectors` and `eigenvalues`.",
         call. = FALSE)
  }

  if (inherits(R, "blm_ld_eigen")) {
    if (!is.null(variants)) {
      stop(
        "`variants` must be NULL for `blm_ld_eigen` input.",
        call. = FALSE
      )
    }
    .validate_blm_ld_eigen_object(R)
    if (!prop_var_missing &&
        (is.na(R$requested_prop_var) ||
         !isTRUE(all.equal(prop_var, R$requested_prop_var)))) {
      stop(
        paste0(
          "`prop_var` cannot change an existing `blm_ld_eigen` object; ",
          "reconstruct it from LD or eigenpairs."
        ),
        call. = FALSE
      )
    }
    if (!negative_eigenvalues_missing) {
      policies <- unique(vapply(
        R$blocks, `[[`, character(1), "negative_eigenvalues"
      ))
      if (length(policies) != 1L || policies != negative_eigenvalues) {
        stop(
          paste0(
            "`negative_eigenvalues` cannot change an existing ",
            "`blm_ld_eigen` object; reconstruct it from LD or eigenpairs."
          ),
          call. = FALSE
        )
      }
    }
    return(R)
  }

  input_regularization_report <- NULL
  if (inherits(R, "blm_ld")) {
    if (!is.null(variants)) {
      stop(
        "`variants` must be NULL for `blm_ld` input.",
        call. = FALSE
      )
    }
    .validate_blm_ld_object(R)
    input_regularization_report <- R$regularization_report
    explicit_ld <- R
  } else {
    explicit_ld <- as_blm_ld(R, variants)
  }

  block_sizes <- vapply(explicit_ld$blocks, `[[`, integer(1), "size")
  block_ends <- cumsum(block_sizes)
  block_starts <- block_ends - block_sizes + 1L
  blocks <- lapply(seq_along(explicit_ld$blocks), function(block_index) {
    source <- explicit_ld$blocks[[block_index]]
    indices <- seq.int(block_starts[[block_index]], block_ends[[block_index]])
    decomposition <- eigen(
      .materialize_ld_block(source), symmetric = TRUE
    )
    .make_blm_ld_eigen_block(
      decomposition$vectors, decomposition$values,
      explicit_ld$variants$ID[indices], source$name, source$parent,
      prop_var, check_eigenvectors = FALSE,
      negative_eigenvalues = negative_eigenvalues
    )
  })
  names(blocks) <- names(explicit_ld$blocks)
  result <- .new_blm_ld_eigen(
    blocks, explicit_ld$variants, explicit_ld$parents, prop_var,
    input_regularization_report
  )
  .validate_blm_ld_eigen_object(result)
  result
}

.as_blm_ld_eigen_from_pairs <- function(
    eigenvectors, eigenvalues, variants, prop_var, check_eigenvectors,
    negative_eigenvalues) {
  list_input <- is.list(eigenvectors) && !is.matrix(eigenvectors)
  if (list_input) {
    if (!length(eigenvectors) || !is.list(eigenvalues) ||
        length(eigenvalues) != length(eigenvectors) ||
        !is.list(variants) || is.data.frame(variants) ||
        length(variants) != length(eigenvectors)) {
      stop(
        paste0(
          "List eigenvector input requires matching nonempty `eigenvalues` ",
          "and `variants` lists."
        ),
        call. = FALSE
      )
    }
    block_names <- names(eigenvectors)
    if (is.null(block_names) || anyNA(block_names) || any(block_names == "") ||
        anyDuplicated(block_names)) {
      stop("Eigenvector lists must have unique, nonempty names.",
           call. = FALSE)
    }
    if (is.null(names(eigenvalues)) ||
        !identical(names(eigenvalues), block_names) ||
        is.null(names(variants)) || !identical(names(variants), block_names)) {
      stop(
        "Eigenvector, `eigenvalues`, and `variants` lists must have identical names.",
        call. = FALSE
      )
    }
  } else {
    eigenvectors <- list(LD = eigenvectors)
    eigenvalues <- list(LD = eigenvalues)
    variants <- list(LD = variants)
    block_names <- "LD"
  }

  blocks <- vector("list", length(eigenvectors))
  tables <- vector("list", length(eigenvectors))
  for (block_index in seq_along(eigenvectors)) {
    name <- block_names[[block_index]]
    vectors <- eigenvectors[[block_index]]
    if (!is.matrix(vectors) || !is.numeric(vectors) ||
        nrow(vectors) < 1L || ncol(vectors) < 1L ||
        ncol(vectors) > nrow(vectors) ||
        !eigen_matrix_is_finite_cpp(vectors)) {
      stop(sprintf(
        "Eigenvectors for block `%s` must be a finite numeric m-by-q matrix.",
        name
      ), call. = FALSE)
    }
    table <- .validate_ld_variants(variants[[block_index]], nrow(vectors), name)
    vector_names <- rownames(vectors)
    if (!is.null(vector_names) &&
        !identical(as.character(vector_names), table$ID)) {
      stop(sprintf(
        "Row names of eigenvectors for block `%s` must match `variants$ID`.",
        name
      ), call. = FALSE)
    }
    blocks[[block_index]] <- .make_blm_ld_eigen_block(
      vectors, eigenvalues[[block_index]], table$ID, name, name,
      prop_var, check_eigenvectors, negative_eigenvalues,
      full_spectrum_supplied = ncol(vectors) == nrow(vectors)
    )
    tables[[block_index]] <- table
  }
  names(blocks) <- block_names
  all_variants <- do.call(rbind, unname(tables))
  rownames(all_variants) <- NULL
  if (anyDuplicated(all_variants$ID)) {
    stop("`variants$ID` must be unique across eigen blocks.", call. = FALSE)
  }
  result <- .new_blm_ld_eigen(
    blocks, all_variants, block_names, prop_var, NULL
  )
  .validate_blm_ld_eigen_object(result)
  result
}

.make_blm_ld_eigen_block <- function(
    eigenvectors, eigenvalues, predictor_names, name, parent, prop_var,
    check_eigenvectors, negative_eigenvalues,
    full_spectrum_supplied = TRUE) {
  if (!is.numeric(eigenvalues) || !is.atomic(eigenvalues) ||
      is.object(eigenvalues) || !is.null(dim(eigenvalues)) ||
      length(eigenvalues) != ncol(eigenvectors) || anyNA(eigenvalues) ||
      any(!is.finite(eigenvalues))) {
    stop(sprintf(
      "`eigenvalues` for block `%s` must match its eigenvector columns.", name
    ), call. = FALSE)
  }
  eigenvalues <- as.numeric(eigenvalues)
  storage.mode(eigenvectors) <- "double"
  order <- order(eigenvalues, decreasing = TRUE)
  eigenvalues <- eigenvalues[order]
  eigenvectors <- eigenvectors[, order, drop = FALSE]
  tolerance <- sqrt(.Machine$double.eps) *
    max(1, max(abs(eigenvalues)))
  minimum_source_eigenvalue <- min(eigenvalues)
  materially_negative <- eigenvalues < -tolerance
  if (any(materially_negative) && negative_eigenvalues == "error") {
    stop(sprintf(
      paste0(
        "LD block `%s` has a materially negative eigenvalue (minimum %.6g; ",
        "tolerance %.6g); regularize it before eigen conversion."
      ),
      name, min(eigenvalues), tolerance
    ), call. = FALSE)
  }
  positive <- eigenvalues > tolerance
  if (!any(positive)) {
    stop(sprintf("LD block `%s` has no positive eigenvalues.", name),
         call. = FALSE)
  }
  discarded_negative_eigenvalues <- sum(materially_negative)
  eigenvalues <- eigenvalues[positive]
  eigenvectors <- eigenvectors[, positive, drop = FALSE]
  if (check_eigenvectors) {
    error <- max(abs(crossprod(eigenvectors) - diag(ncol(eigenvectors))))
    orthogonality_tolerance <- sqrt(.Machine$double.eps) *
      max(1, nrow(eigenvectors))
    if (error > orthogonality_tolerance) {
      stop(sprintf(
        "Eigenvectors for block `%s` must have orthonormal columns.", name
      ), call. = FALSE)
    }
  }

  size <- nrow(eigenvectors)
  trace_basis <- if (negative_eigenvalues == "discard") {
    sum(eigenvalues)
  } else {
    size
  }
  target_trace <- prop_var * trace_basis
  cumulative <- cumsum(eigenvalues)
  trace_tolerance <- sqrt(.Machine$double.eps) * max(1, trace_basis)
  source_eigenspace_complete <- full_spectrum_supplied ||
    (negative_eigenvalues == "error" &&
     abs(sum(eigenvalues) - size) <= trace_tolerance)
  reached <- which(cumulative >= target_trace - trace_tolerance)[1L]
  if (is.na(reached)) {
    reached <- length(eigenvalues)
    warning(sprintf(
      paste0(
        "Supplied positive eigenpairs for block `%s` retain %.6f of its ",
        "correlation trace, below requested `prop_var = %.6f`; all supplied ",
        "eigenpairs were retained."
      ),
      name, cumulative[[reached]] / trace_basis, prop_var
    ), call. = FALSE)
  }
  complete_eigenspace <- source_eigenspace_complete &&
    reached == length(eigenvalues)
  eigenvalues <- eigenvalues[seq_len(reached)]
  eigenvectors <- eigenvectors[, seq_len(reached), drop = FALSE]
  rownames(eigenvectors) <- predictor_names
  colnames(eigenvectors) <- NULL
  retained_trace <- sum(eigenvalues)
  diagonal <- numeric(size)
  for (component in seq_along(eigenvalues)) {
    diagonal <- diagonal +
      eigenvalues[[component]] * eigenvectors[, component]^2
  }
  diagonal_tolerance <- sqrt(.Machine$double.eps) * max(1, size)
  if (negative_eigenvalues == "error" &&
      any(diagonal > 1 + diagonal_tolerance)) {
    stop(sprintf(
      "Eigenpairs for block `%s` are incompatible with a correlation matrix.",
      name
    ), call. = FALSE)
  }
  list(
    name = name,
    parent = parent,
    source_block = name,
    size = as.integer(size),
    rank = as.integer(length(eigenvalues)),
    eigenvectors = eigenvectors,
    eigenvalues = eigenvalues,
    retained_trace = retained_trace,
    trace_basis = trace_basis,
    prop_var = min(1, retained_trace / trace_basis),
    requested_prop_var = prop_var,
    eigenvalue_tolerance = tolerance,
    negative_eigenvalues = negative_eigenvalues,
    discarded_negative_eigenvalues = as.integer(
      discarded_negative_eigenvalues
    ),
    minimum_source_eigenvalue = minimum_source_eigenvalue,
    source_eigenspace_complete = source_eigenspace_complete,
    complete_eigenspace = complete_eigenspace
  )
}

.new_blm_ld_eigen <- function(
    blocks, variants, parents, requested_prop_var, regularization_report) {
  result <- structure(list(
    blocks = blocks,
    variants = variants,
    parents = parents,
    block_table = .ld_eigen_block_table(blocks),
    requested_prop_var = requested_prop_var,
    format_version = .blm_ld_eigen_format_version,
    cross_block_assumption = if (length(blocks) > 1L) "zero" else NULL
  ), class = "blm_ld_eigen")
  if (!is.null(regularization_report)) {
    result$regularization_report <- regularization_report
  }
  result
}

#' Combine eigen LD objects
#'
#' Combines independently constructed [blm_ld_eigen][as_blm_ld_eigen()]
#' objects without reconstructing LD matrices or copying eigen blocks.
#'
#' @param ... One or more `blm_ld_eigen` objects, or one nonempty list of such
#'   objects.
#'
#' @return A combined `blm_ld_eigen` object. Cross-object LD is assumed to be
#' exactly zero.
#'
#' @examples
#' R <- matrix(c(1, 0.2, 0.2, 1), 2)
#' variants1 <- data.frame(
#'   CHR = 1, ID = c("rs1", "rs2"), POS = 1:2,
#'   A1 = c("A", "C"), A0 = c("C", "T")
#' )
#' variants2 <- transform(variants1, CHR = 2, ID = c("rs3", "rs4"))
#' ld1 <- as_blm_ld_eigen(list(chr1 = R), list(chr1 = variants1))
#' ld2 <- as_blm_ld_eigen(list(chr2 = R), list(chr2 = variants2))
#' combined <- combine_blm_ld_eigen(ld1, ld2)
#' @export
combine_blm_ld_eigen <- function(...) {
  objects <- list(...)
  if (length(objects) == 1L && is.list(objects[[1L]]) &&
      !inherits(objects[[1L]], "blm_ld_eigen")) {
    objects <- objects[[1L]]
  }
  if (!length(objects)) {
    stop("At least one `blm_ld_eigen` object is required.", call. = FALSE)
  }
  for (index in seq_along(objects)) {
    tryCatch(
      .validate_blm_ld_eigen_object(objects[[index]]),
      error = function(condition) stop(sprintf(
        "Input eigen LD object %d is invalid: %s",
        index, conditionMessage(condition)
      ), call. = FALSE)
    )
  }
  if (length(objects) == 1L) return(objects[[1L]])
  parents <- unlist(lapply(objects, `[[`, "parents"), use.names = FALSE)
  if (anyDuplicated(parents)) {
    stop("LD parent names must be unique across inputs.", call. = FALSE)
  }
  blocks <- unlist(unname(lapply(objects, `[[`, "blocks")), recursive = FALSE)
  if (is.null(names(blocks)) || anyDuplicated(names(blocks))) {
    stop("LD eigen block names must be unique across inputs.", call. = FALSE)
  }
  variant_columns <- lapply(objects, function(object) names(object$variants))
  if (!all(vapply(
    variant_columns, identical, logical(1), variant_columns[[1L]]
  ))) {
    stop("LD variant metadata must have identical columns across inputs.",
         call. = FALSE)
  }
  variants <- do.call(rbind, lapply(objects, `[[`, "variants"))
  rownames(variants) <- NULL
  if (anyDuplicated(variants$ID)) {
    stop("Variant IDs must be unique across inputs.", call. = FALSE)
  }
  requested <- vapply(objects, `[[`, numeric(1), "requested_prop_var")
  result <- .new_blm_ld_eigen(
    blocks, variants, parents,
    if (all(requested == requested[[1L]])) requested[[1L]] else NA_real_,
    NULL
  )
  has_report <- vapply(
    objects, function(object) !is.null(object$regularization_report), logical(1)
  )
  if (any(has_report)) {
    reports <- lapply(objects, .ld_eigen_complete_regularization_report)
    result$regularization_report <- do.call(rbind, reports)
    rownames(result$regularization_report) <- NULL
  }
  # Every input block was fully validated above and is reused unchanged.
  # Validate only the newly assembled report; the constructor and explicit
  # cross-object checks establish the remaining combined metadata invariants.
  .validate_ld_regularization_report(result)
  result
}

.ld_eigen_complete_regularization_report <- function(object) {
  report <- object$regularization_report
  if (is.null(report)) {
    return(data.frame(
      block = names(object$blocks),
      source_block = names(object$blocks),
      subset_source_block = names(object$blocks),
      parent = vapply(object$blocks, `[[`, character(1), "parent"),
      predictors = as.numeric(vapply(object$blocks, `[[`, integer(1), "size")),
      method = "none", shrink = 0, floor_shrink = 0,
      minimum_eigenvalue_before = NA_real_,
      minimum_eigenvalue_after = NA_real_,
      positive_definite_after = NA,
      stringsAsFactors = FALSE
    ))
  }
  if (!"source_block" %in% names(report)) report$source_block <- report$block
  if (!"subset_source_block" %in% names(report)) {
    report$subset_source_block <- report$block
  }
  if (!"floor_shrink" %in% names(report)) report$floor_shrink <- 0
  report[c(
    "block", "source_block", "subset_source_block", "parent", "predictors",
    "method", "shrink",
    "floor_shrink", "minimum_eigenvalue_before", "minimum_eigenvalue_after",
    "positive_definite_after"
  )]
}

#' @export
print.blm_ld_eigen <- function(x, ...) {
  if (length(list(...))) stop("Additional arguments are not supported.",
                             call. = FALSE)
  total_rank <- sum(vapply(x$blocks, `[[`, integer(1), "rank"))
  cat(sprintf(
    paste0(
      "BayesLinReg eigen LD object: %d variants, %d block%s, ",
      "total rank %d\n"
    ),
    nrow(x$variants), length(x$blocks),
    if (length(x$blocks) == 1L) "" else "s", total_rank
  ))
  invisible(x)
}

.blm_ld_eigen_format_version <- 3L

.ld_eigen_block_table <- function(blocks) {
  sizes <- vapply(blocks, `[[`, integer(1), "size")
  ends <- cumsum(sizes)
  data.frame(
    block = names(blocks),
    parent = vapply(blocks, `[[`, character(1), "parent"),
    predictors = sizes,
    variant_start = ends - sizes + 1L,
    variant_end = ends,
    storage = "truncated_eigen",
    rank = vapply(blocks, `[[`, integer(1), "rank"),
    prop_var = vapply(blocks, `[[`, numeric(1), "prop_var"),
    negative_eigenvalues = vapply(
      blocks, `[[`, character(1), "negative_eigenvalues"
    ),
    discarded_negative_eigenvalues = vapply(
      blocks, `[[`, integer(1), "discarded_negative_eigenvalues"
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.validate_blm_ld_eigen_object <- function(ld) {
  if (!inherits(ld, "blm_ld_eigen") || !is.list(ld) ||
      !identical(ld$format_version, .blm_ld_eigen_format_version)) {
    stop(
      "`ld` must be an object returned by the current `as_blm_ld_eigen()`.",
      call. = FALSE
    )
  }
  if (!is.data.frame(ld$variants) || !is.list(ld$blocks) ||
      !length(ld$blocks) || !is.character(ld$parents) ||
      !length(ld$parents) || anyNA(ld$parents) || any(ld$parents == "") ||
      anyDuplicated(ld$parents) || !is.numeric(ld$requested_prop_var) ||
      length(ld$requested_prop_var) != 1L ||
      (!is.na(ld$requested_prop_var) &&
       (!is.finite(ld$requested_prop_var) || ld$requested_prop_var <= 0 ||
        ld$requested_prop_var > 1))) {
    stop("`ld` has an invalid eigen LD structure.", call. = FALSE)
  }
  required_variants <- c("CHR", "ID", "POS", "A1", "A0")
  if (!all(required_variants %in% names(ld$variants)) ||
      anyNA(ld$variants$ID) || any(ld$variants$ID == "") ||
      anyDuplicated(ld$variants$ID)) {
    stop("`ld` has invalid variant metadata.", call. = FALSE)
  }
  block_names <- names(ld$blocks)
  if (is.null(block_names) || anyNA(block_names) || any(block_names == "") ||
      anyDuplicated(block_names)) {
    stop("`ld` has invalid eigen block names.", call. = FALSE)
  }
  for (block_index in seq_along(ld$blocks)) {
    block <- ld$blocks[[block_index]]
    required <- c(
      "name", "parent", "source_block", "size", "rank", "eigenvectors",
      "eigenvalues", "retained_trace", "prop_var", "requested_prop_var",
      "eigenvalue_tolerance", "trace_basis", "negative_eigenvalues",
      "discarded_negative_eigenvalues", "minimum_source_eigenvalue",
      "source_eigenspace_complete", "complete_eigenspace"
    )
    valid <- is.list(block) && all(required %in% names(block)) &&
      identical(block$name, block_names[[block_index]]) &&
      is.character(block$parent) && length(block$parent) == 1L &&
      block$parent %in% ld$parents && is.character(block$source_block) &&
      length(block$source_block) == 1L && is.integer(block$size) &&
      length(block$size) == 1L && block$size >= 1L &&
      is.integer(block$rank) && length(block$rank) == 1L &&
      block$rank >= 1L && block$rank <= block$size &&
      is.matrix(block$eigenvectors) && is.numeric(block$eigenvectors) &&
      identical(dim(block$eigenvectors), c(block$size, block$rank)) &&
      eigen_matrix_is_finite_cpp(block$eigenvectors) &&
      is.numeric(block$eigenvalues) &&
      length(block$eigenvalues) == block$rank &&
      !anyNA(block$eigenvalues) && all(is.finite(block$eigenvalues)) &&
      all(block$eigenvalues > 0) &&
      all(diff(block$eigenvalues) <= 0) &&
      is.numeric(block$retained_trace) &&
      length(block$retained_trace) == 1L &&
      isTRUE(all.equal(block$retained_trace, sum(block$eigenvalues))) &&
      is.numeric(block$trace_basis) && length(block$trace_basis) == 1L &&
      is.finite(block$trace_basis) && block$trace_basis > 0 &&
      block$retained_trace <= block$trace_basis *
        (1 + sqrt(.Machine$double.eps)) &&
      is.numeric(block$prop_var) && length(block$prop_var) == 1L &&
      is.finite(block$prop_var) && block$prop_var > 0 && block$prop_var <= 1 &&
      isTRUE(all.equal(
        block$prop_var, min(1, block$retained_trace / block$trace_basis)
      )) &&
      is.numeric(block$requested_prop_var) &&
      length(block$requested_prop_var) == 1L &&
      is.finite(block$requested_prop_var) &&
      block$requested_prop_var > 0 && block$requested_prop_var <= 1 &&
      is.numeric(block$eigenvalue_tolerance) &&
      length(block$eigenvalue_tolerance) == 1L &&
      is.finite(block$eigenvalue_tolerance) &&
      block$eigenvalue_tolerance >= 0 &&
      is.character(block$negative_eigenvalues) &&
      length(block$negative_eigenvalues) == 1L &&
      block$negative_eigenvalues %in% c("error", "discard") &&
      is.integer(block$discarded_negative_eigenvalues) &&
      length(block$discarded_negative_eigenvalues) == 1L &&
      block$discarded_negative_eigenvalues >= 0L &&
      (block$negative_eigenvalues == "discard" ||
       block$discarded_negative_eigenvalues == 0L) &&
      is.numeric(block$minimum_source_eigenvalue) &&
      length(block$minimum_source_eigenvalue) == 1L &&
      is.finite(block$minimum_source_eigenvalue) &&
      is.logical(block$source_eigenspace_complete) &&
      length(block$source_eigenspace_complete) == 1L &&
      !is.na(block$source_eigenspace_complete) &&
      is.logical(block$complete_eigenspace) &&
      length(block$complete_eigenspace) == 1L &&
      !is.na(block$complete_eigenspace) &&
      (!block$complete_eigenspace || block$source_eigenspace_complete)
    if (!valid) stop("`ld` contains an invalid eigen block.", call. = FALSE)
  }
  sizes <- vapply(ld$blocks, `[[`, integer(1), "size")
  block_parents <- vapply(ld$blocks, `[[`, character(1), "parent")
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  row_names_match <- all(vapply(seq_along(ld$blocks), function(index) {
    identical(
      rownames(ld$blocks[[index]]$eigenvectors),
      as.character(ld$variants$ID[seq.int(starts[[index]], ends[[index]])])
    )
  }, logical(1)))
  if (sum(sizes) != nrow(ld$variants) ||
      !identical(unique(unname(block_parents)), unname(ld$parents)) ||
      !row_names_match ||
      !identical(ld$block_table, .ld_eigen_block_table(ld$blocks)) ||
      !identical(
        ld$cross_block_assumption,
        if (length(ld$blocks) > 1L) "zero" else NULL
      )) {
    stop("`ld` eigen block metadata are inconsistent.", call. = FALSE)
  }
  .validate_ld_regularization_report(ld)
  invisible(ld)
}
