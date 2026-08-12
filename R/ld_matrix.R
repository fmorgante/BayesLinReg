#' Construct an LD matrix object
#'
#' Constructs the block-aware linkage-disequilibrium representation used by
#' [blm_gwas()]. Exact contiguous block-diagonal structure is detected within
#' every supplied matrix. The matrices contain signed correlations, not
#' squared correlations.
#'
#' @param R A finite symmetric correlation matrix, or a nonempty named list of
#'   such matrices representing exactly zero cross-block LD. Dense base-R
#'   matrices and compressed sparse `Matrix` objects of class `dgCMatrix` or
#'   `dsCMatrix` are supported.
#' @param variants A data frame corresponding to `R`, or a named list of data
#'   frames corresponding to list `R`. Every table must contain `CHR`, `ID`,
#'   `POS`, `A1`, and `A0`, in matrix order. Additional columns are preserved.
#'   `A1` is the allele whose dosage was used to calculate `R`.
#'
#' @return An object of class `blm_ld`. Its internal representation is an
#'   implementation detail; use it as the `ld` argument of [blm_gwas()].
#' @export
as_blm_ld <- function(R, variants) {
  list_input <- is.list(R) && !is.matrix(R) &&
    !inherits(R, "sparseMatrix")
  if (list_input) {
    if (!length(R)) stop("List `R` must be nonempty.", call. = FALSE)
    if (!is.list(variants) || is.data.frame(variants) ||
        length(variants) != length(R)) {
      stop("List `R` requires a matching list `variants`.", call. = FALSE)
    }
    block_names <- names(R)
    variant_names <- names(variants)
    if (is.null(block_names) || anyNA(block_names) || any(block_names == "") ||
        anyDuplicated(block_names)) {
      stop("List `R` must have unique, nonempty names.", call. = FALSE)
    }
    if (is.null(variant_names) || !identical(block_names, variant_names)) {
      stop("`R` and `variants` lists must have identical names.",
           call. = FALSE)
    }
  } else {
    R <- list(LD = R)
    variants <- list(LD = variants)
    block_names <- "LD"
  }

  computational_blocks <- list()
  parent_tables <- vector("list", length(R))
  for (parent_index in seq_along(R)) {
    parent <- block_names[parent_index]
    matrix <- .validate_ld_correlation(R[[parent_index]], parent)
    table <- .validate_ld_variants(
      variants[[parent_index]], ncol(matrix), parent
    )
    matrix_names <- colnames(matrix)
    if (!is.null(matrix_names) && !identical(matrix_names, table$ID)) {
      stop(
        sprintf("Column names of `R[[\"%s\"]]` must match `variants$ID`.",
                parent),
        call. = FALSE
      )
    }
    parent_tables[[parent_index]] <- table
    ranges <- .exact_contiguous_ld_blocks(matrix)
    child_count <- length(ranges)
    children <- lapply(seq_along(ranges), function(child_index) {
      indices <- ranges[[child_index]]
      child_name <- if (child_count == 1L) {
        parent
      } else {
        paste0(parent, ".", child_index)
      }
      .compress_ld_block(
        matrix[indices, indices, drop = FALSE], parent, child_name
      )
    })
    names(children) <- vapply(children, `[[`, character(1), "name")
    computational_blocks <- c(computational_blocks, children)
  }

  computational_names <- make.unique(
    vapply(computational_blocks, `[[`, character(1), "name"), sep = "."
  )
  for (block_index in seq_along(computational_blocks)) {
    computational_blocks[[block_index]]$name <-
      computational_names[[block_index]]
  }
  names(computational_blocks) <- computational_names

  all_variants <- do.call(rbind, unname(parent_tables))
  rownames(all_variants) <- NULL
  if (anyDuplicated(all_variants$ID)) {
    duplicates <- unique(all_variants$ID[duplicated(all_variants$ID)])
    stop(
      sprintf("`variants$ID` must be unique; duplicated ID(s): %s.",
              paste(utils::head(duplicates, 10L), collapse = ", ")),
      call. = FALSE
    )
  }
  block_table <- .ld_block_table(computational_blocks)
  structure(
    list(
      blocks = computational_blocks,
      variants = all_variants,
      parents = block_names,
      block_table = block_table,
      format_version = .blm_ld_format_version,
      cross_block_assumption = if (length(computational_blocks) > 1L) {
        "zero"
      } else {
        NULL
      }
    ),
    class = "blm_ld"
  )
}

#' @export
print.blm_ld <- function(x, ...) {
  if (length(list(...))) {
    stop("Additional arguments are not supported.", call. = FALSE)
  }
  cat(sprintf(
    "BayesLinReg LD object: %d variants, %d parent group%s, %d computational block%s\n",
    nrow(x$variants), length(x$parents), if (length(x$parents) == 1L) "" else "s",
    length(x$blocks), if (length(x$blocks) == 1L) "" else "s"
  ))
  invisible(x)
}

.validate_ld_correlation <- function(matrix, label) {
  matrix <- .validate_gram_block(
    matrix, sprintf("`R[[\"%s\"]]`", label), TRUE
  )
  diagonal <- if (inherits(matrix, "sparseMatrix")) {
    Matrix::diag(matrix)
  } else {
    diag(matrix)
  }
  tolerance <- 1e-6
  if (any(abs(diagonal - 1) > tolerance)) {
    stop(sprintf("LD block `%s` must have a unit diagonal.", label),
         call. = FALSE)
  }
  if (inherits(matrix, "sparseMatrix")) {
    Matrix::diag(matrix) <- 1
  } else {
    diag(matrix) <- 1
  }
  maximum <- if (inherits(matrix, "sparseMatrix")) {
    if (length(matrix@x)) max(abs(matrix@x)) else 0
  } else {
    max(abs(matrix))
  }
  if (maximum > 1 + sqrt(.Machine$double.eps)) {
    stop(sprintf("LD block `%s` contains a correlation outside [-1, 1].",
                 label), call. = FALSE)
  }
  matrix
}

.validate_ld_variants <- function(variants, size, label) {
  if (!is.data.frame(variants)) {
    stop(sprintf("`variants[[\"%s\"]]` must be a data frame.", label),
         call. = FALSE)
  }
  required <- c("CHR", "ID", "POS", "A1", "A0")
  missing <- setdiff(required, names(variants))
  if (length(missing)) {
    stop(sprintf(
      "`variants[[\"%s\"]]` is missing required column(s): %s.",
      label, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  if (nrow(variants) != size) {
    stop(sprintf("`variants[[\"%s\"]]` must have %d rows.", label, size),
         call. = FALSE)
  }
  result <- as.data.frame(variants, stringsAsFactors = FALSE)
  for (name in c("CHR", "ID", "A1", "A0")) {
    result[[name]] <- as.character(result[[name]])
    if (anyNA(result[[name]]) || any(result[[name]] == "")) {
      stop(sprintf("`variants$%s` must not contain missing or empty values.",
                   name), call. = FALSE)
    }
  }
  if (!is.numeric(result$POS) || anyNA(result$POS) ||
      any(!is.finite(result$POS)) || any(result$POS != floor(result$POS)) ||
      any(result$POS < 1)) {
    stop("`variants$POS` must contain positive finite integers.",
         call. = FALSE)
  }
  result$POS <- as.numeric(result$POS)
  result$A1 <- toupper(result$A1)
  result$A0 <- toupper(result$A0)
  if (any(result$A1 == result$A0)) {
    stop("`variants$A1` and `variants$A0` must differ.", call. = FALSE)
  }
  if (anyDuplicated(result$ID)) {
    stop(sprintf("`variants$ID` must be unique within LD block `%s`.", label),
         call. = FALSE)
  }
  result
}

.ld_lower_triplets <- function(matrix) {
  if (inherits(matrix, "sparseMatrix")) {
    lower <- Matrix::forceSymmetric(matrix, uplo = "L")
    entries <- Matrix::summary(lower)
    keep <- entries$i > entries$j & entries$x != 0
    return(list(
      row = as.integer(entries$i[keep]),
      column = as.integer(entries$j[keep]),
      value = as.numeric(entries$x[keep])
    ))
  }
  positions <- which(lower.tri(matrix) & matrix != 0, arr.ind = TRUE)
  if (!nrow(positions)) {
    return(list(row = integer(), column = integer(), value = numeric()))
  }
  list(
    row = as.integer(positions[, 1L]),
    column = as.integer(positions[, 2L]),
    value = as.numeric(matrix[positions])
  )
}

.exact_contiguous_ld_blocks <- function(matrix) {
  entries <- .ld_lower_triplets(matrix)
  .exact_contiguous_ld_triplet_blocks(
    ncol(matrix), entries$row, entries$column, entries$value
  )
}

.exact_contiguous_ld_triplet_blocks <- function(p, row, column, value) {
  if (p == 1L) return(list(1L))
  crossing_delta <- integer(p)
  nonzero <- value != 0
  if (any(nonzero)) {
    starts <- pmin(row[nonzero], column[nonzero])
    ends <- pmax(row[nonzero], column[nonzero])
    crossing_delta <- tabulate(starts, nbins = p) -
      tabulate(ends, nbins = p)
  }
  crossing <- cumsum(crossing_delta)[seq_len(p - 1L)]
  boundaries <- c(which(crossing == 0L), p)
  starts <- c(1L, boundaries[-length(boundaries)] + 1L)
  Map(seq.int, starts, boundaries)
}

.compress_ld_block <- function(matrix, parent, name) {
  entries <- .ld_lower_triplets(matrix)
  .compress_ld_triplets(
    ncol(matrix), entries$row, entries$column, entries$value, parent, name
  )
}

.compress_ld_triplets <- function(p, row, column, value, parent, name) {
  nonzero <- value != 0
  row <- as.integer(row[nonzero])
  column <- as.integer(column[nonzero])
  value <- as.numeric(value[nonzero])
  counts <- tabulate(column, nbins = p)
  indexed_bytes <- 12 * length(value) + 4 * (p + 1)
  last_row <- seq_len(p)
  if (length(value)) {
    split_rows <- split(row, column)
    columns <- as.integer(names(split_rows))
    last_row[columns] <- vapply(split_rows, max, integer(1))
  }
  spans <- pmax(0L, last_row - seq_len(p))
  interval_values <- sum(as.double(spans))
  interval_bytes <- 8 * interval_values + 4 * (p + 1)

  if (interval_bytes <= indexed_bytes) {
    indptr <- c(0, cumsum(spans))
    data <- numeric(indptr[p + 1L])
    if (length(value)) {
      positions <- indptr[column] + (row - column)
      data[positions] <- value
    }
    row_index <- integer()
    type <- 0L
    storage <- "interval_triangular"
  } else {
    order <- order(column, row)
    data <- value[order]
    row_index <- row[order] - 1L
    counts <- tabulate(column[order], nbins = p)
    indptr <- c(0L, cumsum(counts))
    type <- 1L
    storage <- "indexed_triangular"
  }
  list(
    name = name,
    parent = parent,
    size = as.integer(p),
    type = type,
    storage = storage,
    data = as.numeric(data),
    indptr = as.integer(indptr),
    row_index = as.integer(row_index)
  )
}

.blm_ld_format_version <- 1L

.validate_blm_ld_object <- function(ld) {
  if (!inherits(ld, "blm_ld") || !is.list(ld)) {
    stop("`ld` must be an object returned by `as_blm_ld()`.", call. = FALSE)
  }
  if (!identical(ld$format_version, .blm_ld_format_version)) {
    stop(
      paste0(
        "`ld` uses an unsupported internal format; recreate it with the ",
        "current `as_blm_ld()`."
      ),
      call. = FALSE
    )
  }
  if (!is.data.frame(ld$variants) || !is.list(ld$blocks) ||
      !length(ld$blocks) || !is.character(ld$parents) ||
      !length(ld$parents) || anyNA(ld$parents) || any(ld$parents == "") ||
      anyDuplicated(ld$parents)) {
    stop("`ld` has an invalid internal structure.", call. = FALSE)
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
    stop("`ld` has invalid computational block names.", call. = FALSE)
  }
  for (block_index in seq_along(ld$blocks)) {
    block <- ld$blocks[[block_index]]
    required <- c(
      "name", "parent", "size", "type", "storage", "data", "indptr",
      "row_index"
    )
    if (!is.list(block) || !all(required %in% names(block)) ||
        !identical(block$name, block_names[[block_index]]) ||
        !is.character(block$parent) || length(block$parent) != 1L ||
        is.na(block$parent) || !block$parent %in% ld$parents ||
        !is.character(block$storage) || length(block$storage) != 1L ||
        is.na(block$storage) ||
        !is.integer(block$size) || length(block$size) != 1L ||
        is.na(block$size) || block$size < 1L ||
        !is.integer(block$type) || length(block$type) != 1L ||
        !block$type %in% 0:1 || !is.numeric(block$data) ||
        anyNA(block$data) || any(!is.finite(block$data)) ||
        any(abs(block$data) > 1 + sqrt(.Machine$double.eps)) ||
        !is.integer(block$indptr) ||
        length(block$indptr) != block$size + 1L ||
        block$indptr[[1L]] != 0L ||
        block$indptr[[block$size + 1L]] != length(block$data) ||
        any(diff(block$indptr) < 0L) || !is.integer(block$row_index)) {
      stop("`ld` contains an invalid compressed block.", call. = FALSE)
    }
    counts <- diff(block$indptr)
    if (block$type == 0L) {
      if (!identical(block$storage, "interval_triangular") ||
          length(block$row_index) ||
          any(seq_len(block$size) + counts > block$size)) {
        stop("`ld` contains an invalid interval block.", call. = FALSE)
      }
    } else {
      if (!identical(block$storage, "indexed_triangular") ||
          length(block$row_index) != length(block$data) ||
          any(block$row_index < 0L | block$row_index >= block$size)) {
        stop("`ld` contains an invalid indexed block.", call. = FALSE)
      }
    }
  }
  block_sizes <- vapply(ld$blocks, `[[`, integer(1), "size")
  if (sum(block_sizes) != nrow(ld$variants) ||
      !identical(ld$block_table, .ld_block_table(ld$blocks))) {
    stop("`ld` block metadata are inconsistent.", call. = FALSE)
  }
  invisible(ld)
}

.ld_block_table <- function(blocks) {
  sizes <- vapply(blocks, `[[`, integer(1), "size")
  variant_end <- cumsum(sizes)
  data.frame(
    block = names(blocks),
    parent = vapply(blocks, `[[`, character(1), "parent"),
    predictors = sizes,
    variant_start = variant_end - sizes + 1L,
    variant_end = variant_end,
    storage = vapply(blocks, `[[`, character(1), "storage"),
    stored_values = vapply(
      blocks, function(block) length(block$data), integer(1)
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.ld_block_triplets <- function(block) {
  counts <- diff(block$indptr)
  columns <- rep.int(seq_len(block$size), counts)
  if (!length(columns)) {
    return(list(row = integer(), column = integer(), value = numeric()))
  }
  rows <- if (block$type == 0L) {
    seq_along(block$data) - rep.int(block$indptr[-length(block$indptr)], counts) +
      columns
  } else {
    block$row_index + 1L
  }
  list(row = as.integer(rows), column = columns, value = block$data)
}

.materialize_blm_ld <- function(ld, selected = seq_len(nrow(ld$variants))) {
  selected <- as.integer(selected)
  selected_position <- integer(nrow(ld$variants))
  selected_position[selected] <- seq_along(selected)
  result <- matrix(0, length(selected), length(selected))
  offset <- 0L
  for (block in ld$blocks) {
    size <- block$size
    indices <- seq.int(offset + 1L, offset + size)
    kept <- indices[selected_position[indices] > 0L]
    if (length(kept)) {
      kept_positions <- selected_position[kept]
      result[cbind(kept_positions, kept_positions)] <- 1
    }
    for (column in seq_len(size)) {
      start <- block$indptr[column] + 1L
      end <- block$indptr[column + 1L]
      if (start > end) next
      value_positions <- seq.int(start, end)
      rows <- if (block$type == 0L) {
        column + seq_along(value_positions)
      } else {
        block$row_index[value_positions] + 1L
      }
      values <- block$data[value_positions]
      global_column <- offset + column
      global_rows <- offset + rows
      keep <- selected_position[global_rows] > 0L &
        selected_position[global_column] > 0L
      if (any(keep)) {
        selected_rows <- selected_position[global_rows[keep]]
        selected_column <- selected_position[global_column]
        result[cbind(selected_rows, selected_column)] <- values[keep]
        result[cbind(selected_column, selected_rows)] <- values[keep]
      }
    }
    offset <- offset + size
  }
  dimnames(result) <- list(ld$variants$ID[selected], ld$variants$ID[selected])
  result
}
