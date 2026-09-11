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
