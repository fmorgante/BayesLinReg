.ld_regularization_block_map <- function(source_report, fitted_report) {
  fitted_block <- as.character(fitted_report$block)
  source_row <- match(fitted_block, source_report$block)
  unmatched <- is.na(source_row)
  if (any(unmatched)) {
    immediate_source <- if (
      "subset_source_block" %in% names(fitted_report)
    ) {
      as.character(fitted_report$subset_source_block)
    } else if ("source_block" %in% names(fitted_report)) {
      # Legacy reports only recorded original provenance. This remains
      # sufficient for one filtering stage and for unsplit blocks.
      as.character(fitted_report$source_block)
    } else {
      fitted_block
    }
    source_row[unmatched] <- match(
      immediate_source[unmatched], source_report$block
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
