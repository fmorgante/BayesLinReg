#' Report the compiled dense linear-algebra backend
#'
#' Reports whether the package was compiled with `EIGEN_USE_BLAS`, which lets
#' eligible dense Eigen operations call the BLAS implementation linked to R.
#' The default source build uses Eigen's native kernels. See the package
#' `INSTALL` file for opt-in installation instructions.
#'
#' @return A list with `eigen_blas`, indicating whether Eigen BLAS delegation
#'   was enabled; `dense_kernel`, naming the compiled Eigen backend; and
#'   `blas`, the BLAS library reported by the running R installation.
#' @export
#'
#' @examples
#' blm_build_info()
blm_build_info <- function() {
  information <- blm_build_info_cpp()
  blas <- unname(extSoftVersion()["BLAS"])
  information$blas <- if (length(blas) == 1L && !is.na(blas)) {
    blas
  } else {
    NA_character_
  }
  information
}
