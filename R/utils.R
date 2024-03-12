
#' Get last element of a vector-like object
#'
#' @param x Vector or list or anything that works with [`tail()`].
#'
#' @return The last element of `x`
#' @keywords internal
#'
last <- function(x) {
  tail(x, n = 1)
}
