
#' Get last element of a vector-like object
#'
#' @param x Vector or list or anything that works with [`tail()`].
#'
#' @return The last element of `x`
#' @keywords internal
#'
#' @examples
#' last(1:3)
last <- function(x) {
  tail(x, n = 1)
}
