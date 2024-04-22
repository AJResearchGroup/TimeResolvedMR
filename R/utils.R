
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

#' Derives cumulative effects where noise has been added to break ties
#'
#' @param x The parameter of the function (usually time). The fractional part
#'   of the values are assumed to be noise.
#' @param y The values of the function (usually cumulative effect)
#'
#' @return A numerically-derived function
#' @keywords internal
derive_without_noise <- function(x,y) {
  int_x <- round(x)
  unique_x <- unique(int_x)
  collapsed <- stats::approx(
    x = int_x, y = y, xout = unique_x, ties = last, rule = 2, method = "constant"
  )$y
  derivative <- c(0, diff(collapsed) / diff(unique_x))
  stats::approxfun(
    x = unique_x, y = derivative, method = "constant", rule = 2, f = 1
  )
}
