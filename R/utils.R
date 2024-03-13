
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
