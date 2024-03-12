#' Integrate time-dependent estimates based on midpoint rule and compute 95% CIs
#'
#' @param age The ages at which to integrate
#' @param wald_ratio The calculated Wald ratios
#' @param exposure_effect The genetic effect on the exposure at every age
#' @param outcome_variance The variance of the genetic effect on the outcome
#'   at every age
#'
#' @return A data.frame with age, effect (Gamma), variance and 95% CIs
#' @keywords internal
integrate_midpoint <- function(age, wald_ratio, exposure_effect, outcome_variance) {
  Gamma <- c(0, cumsum(diff(age) * wald_ratio))
  Gamma_variance <- outcome_variance / exposure_effect ^ 2
  data.frame(
    age = age, Gamma = Gamma, Variance = Gamma_variance,
    L95 = Gamma - 1.96 * sqrt(Gamma_variance),
    H95 = Gamma + 1.96 * sqrt(Gamma_variance)
  )
}

#' Integrate time-dependent effects using the trapezoidal rule and compute 95%
#' CIs
#'
#' @param age The ages at which to integrate
#' @param wald_ratio The calculated Wald ratios
#' @param exposure_effect The genetic effect on the exposure at every age
#' @param outcome_variance The variance of the genetic effect on the outcome
#'   at every age
#'
#' @return A data.frame with age, effect (Gamma), variance and 95% CIs
#' @keywords internal
integrate_trapezoid <- function(age, wald_ratio, exposure_effect, outcome_variance) {
  # Add (0,0) if it isn't there yet
  if (age[1] != 0) {
    age <- c(0, age)
    wald_ratio <- c(0, wald_ratio)
  }
  Gamma <- pracma::cumtrapz(age, wald_ratio)
  Gamma_variance <- calculate_mr_variance_trapz(
    effect = exposure_effect,
    variance = outcome_variance
  )
  data.frame(
    age = age, Gamma = Gamma, Variance = Gamma_variance,
    L95 = Gamma - 1.96 * sqrt(Gamma_variance),
    H95 = Gamma + 1.96 * sqrt(Gamma_variance)
  )
}
#' Calculate time-dependent effect variance of MR effects estimated by
#' integration using the trapezoidal rules
#'
#' @param effect Numeric vector of total genetic effects on exposure for each
#'   time point
#' @param variance Numeric vector of variance of genetic effect on outcome for
#'   each time point.
#'
#' @return A numeric vector of variance at each time point
#' @keywords internal
calculate_mr_variance_trapz <- function(effect, variance) {
  stopifnot(
    is.numeric(effect),
    is.numeric(variance),
    length(effect) == length(variance),
    length(effect) > 0
  )
  n <- length(effect)
  var1 <- variance[1] / (16 * effect[1] ^ 2)
  if (n < 2) return(var1)

  var2 <- var1 + variance[2] / (4 * effect[2] ^ 2)
  if (n < 3) return (c(var1,var2))

  # At this point, we know that n > 2, so the following is safe:
  point_variances <- c(
    var1,
    vapply(2:n, \(k){ variance[k] / (4 * effect[k] ^2) }, numeric(1))
  )

  c(
    var1,
    var2,
    vapply(3:n, \(k) { var1 + point_variances[k-1] + point_variances[k] }, numeric(1))
  )
}
