#' Numerically integrate effect estimates by age
#'
#' @param age Numeric vector containing ages
#' @param beta Numeric vector containing effect estimates
#' @param method The integration method to use. Either `trapezoidal` (default)
#'   or `midpoint`.
#'
#' @return A numeric matrix with columns age and cumulative integral
#' @keywords internal
integrate_estimates <- function(age, beta, method = c("trapezoidal", "midpoint")){
  stopifnot(
    is.numeric(age),
    is.numeric(beta)
  )
  method <- match.arg(method)
  effects <- switch(method,
    trapezoidal = pracma::cumtrapz(age, beta),
    midpoint =  c(0, cumsum(diff(age) * zoo::rollmean(beta, 2))),
    "You should not be here"
  )
  effects
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

#' Perform time-dependent MR analysis
#'
#' @param age_seq Ages at which to estimate exposure-outcome effect
#' @param exposure_model Time-dependent model
#'   ([`GlmTimeDependentModel`][GlmTimeDependentModel()] or
#'   [`LoessTimeDependentModel`][LoessTimeDependentModel()]) of genetic effects
#'   on exposure
#' @param outcome_model Time-dependent model
#'   ([`AalenTimeDependentModel`][AalenTimeDependentModel()]) of genetic effects
#'   on outcome
#' @param method Integration method for calculating time-dependente
#'   cumulative effects. `trapezoidal` (default) or `midpoint`.
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{age}{Time point of estimated effect}
#'     \item{Gamma}{Cumulative effect}
#'     \item{var}{Variance of cumulative effect}
#'     \item{L95}{Lower end of 95% CI}
#'     \item{H95}{Higher end of 95% CI}
#'   }
#' @export
#' @seealso [calculate_genetic_effects()], [time_dependent_aalen()],
#'   [time_dependent_glm()], [time_dependent_loess()],
#'   [`GlmTimeDependentModel`][GlmTimeDependentModel()],
#'   [`LoessTimeDependentModel`][LoessTimeDependentModel()],
#'   [`AalenTimeDependentModel`][AalenTimeDependentModel()]
#' @examples
#' # Generate random example data
#' library(stats)
#' pgs <- rnorm(1000)
#' exposure <- pgs * 3
#' exposure_age <- rep(40:70, length.out = 1000)
#' covariates <- data.frame(covar = rnorm(1000))
#' outcome <- sample(c(TRUE,FALSE), 1000, replace=TRUE)
#' outcome_age <- rnorm(1000, mean=60, sd = 3) |> pmax(40)
#'
#' # Estimate exposure and outcome effects
#' models <- calculate_genetic_effects(pgs, exposure, outcome, exposure_age, outcome_age,
#'   covariates)
#' time_dependent_MR(40:70, models$exposure, models$outcome)
time_dependent_MR <- function(age_seq, exposure_model, outcome_model,
                              method = c("trapezoidal", "midpoint")) {
  # Checking arguments
  method <- match.arg(method)
  if (!is(exposure_model, "TimeDependentModel"))
    rlang::abort(
      paste0("exposure_model should be a TimeDependentModel. Is ", class(exposure_model))
    )
  if (!is(outcome_model, "AalenTimeDependentModel"))
    rlang::abort(
      paste0("outcome_model should be a AalenTimeDependentModel. Is ", class(outcome_model))
    )
  if (any(is.na(age_seq)) || is.unsorted(age_seq)) {
    age_str <- paste0(age_seq, collapse = ",")
    rlang::abort(paste0("age_seq is an invalid range: ", age_seq))
  }

  midpoints <-  zoo::rollmean(age_seq, 2)
  total_exposure_effects <- totalEffect(exposure_model, age_seq)
  midpoint_exposure_effects <- totalEffect(exposure_model, midpoints)
  outcome_variances <- variance(outcome_model, age_seq)
  outcome_effects <- totalEffect(outcome_model, age_seq)
  dBeta <- diff(outcome_effects) / diff(age_seq)
  wald_ratio <- dBeta / midpoint_exposure_effects
  Gamma <- c(0, integrate_estimates(midpoints, wald_ratio))
  Gamma_variance <- switch(method,
    trapezoidal = calculate_mr_variance_trapz(
      effect = total_exposure_effects,
      variance = outcome_variances
    ),
    midpoint =  outcome_variances / total_exposure_effects ^ 2,
    "Method is neither trapezoidal nor midpoint. There should've been an error earlier."
  )
  data.frame(
    age = c(0, midpoints),
    Gamma = Gamma,
    var = Gamma_variance,
    L95 = Gamma - 1.96 * sqrt(Gamma_variance),
    H95 = Gamma + 1.96 * sqrt(Gamma_variance)
  )

}
