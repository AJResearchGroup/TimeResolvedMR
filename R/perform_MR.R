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

  midpoints <- zoo::rollmean(age_seq, 2)

  total_exposure_effects <- totalEffect(exposure_model, age_seq)
  midpoint_exposure_effects <- totalEffect(exposure_model, midpoints)

  outcome_variances <- variance(outcome_model, age_seq)
  outcome_effects <- totalEffect(outcome_model, age_seq)
  dBeta <- diff(outcome_effects) / diff(age_seq)

  wald_ratio <- dBeta / midpoint_exposure_effects

  switch(method,
    trapezoidal = integrate_trapezoid(
      age = midpoints, wald_ratio = wald_ratio,
      outcome_variance = outcome_variances,
      exposure_effect = total_exposure_effects
    ),
    midpoint = integrate_midpoint(
      age = age_seq, wald_ratio = wald_ratio,
      outcome_variance = outcome_variances,
      exposure_effect = total_exposure_effects
    ),
    "Method is neither trapezoidal nor midpoint. There should've been an error earlier."
  )
}
