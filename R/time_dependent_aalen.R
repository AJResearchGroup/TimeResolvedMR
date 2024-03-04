#' Calculate time-dependent effect of polygenic scores on hazard
#'
#' @param pgs Numeric vector of polygenic scores
#' @param event Event status, coded as specified by [`survival::Surv()`]
#' @param event_age Numeric vector of age at which event occured
#' @param covariates Data frame of covariates
#'
#' @return An instance of [AalenTimeDependentModel()]
#'
#' @seealso
#'   [AalenTimeDependentModel()]
#'   [time_dependent_glm()],
#'   [time_dependent_loess()],
#'   [timereg::aalen()]
#' @export
#'
#' @examples
time_dependent_aalen <- function(pgs, event, event_age, covariates) {
  # PGS effect on outcome
  outcome_model <- timereg::aalen(
    survival::Surv(event_age, event) ~ pgs + .,
    data = covariates
  )
  AalenTimeDependentModel(
    model = outcome_model,
    cumulative_effects = model$cum[, "pgs"]
  )
}
