setOldClass("aalen")

# This is not assigned to the class name because I want a custom constructor
#' @include TimeDependentModel-class.R
#' @include utils.R
setClass(
  "AalenTimeDependentModel",
  contains = "TimeDependentModel",
  slots = c(
    model = "aalen",
    effectFunction = "function",
    varianceFunction = "function"
  )
)

#' Models the relationship between a polygenic score and a phenotype calculated
#' through Aalen's additive hazards model
#'
#' @param model The underlying model
#' @param cumulative_effects Numeric vector of cumulative effects
#' @param cumulative_variances Numeric vector of cumulative variances
#'
#' @slot model The Aalen model of the relationship
#' @slot effectFunction Function which can be used to retrieve the effect at a
#'   certain time point
#' @slot varianceFunction Function to retrieve the variance of the effect at a
#'   certain time point
#'
#' @seealso [LoessTimeDependentModel()], [GlmTimeDependentModel()],
#'   [timereg::aalen()]
#' @docType class
#' @export
AalenTimeDependentModel <- function(model, cumulative_effects, cumulative_variances) {
  if (!is(model, "aalen"))
    rlang::abort("model must be an Aalen model")
  new(
    "AalenTimeDependentModel",
    model = model,
    effectFunction = derive_without_noise(x=model$cum[, "time"], y = cumulative_effects),
    varianceFunction = derive_without_noise(x=model$cum[, "time"], y = cumulative_variances)
  )
}

#' @rdname totalEffect-methods
#' @aliases totalEffect,GlmTimeDependentModel-method
setMethod("totalEffect", "AalenTimeDependentModel", function(this, age){
  this@effectFunction(age)
})

#' Calculate hazard variance at a certain age
#'
#' @param this Model to use for calcularion
#' @param age Numeric vector of time points to retrieve the variances for.
#' @export
#' @docType methods
#' @rdname variance-methods
setGeneric("variance", function(this, age) standardGeneric("variance"))
#' @rdname variance-methods
#' @aliases variance,AalenTimeDependentModel-method
setMethod("variance", "AalenTimeDependentModel", function(this, age) {
  this@varianceFunction(age)
})
