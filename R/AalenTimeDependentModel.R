class_name <- "AalenTimeDependentModel"

# This is not assigned to the class name because I want a custom constructor
setClass(
  class_name,
  contains = "TimeDependentModel",
  slots = c(
    model = "aalen",
    effectFunction = "function"
  )
)

#' Models the relationship between a polygenic score and a phenotype calculated
#' through Aalen's additive hazards model
#'
#' @slot model The Aalen model of the relationship
#' @slot effectFunction Function which can be used to retrieve the effect at a
#'   certain time point
#'
#' @seealso [LoessTimeDependentModel()], [GlmTimeDependentModel()],
#'   [timereg::aalen()]
#' @export
AalenTimeDependentModel <- function(model, cumulative_effects) {
  if (!is(model, "aalen"))
    rlang::abort("model must be an Aalen model")
  midpoints <- zoo::rollmean(model$cum[, "time"], 2)
  point_effects <- diff(cumulative_effects) / diff(model$cum[, "time"])
  new(
    class_name,
    model = model,
    effectFunction = approxfun(x=midpoints, y=point_effects, method = "constant", rule = 2)
  )
}

setMethod("totalEffect", class_name, function(this, age){
  this@effectFunction(age)
})
