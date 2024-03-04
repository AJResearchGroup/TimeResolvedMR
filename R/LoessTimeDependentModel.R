class_name <- "LoessTimeDependentModel"

#' Models the relationship between age and time-dependent genetic effects as
#' estimated through a LOESS model
#'
#' @slot model The LOESS model of the the PGS ~ age relationship
#' @slot stratifiedModels Age-stratified [`glm`][stats::glm()] objects used to
#'   estimate genetic effects.
#'
#' @seealso [GlmTimeDependentModel()], [time_dependent_loess()], [stats::loess()]
#'
#' @export
LoessTimeDependentModel <- setClass(
  class_name,
  contains = "TimeDependentModel",
  slots = c(
    model = "loess",
    stratifiedModels = "list"
  )
)

setMethod("totalEffect", class_name, function(this, age){
  predict(model(this), age, se = TRUE)$fit
})

setGeneric("stratifiedModels", class_name, function(this) standardGeneric("stratifiedModels"))
setMethod("stratifiedModels", class_name, function(this) {
  this@stratifiedModels
})
