setOldClass("loess")

#' Models the relationship between age and time-dependent genetic effects as
#' estimated through a LOESS model
#'
#' @slot model The LOESS model of the the PGS ~ age relationship
#' @slot stratifiedModels Age-stratified [`glm`][stats::glm()] objects used to
#'   estimate genetic effects.
#'
#' @seealso [GlmTimeDependentModel()], [time_dependent_loess()], [stats::loess()]
#' @include TimeDependentModel-class.R
#' @export
LoessTimeDependentModel <- setClass(
  "LoessTimeDependentModel",
  contains = "TimeDependentModel",
  slots = c(
    model = "loess",
    stratifiedModels = "list"
  )
)

setMethod("totalEffect", "LoessTimeDependentModel", function(this, age){
  predict(model(this), age, se = TRUE)$fit
})

#' Named list of age-stratified models used to estimate effect strength for the
#' LOESS model
#' @keywords internal
#' @rdname LoessTimeDependentModel-methods
#' @export
setGeneric("stratifiedModels", function(this) standardGeneric("stratifiedModels"))
setMethod("stratifiedModels", "LoessTimeDependentModel", function(this) {
  this@stratifiedModels
})
